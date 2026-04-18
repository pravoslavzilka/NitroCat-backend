import torch, pickle, json, ast, numpy as np, sys, os
sys.path.insert(0, '/content/clipzyme')
os.environ['MPLBACKEND'] = 'Agg'

from clipzyme import CLIPZyme
from clipzyme.utils.pyg import from_mapped_smiles
from clipzyme.utils.wln_processing import get_bond_changes
from torch_geometric.data import Batch

ORIGINAL_CKPT   = 'files/clipzyme_model.ckpt'
FINETUNED_CKPT1 = 'files/models/model_random_split.ckpt'
FINETUNED_CKPT2 = 'files/models/last.ckpt'          
SCREENING_SET   = 'files/cyp_screening_set.p'
MY_DATA         = 'files/disjoint/my_data.json'

# ── Load screening set ────────────────────────────────────────────
screening  = pickle.load(open(SCREENING_SET, 'rb'))
screen_emb = screening['hiddens'].float()
screen_ids = list(screening['uniprots'])
print(f"Screening set: {len(screen_ids)} enzymes")

# ── Helper to load a fine-tuned checkpoint ────────────────────────
def load_finetuned(ckpt_path: str, label: str):
    print(f"Loading {label}...")
    model = CLIPZyme(checkpoint_path=ORIGINAL_CKPT)
    ckpt  = torch.load(ckpt_path, map_location='cpu')
    model.load_state_dict(ckpt['state_dict'], strict=False)
    model = model.eval()
    print(f"  {label}: epoch {ckpt.get('epoch')}, step {ckpt.get('global_step')}")
    return model

# ── Load models ───────────────────────────────────────────────────
print("Loading original model...")
model_orig = CLIPZyme(checkpoint_path=ORIGINAL_CKPT).eval()

model_ft1 = load_finetuned(FINETUNED_CKPT1, "Fine-tuned 1 (random split)")
model_ft2 = load_finetuned(FINETUNED_CKPT2, "Fine-tuned 2")

# ── Load test data ────────────────────────────────────────────────
data      = json.load(open(MY_DATA))
test_data = [r for r in data if r.get('split') == 'test']
print(f"Test samples: {len(test_data)}")

# ── Encode reaction directly via WLN ─────────────────────────────
def encode_reaction_direct(model, r):
    """Encode a single reaction using the model's reaction encoder directly."""
    reactants = sorted(r['mapped_reactants'])
    products  = sorted(r['mapped_products'])

    reactant_graph, atom_map2new_index = from_mapped_smiles(
        '.'.join(reactants), encode_no_edge=True
    )
    product_graph, _ = from_mapped_smiles(
        '.'.join(products), encode_no_edge=True
    )

    bond_changes = []
    for u, v, btype in r.get('bond_changes', []):
        u_idx = atom_map2new_index.get(int(u))
        v_idx = atom_map2new_index.get(int(v))
        if u_idx is not None and v_idx is not None:
            bond_changes.append((u_idx, v_idx, btype))
    if not bond_changes:
        bond_changes = [(0, 1, 1.0)]
    reactant_graph.bond_changes = bond_changes

    reactant_graph.batch = torch.zeros(reactant_graph.x.shape[0], dtype=torch.long)
    product_graph.batch  = torch.zeros(product_graph.x.shape[0],  dtype=torch.long)

    batch = {
        'reactants': reactant_graph,
        'products':  product_graph,
    }

    with torch.no_grad():
        emb = model.model.encode_reaction(batch)
    return emb.squeeze(0)

# ── Evaluate ──────────────────────────────────────────────────────
def compute_metrics(model, test_data, screen_emb, screen_ids, label):
    screen_norm   = screen_emb / screen_emb.norm(dim=1, keepdim=True)
    screen_id_map = {uid: i for i, uid in enumerate(screen_ids)}

    ranks     = []
    not_found = 0
    failed    = 0

    for i, r in enumerate(test_data):
        uid = ast.literal_eval(r['protein_refs'])[0]
        if uid not in screen_id_map:
            not_found += 1
            continue
        try:
            rxn_emb    = encode_reaction_direct(model, r)
            rxn_norm   = rxn_emb / rxn_emb.norm()
            scores     = (screen_norm @ rxn_norm).numpy()
            sorted_idx = np.argsort(-scores)
            rank       = int(np.where(sorted_idx == screen_id_map[uid])[0][0]) + 1
            ranks.append(rank)
        except Exception as e:
            failed += 1

        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{len(test_data)}", end='\r')

    print(f"\n=== {label} ===")
    print(f"  Evaluated:    {len(ranks)}")
    print(f"  Not in set:   {not_found}")
    print(f"  Failed:       {failed}")
    print(f"  Mean rank:    {np.mean(ranks):.1f}")
    print(f"  Median rank:  {np.median(ranks):.1f}")
    print(f"  Top-1  hit %: {np.mean([r==1  for r in ranks])*100:.1f}%")
    print(f"  Top-5  hit %: {np.mean([r<=5  for r in ranks])*100:.1f}%")
    print(f"  Top-10 hit %: {np.mean([r<=10 for r in ranks])*100:.1f}%")
    print(f"  Top-20 hit %: {np.mean([r<=20 for r in ranks])*100:.1f}%")

    return {
        'label':       label,
        'mean_rank':   float(np.mean(ranks)),
        'median_rank': float(np.median(ranks)),
        'top1':        float(np.mean([r==1  for r in ranks])*100),
        'top5':        float(np.mean([r<=5  for r in ranks])*100),
        'top10':       float(np.mean([r<=10 for r in ranks])*100),
        'top20':       float(np.mean([r<=20 for r in ranks])*100),
        'ranks':       ranks,
    }

# ── Run all three ─────────────────────────────────────────────────
models = [
    (model_orig, "ORIGINAL"),
    (model_ft1,  "FINE-TUNED 1 (random split)"),
    (model_ft2,  "FINE-TUNED 2"),
]

all_metrics = [
    compute_metrics(model, test_data, screen_emb, screen_ids, label)
    for model, label in models
]

# ── Comparison table ──────────────────────────────────────────────
METRICS  = ['mean_rank', 'median_rank', 'top1', 'top5', 'top10', 'top20']
col_w    = 18
label_w  = 20
total_w  = label_w + col_w * len(all_metrics) + 14

print("\n" + "=" * total_w)
print("COMPARISON")
print("=" * total_w)

# Header
header = f"  {'Metric':<{label_w}}"
for m in all_metrics:
    header += f"{m['label'][:col_w-1]:>{col_w}}"
print(header)
print("-" * total_w)

# Rows
for key in METRICS:
    row    = f"  {key:<{label_w}}"
    values = [m[key] for m in all_metrics]

    for val in values:
        is_best = (val == min(values)) if 'rank' in key else (val == max(values))
        marker  = "*" if is_best else " "
        row    += f"{val:>{col_w-1}.1f}{marker}"

    # Delta: each fine-tuned vs original
    orig = all_metrics[0][key]
    deltas = [m[key] - orig for m in all_metrics[1:]]
    delta_str = "  " + "  ".join(
        f"{d:>+6.1f}{'↑' if (d<0 if 'rank' in key else d>0) else ('↓' if d!=0 else '→')}"
        for d in deltas
    )
    print(row + delta_str)

print("=" * total_w)
print("* = best value   |   Deltas vs ORIGINAL: FT1, FT2")