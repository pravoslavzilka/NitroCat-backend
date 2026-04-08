import torch, pickle, json, ast, numpy as np, sys, os
sys.path.insert(0, '/content/clipzyme')
os.environ['MPLBACKEND'] = 'Agg'

from clipzyme import CLIPZyme
from clipzyme.utils.pyg import from_mapped_smiles
from clipzyme.utils.wln_processing import get_bond_changes
from torch_geometric.data import Batch

ORIGINAL_CKPT  = 'files/clipzyme_model.ckpt'
FINETUNED_CKPT = 'files/models/model_random_split.ckpt'
SCREENING_SET  = 'files/cyp_screening_set.p'
MY_DATA        = 'files/my_data.json'

# ── Load screening set ────────────────────────────────────────────
screening  = pickle.load(open(SCREENING_SET, 'rb'))
screen_emb = screening['hiddens'].float()
screen_ids = list(screening['uniprots'])
print(f"Screening set: {len(screen_ids)} enzymes")

# ── Load models ───────────────────────────────────────────────────
print("Loading original model...")
model_orig = CLIPZyme(checkpoint_path=ORIGINAL_CKPT).eval()

print("Loading fine-tuned model...")
model_ft = CLIPZyme(checkpoint_path=ORIGINAL_CKPT)
ft_ckpt  = torch.load(FINETUNED_CKPT, map_location='cpu')
model_ft.load_state_dict(ft_ckpt['state_dict'], strict=False)
model_ft = model_ft.eval()
print(f"Fine-tuned: epoch {ft_ckpt.get('epoch')}, step {ft_ckpt.get('global_step')}")

# ── Load test data ────────────────────────────────────────────────
data      = json.load(open(MY_DATA))
test_data = [r for r in data if r.get('split') == 'test']
print(f"Test samples: {len(test_data)}")

# ── Encode reaction directly via WLN ─────────────────────────────
def encode_reaction_direct(model, r):
    """Encode a single reaction using the model's reaction encoder directly."""
    reactants = sorted(r['mapped_reactants'])
    products  = sorted(r['mapped_products'])

    # Build molecule graphs
    reactant_graph, atom_map2new_index = from_mapped_smiles(
        '.'.join(reactants), encode_no_edge=True
    )
    product_graph, _ = from_mapped_smiles(
        '.'.join(products), encode_no_edge=True
    )

    # Bond changes
    bond_changes = []
    for u, v, btype in r.get('bond_changes', []):
        u_idx = atom_map2new_index.get(int(u))
        v_idx = atom_map2new_index.get(int(v))
        if u_idx is not None and v_idx is not None:
            bond_changes.append((u_idx, v_idx, btype))
    if not bond_changes:
        bond_changes = [(0, 1, 1.0)]
    reactant_graph.bond_changes = bond_changes

    # Add batch index
    reactant_graph.batch = torch.zeros(reactant_graph.x.shape[0], dtype=torch.long)
    product_graph.batch  = torch.zeros(product_graph.x.shape[0],  dtype=torch.long)

    batch = {
        'reactants': reactant_graph,
        'products':  product_graph,
    }

    with torch.no_grad():
        emb = model.model.encode_reaction(batch)  # [1, 1280]
    return emb.squeeze(0)

# ── Evaluate ──────────────────────────────────────────────────────
def compute_metrics(model, test_data, screen_emb, screen_ids, label):
    screen_norm   = screen_emb / screen_emb.norm(dim=1, keepdim=True)
    screen_id_map = {uid: i for i, uid in enumerate(screen_ids)}

    ranks     = []
    in_top20  = []
    not_found = 0
    failed    = 0

    for i, r in enumerate(test_data):
        uid = ast.literal_eval(r['protein_refs'])[0]
        if uid not in screen_id_map:
            not_found += 1
            continue
        try:
            rxn_emb  = encode_reaction_direct(model, r)
            rxn_norm = rxn_emb / rxn_emb.norm()
            scores   = (screen_norm @ rxn_norm).numpy()
            sorted_idx = np.argsort(-scores)
            rank = int(np.where(sorted_idx == screen_id_map[uid])[0][0]) + 1
            ranks.append(rank)
            in_top20.append(rank <= 20)
        except Exception as e:
            failed += 1

        if (i+1) % 50 == 0:
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
        'mean_rank':  np.mean(ranks),
        'median_rank': np.median(ranks),
        'top1':  np.mean([r==1  for r in ranks])*100,
        'top5':  np.mean([r<=5  for r in ranks])*100,
        'top10': np.mean([r<=10 for r in ranks])*100,
        'top20': np.mean([r<=20 for r in ranks])*100,
        'ranks': ranks,
    }

metrics_orig = compute_metrics(model_orig, test_data, screen_emb, screen_ids, "ORIGINAL")
metrics_ft   = compute_metrics(model_ft,   test_data, screen_emb, screen_ids, "FINE-TUNED")

print("\n=== COMPARISON ===")
print(f"{'Metric':<20} {'Original':>12} {'Fine-tuned':>12} {'Delta':>10}")
print("-" * 56)
for key in ['mean_rank', 'median_rank', 'top1', 'top5', 'top10', 'top20']:
    orig_val = metrics_orig[key]
    ft_val   = metrics_ft[key]
    delta    = ft_val - orig_val
    better   = delta < 0 if 'rank' in key else delta > 0
    arrow    = '↑' if better else '↓'
    print(f"  {key:<18} {orig_val:>12.1f} {ft_val:>12.1f} {delta:>+9.1f} {arrow}")