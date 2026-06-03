import torch, pickle, json, ast, numpy as np, gc, os, traceback
from collections import Counter
os.environ['MPLBACKEND'] = 'Agg'

from clipzyme import CLIPZyme
from clipzyme.utils.pyg import from_mapped_smiles

# ── Checkpoints ───────────────────────────────────────────────────
ORIGINAL_CKPT     = 'files/clipzyme_model.ckpt'
CKPT_RANDOM_BS8   = 'files/models/disjoint_v2_random.ckpt'
CKPT_COBACO_BS8   = 'files/models/disjoint_v2_cobaco.ckpt'

# Cleaned dataset (43 records with empty mapped_reactants/products are dropped)
MY_DATA = 'data_processing/my_data_processed.json'

# Each model paired with the screening set encoded by its own protein encoder.
# Mixing model and screening across pairs invalidates the cosine similarity.
MODELS = [
    ('Original (pretrained)', ORIGINAL_CKPT,   'files/screening/clipzyme_default.p'),
    ('Disjoint random bs8',   CKPT_RANDOM_BS8, 'files/screening/disjoint_v2_random.p'),
    ('Disjoint CoBaCo bs8',   CKPT_COBACO_BS8, 'files/screening/disjoint_v2_cobaco.p'),
]


def load_screening(path):
    s = pickle.load(open(path, 'rb'))
    return s['hiddens'].float(), list(s['uniprots'])


def load_model(label, ckpt_path):
    print(f"Loading {label} ...")
    if ckpt_path == ORIGINAL_CKPT:
        model = CLIPZyme(checkpoint_path=ckpt_path).eval()
    else:
        model = CLIPZyme(checkpoint_path=ORIGINAL_CKPT)
        ckpt  = torch.load(ckpt_path, map_location='cpu')
        model.load_state_dict(ckpt['state_dict'], strict=False)
        model = model.eval()
        print(f"  epoch={ckpt.get('epoch')}, step={ckpt.get('global_step')}")
    return model


def encode_reaction(model, r):
    rg, atom_map = from_mapped_smiles(
        '.'.join(sorted(r['mapped_reactants'])), encode_no_edge=True
    )
    pg, _ = from_mapped_smiles(
        '.'.join(sorted(r['mapped_products'])), encode_no_edge=True
    )
    bond_changes = [
        (atom_map[int(u)], atom_map[int(v)], t)
        for u, v, t in r.get('bond_changes', [])
        if atom_map.get(int(u)) is not None and atom_map.get(int(v)) is not None
    ] or [(0, 1, 1.0)]
    rg.bond_changes = bond_changes
    rg.batch = torch.zeros(rg.x.shape[0], dtype=torch.long)
    pg.batch = torch.zeros(pg.x.shape[0], dtype=torch.long)
    with torch.no_grad():
        emb = model.model.encode_reaction({'reactants': rg, 'products': pg})
    return emb.squeeze(0)


def parse_uid(refs_field):
    """protein_refs may be a stringified list, a real list, or a bare string."""
    if isinstance(refs_field, list):
        refs = refs_field
    elif isinstance(refs_field, str):
        s = refs_field.strip()
        if s.startswith('[') and s.endswith(']'):
            try:
                refs = ast.literal_eval(s)
            except Exception:
                refs = []
        else:
            refs = [s] if s else []
    else:
        refs = []
    return refs[0] if refs else None


def compute_metrics(model, test_data, label, screening_path):
    screen_emb, screen_ids = load_screening(screening_path)
    print(f"  Screening: {screening_path} ({len(screen_ids)} enzymes)")
    screen_norm   = screen_emb / screen_emb.norm(dim=1, keepdim=True)
    screen_id_map = {uid: i for i, uid in enumerate(screen_ids)}

    ranks = []
    skipped_no_uid     = 0  # protein_refs missing or unparseable
    skipped_not_in_set = 0  # uid not in screening
    skipped_empty_smi  = 0  # empty mapped_reactants/products
    failed_encode      = 0  # encoding/forward error
    fail_reasons       = Counter()

    for i, r in enumerate(test_data):
        uid = parse_uid(r.get('protein_refs'))
        if not uid:
            skipped_no_uid += 1
            continue
        if uid not in screen_id_map:
            skipped_not_in_set += 1
            continue
        if not r.get('mapped_reactants') or not r.get('mapped_products'):
            skipped_empty_smi += 1
            continue
        try:
            rxn_emb    = encode_reaction(model, r)
            rxn_norm   = rxn_emb / rxn_emb.norm()
            scores     = (screen_norm @ rxn_norm).numpy()
            sorted_idx = np.argsort(-scores)
            rank       = int(np.where(sorted_idx == screen_id_map[uid])[0][0]) + 1
            ranks.append(rank)
        except Exception as e:
            failed_encode += 1
            reason = f"{type(e).__name__}: {str(e)[:80]}"
            fail_reasons[reason] += 1
            if failed_encode <= 3:
                print(f"\n  [first failures] rxnid={r.get('rxnid')}: {reason}")

        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{len(test_data)}", end='\r')

    total = len(test_data)
    skipped = skipped_no_uid + skipped_not_in_set + skipped_empty_smi + failed_encode
    print(f"\n=== {label} ===")
    print(f"  Total test:      {total}")
    print(f"  Evaluated:       {len(ranks)} ({100*len(ranks)/total:.1f}%)")
    print(f"  Skipped total:   {skipped}")
    print(f"    no uniprot id: {skipped_no_uid}")
    print(f"    uid not in screening set: {skipped_not_in_set}")
    print(f"    empty SMILES:  {skipped_empty_smi}")
    print(f"    encode failed: {failed_encode}")
    if fail_reasons:
        print("  Top failure reasons:")
        for reason, n in fail_reasons.most_common(5):
            print(f"    [{n}x] {reason}")

    if ranks:
        print(f"  Mean rank:    {np.mean(ranks):.1f}")
        print(f"  Median rank:  {np.median(ranks):.1f}")
        print(f"  Top-1  hit %: {np.mean([r==1  for r in ranks])*100:.1f}%")
        print(f"  Top-5  hit %: {np.mean([r<=5  for r in ranks])*100:.1f}%")
        print(f"  Top-10 hit %: {np.mean([r<=10 for r in ranks])*100:.1f}%")
        print(f"  Top-20 hit %: {np.mean([r<=20 for r in ranks])*100:.1f}%")
        print(f"  Top-90 hit %: {np.mean([r<=90 for r in ranks])*100:.1f}%")

    return {
        'label':       label,
        'evaluated':   len(ranks),
        'mean_rank':   float(np.mean(ranks))   if ranks else float('nan'),
        'median_rank': float(np.median(ranks)) if ranks else float('nan'),
        'top1':  float(np.mean([r==1  for r in ranks])*100) if ranks else float('nan'),
        'top5':  float(np.mean([r<=5  for r in ranks])*100) if ranks else float('nan'),
        'top10': float(np.mean([r<=10 for r in ranks])*100) if ranks else float('nan'),
        'top20': float(np.mean([r<=20 for r in ranks])*100) if ranks else float('nan'),
        'top90': float(np.mean([r<=90 for r in ranks])*100) if ranks else float('nan'),
        'ranks': ranks,
    }


# ── Load test data ────────────────────────────────────────────────
data      = json.load(open(MY_DATA))
test_data = [r for r in data
             if r.get('split') == 'test' and not r.get('skip', False)]
print(f"Dataset: {MY_DATA}")
print(f"Test samples: {len(test_data)}")

# ── Evaluate each model independently ────────────────────────────
all_metrics = []
for label, ckpt_path, screening_path in MODELS:
    if not os.path.exists(ckpt_path):
        print(f"!! Skipping {label} — checkpoint not found: {ckpt_path}")
        continue
    if not os.path.exists(screening_path):
        print(f"!! Skipping {label} — screening set not found: {screening_path}")
        continue
    model   = load_model(label, ckpt_path)
    metrics = compute_metrics(model, test_data, label, screening_path)
    all_metrics.append(metrics)
    del model
    gc.collect()

# ── Comparison table ──────────────────────────────────────────────
METRICS = ['mean_rank', 'median_rank', 'top1', 'top5', 'top10', 'top20', 'top90']
col_w   = 16
label_w = 22
total_w = label_w + col_w * len(all_metrics) + 14

print("\n" + "=" * total_w)
print("COMPARISON")
print("=" * total_w)
header = f"  {'Metric':<{label_w}}"
for m in all_metrics:
    header += f"{m['label'][:col_w-1]:>{col_w}}"
print(header)
print("-" * total_w)

for key in METRICS:
    row    = f"  {key:<{label_w}}"
    values = [m[key] for m in all_metrics]
    for val in values:
        is_best = (val == min(v for v in values if not np.isnan(v))) if 'rank' in key \
             else (val == max(v for v in values if not np.isnan(v)))
        row += f"{val:>{col_w-1}.1f}{'*' if is_best else ' '}"
    orig      = all_metrics[0][key] if all_metrics else float('nan')
    deltas    = [m[key] - orig for m in all_metrics[1:]]
    delta_str = "  " + "  ".join(
        f"{d:>+5.1f}{'↑' if (d<0 if 'rank' in key else d>0) else ('↓' if d!=0 else '→')}"
        for d in deltas
    )
    print(row + delta_str)

print("=" * total_w)
print("* = best   |   Deltas vs Original")
