import os, json, pickle, ast, time, gc, ast as _ast
import numpy as np
import requests
import torch

os.environ['MPLBACKEND'] = 'Agg'

# Reuse the exact pipeline main.py / query.py uses: model, screening set,
# RXNMapper, and the balancing logic. Importing query.py loads the model.
from app.query import model, candidates, candidate_ids, mapper, _balance

# ── Config ────────────────────────────────────────────────────────────────────
TEST_DATA      = 'data_processing/my_data.json'
EC2UNIPROT     = 'data_processing/ec2uniprot.p'
BRENDA_JSON    = 'files/brenda_2026_1.json'
CLAIRE_URL       = 'http://127.0.0.1:4550/predict'
CLAIRE_BATCH_URL = 'http://127.0.0.1:4550/predict_batch'
CLAIRE_TIMEOUT   = 600
CLAIRE_BATCH_SZ  = 200
CLAIRE_TOPK      = int(os.environ.get('CLAIRE_TOPK', '5'))


# ── Helpers ───────────────────────────────────────────────────────────────────
def parse_protein_refs(refs):
    if isinstance(refs, list):
        return refs
    if isinstance(refs, str) and refs.strip():
        try:
            return _ast.literal_eval(refs)
        except Exception:
            return []
    return []


def ec_prefix2(ec):
    parts = ec.split('.')
    return '.'.join(parts[:2]) if len(parts) >= 2 else None


def build_uniprot2ec():
    """uniprot → set of ECs from BRENDA + ec2uniprot.p + my_data.json."""
    u2e = {}
    print("  Loading BRENDA ...")
    t = time.time()
    b = json.load(open(BRENDA_JSON))
    for ec, e in b['data'].items():
        for _, rec in e.get('protein', {}).items():
            if isinstance(rec, dict):
                for acc in rec.get('accessions') or []:
                    u2e.setdefault(acc, set()).add(ec)
    print(f"    {time.time()-t:.1f}s — {len(u2e)} uniprots from BRENDA")

    ec2u = pickle.load(open(EC2UNIPROT, 'rb'))
    for ec, uids in ec2u.items():
        for uid in uids:
            u2e.setdefault(uid, set()).add(ec)

    for r in json.load(open(TEST_DATA)):
        ec = r.get('ec')
        if not ec:
            continue
        for uid in parse_protein_refs(r.get('protein_refs')):
            u2e.setdefault(uid, set()).add(ec)

    return u2e


def encode_reaction_via_pipeline(substrate, product):
    """Mirror query.py's pipeline: balance → atom-map → extract_reaction_features.
    Returns (rxn_vec, comments) or (None, error_str)."""
    balanced, _ = _balance(substrate, product)
    try:
        out = mapper.get_attention_guided_atom_maps([balanced])[0]
    except Exception as e:
        return None, f"rxnmapper failed: {e}"
    mapped, conf = out['mapped_rxn'], float(out['confidence'])
    if ':1]' not in mapped or conf < 0.1:
        return None, f"low/no mapping (conf={conf:.2f})"
    try:
        with torch.no_grad():
            v = model.extract_reaction_features(reaction=mapped)
    except Exception as e:
        return None, f"encode failed: {e}"
    return v, None


def claire_ec2_batch(rxns, topk=CLAIRE_TOPK):
    """Return a list of EC2-prefix sets (one set per reaction). Empty set on failure."""
    out = [set() for _ in rxns]
    for i in range(0, len(rxns), CLAIRE_BATCH_SZ):
        chunk = rxns[i:i + CLAIRE_BATCH_SZ]
        try:
            r = requests.post(CLAIRE_BATCH_URL,
                              json={'smiles': chunk, 'topk': topk},
                              timeout=CLAIRE_TIMEOUT)
            r.raise_for_status()
            for j, preds in enumerate(r.json()):
                out[i + j] = {p['ec2'] for p in preds if p.get('ec2')}
        except Exception as e:
            print(f"  CLAIRE batch [{i}:{i+len(chunk)}] failed: {e}")
    return out


def rank_in(scores_subset, ids_subset, true_uid):
    if true_uid not in ids_subset:
        return None
    order = np.argsort(-scores_subset)
    pos   = ids_subset.index(true_uid)
    return int(np.where(order == pos)[0][0]) + 1


def summarise(label, ranks, not_in_set, failed, no_ec=0):
    print(f"\n=== {label} ===")
    print(f"  Evaluated: {len(ranks)} | Not in set: {not_in_set} | Failed: {failed} | No EC: {no_ec}")
    if not ranks:
        return {'label': label}
    arr = np.asarray(ranks)
    out = {
        'label':       label,
        'n':           len(ranks),
        'mean_rank':   float(arr.mean()),
        'median_rank': float(np.median(arr)),
        'top1':        float((arr == 1).mean() * 100),
        'top5':        float((arr <= 5).mean() * 100),
        'top10':       float((arr <= 10).mean() * 100),
        'top20':       float((arr <= 20).mean() * 100),
        'top50':       float((arr <= 50).mean() * 100),
    }
    print(f"  Mean rank:    {out['mean_rank']:.1f}")
    print(f"  Median rank:  {out['median_rank']:.1f}")
    print(f"  Top-1  hit %: {out['top1']:.1f}%")
    print(f"  Top-5  hit %: {out['top5']:.1f}%")
    print(f"  Top-10 hit %: {out['top10']:.1f}%")
    print(f"  Top-20 hit %: {out['top20']:.1f}%")
    print(f"  Top-50 hit %: {out['top50']:.1f}%")
    return out


# ── Setup ─────────────────────────────────────────────────────────────────────
print(f"Screening set: {len(candidate_ids):,} enzymes")
candidate_idset = set(candidate_ids)
candidate_pos   = {uid: i for i, uid in enumerate(candidate_ids)}
cand_t = candidates.float()  # already normalised by extract_reaction_features

print("Building uniprot → EC index ...")
u2e = build_uniprot2ec()
print(f"  total uniprots: {len(u2e):,}")

# Per-prefix index into the screening tensor
prefix2idx = {}
for i, uid in enumerate(candidate_ids):
    for ec in u2e.get(uid, ()):
        p = ec_prefix2(ec)
        if p:
            prefix2idx.setdefault(p, []).append(i)
prefix2idx = {p: np.asarray(sorted(set(idx))) for p, idx in prefix2idx.items()}
print(f"  EC2 prefixes covered in screening: {len(prefix2idx)} "
      f"({sum(len(v) for v in prefix2idx.values()):,} (uid,prefix) entries)")

test_data = [r for r in json.load(open(TEST_DATA)) if r.get('split') == 'test']
limit = int(os.environ.get('BENCHCLAIRE_LIMIT', '0'))
if limit:
    test_data = test_data[:limit]
    print(f"Test samples (LIMIT={limit}): {len(test_data)}")
else:
    print(f"Test samples: {len(test_data)}")

# Pre-batch all CLAIRE predictions in one go (≈80× faster than per-call).
print(f"Pre-computing CLAIRE EC2 top-{CLAIRE_TOPK} predictions for {len(test_data)} reactions ...")
t = time.time()
all_rxns = [
    '.'.join(r.get('reactants', [])) + '>>' + '.'.join(r.get('products', []))
    for r in test_data
]
ec2_pred = claire_ec2_batch(all_rxns, topk=CLAIRE_TOPK)
have = sum(1 for s in ec2_pred if s)
print(f"  done in {time.time()-t:.1f}s "
      f"({have:,}/{len(ec2_pred)} got predictions, "
      f"avg {sum(len(s) for s in ec2_pred)/max(have,1):.2f} families/query)")


# ── Evaluate ──────────────────────────────────────────────────────────────────
ranks_full,     not_in_full,     failed_full     = [], 0, 0
ranks_filtered, not_in_filtered, failed_filtered, no_ec = [], 0, 0, 0

t_start = time.time()
for i, r in enumerate(test_data):
    uids = parse_protein_refs(r.get('protein_refs'))
    true_uid = uids[0] if uids else None
    if not true_uid or true_uid not in candidate_idset:
        not_in_full += 1
        not_in_filtered += 1
        continue

    sub  = '.'.join(r.get('reactants', []))
    prod = '.'.join(r.get('products', []))
    if not sub or not prod:
        failed_full += 1
        failed_filtered += 1
        continue

    rxn_vec, err = encode_reaction_via_pipeline(sub, prod)
    if rxn_vec is None:
        failed_full += 1
        failed_filtered += 1
        continue

    scores = (cand_t @ rxn_vec.T).squeeze().cpu().numpy()

    # Option 1 — full screening set
    order = np.argsort(-scores)
    pos   = candidate_pos[true_uid]
    ranks_full.append(int(np.where(order == pos)[0][0]) + 1)

    # Option 2 — filter by union of CLAIRE-predicted EC2 prefixes (top-K)
    ec2_set = ec2_pred[i]
    if not ec2_set:
        no_ec += 1
    else:
        # Union of candidate indices across all predicted families
        idx_lists = [prefix2idx[p] for p in ec2_set if p in prefix2idx]
        if not idx_lists:
            failed_filtered += 1
        else:
            idx        = np.unique(np.concatenate(idx_lists))
            sub_scores = scores[idx]
            sub_ids    = [candidate_ids[j] for j in idx]
            sr = rank_in(sub_scores, sub_ids, true_uid)
            if sr is None:
                not_in_filtered += 1
            else:
                ranks_filtered.append(sr)

    if (i + 1) % 25 == 0:
        elapsed = time.time() - t_start
        rate    = (i + 1) / elapsed
        eta     = (len(test_data) - (i + 1)) / rate
        print(f"  {i+1}/{len(test_data)}  "
              f"({rate:.1f}/s, ETA {eta/60:.1f}min)  "
              f"full ok={len(ranks_full)} fail={failed_full}  "
              f"filt ok={len(ranks_filtered)} fail={failed_filtered}")

# ── Report ────────────────────────────────────────────────────────────────────
m_full = summarise('CLIPZyme — full screening set',     ranks_full,     not_in_full,     failed_full)
m_filt = summarise('CLIPZyme + CLAIRE EC2 filter',      ranks_filtered, not_in_filtered, failed_filtered, no_ec)

METRICS = ['n', 'mean_rank', 'median_rank', 'top1', 'top5', 'top10', 'top20', 'top50']
if all(k in m_full for k in METRICS) and all(k in m_filt for k in METRICS):
    print("\n" + "=" * 70)
    print(f"  {'Metric':<14}{'Full':>16}{'CLAIRE-filtered':>20}{'Δ':>14}")
    print("-" * 70)
    for k in METRICS:
        a, b = m_full[k], m_filt[k]
        d = b - a
        arrow = '↑' if (d > 0 and 'rank' not in k) or (d < 0 and 'rank' in k) else (
                '↓' if d != 0 else '→')
        print(f"  {k:<14}{a:>16.2f}{b:>20.2f}{d:>+12.2f} {arrow}")
    print("=" * 70)

del model
gc.collect()
