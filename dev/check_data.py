# validate_for_clipzyme.py
# Run locally: python validate_for_clipzyme.py files/disjoint/my_data_ready.json

import json, re, ast, sys
from collections import Counter

FILE = sys.argv[1] if len(sys.argv) > 1 else 'files/disjoint/my_data_final.json'
data = json.load(open(FILE))
print(f"Validating {len(FILE)} — {len(data)} records\n")

def has_atom_map(smi):
    return bool(re.search(r':\d+', smi))

issues     = []   # (record_index, uid, issue_description)
will_skip  = []   # records CLIPZyme will silently skip
will_crash = []   # records that will cause a crash

for i, r in enumerate(data):
    record_issues = []

    # ── Get UID ───────────────────────────────────────────────────
    try:
        pr = r['protein_refs']
        uid = ast.literal_eval(pr)[0] if isinstance(pr, str) else pr[0]
    except Exception as e:
        uid = '???'
        record_issues.append(f"protein_refs unparseable: {e}")

    # ── Fields that cause skip_sample() to return True ────────────
    # 1. EC missing or contains '-'
    ec = r.get('ec')
    if not ec or ec is None:
        record_issues.append("ec is None → skip_sample")
    elif '-' in str(ec):
        record_issues.append(f"ec='{ec}' contains '-' → skip_sample")

    # 2. protein_refs not a string
    if not isinstance(r.get('protein_refs'), str):
        record_issues.append(f"protein_refs is {type(r.get('protein_refs')).__name__} not str → skip_sample")

    # 3. protein_db invalid
    if r.get('protein_db') not in ['uniprot', 'swissprot']:
        record_issues.append(f"protein_db='{r.get('protein_db')}' invalid → skip_sample")

    # 4. split invalid
    if r.get('split') not in ['train', 'dev', 'test']:
        record_issues.append(f"split='{r.get('split')}' invalid → skip_sample")

    # ── Fields that cause __getitem__ crash ───────────────────────
    # 5. Unmapped SMILES → "SMILES must contain atom map numbers"
    mapped_r = r.get('mapped_reactants', [])
    mapped_p = r.get('mapped_products', [])
    if any(not has_atom_map(s) for s in mapped_r + mapped_p):
        record_issues.append("unmapped SMILES → getitem crash (skipped by collate)")

    # 6. Empty bond_changes
    if not r.get('bond_changes'):
        record_issues.append("empty bond_changes → getitem crash")

    # 7. rule_id missing
    if r.get('rule_id') is None:
        record_issues.append("rule_id missing → all samples same split")

    # 8. mapped_reactants/products empty
    if not mapped_r:
        record_issues.append("mapped_reactants empty → crash")
    if not mapped_p:
        record_issues.append("mapped_products empty → crash")

    # ── Classify ──────────────────────────────────────────────────
    if record_issues:
        issues.append((i, uid, record_issues))
        skip_triggers = ['skip_sample', 'protein_refs', 'protein_db', 'split', 'ec']
        if any(any(t in issue for t in skip_triggers) for issue in record_issues):
            will_skip.append(i)
        else:
            will_crash.append(i)

# ── Summary ───────────────────────────────────────────────────────
print(f"{'='*60}")
print(f"SUMMARY")
print(f"{'='*60}")
print(f"Total records:          {len(data)}")
print(f"Clean records:          {len(data) - len(set(will_skip + will_crash))}")
print(f"Will be skipped:        {len(set(will_skip))}  ← skip_sample() returns True")
print(f"Will crash getitem:     {len(set(will_crash))}  ← returns None, discarded by collate")
print(f"Effective training set: {len(data) - len(issues)} records")

# ── Issue breakdown ───────────────────────────────────────────────
all_issue_types = Counter()
for _, _, record_issues in issues:
    for issue in record_issues:
        # Extract just the issue type
        key = issue.split('→')[0].strip().split('=')[0].strip()
        all_issue_types[key] += 1

print(f"\n── Issue breakdown ──────────────────────────────────────────")
for issue_type, count in all_issue_types.most_common():
    print(f"  {count:>5}x  {issue_type}")

# ── Show first 10 problematic records ────────────────────────────
print(f"\n── First 10 problematic records ─────────────────────────────")
for i, uid, record_issues in issues[:10]:
    print(f"\n  [{i}] uid={uid} split={data[i].get('split')}")
    for issue in record_issues:
        print(f"       ❌ {issue}")

# ── Per-split breakdown ───────────────────────────────────────────
print(f"\n── Per-split effective records ──────────────────────────────")
problem_idx = set(i for i, _, _ in issues)
for split in ['train', 'dev', 'test']:
    total   = sum(1 for r in data if r.get('split') == split)
    problem = sum(1 for i, r in enumerate(data)
                  if r.get('split') == split and i in problem_idx)
    clean   = total - problem
    print(f"  {split:<6}: {clean}/{total} clean ({clean/total*100:.1f}% usable)")

# ── Final verdict ─────────────────────────────────────────────────
print(f"\n{'='*60}")
usable = len(data) - len(issues)
if usable == len(data):
    print("✅ ALL RECORDS VALID — ready for training")
elif usable >= len(data) * 0.8:
    print(f"⚠️  {usable}/{len(data)} records usable ({usable/len(data)*100:.1f}%)")
    print("   Acceptable — proceed with training")
else:
    print(f"❌ ONLY {usable}/{len(data)} records usable ({usable/len(data)*100:.1f}%)")
    print("   Fix issues before training")
print(f"{'='*60}")