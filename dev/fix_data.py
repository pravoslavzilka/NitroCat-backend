# remap_safe.py — saves after every record, nothing can be lost
import json, re, ast, os
from rxnmapper import RXNMapper
from clipzyme.utils.wln_processing import get_bond_changes
from collections import Counter

INPUT  = 'files/disjoint/my_data_ready.json'
OUTPUT = 'files/disjoint/my_data_final.json'

# Copy input to output first so we always have a valid file
import shutil
if not os.path.exists(OUTPUT):
    shutil.copy(INPUT, OUTPUT)
    print(f"Initialized output from input")

data = json.load(open(OUTPUT))
print(f"Loaded {len(data)} records from {OUTPUT}")

def has_atom_map(smi):
    return bool(re.search(r':\d+', smi))

def strip_atom_map(smi):
    return re.sub(r':\d+', '', smi)

unmapped_idx = [
    i for i, r in enumerate(data)
    if any(not has_atom_map(s)
           for s in r.get('mapped_reactants', []) + r.get('mapped_products', []))
    and not r.get('_atom_map_failed')
]
print(f"Records to fix: {len(unmapped_idx)}")

if len(unmapped_idx) == 0:
    print("✅ Nothing to fix!")
else:
    rxn_mapper = RXNMapper()
    fixed  = 0
    failed = 0

    for j, i in enumerate(unmapped_idx):
        r     = data[i]
        raw_r = r.get('reactants') or [strip_atom_map(s) for s in r['mapped_reactants']]
        raw_p = r.get('products')  or [strip_atom_map(s) for s in r['mapped_products']]
        rxn   = "{}>>{}".format('.'.join(raw_r), '.'.join(raw_p))

        try:
            result = rxn_mapper.get_attention_guided_atom_maps([rxn])[0]
            mr, mp = result['mapped_rxn'].split('>>')
            r['mapped_reactants'] = mr.split('.')
            r['mapped_products']  = mp.split('.')
            r.pop('_atom_map_failed', None)

            # Recompute bond changes immediately
            try:
                rxn_str = "{}>>{}".format(
                    '.'.join(sorted(r['mapped_reactants'])),
                    '.'.join(sorted(r['mapped_products']))
                )
                bc = list(get_bond_changes(rxn_str))
                r['bond_changes'] = bc if bc else [[0, 1, 1.0]]
            except:
                r['bond_changes'] = [[0, 1, 1.0]]

            fixed += 1

        except Exception as e:
            r['_atom_map_failed'] = True
            failed += 1

        # ← Save after EVERY record — nothing lost if script crashes
        data[i] = r
        json.dump(data, open(OUTPUT, 'w'), indent=2)

        if (j + 1) % 100 == 0:
            fully = sum(1 for r in data
                       if all(has_atom_map(s)
                              for s in r['mapped_reactants'] + r['mapped_products']))
            print(f"  [{j+1}/{len(unmapped_idx)}] fixed={fixed} failed={failed} "
                  f"fully_mapped={fully}/{len(data)}")

print(f"\n✅ Done — fixed={fixed} failed={failed}")

# Final check
data = json.load(open(OUTPUT))
fully = sum(1 for r in data
            if all(has_atom_map(s)
                   for s in r['mapped_reactants'] + r['mapped_products']))
print(f"Final: {fully}/{len(data)} fully mapped ({fully/len(data)*100:.1f}%)")
print(f"Split: {dict(Counter(r['split'] for r in data))}")