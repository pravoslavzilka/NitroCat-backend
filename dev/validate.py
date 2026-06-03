import json
import ast
import re
import pickle
from pathlib import Path
from collections import Counter, defaultdict

DIR = Path(__file__).parent

data = json.load(open(DIR / "my_data_processed.json"))
ec2uniprot = pickle.load(open(DIR / "ec2uniprot.p", "rb"))
uniprot2sequence = pickle.load(open(DIR / "uniprot2sequence.p", "rb"))

missing_ec = []
missing_seq = []

for item in data:
    rxnid = item["rxnid"]
    ec = item.get("ec", "")
    refs = ast.literal_eval(item["protein_refs"])

    if ec not in ec2uniprot:
        missing_ec.append((rxnid, ec))

    for uid in refs:
        if uid not in uniprot2sequence:
            missing_seq.append((rxnid, uid))

print(f"Total reactions checked: {len(data)}")
print(f"EC numbers missing from ec2uniprot.p:    {len(missing_ec)}")
print(f"UniProt IDs missing from uniprot2sequence.p: {len(missing_seq)}")

if missing_ec:
    print("\nMissing EC numbers (rxnid, ec):")
    for entry in missing_ec[:20]:
        print(f"  {entry}")
    if len(missing_ec) > 20:
        print(f"  ... and {len(missing_ec) - 20} more")

if missing_seq:
    print("\nMissing sequences (rxnid, uniprot_id):")
    for entry in missing_seq[:20]:
        print(f"  {entry}")
    if len(missing_seq) > 20:
        print(f"  ... and {len(missing_seq) - 20} more")

# --- Atom mapping validation ---
# CLIPZyme requires: unique map numbers on each side, and every heavy atom
# in mapped_reactants/mapped_products must carry a map number.

def extract_map_numbers(smiles_list):
    nums = []
    for smi in smiles_list:
        nums.extend(int(n) for n in re.findall(r':(\d+)', smi))
    return nums

def has_unmapped_heavy_atoms(smiles_list):
    # Unbracketed atoms (organic subset) are always unmapped in mapped SMILES.
    # Bracketed atoms without :N are also unmapped.
    for smi in smiles_list:
        if re.search(r'(?<!\[)(?<![A-Za-z])(B(?!r)|C(?!l)|N|O|P|S(?!i)|F|Cl|Br|I)(?!\])', smi):
            return True
        if re.search(r'\[[A-Z][a-z]?\]|\[[A-Z][a-z]?H\d*\]', smi):
            return True
    return False

bad_mapping = []

for item in data:
    rxnid = item["rxnid"]
    r_maps = extract_map_numbers(item.get("mapped_reactants", []))
    p_maps = extract_map_numbers(item.get("mapped_products", []))

    issues = []

    # Duplicate map numbers within each side
    for n, c in Counter(r_maps).items():
        if c > 1:
            issues.append(f"duplicate map num {n} in reactants ({c}x)")
    for n, c in Counter(p_maps).items():
        if c > 1:
            issues.append(f"duplicate map num {n} in products ({c}x)")

    # Warn if any SMILES looks like it has unbracketed (unmapped) heavy atoms
    if has_unmapped_heavy_atoms(item.get("mapped_reactants", [])):
        issues.append("possible unmapped heavy atom in mapped_reactants")
    if has_unmapped_heavy_atoms(item.get("mapped_products", [])):
        issues.append("possible unmapped heavy atom in mapped_products")

    if issues:
        bad_mapping.append((rxnid, issues))

print(f"\nAtom-mapping issues (CLIPZyme): {len(bad_mapping)}")
if bad_mapping:
    for rxnid, issues in bad_mapping[:20]:
        print(f"  rxnid {rxnid}:")
        for iss in issues:
            print(f"    - {iss}")
    if len(bad_mapping) > 20:
        print(f"  ... and {len(bad_mapping) - 20} more")

# --- Disjoint split validation ---
# No reaction (identified by orig_rxn_text) may appear in more than one split.

rxn_to_splits = defaultdict(set)
for item in data:
    rxn_to_splits[item["orig_rxn_text"]].add(item.get("split", "unknown"))

overlap = {rxn: splits for rxn, splits in rxn_to_splits.items() if len(splits) > 1}

print(f"\nReactions appearing in multiple splits: {len(overlap)}")
if overlap:
    for rxn, splits in list(overlap.items())[:20]:
        print(f"  splits {splits}: {rxn[:80]}{'...' if len(rxn) > 80 else ''}")
    if len(overlap) > 20:
        print(f"  ... and {len(overlap) - 20} more")

if not missing_ec and not missing_seq and not bad_mapping and not overlap:
    print("\nAll checks passed.")
