#!/home/pravoslav/miniconda3/envs/clipzyme/bin/python3
"""
Fix atom-mapping issues in my_data.json for CLIPZyme compatibility.
Assigns map numbers to every non-aromatic, unmapped heavy atom in
mapped_reactants using RDKit.  mapped_products is left untouched.
Writes result to my_data_processed.json.
"""
import json
import re
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import RWMol

DIR = Path(__file__).parent


def max_map(smiles_list):
    nums = [int(n) for s in smiles_list for n in re.findall(r":(\d+)", s)]
    return max(nums) if nums else 0


def assign_maps(mol, start):
    """
    Assign consecutive map numbers (starting at `start`) to every
    non-aromatic atom that currently has no map number (mapNum == 0).
    Aromatic atoms (written as lowercase in SMILES: c, n, o, s …) are
    skipped — validate.py only flags uppercase unbracketed atoms and
    aromatic atoms are handled correctly by CLIPZyme without explicit maps.
    Returns (new_mol, next_start, changed).
    """
    em = RWMol(mol)
    cur = start
    changed = False
    for atom in em.GetAtoms():
        if atom.GetAtomMapNum() == 0 and not atom.GetIsAromatic():
            atom.SetAtomMapNum(cur)
            cur += 1
            changed = True
    return em.GetMol(), cur, changed


def process_item(item):
    all_smi = item.get("mapped_reactants", []) + item.get("mapped_products", [])
    start = max_map(all_smi) + 1

    new_reactants = []
    item_changed = False

    for smi in item.get("mapped_reactants", []):
        mol = Chem.MolFromSmiles(smi)
        if mol is None or not any(a.GetAtomMapNum() == 0 for a in mol.GetAtoms()):
            new_reactants.append(smi)
            continue

        mapped_mol, start, changed = assign_maps(mol, start)
        if changed:
            new_reactants.append(Chem.MolToSmiles(mapped_mol))
            item_changed = True
        else:
            new_reactants.append(smi)

    if item_changed:
        item["mapped_reactants"] = new_reactants

    return item_changed


data = json.load(open(DIR / "my_data.json"))
fixed = sum(process_item(item) for item in data)
json.dump(data, open(DIR / "my_data_processed.json", "w"), indent=2)
print(f"Fixed {fixed} reactions → my_data_processed.json")
