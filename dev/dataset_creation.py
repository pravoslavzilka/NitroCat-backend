import json
import pickle
import re
import torch
import torch.nn.functional as F
from collections import Counter

from clipzyme import CLIPZyme
from clipzyme.utils.pyg import from_mapped_smiles
from clipzyme.utils.wln_processing import get_bond_changes
from rdkit import Chem
from rxnmapper import RXNMapper

CKPT          = "files/clipzyme_model.ckpt"
SCREENING_SET = "files/cyp_screening_set.p"


# ── Internal helpers ──────────────────────────────────────────────────────────

def _heavy_atom_counts(smiles: str) -> Counter:
    counts = Counter()
    for part in smiles.split("."):
        mol = Chem.MolFromSmiles(part.strip())
        if mol:
            for a in mol.GetAtoms():
                if a.GetSymbol() != "H":
                    counts[a.GetSymbol()] += 1
    return counts


def _balance_reaction(reaction: str) -> str:
    reactants, products = reaction.split(">>")
    r_counts = _heavy_atom_counts(reactants)
    p_counts = _heavy_atom_counts(products)

    co_reactants = [("O", {"O": 1}), ("O=O", {"O": 2}), ("C", {"C": 1})]
    missing = {el: p_counts[el] - r_counts.get(el, 0)
               for el in p_counts if p_counts[el] > r_counts.get(el, 0)}
    added_r, remaining = [], missing.copy()
    for smi, counts in co_reactants:
        while remaining and any(el in remaining for el in counts) and all(
            counts.get(el, 0) >= remaining.get(el, 0) for el in remaining if el in counts
        ):
            added_r.append(smi)
            for el in list(remaining):
                remaining[el] -= counts.get(el, 0)
                if remaining[el] <= 0:
                    del remaining[el]

    byproducts = [("O=C=O", {"C": 1, "O": 2}), ("Br", {"Br": 1}),
                  ("Cl", {"Cl": 1}), ("O", {"O": 1})]
    full_r = reactants + ("." + ".".join(added_r) if added_r else "")
    r2 = _heavy_atom_counts(full_r)
    excess = {el: r2[el] - p_counts.get(el, 0)
              for el in r2 if r2[el] > p_counts.get(el, 0)}
    added_p, exc_rem = [], excess.copy()
    for smi, counts in byproducts:
        while exc_rem and all(
            el in exc_rem and exc_rem[el] >= n for el, n in counts.items()
        ):
            added_p.append(smi)
            for el, n in counts.items():
                exc_rem[el] -= n
                if exc_rem[el] <= 0:
                    del exc_rem[el]

    full_p = products + ("." + ".".join(added_p) if added_p else "")
    return f"{full_r}>>{full_p}"


def _fix_mapped(mapped: str) -> str:
    """
    Two-pass cleanup of a mapped reaction SMILES:
    1. Assign fresh map numbers to any atoms RXNMapper left unmapped.
    2. Remove atoms whose map number appears on only one side (spectators
       CLIPZyme's WLN cannot handle because both graphs must be the same size).
    """
    existing = {int(x) for x in re.findall(r":(\d+)", mapped)}
    next_num = max(existing) + 1 if existing else 1
    sides = []
    for smi in mapped.split(">>"):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            sides.append(smi)
            continue
        for atom in mol.GetAtoms():
            if not atom.HasProp("molAtomMapNumber"):
                atom.SetAtomMapNum(next_num)
                next_num += 1
        sides.append(Chem.MolToSmiles(mol))
    mapped = ">>".join(sides)

    r_smi, p_smi = mapped.split(">>")
    only_r = {int(x) for x in re.findall(r":(\d+)", r_smi)} - \
             {int(x) for x in re.findall(r":(\d+)", p_smi)}
    only_p = {int(x) for x in re.findall(r":(\d+)", p_smi)} - \
             {int(x) for x in re.findall(r":(\d+)", r_smi)}

    if only_r or only_p:
        def strip_atoms(smi, to_remove):
            mol = Chem.RWMol(Chem.MolFromSmiles(smi))
            for idx in sorted(
                [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum() in to_remove],
                reverse=True,
            ):
                mol.RemoveAtom(idx)
            try:
                Chem.SanitizeMol(mol)
            except Exception:
                # Removing atoms can break aromatic rings leaving stale aromaticity flags;
                # clear them and re-perceive so kekulization succeeds.
                for atom in mol.GetAtoms():
                    atom.SetIsAromatic(False)
                for bond in mol.GetBonds():
                    if bond.GetBondType() == Chem.BondType.AROMATIC:
                        bond.SetBondType(Chem.BondType.SINGLE)
                try:
                    Chem.SanitizeMol(mol)
                except Exception:
                    return smi
            return Chem.MolToSmiles(mol)

        r_smi = strip_atoms(r_smi, only_r) if only_r else r_smi
        p_smi = strip_atoms(p_smi, only_p) if only_p else p_smi
        mapped = f"{r_smi}>>{p_smi}"

    return mapped


def _encode_reaction(model, mapper, reaction: str) -> torch.Tensor:
    balanced = _balance_reaction(reaction)
    out        = mapper.get_attention_guided_atom_maps([balanced])[0]
    mapped     = _fix_mapped(out["mapped_rxn"])
    confidence = float(out["confidence"])
    print(f"  Mapped (confidence {confidence:.4f}): {mapped}")

    try:
        with torch.no_grad():
            return model.extract_reaction_features(reaction=mapped).squeeze()
    except KeyError:
        pass

    reactants_smi, products_smi = mapped.split(">>")
    reactant_graph, atom_map2new_index = from_mapped_smiles(reactants_smi, encode_no_edge=True)
    product_graph,  _                  = from_mapped_smiles(products_smi,  encode_no_edge=True)

    bond_changes = []
    for u, v, btype in get_bond_changes(mapped):
        u_idx = atom_map2new_index.get(int(u))
        v_idx = atom_map2new_index.get(int(v))
        if u_idx is not None and v_idx is not None:
            bond_changes.append((u_idx, v_idx, btype))
    if not bond_changes:
        bond_changes = [(0, 1, 1.0)]

    reactant_graph.bond_changes = bond_changes
    reactant_graph.batch = torch.zeros(reactant_graph.x.shape[0], dtype=torch.long)
    product_graph.batch  = torch.zeros(product_graph.x.shape[0],  dtype=torch.long)

    with torch.no_grad():
        return model.model.encode_reaction(
            {"reactants": reactant_graph, "products": product_graph}
        ).squeeze()


# ── Public API ────────────────────────────────────────────────────────────────

def load_models(ckpt: str = CKPT, screening_set: str = SCREENING_SET) -> dict:
    """
    Load and return all heavy resources (CLIPZyme, RXNMapper, screening set).
    Call once and pass the result to compute_similarities() to avoid reloading.

    Returns
    -------
    {
        "model":        CLIPZyme,
        "mapper":       RXNMapper,
        "enzyme_embs":  Tensor (N, 1280),
        "enzyme_ids":   list[str],
    }
    """
    print("Loading CLIPZyme model...")
    model = CLIPZyme(checkpoint_path=ckpt).eval()
    mapper = RXNMapper()

    print("Loading screening set...")
    with open(screening_set, "rb") as f:
        screening = pickle.load(f)

    return {
        "model":       model,
        "mapper":      mapper,
        "enzyme_embs": screening["hiddens"].float(),
        "enzyme_ids":  list(screening["uniprots"]),
    }


def compute_similarities(comparisons: list[dict], resources: dict) -> list[dict]:
    """
    Compute cosine similarities between reactions and enzymes.

    Parameters
    ----------
    comparisons : list of dicts, each with:
        - "reaction" : str   reaction SMILES  (reactants>>products, unmapped)
        - "enzymes"  : list[str]  UniProt IDs to compare against

    resources : dict returned by load_models()

    Returns
    -------
    list of dicts, each with:
        - "reaction"      : str
        - "similarities"  : list of {"uniprot": str, "similarity": float}
    """
    model       = resources["model"]
    mapper      = resources["mapper"]
    enzyme_embs = resources["enzyme_embs"]
    enzyme_ids  = resources["enzyme_ids"]

    results = []
    for entry in comparisons:
        reaction = entry["reaction"]
        print(f"\nEncoding reaction: {reaction}")

        rxn_vec  = _encode_reaction(model, mapper, reaction)
        rxn_norm = F.normalize(rxn_vec.unsqueeze(0), dim=1)

        print(f"\n{'UniProt':<15} {'Similarity':>10}")
        print("-" * 27)

        similarities = []
        for uniprot in entry["enzymes"]:
            if uniprot not in enzyme_ids:
                print(f"{uniprot:<15} {'NOT FOUND':>10}")
                continue
            idx         = enzyme_ids.index(uniprot)
            enzyme_norm = F.normalize(enzyme_embs[idx].unsqueeze(0), dim=1)
            similarity  = (rxn_norm @ enzyme_norm.T).item()
            similarities.append({"uniprot": uniprot, "similarity": round(similarity, 6)})
            print(f"{uniprot:<15} {similarity:>10.6f}")

        results.append({"reaction": reaction, "similarities": similarities})

    return results


# ── Standalone use ────────────────────────────────────────────────────────────

COMPARISONS = [
    {
        "reaction": "CC(=O)Oc1ccccc1C(=O)O>>CC(=O)O.OC(=O)c1ccccc1O",
        "enzymes":  ["P08684", "Q9NR61", "P33261"],
    },
    {
        "reaction": "CC1=C2CC[C@@]3(C)CC/C=C(/C)[C@H]3C[C@@H]2CC1(C)C>>OC1(C)CC[C@@]2(C)[C@@H]3CC=C(C)C[C@H]3CC[C@@]2(C)C1",
        "enzymes":  ["P08684", "P33261"],
    },
]

if __name__ == "__main__":
    resources = load_models()
    results   = compute_similarities(COMPARISONS, resources)
    print("\n" + "=" * 50)
    print(json.dumps(results, indent=2))
