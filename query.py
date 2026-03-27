import pickle
from collections import Counter

from clipzyme import CLIPZyme
from rxnmapper import RXNMapper
from rdkit import Chem

# ── Load once ─────────────────────────────────────────────────────────────────

screening     = pickle.load(open("files/clipzyme_screening_set.p", "rb"))
candidates    = screening["hiddens"]
candidate_ids = screening["uniprots"]

model  = CLIPZyme(checkpoint_path="files/clipzyme_model.ckpt").eval()
mapper = RXNMapper()

# ── Helpers ───────────────────────────────────────────────────────────────────

def _validate(smiles: str, label: str):
    if not smiles or not smiles.strip():
        return f"{label} is empty."
    for part in smiles.strip().split("."):
        if not part.strip():
            return f"{label} has an empty component (double dot or trailing dot)."
        if Chem.MolFromSmiles(part.strip()) is None:
            return f"{label} contains invalid SMILES: '{part}'."
    return None

def _canonical(smiles: str) -> str:
    parts = [Chem.MolToSmiles(Chem.MolFromSmiles(p.strip()))
             for p in smiles.split(".") if Chem.MolFromSmiles(p.strip())]
    return ".".join(sorted(parts))

def _balance(substrate: str, product: str):
    comments = []

    def counts(smi):
        total = Counter()
        for p in smi.split("."):
            mol = Chem.MolFromSmiles(p.strip())
            if mol:
                for a in mol.GetAtoms():
                    total[a.GetSymbol()] += 1
        return total

    missing = {el: n - counts(substrate).get(el, 0)
               for el, n in counts(product).items()
               if n > counts(substrate).get(el, 0)}

    if not missing:
        return f"{substrate}>>{product}", comments

    co_reactants = [("O", "water"), ("O=O", "molecular oxygen"), ("[H][H]", "hydrogen")]
    added, remaining = [], missing.copy()

    for smi, name in co_reactants:
        if not remaining:
            break
        mol = Chem.MolFromSmiles(smi)
        co  = Counter(a.GetSymbol() for a in mol.GetAtoms())
        while remaining and any(co.get(el, 0) >= remaining.get(el, 0) for el in remaining):
            added.append(smi)
            for el in list(remaining):
                remaining[el] -= co.get(el, 0)
                if remaining[el] <= 0:
                    del remaining[el]

    if remaining:
        added.append("O")
        comments.append(f"Could not fully balance {list(remaining.keys())} — water added as fallback.")

    if added:
        names = list(dict.fromkeys(n for s, n in co_reactants if s in added))
        comments.append(f"Auto-added co-reactant(s): {', '.join(names)}.")

    return f"{substrate + '.' + '.'.join(added)}>>{product}", comments

# ── Main ──────────────────────────────────────────────────────────────────────

def query_enzymes(substrate_smiles: str, product_smiles: str, top_k: int = 20) -> dict:
    """
    Find enzymes likely to catalyse a reaction.

    Returns
    -------
    {
        "status":   "success" | "error",
        "result":   [{"rank": 1, "uniprot": "Q8KLT9", "score": 0.5181}, ...]
                    OR  "error description string",
        "comments": ["note 1", "note 2", ...]
    }
    """
    comments = []

    err = _validate(substrate_smiles, "Substrate")
    if err:
        return {"status": "error", "result": err, "comments": comments}

    err = _validate(product_smiles, "Product")
    if err:
        return {"status": "error", "result": err, "comments": comments}

    if _canonical(substrate_smiles) == _canonical(product_smiles):
        return {"status": "error",
                "result": "Substrate and product are identical — no reaction to encode.",
                "comments": comments}

    balanced, balance_comments = _balance(substrate_smiles, product_smiles)
    comments.extend(balance_comments)

    try:
        out        = mapper.get_attention_guided_atom_maps([balanced])[0]
        mapped     = out["mapped_rxn"]
        confidence = float(out["confidence"])
        comments.append(f"Atom mapping confidence: {confidence:.2f}.")
    except Exception as e:
        return {"status": "error", "result": f"Atom mapping failed: {e}", "comments": comments}

    if ":1]" not in mapped:
        return {"status": "error",
                "result": "RXNMapper returned no atom numbers. Check your SMILES.",
                "comments": comments}

    if confidence < 0.1:
        return {"status": "error",
                "result": f"Atom mapping confidence too low ({confidence:.2f}).",
                "comments": comments}

    if confidence < 0.5:
        comments.append(f"Warning: moderate mapping confidence ({confidence:.2f}). Results may be less reliable.")

    try:
        rxn_vec = model.extract_reaction_features(reaction=mapped)
        scores  = (candidates @ rxn_vec.T).squeeze()
        top_idx = scores.argsort(descending=True)[:top_k]
    except Exception as e:
        return {"status": "error", "result": f"CLIPZyme encoding failed: {e}", "comments": comments}

    comments.append(f"Ranked {len(candidate_ids):,} enzyme candidates.")

    return {
        "status": "success",
        "result": [
            {"rank": i + 1, "uniprot": candidate_ids[j], "score": round(float(scores[j]), 4)}
            for i, j in enumerate(top_idx)
        ],
        "comments": comments,
    }


# ── Pretty print helper ───────────────────────────────────────────────────────

def print_result(output: dict):
    print(f"\nStatus:  {output['status'].upper()}")

    if output["status"] == "success":
        print("Result:")
        for e in output["result"]:
            print(f"  #{e['rank']:2d}  {e['uniprot']:15s}  score={e['score']:.4f}")
    else:
        print(f"Error:   {output['result']}")

    if output["comments"]:
        print("Comments:")
        for c in output["comments"]:
            print(f"  • {c}")


# ── Run ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    # 7. Invalid SMILES — gibberish
    print_result(query_enzymes("not_a_molecule", "CC=O"))

    # 8. Invalid SMILES — broken ring closure
    print_result(query_enzymes("C1CC", "CCC"))

    # 9. Empty substrate
    print_result(query_enzymes("", "CC=O"))

    # 10. Substrate equals product
    print_result(query_enzymes("CC(O)C", "CC(O)C"))

    # 11. Double dot in SMILES
    print_result(query_enzymes("CC(O)C..O", "CC(=O)C"))

    # 12. Valid SMILES but no chemical change possible to map
    print_result(query_enzymes("C", "C"))