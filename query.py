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

def _atom_counts(smi: str) -> Counter:
    """Count heavy atoms only — RDKit handles H implicitly so we skip it."""
    total = Counter()
    for part in smi.split("."):
        mol = Chem.MolFromSmiles(part.strip())
        if mol:
            for a in mol.GetAtoms():
                if a.GetSymbol() != "H":   # ignore implicit/explicit H
                    total[a.GetSymbol()] += 1
    return total

def _balance(substrate: str, product: str):
    """
    Two-sided balancing using heavy atoms only (H is implicit in RDKit).

    Pass 1 — add co-reactants to substrate:
        atoms in product but not substrate (e.g. hydroxylation needs O)
    Pass 2 — add byproducts to product:
        atoms in substrate but not product (e.g. elimination loses Br)
    """
    comments = []

    s_counts = _atom_counts(substrate)
    p_counts = _atom_counts(product)

    # ── Pass 1: co-reactants ──────────────────────────────────────────────────
    # heavy-atom counts for each co-reactant (H not counted)
    co_reactants = [
        ("O",    {"O": 1}, "water"),
        ("O=O",  {"O": 2}, "molecular oxygen"),
        ("C",    {"C": 1}, "methyl donor (SAM placeholder)"),
        ("CC",   {"C": 2}, "ethyl donor"),
        ("CCC",  {"C": 3}, "propyl donor"),
    ]

    missing = {
        el: p_counts[el] - s_counts.get(el, 0)
        for el in p_counts if p_counts[el] > s_counts.get(el, 0)
    }

    added_reactants, remaining = [], missing.copy()

    for smi, co_counts, name in co_reactants:
        if not remaining:
            break
        while remaining and all(
            co_counts.get(el, 0) >= remaining.get(el, 0)
            for el in remaining if el in co_counts
        ) and any(el in remaining for el in co_counts):
            added_reactants.append(smi)
            for el in list(remaining):
                remaining[el] -= co_counts.get(el, 0)
                if remaining[el] <= 0:
                    del remaining[el]

    if remaining:
        # Don't add fallback water — just warn and proceed with unbalanced reaction
        comments.append(
            f"Could not balance missing atoms {list(remaining.keys())} "
            f"— may need cofactor (e.g. SAM, NADPH, coenzyme A). Proceeding as-is."
        )

    if added_reactants:
        names = list(dict.fromkeys(n for s, c, n in co_reactants if s in added_reactants))
        comments.append(f"Auto-added co-reactant(s): {', '.join(names)}.")

    full_substrate = substrate + ("." + ".".join(added_reactants) if added_reactants else "")
    s_counts_new   = _atom_counts(full_substrate)

    # ── Pass 2: byproducts ────────────────────────────────────────────────────
    # heavy-atom definitions only — H is omitted since it is never counted
    # CO2 MUST come before water so decarboxylation doesn't consume O prematurely
    byproducts = [
        ("O=C=O", {"C":  1, "O": 2}, "CO2"),
        ("Br",    {"Br": 1},         "HBr"),
        ("Cl",    {"Cl": 1},         "HCl"),
        ("F",     {"F":  1},         "HF"),
        ("N",     {"N":  1},         "ammonia"),
        ("S",     {"S":  1},         "H2S"),
        ("O",     {"O":  1},         "water"),
    ]

    excess = {
        el: s_counts_new[el] - p_counts.get(el, 0)
        for el in s_counts_new if s_counts_new[el] > p_counts.get(el, 0)
    }

    added_products, excess_remaining = [], excess.copy()

    for smi, bp_counts, name in byproducts:
        if not excess_remaining:
            break
        # ALL heavy atoms of byproduct must be present in excess
        while excess_remaining and all(
            el in excess_remaining and excess_remaining[el] >= n
            for el, n in bp_counts.items()
        ):
            added_products.append(smi)
            for el, n in bp_counts.items():
                excess_remaining[el] -= n
                if excess_remaining[el] <= 0:
                    del excess_remaining[el]

    if excess_remaining:
        # Don't add random water — just warn and proceed as-is
        comments.append(
            f"Could not balance excess heavy atoms {dict(excess_remaining)} "
            f"— may involve uncommon cofactors (SAM, NADPH, CoA). Proceeding as-is."
        )

    if added_products:
        bp_names = []
        for smi, _, name in byproducts:
            if smi in added_products and name not in bp_names:
                bp_names.append(name)
        comments.append(f"Auto-added byproduct(s): {', '.join(bp_names)}.")

    full_product = product + ("." + ".".join(added_products) if added_products else "")
    return f"{full_substrate}>>{full_product}", comments


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
    cases = [
        # (label, substrate, product)
        ("Elimination — HBr leaves",
            "CC1=C2CC[C@@]3(C)CC/C=C(/C)[C@H]3C[C@@H]2CC1(C)C",     "OC1(C)CC[C@@]2(C)[C@@H]3CC=C(C)C[C@H]3CC[C@@]2(C)C1"),

    ]

    for label, sub, prod in cases:
        print(f"\n=== {label} ===")
        print(f"  {sub} >> {prod}")
        print_result(query_enzymes(sub, prod))