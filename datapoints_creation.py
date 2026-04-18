import ast
import csv
import json

import torch
import torch.nn.functional as F

from dataset_creation import load_models, _encode_reaction

DATA_FILE = "files/disjoint/my_data_updated.json"
OUT_CSV   = "datapoints_similarities.csv"


def build_comparisons(data: list[dict]) -> list[dict]:
    comparisons = []
    for entry in data:
        reaction = entry["orig_rxn_text"]
        try:
            enzymes = ast.literal_eval(entry["protein_refs"])
        except (ValueError, SyntaxError):
            print(f"  Skipping rxnid {entry['rxnid']}: bad protein_refs {entry['protein_refs']!r}")
            continue
        if not enzymes:
            continue
        comparisons.append({"reaction": reaction, "enzymes": enzymes, "rxnid": entry["rxnid"]})
    return comparisons


def build_enzyme_order(comparisons: list[dict]) -> list[str]:
    """
    Build the column order so that the diagonal holds the dataset pairs:
    enzyme_i (first enzyme of reaction_i) is at column i.
    Remaining enzymes are appended in encounter order.
    """
    diagonal_enzymes = [c["enzymes"][0] for c in comparisons]
    seen = set(diagonal_enzymes)
    extra = []
    for c in comparisons:
        for e in c["enzymes"][1:]:
            if e not in seen:
                seen.add(e)
                extra.append(e)
    return diagonal_enzymes + extra


if __name__ == "__main__":
    print(f"Loading data from {DATA_FILE}...")
    with open(DATA_FILE) as f:
        data = json.load(f)
    print(f"  {len(data)} entries loaded")

    comparisons = build_comparisons(data)
    print(f"  {len(comparisons)} valid comparisons built")

    resources = load_models()
    model     = resources["model"]
    mapper    = resources["mapper"]
    emb_map   = resources["enzyme_embs"]
    emb_ids   = resources["enzyme_ids"]

    enzyme_cols = build_enzyme_order(comparisons)
    print(f"  {len(enzyme_cols)} unique enzymes (columns)")

    # Pre-normalise all needed enzyme embeddings
    enzyme_norm = {}
    for uniprot in enzyme_cols:
        if uniprot in emb_ids:
            idx = emb_ids.index(uniprot)
            enzyme_norm[uniprot] = F.normalize(emb_map[idx].unsqueeze(0), dim=1)

    reactions = [c["reaction"] for c in comparisons]

    # Encode all reactions (skipping failures)
    rxn_vecs: dict[str, torch.Tensor | None] = {}
    for rxn in reactions:
        if rxn in rxn_vecs:
            continue
        print(f"\nEncoding: {rxn}")
        try:
            vec = _encode_reaction(model, mapper, rxn)
            rxn_vecs[rxn] = F.normalize(vec.unsqueeze(0), dim=1)
        except Exception as exc:
            print(f"  ERROR encoding reaction: {exc}")
            rxn_vecs[rxn] = None

    # Write CSV: first row = header (enzymes), first column = reaction SMILES
    with open(OUT_CSV, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["reaction"] + enzyme_cols)
        for rxn in reactions:
            norm = rxn_vecs.get(rxn)
            if norm is None:
                row = [rxn] + [""] * len(enzyme_cols)
            else:
                row = [rxn]
                for uniprot in enzyme_cols:
                    if uniprot not in enzyme_norm:
                        row.append("")
                    else:
                        sim = (norm @ enzyme_norm[uniprot].T).item()
                        row.append(round(sim, 6))
            writer.writerow(row)

    print(f"\nSaved {len(reactions)} x {len(enzyme_cols)} matrix to {OUT_CSV}")
    print("Diagonal entries are the dataset reaction–enzyme pairs.")
