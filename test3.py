
import json
import pickle
import ast

# ── Load data ─────────────────────────────────────────────────────────────────

with open("files/disjoint/my_data.json") as f:
    my_data = json.load(f)

with open("files/disjoint/uniprot2sequence.p", "rb") as f:
    uniprot2seq = pickle.load(f)

# Reverse map: sequence → uniprot_id
seq2uniprot = {seq.strip(): uid for uid, seq in uniprot2seq.items()}

print(f"Reactions:  {len(my_data)}")
print(f"Proteins:   {len(uniprot2seq)}")


# ── Build index: uniprot_id → list of reactions ───────────────────────────────

def build_reaction_index(my_data: list) -> dict:
    """Index reactions by UniProt ID for O(1) lookups."""
    index = {}
    for rxn in my_data:
        # protein_refs is stored as a string like "['P0DO79', 'Q9Y6K9']"
        try:
            refs = ast.literal_eval(rxn.get("protein_refs", "[]"))
        except (ValueError, SyntaxError):
            refs = []

        for uniprot_id in refs:
            if uniprot_id not in index:
                index[uniprot_id] = []
            index[uniprot_id].append(rxn)

    return index


reaction_index = build_reaction_index(my_data)
print(f"UniProt IDs with reactions: {len(reaction_index)}")


# ── Lookup functions ──────────────────────────────────────────────────────────

def find_reactions_by_uniprot(uniprot_id: str) -> dict:
    """Find all reactions for a UniProt ID."""
    reactions = reaction_index.get(uniprot_id, [])
    return {
        "status":      "found" if reactions else "not_found",
        "uniprot_id":  uniprot_id,
        "n_reactions": len(reactions),
        "reactions":   reactions,
    }


def find_reactions_by_sequence(sequence: str) -> dict:
    """Find all reactions for an enzyme given its amino acid sequence."""
    uniprot_id = seq2uniprot.get(sequence.strip())
    if not uniprot_id:
        return {
            "status":      "not_found",
            "uniprot_id":  None,
            "n_reactions": 0,
            "reactions":   [],
            "message":     "Sequence not in uniprot2sequence.p",
        }
    return find_reactions_by_uniprot(uniprot_id)


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    # Your enzyme sequence
    my_sequence = "MDPIMYYSLLTLAFIITIKLVLQIQSRRLKNLPPGPPTLPIIGNLHHLKPPIHRTFSTLSQKYGDIISLWFGSRLVVVVSSPSLVQECFTKNDVVLANRPRFLTGKYIFYNYSTLGSASYGDHWRNLRRITTLDVLSNNRLNSFIDIRRDEAMRLVQKLGHDTKASDSEGFVKVELRSRLTEMTFNEMMRMISGKRYYGEDIDVSDVEEASQFREIISEMLSLLGANNKGDFLPLLRLFDFEDLEKRLKRIAKRADAFLQGLIEEHRVGKHSADTMIEHLLKMQESQPEYYSDLMIKGLIQAMLLAGTDTSAVTIEWVMAELLNHPEALKKVKDELDTQIGKDRLVNEQDTPKLPYLQNVISEALRLHPPAPLLLPHSSSEAFTLGGYNIPRDTILLTNAWLIHRDPKLWSDAASFKPERFEKEGEVNKLIAFGLGRRACPGLGLAQRTVGYTVGLLIQCFEWKRESEEKLDMMEDKGVTMPKRIPLEALCKARPIVNDVMK"

    result = find_reactions_by_sequence(my_sequence)

    print(f"UniProt:   {result['uniprot_id']}")
    print(f"Reactions: {result['n_reactions']}")
    for rxn in result["reactions"]:
        print(f"\n  rxnid:     {rxn['rxnid']}")
        print(f"  ec:        {rxn['ec']}")
        print(f"  reactants: {rxn['reactants']}")
        print(f"  products:  {rxn['products']}")