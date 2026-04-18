import json
import pickle
import ast
from collections import defaultdict

with open("files/disjoint/my_data.json") as f:
    my_data = json.load(f)

with open("files/disjoint/ec2uniprot.p", "rb") as f:
    ec2uniprot = pickle.load(f)

# Build reverse map: uniprot_id → list of EC numbers
uniprot2ec = defaultdict(list)
for ec, uniprots in ec2uniprot.items():
    for uid in uniprots:
        uniprot2ec[uid].append(ec)


def resolve_ec(current_ec: str, protein_refs: list) -> tuple[str, str]:
    """
    Resolve the best EC number for a reaction.
    Returns (new_ec, reason).

    Rules:
    1. No mapping found        → keep current EC
    2. One exact match         → use it
    3. Multiple ECs            → prefer the one matching current EC if present,
                                 otherwise pick most specific (fewest wildcards)
    """
    ecs = set()
    for uid in protein_refs:
        ecs.update(uniprot2ec.get(uid, []))

    if not ecs:
        return current_ec, "no_mapping"

    if len(ecs) == 1:
        ec = list(ecs)[0]
        reason = "same" if ec == current_ec else "updated"
        return ec, reason

    # Multiple ECs — prefer current if it's in the set
    if current_ec in ecs:
        return current_ec, "kept_current_from_multiple"

    # Otherwise pick most specific: fewest wildcards and dashes
    def specificity(ec: str) -> int:
        return ec.count("-") + ec.count(".")  * -1  # fewer dashes = more specific

    best = sorted(ecs, key=lambda e: (e.count("-"), e.count(".") * -1))[0]
    return best, f"picked_from_multiple({len(ecs)})"


# ── Apply updates ─────────────────────────────────────────────────────────────

stats = defaultdict(int)
updated_data = []

for rxn in my_data:
    try:
        refs = ast.literal_eval(rxn.get("protein_refs", "[]"))
    except (ValueError, SyntaxError):
        refs = []

    new_ec, reason = resolve_ec(rxn["ec"], refs)

    stats[reason] += 1

    updated_rxn = {**rxn, "ec": new_ec, "ec_original": rxn["ec"]}
    updated_data.append(updated_rxn)


# ── Save ──────────────────────────────────────────────────────────────────────

with open("files/disjoint/my_data_updated.json", "w") as f:
    json.dump(updated_data, f, indent=2)

# ── Report ────────────────────────────────────────────────────────────────────

print(f"\n=== Update Summary ===")
print(f"Total reactions:                    {len(updated_data)}")
print(f"No mapping — kept current:          {stats['no_mapping']}")
print(f"Same EC confirmed:                  {stats['same']}")
print(f"Updated to new EC:                  {stats['updated']}")
print(f"Kept current from multiple options: {stats['kept_current_from_multiple']}")
for key, count in stats.items():
    if key.startswith("picked_from_multiple"):
        print(f"Picked best from multiple:          {count}")

print(f"\nSaved → my_data_updated.json")
print(f"Original EC preserved in 'ec_original' field")

# Show sample of changes
print(f"\n=== Sample changes ===")
changes = [r for r in updated_data if r["ec"] != r["ec_original"]][:5]
for r in changes:
    refs = ast.literal_eval(r.get("protein_refs", "[]"))
    print(f"  rxnid={r['rxnid']}  {r['ec_original']} → {r['ec']}  (refs={refs})")