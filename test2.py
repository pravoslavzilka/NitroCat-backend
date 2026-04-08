import json
import re
import requests

with open("brenda_2026_1.json") as f:
    brenda = json.load(f)

def get_enzyme_data(ec: str, uniprot_id: str) -> dict:
    ec_entry = brenda["data"][ec]

    # Step 1 — find the BRENDA protein ID for this UniProt accession
    protein_id = None
    for pid, protein in ec_entry["protein"].items():
        if uniprot_id in protein.get("accessions", []):
            protein_id = pid
            break

    if not protein_id:
        raise ValueError(f"UniProt ID {uniprot_id} not found in EC {ec}")

    print(f"Found protein_id: {protein_id}, organism: {ec_entry['protein'][protein_id]['organism']}")

    # Step 2 — collect all data fields that mention this protein_id
    data_fields = [
        "km_value", "turnover_number", "kcat_km_value",
        "ki_value", "ic50_value",
        "ph_optimum", "ph_range", "ph_stability",
        "temperature_optimum", "temperature_range", "temperature_stability",
        "specific_activity", "molecular_weight", "pi_value",
        "inhibitor", "cofactor", "activating_compound",
        "substrates_products", "natural_substrates_products",
        "subunits", "localization", "source_tissue",
    ]

    result = {"protein_id": protein_id, "ec": ec, "uniprot": uniprot_id}

    for field in data_fields:
        entries = ec_entry.get(field, [])
        # Keep only entries linked to our protein
        matched = [e for e in entries if protein_id in e.get("proteins", [])]
        if matched:
            result[field] = matched

    # Step 3 — resolve literature references for context
    all_ref_ids = set()
    for entries in result.values():
        if isinstance(entries, list):
            for e in entries:
                all_ref_ids.update(e.get("references", []))

    result["literature"] = {
        rid: ec_entry["reference"][rid]
        for rid in all_ref_ids
        if rid in ec_entry["reference"]
    }

    return result






def summarize_enzyme(ec_entry: dict, protein_id: str):
    """Show everything available for a protein before filtering."""

    numeric_fields = [
        "km_value", "turnover_number", "kcat_km_value", "ki_value",
        "ph_optimum", "ph_range", "temperature_optimum", "temperature_range",
        "specific_activity", "molecular_weight",
    ]
    text_fields = [
        "inhibitor", "cofactor", "activating_compound",
        "substrates_products", "subunits",
    ]

    print(f"\n{'='*60}")
    print(f"SUMMARY FOR PROTEIN ID: {protein_id}")
    print(f"{'='*60}")

    for field in numeric_fields:
        entries = [e for e in ec_entry.get(field, []) if protein_id in e.get("proteins", [])]
        if not entries:
            continue

        print(f"\n--- {field.upper()} ({len(entries)} entries) ---")

        # extract all unique substrates
        substrates = set()
        comments = set()
        for e in entries:
            match = re.match(r"[\d.eE+\-]+\s*\{(.+)\}", e["value"].strip())
            if match:
                substrates.add(match.group(1))
            if e.get("comment"):
                comments.add(e["comment"])

        if substrates:
            print(f"  Substrates: {sorted(substrates)}")
        if comments:
            print(f"  Comments:   {sorted(comments)}")

        # show value range
        values = []
        for e in entries:
            try:
                values.append(float(e["value"].split()[0]))
            except:
                pass
        if values:
            print(f"  Value range: {min(values)} — {max(values)}")

    for field in text_fields:
        entries = [e for e in ec_entry.get(field, []) if protein_id in e.get("proteins", [])]
        if not entries:
            continue

        print(f"\n--- {field.upper()} ({len(entries)} entries) ---")
        for e in entries:
            print(f"  {e['value']}  [{e.get('comment', '')}]")


def find_protein_id(ec: str, uniprot_id: str) -> tuple[str, dict]:
    """Find BRENDA protein_id from a UniProt accession."""
    ec_entry = brenda["data"].get(ec)
    if not ec_entry:
        raise ValueError(f"EC {ec} not found in BRENDA")

    for pid, protein in ec_entry.get("protein", {}).items():
        if uniprot_id in protein.get("accessions", []):
            return pid, ec_entry

    raise ValueError(f"UniProt ID {uniprot_id} not found under EC {ec}")




protein_id, ec_entry = find_protein_id("1.1.1.1", "P00326")

# Print ONLY the protein_id
print(f"Protein ID: {protein_id}")

# And just the protein info, not the whole ec_entry
print(f"Organism: {ec_entry['protein'][protein_id]['organism']}")
print(f"Accessions: {ec_entry['protein'][protein_id]['accessions']}")


summarize_enzyme(ec_entry, protein_id)

