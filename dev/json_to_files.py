import json
import ast
import pickle
import requests
from collections import defaultdict
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
DATA_PATH = SCRIPT_DIR / "my_data.json"
OUT_DIR = SCRIPT_DIR


def fetch_sequences_batch(uniprot_ids, batch_size=100):
    sequences = {}
    ids = list(uniprot_ids)
    url = "https://rest.uniprot.org/uniprotkb/search"
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i + batch_size]
        query = " OR ".join(f"accession:{uid}" for uid in batch)
        cursor = None
        while True:
            params = {
                "query": query,
                "fields": "accession,sequence",
                "format": "json",
                "size": batch_size,
            }
            if cursor:
                params["cursor"] = cursor
            resp = requests.get(url, params=params, timeout=60)
            resp.raise_for_status()
            data = resp.json()
            for entry in data.get("results", []):
                acc = entry["primaryAccession"]
                seq = entry.get("sequence", {}).get("value", "")
                if seq:
                    sequences[acc] = seq
            next_link = resp.links.get("next", {}).get("url")
            if not next_link:
                break
            cursor = next_link.split("cursor=")[-1].split("&")[0]
        print(f"  Fetched {len(sequences)} sequences so far ({i + len(batch)}/{len(ids)})")
    return sequences


def main():
    with open(DATA_PATH) as f:
        data = json.load(f)

    ec2uniprot = defaultdict(list)
    all_uniprots = set()

    for item in data:
        ec = item.get("ec", "")
        refs = ast.literal_eval(item["protein_refs"])
        for uid in refs:
            if uid not in ec2uniprot[ec]:
                ec2uniprot[ec].append(uid)
            all_uniprots.add(uid)

    print(f"Unique EC numbers: {len(ec2uniprot)}")
    print(f"Unique UniProt IDs: {len(all_uniprots)}")

    print("Fetching sequences from UniProt...")
    uniprot2sequence = fetch_sequences_batch(all_uniprots)
    print(f"Retrieved sequences for {len(uniprot2sequence)} / {len(all_uniprots)} proteins")

    with open(f"{OUT_DIR}/ec2uniprot.p", "wb") as f:
        pickle.dump(dict(ec2uniprot), f)
    print(f"Saved {OUT_DIR}/ec2uniprot.p")

    with open(f"{OUT_DIR}/uniprot2sequence.p", "wb") as f:
        pickle.dump(uniprot2sequence, f)
    print(f"Saved {OUT_DIR}/uniprot2sequence.p")


if __name__ == "__main__":
    main()
