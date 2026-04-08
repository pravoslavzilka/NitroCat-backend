import json
import pickle
import ast
import requests
import time
from collections import defaultdict

# Get all unique UniProt IDs from your training data
data = json.load(open('files/my_data.json'))

uniprot_ids = set()
for record in data:
    refs = record.get('protein_refs', '[]')
    uids = ast.literal_eval(refs) if isinstance(refs, str) else refs
    uniprot_ids.update(uids)

print(f"Unique UniProt IDs: {len(uniprot_ids)}")

# Fetch EC numbers from UniProt REST API
def get_ec(uniprot_id):
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        resp = requests.get(url, timeout=10)
        if resp.status_code != 200:
            return []
        data = resp.json()
        ecs = []
        for db in data.get('dbReferences', []):
            if db.get('type') == 'EC':
                ecs.append(db.get('id', ''))
        # Also check proteinDescription
        for rec in data.get('proteinDescription', {}).get('recommendedName', {}).get('ecNumbers', []):
            ecs.append(rec.get('value', ''))
        return list(set(filter(None, ecs)))
    except:
        return []

# Fetch EC for all your enzymes
uid2ec = {}
for i, uid in enumerate(uniprot_ids):
    ecs = get_ec(uid)
    uid2ec[uid] = ecs
    if i % 50 == 0:
        print(f"  {i}/{len(uniprot_ids)} — {uid}: {ecs}")
    time.sleep(0.1)  # rate limit

print(f"\nFetched EC numbers for {sum(1 for v in uid2ec.values() if v)} enzymes")

# Build ec2uniprot mapping
cyp_ec2uniprot = defaultdict(set)
for uid, ecs in uid2ec.items():
    for ec in ecs:
        cyp_ec2uniprot[ec].add(uid)

print(f"Unique EC numbers found: {len(cyp_ec2uniprot)}")

# Merge with existing ec2uniprot.p
ec2uniprot = pickle.load(open('files/ec2uniprot.p', 'rb'))

for ec, uids in cyp_ec2uniprot.items():
    if ec in ec2uniprot:
        ec2uniprot[ec] = list(set(ec2uniprot[ec]) | uids)
    else:
        ec2uniprot[ec] = list(uids)

pickle.dump(ec2uniprot, open('files/ec2uniprot.p', 'wb'))
print(f"Saved → files/ec2uniprot.p ({len(ec2uniprot)} total EC numbers)")