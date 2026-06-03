import json
import random
from collections import defaultdict
from pathlib import Path

SEED       = 42
TRAIN_FRAC = 0.80
DEV_FRAC   = 0.10
# test gets the remainder

DIR       = Path(__file__).parent
DATA_PATH = DIR / "my_data.json"

data = json.load(open(DATA_PATH))
print(f"Loaded {len(data)} records")

# Canonical reaction key: sorted reactants >> sorted products
def rxn_key(r):
    return ">>".join([
        ".".join(sorted(r.get("reactants", []))),
        ".".join(sorted(r.get("products", [])))
    ])

# Group records by reaction
rxn_to_records = defaultdict(list)
for rec in data:
    rxn_to_records[rxn_key(rec)].append(rec)

unique_rxns = list(rxn_to_records.keys())
print(f"Unique reactions: {len(unique_rxns)}")

# Shuffle reactions and split
random.seed(SEED)
random.shuffle(unique_rxns)

n      = len(unique_rxns)
n_train = int(n * TRAIN_FRAC)
n_dev   = int(n * DEV_FRAC)

train_rxns = set(unique_rxns[:n_train])
dev_rxns   = set(unique_rxns[n_train:n_train + n_dev])
test_rxns  = set(unique_rxns[n_train + n_dev:])

print(f"Split reactions — train: {len(train_rxns)}, dev: {len(dev_rxns)}, test: {len(test_rxns)}")

# Assign split to all records
for rec in data:
    key = rxn_key(rec)
    if key in train_rxns:
        rec["split"] = "train"
    elif key in dev_rxns:
        rec["split"] = "dev"
    else:
        rec["split"] = "test"

# Summary
split_counts = {}
for rec in data:
    split_counts[rec["split"]] = split_counts.get(rec["split"], 0) + 1

print(f"Record counts  — {split_counts}")

# Verify disjointness
train_keys = {rxn_key(r) for r in data if r["split"] == "train"}
dev_keys   = {rxn_key(r) for r in data if r["split"] == "dev"}
test_keys  = {rxn_key(r) for r in data if r["split"] == "test"}

assert len(train_keys & dev_keys)  == 0, "LEAK: train/dev share reactions!"
assert len(train_keys & test_keys) == 0, "LEAK: train/test share reactions!"
assert len(dev_keys   & test_keys) == 0, "LEAK: dev/test share reactions!"
print("Disjointness verified — no reaction leaks between splits")

json.dump(data, open(DATA_PATH, "w"), indent=2)
print(f"Saved → {DATA_PATH}")
