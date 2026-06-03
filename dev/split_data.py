import json
import random
from pathlib import Path

SEED = 42
TRAIN_RATIO = 0.80
VAL_RATIO   = 0.10

DIR         = Path(__file__).parent
SOURCE_PATH = DIR / "my_data.json"
REF_PATH    = DIR.parent / "files" / "my_data.json"
OUT_PATH    = DIR / "my_data.json"

data = json.load(open(SOURCE_PATH))
ref  = json.load(open(REF_PATH))

inherited = {r["rxnid"]: r["split"] for r in ref}

new_records = [r for r in data if r["rxnid"] not in inherited]

random.seed(SEED)
random.shuffle(new_records)

n        = len(new_records)
n_train  = int(n * TRAIN_RATIO)
n_val    = int(n * VAL_RATIO)

for i, r in enumerate(new_records):
    if i < n_train:
        r["split"] = "train"
    elif i < n_train + n_val:
        r["split"] = "val"
    else:
        r["split"] = "test"

inherited_count = 0
new_count       = 0
for r in data:
    if r["rxnid"] in inherited:
        r["split"] = inherited[r["rxnid"]]
        inherited_count += 1
    else:
        new_count += 1

counts = {}
for r in data:
    counts[r["split"]] = counts.get(r["split"], 0) + 1

print(f"Total records:      {len(data)}")
print(f"Inherited splits:   {inherited_count}")
print(f"Newly assigned:     {new_count}")
print(f"Split counts:       {counts}")

json.dump(data, open(OUT_PATH, "w"), indent=2)
print(f"Saved → {OUT_PATH}")
