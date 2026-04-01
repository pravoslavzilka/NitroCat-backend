"""
prepare_training_data.py
========================
Converts cyp_clipzyme_minimal.csv to CLIPZyme training format.

Input columns:  reaction, sequence, protein_id
Output columns: enzyme_id, smiles, split

Usage:
    python prepare_training_data.py
    python prepare_training_data.py --input my_file.csv --output files/my_data.csv
"""

import os
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",  "-i", default="files/cyp_clipzyme_minimal.csv")
    parser.add_argument("--output", "-o", default="files/my_data.csv")
    parser.add_argument("--train",  type=float, default=0.80)
    parser.add_argument("--val",    type=float, default=0.10)
    # test = remainder
    args = parser.parse_args()

    # ── Load ──────────────────────────────────────────────────────────────────
    print(f"Loading {args.input}...")
    df = pd.read_csv(args.input)
    print(f"  {len(df):,} rows  |  {df['protein_id'].nunique()} unique enzymes")

    # ── Convert columns ───────────────────────────────────────────────────────
    df = df.rename(columns={
        "reaction":   "smiles",
        "protein_id": "enzyme_id",
    })[["enzyme_id", "smiles"]]   # drop sequence — CLIPZyme uses CIF structures

    # ── Clean ─────────────────────────────────────────────────────────────────
    before = len(df)
    df = df.dropna().drop_duplicates(subset=["enzyme_id", "smiles"])
    if len(df) < before:
        print(f"  Removed {before - len(df)} null/duplicate rows")

    # ── Split 80/10/10 ────────────────────────────────────────────────────────
    df = df.sample(frac=1, random_state=42).reset_index(drop=True)
    n  = len(df)
    t  = int(n * args.train)
    v  = int(n * (args.train + args.val))

    df["split"] = "train"
    df.loc[t:v, "split"] = "val"
    df.loc[v:,  "split"] = "test"

    counts = df["split"].value_counts()
    print(f"  train: {counts.get('train', 0):,}")
    print(f"  val:   {counts.get('val',   0):,}")
    print(f"  test:  {counts.get('test',  0):,}")

    # ── Save ──────────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"\nSaved → {args.output}")
    print("\nFirst 3 rows:")
    print(df.head(3).to_string(index=False))

if __name__ == "__main__":
    main()