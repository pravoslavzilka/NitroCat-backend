#!/usr/bin/env python3
"""
search.py — RHEA reaction similarity search.

Programmatic use:
    from search import RheaSearcher

    searcher = RheaSearcher()          # loads caches once (~1 s)
    results  = searcher.find("CC(=O)O.[H][H]>>CC(O)")
    results  = searcher.find("...", top_k=10, fetch_pubs=False)

Each result dict:
    rhea_id, score, score_morgan, score_drfp,
    ec, confidence, url, publications

CLI use:
    python search.py "reactants>>products"
    python search.py --top 10 --no-pubs "reactants>>products"
    python search.py --reestimate   "reactants>>products"
"""
import argparse
import json
import random
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from drfp import DrfpEncoder
from rdkit import Chem

BASE_DIR  = Path(__file__).parent
CACHE_DIR = BASE_DIR / "rhea_cache"

sys.path.insert(0, str(BASE_DIR))
from simi_core import (
    _encode_reaction, load_fingerprint_cache, load_ec_cache,
    fetch_publications_parallel, load_uniprot_cache, build_uniprot_cache,
    DRFP_FPS_CACHE, DRFP_IDS_CACHE, _cache_ready,
)

ALPHA_CACHE = CACHE_DIR / "alpha.json"
CYP_DATA        = BASE_DIR  / "files" / "my_data_updated.json"

ALPHA_STEPS     = 21
TUNE_PER_SOURCE = 150
SEED            = 42


# ── low-level helpers (module-level so they stay picklable) ──────────────────

def strip_atom_map(rxn_smi: str):
    parts = rxn_smi.split(">>")
    if len(parts) != 2:
        return None
    out = []
    for side in parts:
        mols = []
        for smi in side.split("."):
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                return None
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)
            mols.append(Chem.MolToSmiles(mol))
        out.append(".".join(mols))
    return ">>".join(out)


def _encode_drfp(smi: str):
    try:
        return np.array(DrfpEncoder.encode([smi])[0], dtype=np.uint8)
    except Exception:
        return None


def _tanimoto(matrix: np.ndarray, q: np.ndarray) -> np.ndarray:
    inter = np.sum(matrix & q, axis=1)
    union = np.sum(matrix | q, axis=1)
    return np.divide(inter, union, out=np.zeros(len(matrix), dtype=np.float32),
                     where=union != 0)


def _mrr(sims, q_idx, pos):
    row = sims.copy()
    if q_idx >= 0:
        row[q_idx] = -1.0
    for rank, idx in enumerate(np.argsort(row)[::-1], 1):
        if idx in pos:
            return 1.0 / rank
    return 0.0


# ── RheaSearcher ─────────────────────────────────────────────────────────────

class RheaSearcher:
    """
    Load RHEA fingerprint matrices and EC metadata once, then call
    .find(smiles) repeatedly without reloading.

    Parameters
    ----------
    reestimate : bool
        If True, re-run the α grid search even if a cached value exists.
    verbose : bool
        Print loading / tuning progress.
    """

    def __init__(self, reestimate: bool = False, verbose: bool = True):
        self._verbose = verbose
        self._load_matrices()
        self._alpha = self._get_alpha(reestimate)

    # ── public API ────────────────────────────────────────────────────────────

    def find(
        self,
        smiles: str,
        top_k: int = 5,
        fetch_pubs: bool = True,
    ) -> list[dict]:
        """
        Find the top_k most similar RHEA reactions for a query SMILES.

        Parameters
        ----------
        smiles     : reaction SMILES  (atom mapping stripped automatically)
        top_k      : number of hits to return
        fetch_pubs : fetch up to 2 publications per hit from Europe PMC

        Returns
        -------
        List of dicts, sorted by descending score:
            rhea_id       str    e.g. "33569"
            score         float  combined similarity  [0, 1]
            score_morgan  float  Morgan-only Tanimoto
            score_drfp    float  DRFP-only Tanimoto
            ec            list   EC numbers for this reaction
            confidence    float  EC consensus among top-10 hits  [0, 1]
            url           str    https://www.rhea-db.org/rhea/{id}
            publications  list   [{pmid, title, authors, journal, year}]
        """
        clean = strip_atom_map(smiles) or smiles
        qm    = _encode_reaction(clean)
        qd    = _encode_drfp(clean)

        if qm is None:
            raise ValueError(f"Could not encode Morgan FP for: {smiles}")
        if qd is None:
            raise ValueError(f"Could not encode DRFP for: {smiles}")

        sm       = _tanimoto(self._morgan_fps, qm)
        sd       = _tanimoto(self._drfp_fps,   qd)
        combined = self._alpha * sm + (1 - self._alpha) * sd

        top_idxs      = np.argsort(combined)[::-1][:max(top_k, 10)]
        confidence    = self._ec_consensus(top_idxs)
        top_idxs      = top_idxs[:top_k]
        top_rhea_ids  = [self._common_ids[i] for i in top_idxs]
        top_master_ids = [self._dir_map.get(rid, rid) for rid in top_rhea_ids]

        pub_map = (fetch_publications_parallel(top_master_ids, self._uniprot_map, page_size=2)
                   if fetch_pubs else {})

        results = []
        for rank, i in enumerate(top_idxs, 1):
            rid       = self._common_ids[i]
            master_id = self._dir_map.get(rid, rid)
            results.append({
                "rank":         rank,
                "rhea_id":      rid,
                "reaction":     self._smi_map.get(rid, ""),
                "score":        float(combined[i]),
                "score_morgan": float(sm[i]),
                "score_drfp":   float(sd[i]),
                "ec":           self._ec_map_full.get(master_id, []),
                "confidence":   confidence,
                "url":          f"https://www.rhea-db.org/rhea/{rid}",
                "publications": pub_map.get(master_id, []),
            })
        return results

    @property
    def alpha(self) -> float:
        """Optimal Morgan weight (α). DRFP weight = 1 − α."""
        return self._alpha

    # ── internal loading ──────────────────────────────────────────────────────

    def _log(self, msg):
        if self._verbose:
            print(msg)

    def _load_matrices(self):
        if not _cache_ready.is_set():
            self._log("Waiting for cache setup to complete …")
            _cache_ready.wait()
        self._log("Loading RHEA caches …")
        morgan_full, morgan_ids, morgan_smiles = load_fingerprint_cache()
        drfp_full = np.load(DRFP_FPS_CACHE)
        with open(DRFP_IDS_CACHE) as f:
            drfp_ids = json.load(f)

        drfp_set      = set(drfp_ids)
        common        = [rid for rid in morgan_ids if rid in drfp_set]
        m_idx         = {rid: i for i, rid in enumerate(morgan_ids)}
        d_idx         = {rid: i for i, rid in enumerate(drfp_ids)}

        self._common_ids  = common
        self._morgan_fps  = morgan_full[[m_idx[r] for r in common]]
        self._drfp_fps    = drfp_full[[d_idx[r]   for r in common]]
        self._smi_map     = dict(zip(morgan_ids, morgan_smiles))

        self._ec_map, self._dir_map = load_ec_cache()
        try:
            self._uniprot_map, ec_supplement = load_uniprot_cache()
        except FileNotFoundError:
            self._uniprot_map, ec_supplement = build_uniprot_cache(self._ec_map)

        # Merge EC sources: RHEA ec_map takes priority, UniProt supplement fills gaps
        self._ec_map_full = {**ec_supplement, **self._ec_map}

        master_of = {rid: self._dir_map.get(rid, rid) for rid in common}
        self._rhea_ec_list = [self._ec_map_full.get(master_of[rid], []) for rid in common]

        # EC prefix index (used for weight tuning)
        self._prefix_idx = {t: defaultdict(set) for t in range(1, 5)}
        for i, ecs in enumerate(self._rhea_ec_list):
            for ec in ecs:
                parts = ec.split(".")
                if len(parts) == 4:
                    for t in range(1, 5):
                        self._prefix_idx[t][".".join(parts[:t])].add(i)

        ec_covered = sum(1 for ec in self._rhea_ec_list if ec)
        self._log(f"  DB: {len(common)} reactions loaded  "
                  f"({ec_covered} with EC, {len(common)-ec_covered} without)")

    def _get_alpha(self, reestimate: bool) -> float:
        if not reestimate and ALPHA_CACHE.exists():
            with open(ALPHA_CACHE) as f:
                cached = json.load(f)
            alpha = cached["alpha"]
            self._log(f"  α = {alpha:.2f} (loaded from cache, "
                      f"tuned on {cached.get('n_samples','?')} samples)")
            return alpha

        self._log("Estimating optimal α …")
        return self._estimate_alpha()

    def _estimate_alpha(self) -> float:
        random.seed(SEED)
        samples = []

        # source 1: RHEA self-retrieval
        annotated = [(i, cid) for i, cid in enumerate(self._common_ids)
                     if self._rhea_ec_list[i]]
        for db_idx, cid in random.sample(annotated, min(TUNE_PER_SOURCE, len(annotated))):
            smi = self._smi_map.get(cid)
            if not smi:
                continue
            qm, qd = _encode_reaction(smi), _encode_drfp(smi)
            if qm is None or qd is None:
                continue
            pos = set()
            for ec in self._rhea_ec_list[db_idx]:
                parts = ec.split(".")
                if len(parts) == 4:
                    pos |= self._prefix_idx[4].get(".".join(parts), set())
            pos.discard(db_idx)
            if pos:
                samples.append((qm, qd, db_idx, pos))

        # source 2: CYP validation data
        if CYP_DATA.exists():
            with open(CYP_DATA) as f:
                cyp = json.load(f)
            random.shuffle(cyp)
            added = 0
            for rec in cyp:
                if added >= TUNE_PER_SOURCE:
                    break
                smi = strip_atom_map(rec.get("orig_rxn_text", ""))
                ec  = rec.get("ec", "")
                if not smi or not ec or "-" in ec:
                    continue
                qm, qd = _encode_reaction(smi), _encode_drfp(smi)
                if qm is None or qd is None:
                    continue
                pos = self._prefix_idx[4].get(ec, set())
                if pos:
                    samples.append((qm, qd, -1, pos))
                    added += 1

        self._log(f"  Tuning on {len(samples)} samples …")

        # pre-compute similarity vectors once, then grid-search α
        precomp = [(_tanimoto(self._morgan_fps, qm),
                    _tanimoto(self._drfp_fps,   qd),
                    qi, pos)
                   for qm, qd, qi, pos in samples]

        best_a, best_mrr = 0.0, -1.0
        for a in np.linspace(0.0, 1.0, ALPHA_STEPS):
            mrr = float(np.mean([_mrr(a * sm + (1 - a) * sd, qi, pos)
                                 for sm, sd, qi, pos in precomp]))
            if mrr > best_mrr:
                best_mrr, best_a = mrr, float(a)

        self._log(f"  Best α = {best_a:.2f}  "
                  f"(Morgan {best_a:.0%}, DRFP {1-best_a:.0%})")
        with open(ALPHA_CACHE, "w") as f:
            json.dump({"alpha": best_a, "n_samples": len(samples)}, f)
        return best_a

    def _ec_consensus(self, top_idxs) -> float:
        top1_ecs = self._rhea_ec_list[top_idxs[0]]
        if not top1_ecs:
            return 0.0
        pred  = top1_ecs[0].split(".")[0]
        votes = []
        for idx in top_idxs[:10]:
            for ec in self._rhea_ec_list[idx]:
                votes.append(ec.split(".")[0])
                break
        return votes.count(pred) / len(votes) if votes else 0.0


# ── CLI ───────────────────────────────────────────────────────────────────────

def _print_results(results, query_smi, alpha):
    conf       = results[0]["confidence"]
    conf_label = "✓ high confidence" if conf >= 0.78 else "⚠ low confidence"
    print(f"\nQuery : {query_smi}")
    print(f"Method: {alpha:.0%}·Morgan + {1-alpha:.0%}·DRFP  "
          f"| Confidence: {conf:.0%}  {conf_label}\n{'─'*72}")
    for r in results:
        ec_str = ", ".join(r["ec"]) if r["ec"] else "—"
        print(f"  #{r['rank']}  RHEA:{r['rhea_id']:<8} "
              f"score={r['score']:.4f}  (M={r['score_morgan']:.3f}  D={r['score_drfp']:.3f})")
        print(f"      EC : {ec_str}")
        print(f"      URL: {r['url']}")
        for p in r["publications"]:
            title = p["title"][:70] + "…" if len(p["title"]) > 70 else p["title"]
            print(f"      [{p['year']}] {title}")
            print(f"             PMID:{p['pmid']}  {p['journal']}")
        print()
    print(f"{'─'*72}")


def main():
    QUERIES = [
        "CC(=O)O.[H][H]>>CC(O)",
    ]

    searcher = RheaSearcher()
    for smiles in QUERIES:
        try:
            t0      = time.perf_counter()
            results = searcher.find(smiles, top_k=5, fetch_pubs=True)
            elapsed = (time.perf_counter() - t0) * 1000
            _print_results(results, smiles, searcher.alpha)
            print(f"  Search completed in {elapsed:.0f} ms\n")
        except ValueError as e:
            print(f"  Error: {e}\n")


if __name__ == "__main__":
    main()
