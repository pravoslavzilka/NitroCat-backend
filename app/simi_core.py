import requests
import numpy as np
import json
import threading
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import os

BASE_DIR  = Path(__file__).resolve().parent.parent
CACHE_DIR = BASE_DIR / "rhea_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

FPS_PATH       = CACHE_DIR / "rhea_fps.npy"
META_PATH      = CACHE_DIR / "rhea_meta.json"
EC_PATH        = CACHE_DIR / "rhea_ec.json"
DIR_PATH       = CACHE_DIR / "rhea_dir.json"
DRFP_FPS_CACHE = CACHE_DIR / "bench_fps_drfp.npy"
DRFP_IDS_CACHE = CACHE_DIR / "bench_ids_drfp.json"

RHEA_URL          = "https://ftp.expasy.org/databases/rhea/tsv/rhea-reaction-smiles.tsv"
RHEA_EC_URL       = "https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv"
RHEA_DIR_URL      = "https://ftp.expasy.org/databases/rhea/tsv/rhea-directions.tsv"
RHEA_UNIPROT_URL  = "https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot_sprot.tsv"
UNIPROT_API_URL   = "https://rest.uniprot.org/uniprotkb"
UNIPROT_PATH      = CACHE_DIR / "rhea_uniprot.json"

NUM_WORKERS   = os.cpu_count()
MORGAN_RADIUS = 2
MORGAN_NBITS  = 2048


# ---------------------------------------------------------------------------
# Fingerprint encoding  (module-level so ProcessPoolExecutor can pickle them)
# ---------------------------------------------------------------------------

_morgan_gen = None

def _get_morgan_gen():
    global _morgan_gen
    if _morgan_gen is None:
        from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
        _morgan_gen = GetMorganGenerator(radius=MORGAN_RADIUS, fpSize=MORGAN_NBITS)
    return _morgan_gen


def _side_fp(smiles_part):
    """OR-fold Morgan FPs for all dot-separated molecules on one reaction side."""
    gen = _get_morgan_gen()
    acc = None
    for smi in smiles_part.split("."):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        Chem.SetBondStereoFromDirections(mol)
        arr = gen.GetFingerprintAsNumPy(mol).astype(np.uint8)
        acc = arr if acc is None else (acc | arr)
    return acc


def _encode_reaction(smi):
    """
    Encode a reaction SMILES as a (3*MORGAN_NBITS,) bool array:
        [reactants_fp | products_fp | diff_fp]
    Returns None if the SMILES is invalid.
    """
    try:
        parts = smi.split(">>")
        if len(parts) != 2:
            return None
        r_fp = _side_fp(parts[0])
        p_fp = _side_fp(parts[1])
        if r_fp is None or p_fp is None:
            return None
        diff_fp = r_fp ^ p_fp
        return np.concatenate([r_fp, p_fp, diff_fp]).astype(bool)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Fingerprint cache
# ---------------------------------------------------------------------------

def build_fingerprint_cache():
    print("Downloading RHEA reaction SMILES...")
    resp = requests.get(RHEA_URL)
    resp.raise_for_status()
    lines = resp.text.strip().split("\n")[1:]

    raw_ids, raw_smiles = [], []
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        raw_ids.append(parts[0])
        raw_smiles.append(parts[1])

    print(f"Downloaded {len(raw_ids)} reactions. "
          f"Encoding Morgan FPs (parallel, {NUM_WORKERS} workers)...")

    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as ex:
        encoded = list(ex.map(_encode_reaction, raw_smiles, chunksize=200))

    final_ids    = [rid for rid, enc in zip(raw_ids, encoded)    if enc is not None]
    final_smiles = [smi for smi, enc in zip(raw_smiles, encoded) if enc is not None]
    fps_array    = np.array([enc for enc in encoded if enc is not None], dtype=bool)

    print(f"Encoded {len(final_ids)} / {len(raw_ids)} reactions.")
    np.save(FPS_PATH, fps_array)
    with open(META_PATH, "w") as f:
        json.dump({"ids": final_ids, "smiles": final_smiles}, f)
    print(f"Fingerprint cache saved → {CACHE_DIR}/")
    return fps_array, final_ids, final_smiles


def load_fingerprint_cache():
    fps = np.load(FPS_PATH)
    with open(META_PATH) as f:
        meta = json.load(f)
    return fps, meta["ids"], meta["smiles"]


def build_drfp_cache():
    from drfp import DrfpEncoder
    _, rhea_ids, rhea_smiles = load_fingerprint_cache()
    print(f"Encoding DRFP for {len(rhea_ids)} reactions...")
    fps, ids = [], []
    for rid, smi in zip(rhea_ids, rhea_smiles):
        try:
            fp = np.array(DrfpEncoder.encode([smi])[0], dtype=np.uint8)
            fps.append(fp)
            ids.append(rid)
        except Exception:
            pass
    fps_array = np.array(fps, dtype=np.uint8)
    np.save(DRFP_FPS_CACHE, fps_array)
    with open(DRFP_IDS_CACHE, "w") as f:
        json.dump(ids, f)
    print(f"DRFP cache saved → {len(ids)} reactions")
    return fps_array, ids


# ---------------------------------------------------------------------------
# EC / direction cache
# ---------------------------------------------------------------------------

def build_ec_cache():
    print("Downloading RHEA EC mappings...")
    resp = requests.get(RHEA_EC_URL)
    resp.raise_for_status()
    ec_map = {}
    for line in resp.text.strip().split("\n")[1:]:
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        master_id, ec = parts[2], parts[3]
        ec_map.setdefault(master_id, [])
        if ec not in ec_map[master_id]:
            ec_map[master_id].append(ec)

    print("Downloading RHEA direction mappings...")
    resp2 = requests.get(RHEA_DIR_URL)
    resp2.raise_for_status()
    dir_map = {}
    for line in resp2.text.strip().split("\n")[1:]:
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        master, lr, rl, bi = parts[0], parts[1], parts[2], parts[3]
        for rid in (master, lr, rl, bi):
            if rid:
                dir_map[rid] = master

    with open(EC_PATH, "w") as f:
        json.dump(ec_map, f)
    with open(DIR_PATH, "w") as f:
        json.dump(dir_map, f)
    print(f"EC/direction cache saved → {CACHE_DIR}/")
    return ec_map, dir_map


def load_ec_cache():
    with open(EC_PATH) as f:
        ec_map = json.load(f)
    with open(DIR_PATH) as f:
        dir_map = json.load(f)
    return ec_map, dir_map


# ---------------------------------------------------------------------------
# UniProt mapping cache  (master RHEA ID → [Swiss-Prot accessions])
# Also builds an EC supplement for reactions not covered by rhea2ec.tsv
# ---------------------------------------------------------------------------

def _batch_fetch_ec(accessions: list[str]) -> dict[str, list[str]]:
    """Return {accession: [ec_numbers]} for a batch of UniProt accessions."""
    if not accessions:
        return {}
    query = " OR ".join(f"accession:{a}" for a in accessions)
    try:
        resp = requests.get(
            f"{UNIPROT_API_URL}/stream",
            params={"query": query, "format": "tsv", "fields": "accession,ec"},
            timeout=60,
        )
        resp.raise_for_status()
        result = {}
        for line in resp.text.strip().split("\n")[1:]:
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            acc, ec_raw = parts[0], parts[1]
            ecs = [e.strip() for e in ec_raw.split(";")
                   if e.strip() and "-" not in e.strip().split(".")[-1]]
            if ecs:
                result[acc] = ecs
        return result
    except Exception:
        return {}


def build_uniprot_cache(ec_map: dict | None = None):
    """
    Download rhea2uniprot_sprot.tsv and build:
      - uniprot_map: master_id → [accessions]
      - ec_supplement: master_id → [EC numbers] for reactions missing from rhea2ec.tsv

    Pass ec_map (from load_ec_cache) to skip reactions already covered.
    """
    print("Downloading RHEA → UniProt Swiss-Prot mapping...")
    resp = requests.get(RHEA_UNIPROT_URL)
    resp.raise_for_status()
    uniprot_map: dict[str, list[str]] = {}
    for line in resp.text.strip().split("\n")[1:]:
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        master_id, acc = parts[2], parts[3]
        if master_id not in uniprot_map:
            uniprot_map[master_id] = []
        if acc not in uniprot_map[master_id]:
            uniprot_map[master_id].append(acc)

    # Build EC supplement: for master IDs with no RHEA EC, fetch EC from UniProt
    print("Building EC supplement from UniProt (batch queries)...")
    ec_map = ec_map or {}
    # Collect one representative accession per uncovered master ID
    acc_to_masters: dict[str, list[str]] = {}
    for master_id, accs in uniprot_map.items():
        if not ec_map.get(master_id) and accs:
            rep = accs[0]
            acc_to_masters.setdefault(rep, []).append(master_id)

    all_rep_accs = list(acc_to_masters.keys())
    ec_by_acc: dict[str, list[str]] = {}
    batch_size = 50
    n_batches = (len(all_rep_accs) - 1) // batch_size + 1
    for i in range(0, len(all_rep_accs), batch_size):
        batch = all_rep_accs[i : i + batch_size]
        ec_by_acc.update(_batch_fetch_ec(batch))
        if (i // batch_size + 1) % 10 == 0 or i + batch_size >= len(all_rep_accs):
            print(f"  EC batches: {i // batch_size + 1}/{n_batches}")

    ec_supplement: dict[str, list[str]] = {}
    for acc, master_ids in acc_to_masters.items():
        ecs = ec_by_acc.get(acc, [])
        if ecs:
            for mid in master_ids:
                ec_supplement[mid] = ecs

    data = {"uniprot": uniprot_map, "ec_supplement": ec_supplement}
    with open(UNIPROT_PATH, "w") as f:
        json.dump(data, f)
    covered = len(ec_map) + len(ec_supplement)
    print(f"UniProt cache saved → {UNIPROT_PATH}  "
          f"(EC supplement: +{len(ec_supplement)} reactions)")
    return uniprot_map, ec_supplement


def load_uniprot_cache() -> tuple[dict, dict]:
    """Return (uniprot_map, ec_supplement). Handles both old and new cache formats."""
    with open(UNIPROT_PATH) as f:
        data = json.load(f)
    if isinstance(data, dict) and "uniprot" in data:
        return data["uniprot"], data.get("ec_supplement", {})
    # Legacy format: bare uniprot_map, no EC supplement
    return data, {}


# ---------------------------------------------------------------------------
# Publications via UniProt (accurate, enzyme-linked papers)
# ---------------------------------------------------------------------------

def _fetch_publications(master_id, uniprot_accs, page_size=5):
    """Fetch publications for a RHEA reaction via its linked UniProt proteins."""
    pubs = []
    seen_pmids = set()
    for acc in uniprot_accs[:3]:
        try:
            resp = requests.get(
                f"{UNIPROT_API_URL}/{acc}",
                params={"format": "json"},
                timeout=15,
            )
            resp.raise_for_status()
            data = resp.json()
            for ref in data.get("references", []):
                cit = ref.get("citation", {})
                if cit.get("citationType") != "journal article":
                    continue
                pmid = cit.get("id", "")
                if not pmid or pmid.startswith("CI-") or pmid in seen_pmids:
                    continue
                seen_pmids.add(pmid)
                pubs.append({
                    "pmid":    pmid,
                    "title":   cit.get("title", ""),
                    "authors": ", ".join(cit.get("authors", [])[:3]),
                    "journal": cit.get("journal", ""),
                    "year":    str(cit.get("publicationDate", ""))[:4],
                })
                if len(pubs) >= page_size:
                    break
        except Exception:
            pass
        if len(pubs) >= page_size:
            break
    return master_id, pubs


def fetch_publications_parallel(master_ids, uniprot_map, page_size=5):
    """Fetch publications for multiple RHEA master IDs in parallel via UniProt."""
    pub_map = {}
    args = [(mid, uniprot_map.get(mid, []), page_size) for mid in master_ids]
    with ThreadPoolExecutor(max_workers=min(len(master_ids), 10)) as ex:
        futures = {ex.submit(_fetch_publications, *a): a[0] for a in args}
        for fut in futures:
            mid, pubs = fut.result()
            pub_map[mid] = pubs
    return pub_map


# ---------------------------------------------------------------------------
# Similarity search
# ---------------------------------------------------------------------------

def find_most_similar(query_smiles, rhea_fps, rhea_ids, rhea_smiles,
                      ec_map, dir_map, top_k=5, pub_page_size=5):
    q_fp = _encode_reaction(query_smiles)
    if q_fp is None:
        raise ValueError(f"Could not encode query SMILES: {query_smiles}")

    # Fully vectorized Tanimoto: (N,D) & (D,) broadcasts to (N,D)
    intersection = np.sum(rhea_fps & q_fp, axis=1)
    union        = np.sum(rhea_fps | q_fp, axis=1)
    sims         = np.divide(intersection, union,
                             out=np.zeros(len(rhea_fps), dtype=float),
                             where=union != 0)

    top_indices  = np.argsort(sims)[::-1][:top_k]
    top_rhea_ids = [rhea_ids[i] for i in top_indices]

    print(f"Fetching publications for top-{top_k} reactions...")
    pub_map = fetch_publications_parallel(top_rhea_ids, page_size=pub_page_size)

    results = []
    for i in top_indices:
        rid       = rhea_ids[i]
        master_id = dir_map.get(rid, rid)
        results.append({
            "rhea_id":      rid,
            "smiles":       rhea_smiles[i],
            "score":        float(sims[i]),
            "url":          f"https://www.rhea-db.org/rhea/{rid}",
            "ec_numbers":   ec_map.get(master_id, []),
            "publications": pub_map.get(rid, []),
        })
    return results


# ---------------------------------------------------------------------------
# First-run setup
# ---------------------------------------------------------------------------

def setup_caches():
    """Build any missing cache files. Safe to call repeatedly (skips existing files)."""
    # 1. Morgan fingerprint cache — must come before DRFP (provides SMILES)
    if not FPS_PATH.exists() or not META_PATH.exists():
        print("Morgan fingerprint cache missing — downloading and building...")
        build_fingerprint_cache()

    # 2. DRFP fingerprint cache — built from Morgan SMILES
    if not DRFP_FPS_CACHE.exists() or not DRFP_IDS_CACHE.exists():
        print("DRFP fingerprint cache missing — building...")
        build_drfp_cache()

    # 3. EC / direction cache
    if not EC_PATH.exists() or not DIR_PATH.exists():
        print("EC/direction cache missing — downloading...")
        build_ec_cache()

    # 4. UniProt + EC supplement cache — needs EC cache for gap detection
    if not UNIPROT_PATH.exists():
        print("UniProt cache missing — downloading and building...")
        ec_map, _ = load_ec_cache()
        build_uniprot_cache(ec_map)


_cache_ready = threading.Event()

def _background_setup():
    try:
        setup_caches()
    finally:
        _cache_ready.set()

threading.Thread(target=_background_setup, daemon=True, name="rhea-cache-setup").start()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if FPS_PATH.exists() and META_PATH.exists():
        print("Loading cached fingerprints...")
        rhea_fps, rhea_ids, rhea_smiles = load_fingerprint_cache()
        print(f"Loaded {len(rhea_ids)} reactions from cache.")
    else:
        rhea_fps, rhea_ids, rhea_smiles = build_fingerprint_cache()

    if EC_PATH.exists() and DIR_PATH.exists():
        print("Loading cached EC mappings...")
        ec_map, dir_map = load_ec_cache()
    else:
        ec_map, dir_map = build_ec_cache()

    query = "CC(=O)CCC(C(=O)OC)C(=O)OC.N>>C[C@H](CCC(C(=O)OC)C(=O)OC)N"
    results = find_most_similar(query, rhea_fps, rhea_ids, rhea_smiles,
                                ec_map, dir_map, top_k=5, pub_page_size=3)

    print(f"\nTop matches for: {query}\n")
    for r in results:
        print(f"  RHEA:{r['rhea_id']}  score={r['score']:.4f}")
        print(f"  SMILES:  {r['smiles']}")
        print(f"  URL:     {r['url']}")
        ec_str = ", ".join(r['ec_numbers']) if r['ec_numbers'] else "—"
        print(f"  EC:      {ec_str}")
        if r['publications']:
            print(f"  Publications ({len(r['publications'])}):")
            for p in r['publications']:
                print(f"    [{p['year']}] {p['title']}")
                print(f"           {p['authors']}")
                print(f"           {p['journal']}"
                      + (f"  PMID:{p['pmid']}" if p['pmid'] else ""))
        else:
            print(f"  Publications: —")
        print()
