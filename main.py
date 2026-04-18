import os
import requests
from tqdm import tqdm
import zipfile
import tarfile 

import os
import sys
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

logger.info("=== Starting up ===")
logger.info(f"Python: {sys.version}")
logger.info(f"Files present: {os.listdir('files') if os.path.exists('files') else 'files/ missing'}")

# Download files if missing
logger.info("Checking model files...")
# ... your download code ...

logger.info("Loading screening set...")
# ... your pickle.load ...

logger.info("Loading CLIPZyme model...")
# ... your CLIPZyme load ...

logger.info("=== Ready ===")




FILES = {
    "files/clipzyme_model.ckpt": "https://zenodo.org/records/15161343/files/clipzyme_model.zip?download=1",
    "files/clipzyme_screening_set.p": "https://zenodo.org/records/15161343/files/clipzyme_data.zip?download=1",
    "files/brenda_2026_1.json":"https://zenodo.org/records/19475027/files/brenda_2026_1.json.tar.gz?download=1",
}



def download_if_missing(path: str, url: str):
    if os.path.exists(path):
        return
    os.makedirs("files", exist_ok=True)

    # Determine archive type from URL
    is_targz = ".tar.gz" in url or ".tgz" in url
    is_zip   = ".zip" in url

    archive_path = path + (".tar.gz" if is_targz else ".zip")

    print(f"Downloading {url}...")
    resp = requests.get(url, stream=True)
    with open(archive_path, "wb") as f:
        for chunk in resp.iter_content(chunk_size=1024 * 1024):
            f.write(chunk)

    if is_targz:
        print("Extracting .tar.gz...")
        with tarfile.open(archive_path, "r:gz") as t:
            t.extractall("files/")
        os.remove(archive_path)

    elif is_zip:
        print("Extracting .zip...")
        with zipfile.ZipFile(archive_path, "r") as z:
            z.extractall("files/")
        os.remove(archive_path)

    else:
        print(f"Unknown archive format for {url}, skipping extraction")

    print(f"Done: {path}")

# Download before loading model
for path, url in FILES.items():
    download_if_missing(path, url)

# Now import the rest — model loads after files are present
from query import query_enzymes
from uniprot import enrich_results



from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional

from query import query_enzymes
from enrich import enrich_results
from simi_search import RheaSearcher

# ── App ───────────────────────────────────────────────────────────────────────

app = FastAPI(
    title="Enzyme Screening API",
    description="Given substrate and product SMILES, returns ranked enzyme candidates with kinetic data.",
    version="1.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── Request / Response models ─────────────────────────────────────────────────

class ScreenRequest(BaseModel):
    substrate_smiles: str
    product_smiles:   str
    top_k:            Optional[int] = 20
    enrich:           Optional[bool] = True

    class Config:
        json_schema_extra = {
            "example": {
                "substrate_smiles": "CC(O)C",
                "product_smiles":   "CC(=O)C",
                "top_k":            20,
                "enrich":           True,
            }
        }


class ReactionRequest(BaseModel):
    substrate_smiles: str
    product_smiles:   str
    top_k:            Optional[int] = 5
    pub_fetch:           Optional[bool] = True

    class Config:
        json_schema_extra = {
            "example": {
                "substrate_smiles": "CC(O)C",
                "product_smiles":   "CC(=O)C",
                "top_k":            5,
                "pub_fetch":           True,
            }
        }


# ── Endpoints ─────────────────────────────────────────────────────────────────

@app.get("/health")
def health():
    return {"status": "ok"}


@app.get("/health/brenda")
def health_brenda():
    """Check whether brenda_2026_1.json is present, readable and valid."""
    path = "files/brenda_2026_1.json"

    # Check file exists
    if not os.path.exists(path):
        raise HTTPException(status_code=503, detail="brenda_2026_1.json not found — not downloaded yet")

    # Check file size (should be at least 50MB)
    size_mb = os.path.getsize(path) / (1024 * 1024)
    if size_mb < 50:
        raise HTTPException(
            status_code=503,
            detail=f"brenda_2026_1.json looks incomplete — only {size_mb:.1f} MB (expected >50 MB)"
        )

    # Check it is valid JSON and has expected structure
    try:
        import json
        with open(path) as f:
            data = json.load(f)

        if "data" not in data:
            raise HTTPException(status_code=503, detail="brenda_2026_1.json missing 'data' key — file may be corrupt")

        ec_count      = len(data["data"])
        sample_ec     = next(iter(data["data"]))  # first EC number
        protein_count = sum(len(e.get("protein", {})) for e in data["data"].values())

        return {
            "status":        "ok",
            "size_mb":       round(size_mb, 1),
            "ec_numbers":    ec_count,
            "proteins":      protein_count,
            "sample_ec":     sample_ec,
            "release":       data.get("release", "unknown"),
        }

    except json.JSONDecodeError as e:
        raise HTTPException(status_code=503, detail=f"brenda_2026_1.json is not valid JSON: {e}")
    except Exception as e:
        raise HTTPException(status_code=503, detail=f"brenda_2026_1.json check failed: {e}")


@app.post("/screen")
def screen(req: ScreenRequest):
    """
    Find enzymes likely to catalyse the given reaction.

    Returns
    -------
    {
        "status":   "success" | "error",
        "result":   [ enzyme dicts ] | "error description",
        "comments": [ "note 1", ... ]
    }
    """
    # Step 1: CLIPZyme ranking
    output = query_enzymes(
        substrate_smiles = req.substrate_smiles,
        product_smiles   = req.product_smiles,
        top_k            = req.top_k,
    )

    if output["status"] == "error":
        return output

    # Step 2: enrich with UniProt + SABIO-RK
    if req.enrich:
        output = enrich_results(output)

    return output




@app.post("/reaction")
def simi_find(req: ReactionRequest):
    """
    Find enzymes likely to catalyse the given reaction.

    Returns
    -------
    {
        "status":   "success" | "error",
        "result":   [ enzyme dicts ] | "error description",
        "comments": [ "note 1", ... ]
    }
    """

    searcher = RheaSearcher()

    results = searcher.find(
        smiles = req.substrate_smiles + ">>" + req.product_smiles,
        top_k=req.top_k,
        fetch_pubs=req.pub_fetch)

    return results