import os
import requests
from tqdm import tqdm
import zipfile

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
}



def download_if_missing(path: str, url: str):
    if os.path.exists(path):
        return
    os.makedirs("files", exist_ok=True)
    zip_path = path + ".zip"
    
    print(f"Downloading {url}...")
    resp = requests.get(url, stream=True)
    with open(zip_path, "wb") as f:
        for chunk in resp.iter_content(chunk_size=1024 * 1024):
            f.write(chunk)
    
    # Unzip if it is a zip file
    if url.endswith(".zip?download=1") or url.endswith(".zip"):
        print("Unzipping...")
        with zipfile.ZipFile(zip_path, "r") as z:
            z.extractall("files/")
        os.remove(zip_path)
    
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
from uniprot import enrich_results

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

# ── Endpoints ─────────────────────────────────────────────────────────────────

@app.get("/health")
def health():
    return {"status": "ok"}


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