# NitroCat Backend — Enzyme Screening API

A FastAPI service that, given a reaction (substrate + product SMILES), returns ranked
enzyme candidates with kinetic and annotation data. Enzyme ranking is powered by
[CLIPZyme](https://github.com/pgmikhael/CLIPZyme) (and the
[CoBaCo](https://github.com/pravoslavzilka/clipzyme_CoBaCo) fine-tuned variant); results are
enriched from UniProt, SABIO-RK and BRENDA, with a complementary Rhea reaction-similarity
search.

## Layout

```
app/                 # deployed FastAPI package (this is what ships)
├── main.py          #   FastAPI app + endpoints + startup model/data download
├── query.py         #   CLIPZyme screening (loads model + screening set on import)
├── uniprot.py       #   UniProt / SABIO-RK / BRENDA fetchers
├── enrich.py        #   enrich ranked hits with annotation + kinetic data
├── brenda.py        #   BRENDA JSON lookups
├── enzyme_finder.py #   query + enrichment helper
├── simi_core.py     #   Rhea fingerprint cache + similarity primitives
└── simi_search.py   #   Rhea reaction-similarity search (RheaSearcher)
dev/                 # one-off scripts: benchmarking, dataset prep, tests (NOT deployed)
files/               # model checkpoint, screening set, BRENDA JSON (git-ignored, see below)
data_processing/     # local data-prep artifacts (large blobs git-ignored)
Dockerfile           # build + run (uvicorn app.main:app)
```

## API

Base URL: the deployed service root. Interactive docs at `/docs`.

| Method | Path | Body | Description |
|--------|------|------|-------------|
| `GET`  | `/health` | — | Liveness check (`{"status":"ok"}`). |
| `GET`  | `/health/brenda` | — | Verifies `brenda_2026_1.json` is present, complete (>50 MB) and valid; returns EC/protein counts. |
| `POST` | `/screen` | `ScreenRequest` | CLIPZyme enzyme ranking for a reaction, optionally enriched with UniProt/SABIO-RK. |
| `POST` | `/reaction` | `ReactionRequest` | Rhea reaction-similarity search for the reaction. |

**`ScreenRequest`**
```json
{ "substrate_smiles": "CC(O)C", "product_smiles": "CC(=O)C", "top_k": 20, "enrich": true }
```

**`ReactionRequest`**
```json
{ "substrate_smiles": "CC(O)C", "product_smiles": "CC(=O)C", "top_k": 5, "pub_fetch": true }
```

Responses follow `{ "status": "success"|"error", "result": [...], "comments": [...] }`.

## Running locally

Dependencies are heavy (PyTorch CPU, PyTorch-Geometric, CLIPZyme, RDKit, rxnmapper). The
versions are pinned in the [`Dockerfile`](Dockerfile); `requirements.txt` lists the Python
packages. The simplest reproducible path is the container:

```bash
docker build -t nitrocat-backend .
docker run -p 8000:8000 nitrocat-backend
# -> http://localhost:8000/docs
```

To run without Docker, install the pinned deps into a Python 3.10 env, then from the repo
root (so the `app` package and `files/` resolve correctly):

```bash
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

> Run from the **repository root**, not from inside `app/`: the service resolves `files/`
> and `rhea_cache/` relative to the project root.

## Model & data files

On startup `app/main.py` downloads the required artifacts into `files/` if missing:

- `clipzyme_model.ckpt` and `clipzyme_screening_set.p` — from Zenodo
  ([record 15161343](https://doi.org/10.5281/zenodo.11187747)).
- `brenda_2026_1.json` — from Zenodo (record 19475027).

These (and other large blobs under `files/` and `data_processing/`) are **not** stored in
git — see [`.gitignore`](.gitignore). The Rhea similarity cache is built at runtime under
`rhea_cache/`.

## Deployment

Containerised via the [`Dockerfile`](Dockerfile) and deployed on Railway. The image installs
pinned dependencies, copies the `app/` package, and runs `uvicorn app.main:app` on port 8000.
`.railwayignore` keeps `dev/`, `data_processing/`, `files/` and caches out of the upload.

## Related

- **Model:** [`clipzyme_CoBaCo`](https://github.com/pravoslavzilka/clipzyme_CoBaCo) — the
  CoBaCo (Constraint Batch Construction) fork of CLIPZyme used to train the screening model.
