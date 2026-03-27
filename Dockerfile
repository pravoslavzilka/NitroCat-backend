FROM python:3.10.14-slim

WORKDIR /app

# System dependencies
RUN apt-get update && apt-get install -y \
    libxrender1 libxext6 gcc g++ git \
    && rm -rf /var/lib/apt/lists/*

# PyTorch first — must be before torch_cluster/scatter/sparse
RUN pip install torch==2.3.0 --index-url https://download.pytorch.org/whl/cpu

# Torch geometric packages — pre-built wheels for torch 2.3.0 cpu
RUN pip install \
    torch_cluster==1.6.3 \
    torch_scatter==2.1.2 \
    torch_sparse==0.6.18 \
    torch-geometric==2.7.0 \
    -f https://data.pyg.org/whl/torch-2.3.0+cpu.html

# ESM must come before clipzyme
RUN pip install fair-esm==2.0.0

# Core clipzyme dependencies
RUN pip install \
    pytorch-lightning==2.0.0 \
    MolVS==0.1.1 \
    p_tqdm==1.4.2 \
    wandb==0.25.1

# CLIPZyme
RUN pip install clipzyme==0.0.12

# Reaction mapping
RUN pip install rxnmapper==0.4.3

# Chemistry
RUN pip install rdkit==2025.9.6

# API and utilities
RUN pip install \
    fastapi==0.135.2 \
    uvicorn==0.42.0 \
    requests==2.33.0 \
    zeep==4.3.2

# Copy app code
COPY main.py query.py uniprot.py ./

# Empty files dir — populated at startup from Zenodo
RUN mkdir -p files

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]