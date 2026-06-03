FROM python:3.10.14-slim

WORKDIR /app

RUN apt-get update && apt-get install -y \
    libxrender1 libxext6 gcc g++ git \
    && rm -rf /var/lib/apt/lists/*

# PyTorch first
RUN pip install torch==2.3.0 --index-url https://download.pytorch.org/whl/cpu

# Torch geometric pre-built wheels
RUN pip install \
    torch_cluster==1.6.3 \
    torch_scatter==2.1.2 \
    torch_sparse==0.6.18 \
    torch-geometric==2.7.0 \
    -f https://data.pyg.org/whl/torch-2.3.0+cpu.html

# All remaining packages from pip freeze
RUN pip install \
    aiohttp==3.13.3 \
    annotated-types==0.7.0 \
    anyio==4.13.0 \
    attrs==26.1.0 \
    biopython==1.86 \
    biotite==1.2.0 \
    certifi==2026.2.25 \
    charset-normalizer==3.4.6 \
    click==8.3.1 \
    cloudpathlib==0.23.0 \
    decorator==5.2.1 \
    dill==0.4.1 \
    diskcache==5.6.3 \
    einops==0.8.2 \
    fair-esm==2.0.0 \
    fastapi==0.135.2 \
    filelock==3.25.2 \
    frozenlist==1.8.0 \
    fsspec==2026.2.0 \
    GitPython==3.1.46 \
    h11==0.16.0 \
    httpcore==1.0.9 \
    httpx==0.28.1 \
    huggingface_hub==0.36.2 \
    idna==3.11 \
    isodate==0.7.2 \
    Jinja2==3.1.6 \
    joblib==1.5.3 \
    lightning-utilities==0.15.3 \
    lxml==6.0.2 \
    markdown-it-py==4.0.0 \
    MarkupSafe==3.0.3 \
    MolVS==0.1.1 \
    mpmath==1.3.0 \
    msgpack==1.1.2 \
    msgpack-numpy==0.4.8 \
    multidict==6.7.1 \
    multiprocess==0.70.19 \
    networkx==3.4.2 \
    numpy==2.2.6 \
    p_tqdm==1.4.2 \
    packaging==24.2 \
    pandas==2.3.3 \
    pathos==0.3.5 \
    pillow==12.1.1 \
    platformdirs==4.9.4 \
    pooch==1.9.0 \
    pox==0.3.7 \
    ppft==1.7.8 \
    protobuf==6.33.6 \
    psutil==7.2.2 \
    pydantic==2.12.5 \
    pydantic-settings==2.13.1 \
    pyparsing==3.3.2 \
    python-dateutil==2.9.0.post0 \
    python-dotenv==1.2.2 \
    pytorch-lightning==2.0.0 \
    pytz==2026.1.post1 \
    PyYAML==6.0.3 \
    rdkit==2025.9.6 \
    regex==2026.2.28 \
    requests==2.33.0 \
    requests-file==3.0.1 \
    requests-toolbelt==1.0.0 \
    rich==14.3.3 \
    rxn-chem-utils==1.6.0 \
    rxn-utils==2.0.0 \
    rxnmapper==0.4.3 \
    safetensors==0.7.0 \
    scikit-learn==1.7.2 \
    scipy==1.15.3 \
    sentry-sdk==2.56.0 \
    six==1.17.0 \
    starlette==1.0.0 \
    sympy==1.14.0 \
    tenacity==9.1.4 \
    threadpoolctl==3.6.0 \
    tokenizers==0.21.4 \
    torchmetrics==1.9.0 \
    tqdm==4.67.3 \
    transformers==4.48.1 \
    typing_extensions==4.15.0 \
    tzdata==2025.3 \
    urllib3==2.6.3 \
    uvicorn==0.42.0 \
    wandb==0.25.1 \
    wget==3.2 \
    wcwidth==0.6.0 \
    xxhash==3.6.0 \
    yarl==1.23.0 \
    zeep==4.3.2 \
    zstd==1.5.7.3 \
    drfp==0.3.7

# CLIPZyme last — depends on everything above
RUN pip install clipzyme==0.0.12

# Copy app code (FastAPI package)
COPY app/ ./app/

RUN mkdir -p files

EXPOSE 8000

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]