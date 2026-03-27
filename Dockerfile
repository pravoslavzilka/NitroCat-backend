FROM python:3.10.14-slim

WORKDIR /app

RUN apt-get update && apt-get install -y \
    libxrender1 libxext6 gcc g++ git \
    && rm -rf /var/lib/apt/lists/*

RUN pip install torch==2.3.0 --index-url https://download.pytorch.org/whl/cpu

RUN pip install torch_cluster torch_scatter torch_sparse torch_geometric \
    -f https://data.pyg.org/whl/torch-2.3.0+cpu.html

RUN pip install clipzyme rxnmapper

RUN pip install fair-esm==2.0.0

RUN pip install p_tqdm molvs pytorch_lightning==2.0.0 wandb

RUN pip install fastapi uvicorn requests

COPY main.py query.py uniprot.py ./

RUN mkdir -p files

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]