FROM python:3.10-slim

WORKDIR /app

RUN apt-get update && apt-get install -y \
    libxrender1 libxext6 \
    && rm -rf /var/lib/apt/lists/*

RUN pip install torch==2.3.0 --index-url https://download.pytorch.org/whl/cpu

RUN pip install torch_cluster torch_scatter torch_sparse torch_geometric \
    -f https://data.pyg.org/whl/torch-2.3.0+cpu.html

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY main.py query.py uniprot.py ./

# Create empty files directory — will be populated at startup from Zenodo
RUN mkdir -p files

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]