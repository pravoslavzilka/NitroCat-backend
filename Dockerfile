FROM python:3.10-slim

WORKDIR /app

# Install system dependencies RDKit needs
RUN apt-get update && apt-get install -y \
    libxrender1 libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Install PyTorch first — must be before other torch packages
RUN pip install torch==2.3.0 --index-url https://download.pytorch.org/whl/cpu

# Install torch geometric packages against exact torch version
RUN pip install torch_cluster torch_scatter torch_sparse torch_geometric \
    -f https://data.pyg.org/whl/torch-2.3.0+cpu.html

# Install remaining dependencies
COPY requirements.txt .
RUN pip install -r requirements.txt

# Copy app code
COPY main.py query.py uniprot.py ./

# Copy model files — this is the 5GB, happens at build time once
COPY files/ ./files/

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]