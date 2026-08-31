# =========================================================
# BASE IMAGE
# =========================================================

FROM python:3.10-slim

# =========================================================
# WORKING DIRECTORY
# =========================================================

WORKDIR /app

# =========================================================
# SYSTEM DEPENDENCIES
# =========================================================

RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    git \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# =========================================================
# PYTHON DEPENDENCIES
# =========================================================

COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

# =========================================================
# APPLICATION
# =========================================================

COPY app.py .

# =========================================================
# FRONTEND
# =========================================================

COPY static/ ./static/
COPY templates/ ./templates/

# =========================================================
# MACHINE LEARNING MODELS
# =========================================================

COPY label_encoder.pkl .
COPY best_model_LR.sav .

# =========================================================
# BLAST DATABASE
# =========================================================

COPY pathway_db.pdb .
COPY pathway_db.phr .
COPY pathway_db.pin .
COPY pathway_db.pot .
COPY pathway_db.psq .
COPY pathway_db.ptf .
COPY pathway_db.pto .

# =========================================================
# PATHWAY MAPPING
# =========================================================

COPY pathway_map.csv .

# =========================================================
# TEMPORARY DIRECTORIES
# =========================================================

RUN mkdir -p /app/data /app/temp /app/uploads

# =========================================================
# PORT
# =========================================================

EXPOSE 10000

# =========================================================
# HEALTH CHECK
# =========================================================

HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD python -c "import os, urllib.request; port=os.environ.get('PORT', '10000'); urllib.request.urlopen(f'http://127.0.0.1:{port}/health').read()"

# =========================================================
# START APPLICATION
# =========================================================

CMD ["sh", "-c", "gunicorn --bind 0.0.0.0:${PORT:-10000} --workers 1 --timeout 300 app:app"]
