# =========================================================
# BASE IMAGE - Use Python 3.11 (most stable with pandas)
# =========================================================
FROM python:3.11-slim

# =========================================================
# ENVIRONMENT VARIABLES
# =========================================================
ENV PORT=10000 \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1

# =========================================================
# WORKING DIRECTORY
# =========================================================
WORKDIR /app

# =========================================================
# SYSTEM DEPENDENCIES
# =========================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
    ncbi-blast+ \
    git \
    build-essential \
    curl \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# =========================================================
# PYTHON DEPENDENCIES - Use pre-built wheels (no compilation)
# =========================================================
COPY requirements.txt .
RUN pip install --upgrade pip setuptools wheel && \
    pip install --only-binary :all: \
    numpy \
    pandas \
    scikit-learn \
    biopython && \
    pip install --no-binary :none: gunicorn flask -r requirements.txt

# =========================================================
# TEMPORARY DIRECTORIES
# =========================================================
RUN mkdir -p /app/data /app/temp /app/uploads && \
    chmod 755 /app/data /app/temp /app/uploads

# =========================================================
# APPLICATION CODE
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
# PORT EXPOSURE
# =========================================================
EXPOSE 10000

# =========================================================
# HEALTH CHECK
# =========================================================
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -f http://127.0.0.1:${PORT}/health || exit 1

# =========================================================
# NON-ROOT USER
# =========================================================
RUN useradd -m -u 1000 appuser && chown -R appuser:appuser /app
USER appuser

# =========================================================
# START APPLICATION
# =========================================================
CMD exec gunicorn \
    --bind 0.0.0.0:${PORT} \
    --workers 2 \
    --worker-class sync \
    --timeout 300 \
    --access-logfile - \
    --error-logfile - \
    --log-level info \
    app:app
