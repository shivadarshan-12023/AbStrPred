# =========================================================
# BASE IMAGE
# =========================================================
FROM python:3.11-slim

# =========================================================
# ENVIRONMENT
# =========================================================
ENV PYTHONUNBUFFERED=1 \
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
    build-essential \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# =========================================================
# VERIFY PYTHON
# =========================================================
RUN python --version && \
    python -c "import sys; assert sys.version_info[:2] == (3, 11)"

# =========================================================
# UPGRADE PIP
# =========================================================
RUN pip install --upgrade pip setuptools wheel

# =========================================================
# COPY REQUIREMENTS
# =========================================================
COPY requirements.txt .

# =========================================================
# INSTALL PYTHON DEPENDENCIES
# =========================================================
RUN pip install --no-cache-dir \
    --prefer-binary \
    -r requirements.txt

# =========================================================
# VERIFY IMPORTANT PACKAGES
# =========================================================
RUN python -c "\
import numpy; \
import pandas; \
import sklearn; \
import Bio; \
import flask; \
import gunicorn; \
print('All required packages imported successfully')"

# =========================================================
# APPLICATION DIRECTORIES
# =========================================================
RUN mkdir -p /app/data /app/temp /app/uploads

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
# CREATE NON-ROOT USER
# =========================================================
RUN useradd -m -u 1000 appuser && \
    chown -R appuser:appuser /app

USER appuser

# =========================================================
# HEALTH CHECK
# =========================================================
HEALTHCHECK --interval=30s \
    --timeout=10s \
    --start-period=60s \
    --retries=3 \
    CMD curl -f http://127.0.0.1:${PORT}/health || exit 1

# =========================================================
# START GUNICORN
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
