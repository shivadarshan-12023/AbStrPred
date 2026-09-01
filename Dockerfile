
# BASE IMAGE - EXPLICITLY Python 3.11 (avoid 3.14 auto-select)
# =========================================================
FROM python:3.11.9-slim
 
# =========================================================
# ENVIRONMENT VARIABLES
# =========================================================
ENV PORT=10000 \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PYTHONVERSION=3.11
 
# =========================================================
# VERIFY PYTHON VERSION IMMEDIATELY
# =========================================================
RUN python --version && \
    python -c "import sys; assert sys.version_info[:2] == (3, 11), f'Expected Python 3.11, got {sys.version}'"
 
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
    ca-certificates \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
 
# =========================================================
# UPGRADE PIP FIRST (critical for wheel support)
# =========================================================
RUN pip install --upgrade \
    pip==24.0 \
    setuptools==69.0.0 \
    wheel==0.42.0
 
# =========================================================
# INSTALL CRITICAL PACKAGES AS PRE-BUILT BINARIES ONLY
# This must succeed before any other installs
# =========================================================
RUN pip install --upgrade --only-binary :all: \
    numpy==1.26.4 \
    pandas==2.1.4 \
    scipy==1.11.4 \
    scikit-learn==1.3.2
 
# =========================================================
# INSTALL REMAINING DEPENDENCIES
# =========================================================
COPY requirements.txt .
RUN pip install --no-cache-dir \
    --only-binary :all: \
    --prefer-binary \
    -r requirements.txt
 
# =========================================================
# VERIFY ALL IMPORTS WORK
# =========================================================
RUN python -c "\
import sys; \
print(f'Python: {sys.version}'); \
import numpy as np; print(f'✓ numpy {np.__version__}'); \
import pandas as pd; print(f'✓ pandas {pd.__version__}'); \
import sklearn; print(f'✓ sklearn {sklearn.__version__}'); \
import Bio; print(f'✓ biopython {Bio.__version__}'); \
import flask; print(f'✓ flask {flask.__version__}'); \
import gunicorn; print(f'✓ gunicorn {gunicorn.__version__}')"
 
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
# FINAL VERIFICATION
# =========================================================
RUN python -c "print('✓ Application ready to start')"
 
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
 
