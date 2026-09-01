# =========================================================
# BASE IMAGE
# =========================================================
FROM python:3.10-slim

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
# SYSTEM DEPENDENCIES (optimized for build speed)
# =========================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
    ncbi-blast+ \
    git \
    build-essential \
    curl \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# =========================================================
# PYTHON DEPENDENCIES
# =========================================================
COPY requirements.txt .
RUN pip install --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r requirements.txt gunicorn

# =========================================================
# TEMPORARY DIRECTORIES (created early for consistency)
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
# BLAST DATABASE (protein BLAST format)
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
# PORT EXPOSURE (Render uses PORT environment variable)
# =========================================================
EXPOSE 10000

# =========================================================
# HEALTH CHECK (Render-compatible)
# =========================================================
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -f http://127.0.0.1:${PORT}/health || exit 1

# =========================================================
# NON-ROOT USER (optional but recommended for Render)
# =========================================================
RUN useradd -m -u 1000 appuser && chown -R appuser:appuser /app
USER appuser

# =========================================================
# START APPLICATION (Render-compatible)
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
