import os
import io
import csv
import subprocess
import pandas as pd
import joblib
import numpy as np
import torch
import esm
import logging
import shutil
import tempfile
from flask import Flask, request, render_template, send_file, session, redirect, url_for

# ---------- Logging Setup ----------
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[logging.StreamHandler()]  # Only StreamHandler for HF Spaces
)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

MODEL_PATH = os.path.join(BASE_DIR, "best_model_LR.sav")
ENCODER_PATH = os.path.join(BASE_DIR, "label_encoder.pkl")

BLAST_DB = os.path.join(BASE_DIR, "pathway_db")
PATHWAY_MAP = os.path.join(BASE_DIR, "pathway_map.csv")

logger.info(f"BASE_DIR: {BASE_DIR}")
logger.info(f"MODEL_PATH: {MODEL_PATH}")
logger.info(f"ENCODER_PATH: {ENCODER_PATH}")
logger.info(f"BLAST_DB: {BLAST_DB}")
logger.info(f"PATHWAY_MAP: {PATHWAY_MAP}")

BITSCORE_THRESHOLD = 80.0
CONFIDENCE_THRESHOLD = 0.7

# Use temp directory for results in HF Spaces
RESULT_CSV = os.path.join(tempfile.gettempdir(), "combined_results.csv")

# ---------- Check and Setup BLAST ----------
def find_blast():
    """Find BLAST executable in common locations"""
    possible_paths = [
        "/usr/bin/blastp",
        "/usr/local/bin/blastp",
        "blastp",
        shutil.which("blastp")
    ]
    
    for path in possible_paths:
        if path and os.path.exists(path):
            logger.info(f"Found BLAST at: {path}")
            return path
    
    logger.warning("BLAST not found in common locations")
    return None

# Initialize BLAST_PATH properly BEFORE using it
BLAST_PATH = find_blast()

if BLAST_PATH:
    logger.info(f"BLAST executable: {BLAST_PATH}")
else:
    logger.warning("BLAST executable not found! BLAST searches will be disabled.")

# ---------- Load Models ----------
model = None
encoder = None
esm_model = None
batch_converter = None

try:
    if os.path.exists(MODEL_PATH):
        model = joblib.load(MODEL_PATH)
        logger.info("LR model loaded successfully.")
    else:
        logger.warning(f"Model file not found at {MODEL_PATH}")
except Exception as e:
    logger.exception(f"Failed to load LR model: {e}")

try:
    if os.path.exists(ENCODER_PATH):
        encoder = joblib.load(ENCODER_PATH)
        logger.info("Label encoder loaded.")
    else:
        logger.warning(f"Encoder file not found at {ENCODER_PATH}")
except Exception as e:
    logger.exception(f"Failed to load label encoder: {e}")

try:
    logger.info("Loading ESM2 model...")
    esm_model, alphabet = esm.pretrained.load_model_and_alphabet("esm2_t6_8M_UR50D")
    esm_model.eval()
    esm_model = esm_model.to("cpu")
    batch_converter = alphabet.get_batch_converter()
    logger.info("ESM2 model loaded and set to CPU.")
except Exception as e:
    logger.exception(f"Failed to load ESM model: {e}")

# ---------- Functions ----------
def esm2_320_embed(sequence):
    """Generate ESM2 embeddings for a protein sequence"""
    if esm_model is None or batch_converter is None:
        raise RuntimeError("ESM model not available")
    try:
        batch_labels, batch_strs, batch_tokens = batch_converter([("seq1", sequence)])
        with torch.no_grad():
            results = esm_model(batch_tokens, repr_layers=[6], return_contacts=False)
        token_representations = results["representations"][6]
        seq_repr = token_representations[0, 1:-1].detach().cpu().numpy()
        return seq_repr.mean(axis=0)
    except Exception as e:
        logger.error(f"Error generating embedding: {e}")
        raise

def check_blast_availability():
    """Check if BLAST is available and working"""
    if not BLAST_PATH or not os.path.exists(BLAST_PATH):
        logger.error(f"BLAST executable not found at {BLAST_PATH}")
        logger.error("To fix: install ncbi-blast+ in your requirements.txt")
        return False
    
    try:
        result = subprocess.run([BLAST_PATH, "-version"], capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            version_info = result.stdout.split()[1] if result.stdout else "version unknown"
            logger.info(f"BLAST is working correctly: {version_info}")
            return True
        else:
            logger.error(f"BLAST check failed: {result.stderr}")
            return False
    except Exception as e:
        logger.error(f"Error checking BLAST: {e}")
        return False

def check_blast_db():
    """Check if BLAST database exists and is properly formatted"""
    if not BLAST_DB:
        logger.error("BLAST_DB not configured")
        return False
    
    # Correct file extensions for BLAST protein databases
    # .pdb is NOT a BLAST extension (it's Protein Data Bank)
    # Correct extensions: .phr (header), .pin (index), .psq (sequence)
    db_files = {
        f"{BLAST_DB}.phr": "Header file",
        f"{BLAST_DB}.pin": "Index file", 
        f"{BLAST_DB}.psq": "Sequence file"
    }
    
    found_files = []
    missing_files = []
    
    for filepath, description in db_files.items():
        if os.path.exists(filepath):
            file_size = os.path.getsize(filepath)
            found_files.append(f"{os.path.basename(filepath)} ✓ ({description}, {file_size} bytes)")
        else:
            missing_files.append(f"{os.path.basename(filepath)} ✗ ({description})")
    
    if found_files:
        logger.info(f"BLAST database found at {BLAST_DB}:")
        for f in found_files:
            logger.info(f"  - {f}")
    
    if missing_files:
        logger.error(f"BLAST database incomplete at {BLAST_DB}. Missing:")
        for f in missing_files:
            logger.error(f"  - {f}")
        logger.error("Solution: Run 'makeblastdb -in your_sequences.fasta -dbtype prot -out pathway_db'")
        return False
    
    return len(found_files) == len(db_files)

def run_blast_and_get_dataframe(temp_fasta, blast_output):
    """Run BLAST search and return results as DataFrame"""
    
    # CHECK 1: BLAST executable exists
    if not BLAST_PATH or not os.path.exists(BLAST_PATH):
        logger.error(f"BLAST executable not found at {BLAST_PATH}")
        logger.error("FIX: Add 'ncbi-blast==2.14.1' to requirements.txt")
        return pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"])
    
    # CHECK 2: BLAST database exists and is formatted
    if not check_blast_db():
        logger.error(f"BLAST database problem at {BLAST_DB}")
        logger.error("Possible causes:")
        logger.error("  1. Database files don't exist")
        logger.error("  2. Database is not formatted (missing .phr, .pin, .psq files)")
        logger.error("  3. Wrong path configured")
        logger.error("FIX: See /debug_blast endpoint for current paths")
        return pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"])
    
    # CHECK 3: Query file exists and has content
    if not os.path.exists(temp_fasta) or os.path.getsize(temp_fasta) == 0:
        logger.error(f"Query FASTA file is empty or missing: {temp_fasta}")
        return pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"])
    
    try:
        # Build BLAST command
        cmd = [
            BLAST_PATH,
            "-query", temp_fasta,
            "-db", BLAST_DB,
            "-out", blast_output,
            "-outfmt", "6 qseqid sseqid pident length evalue bitscore",
            "-evalue", "1e-5"  # Default e-value threshold
        ]
        
        logger.info(f"Running BLAST command: {' '.join(cmd)}")
        logger.info(f"Query file: {temp_fasta} ({os.path.getsize(temp_fasta)} bytes)")
        logger.info(f"Database: {BLAST_DB}")
        
        # Run BLAST with timeout
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        
        # Log both stdout and stderr
        if proc.stdout:
            logger.info(f"BLAST stdout: {proc.stdout[:500]}")
        if proc.stderr:
            logger.warning(f"BLAST stderr: {proc.stderr[:500]}")
        
        # Check return code
        if proc.returncode != 0:
            logger.error(f"BLAST failed with return code {proc.returncode}")
            logger.error(f"Full stderr: {proc.stderr}")
            return pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"])
        
        # Check output file
        if not os.path.exists(blast_output):
            logger.error(f"BLAST did not create output file: {blast_output}")
            return pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"])
        
        output_size = os.path.getsize(blast_output)
        logger.info(f"BLAST output file size: {output_size} bytes")
        
        if output_size == 0:
            logger.warning("BLAST search completed but found NO HITS")
            logger.info("This could mean:")
            logger.info("  1. Query sequences are too divergent from database")
            logger.info("  2. E-value threshold is too strict")
            logger.info("  3. Database is empty or corrupted")
            return pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"])
        
        # Parse BLAST output
        cols = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]
        df = pd.read_csv(blast_output, sep="\t", names=cols)
        df["bitscore"] = df["bitscore"].astype(float)
        
        logger.info(f"BLAST completed successfully with {len(df)} hit(s)")
        if len(df) > 0:
            logger.info(f"Top hit: {df.iloc[0]['sseqid']} (bitscore: {df.iloc[0]['bitscore']:.2f})")
        
        return df
        
    except subprocess.TimeoutExpired:
        logger.error("BLAST search timed out (>300 seconds)")
        return pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"])
    except pd.errors.ParserError as e:
        logger.error(f"Failed to parse BLAST output file: {e}")
        logger.info(f"Output file location: {blast_output}")
        try:
            with open(blast_output, 'r') as f:
                content = f.read()
                logger.info(f"Raw output content: {content[:500]}")
        except:
            pass
        return pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"])
    except Exception as e:
        logger.exception(f"Unexpected error during BLAST: {e}")
        return pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"])

def run_startup_diagnostics():
    """Run comprehensive diagnostics at startup"""
    logger.info("=" * 70)
    logger.info("STARTUP DIAGNOSTICS - BLAST CONFIGURATION")
    logger.info("=" * 70)
    
    # Check BLAST executable
    logger.info("\n1. BLAST EXECUTABLE:")
    logger.info(f"   BLAST_PATH: {BLAST_PATH}")
    if BLAST_PATH:
        logger.info(f"   Exists: {os.path.exists(BLAST_PATH)}")
        if os.path.exists(BLAST_PATH):
            logger.info(f"   Absolute path: {os.path.abspath(BLAST_PATH)}")
            logger.info(f"   Size: {os.path.getsize(BLAST_PATH)} bytes")
    logger.info(f"   Available: {check_blast_availability()}")
    
    # Check database files
    logger.info("\n2. BLAST DATABASE:")
    logger.info(f"   BLAST_DB: {BLAST_DB}")
    logger.info(f"   Absolute path: {os.path.abspath(BLAST_DB) if BLAST_DB else 'N/A'}")
    
    if BLAST_DB and os.path.exists(os.path.dirname(BLAST_DB)):
        db_dir = os.path.dirname(BLAST_DB)
        logger.info(f"   Directory: {db_dir}")
        logger.info(f"   Files in directory:")
        try:
            files = os.listdir(db_dir)
            for f in sorted(files):
                fpath = os.path.join(db_dir, f)
                if os.path.isfile(fpath):
                    size = os.path.getsize(fpath)
                    logger.info(f"     - {f} ({size} bytes)")
        except Exception as e:
            logger.error(f"   Could not list directory: {e}")
    else:
        logger.error(f"   Database directory does not exist: {os.path.dirname(BLAST_DB)}")
    
    logger.info(f"   DB check result: {check_blast_db()}")
    
    # Check mapping file
    logger.info("\n3. PATHWAY MAPPING FILE:")
    logger.info(f"   PATHWAY_MAP: {PATHWAY_MAP}")
    logger.info(f"   Exists: {os.path.exists(PATHWAY_MAP)}")
    if os.path.exists(PATHWAY_MAP):
        try:
            map_df = pd.read_csv(PATHWAY_MAP)
            logger.info(f"   Rows: {len(map_df)}")
            logger.info(f"   Columns: {list(map_df.columns)}")
        except Exception as e:
            logger.error(f"   Error reading file: {e}")
    
    # Check models
    logger.info("\n4. ML MODELS:")
    logger.info(f"   LR Model: {model is not None}")
    logger.info(f"   Encoder: {encoder is not None}")
    logger.info(f"   ESM Model: {esm_model is not None}")
    
    logger.info("=" * 70)

# ---------- Flask app ----------
app = Flask(__name__)
app.secret_key = os.environ.get('SECRET_KEY', 'dev-secret-key-change-in-production')
latest_predictions = []

# Log startup diagnostics
run_startup_diagnostics()

@app.after_request
def add_cache_control(response):
    """Add cache control headers to prevent caching"""
    response.headers["Cache-Control"] = "no-cache, no-store, must-revalidate"
    response.headers["Pragma"] = "no-cache"
    response.headers["Expires"] = "0"
    return response

@app.route("/", methods=["GET", "POST"])
def index():
    global latest_predictions
    latest_predictions = []
    predictions = []

    if request.method == "POST":
        sequences = []
        headers = []

        try:
            # --- Parse uploaded file ---
            uploaded_file = request.files.get("fasta_file")
            if uploaded_file and uploaded_file.filename != "":
                seq = ""
                header = ""
                for line in uploaded_file:
                    line = line.decode().strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if seq:
                            sequences.append(seq)
                            headers.append(header)
                            seq = ""
                        header = line[1:]
                    else:
                        seq += line
                if seq:
                    sequences.append(seq)
                    headers.append(header)

            # --- Parse textarea ---
            sequence_text = request.form.get("sequence_text", "").strip()
            if sequence_text:
                seq = ""
                header = ""
                for line in sequence_text.splitlines():
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        if seq:
                            sequences.append(seq)
                            headers.append(header)
                            seq = ""
                        header = line[1:]
                    else:
                        seq += line
                if seq:
                    sequences.append(seq)
                    headers.append(header if header else "sequence_from_text")

            if not sequences:
                return "⚠️ No valid sequences provided.", 400

            logger.info(f"Processing {len(sequences)} sequences")

            # --- Gene family prediction with confidence check ---
            try:
                features = [esm2_320_embed(seq) for seq in sequences]
                features = np.array(features)
                gene_preds = []
                
                if model is None:
                    logger.warning("Model not loaded, using 'Unknown' for all predictions")
                    gene_preds = ["Unknown"] * len(sequences)
                else:
                    probs = model.predict_proba(features)
                    max_probs = probs.max(axis=1)
                    pred_indices = probs.argmax(axis=1)
                    
                    for idx, prob in zip(pred_indices, max_probs):
                        if prob >= CONFIDENCE_THRESHOLD:
                            gene = encoder.inverse_transform([idx])[0] if encoder else str(idx)
                        else:
                            gene = "Unknown"
                        gene_preds.append(gene)
                    logger.info(f"Gene family predictions: {len([p for p in gene_preds if p != 'Unknown'])} confident")
            except Exception as e:
                logger.exception("Gene family prediction failed")
                gene_preds = ["Unknown"] * len(sequences)

            # --- Save temp fasta with BLAST-safe headers ---
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
                temp_fasta = tmp.name
                header_map = {}
                for h, s in zip(headers, sequences):
                    safe_h = h.replace(" ", "_")
                    header_map[safe_h] = h
                    tmp.write(f">{safe_h}\n{s}\n")

            with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp:
                blast_output = tmp.name

            # --- Run BLAST ---
            blast_df = run_blast_and_get_dataframe(temp_fasta, blast_output)

            # --- Load mapping file ---
            map_df = pd.DataFrame()
            if os.path.exists(PATHWAY_MAP):
                try:
                    map_df = pd.read_csv(PATHWAY_MAP, dtype=str)
                    logger.info(f"Mapping file loaded with {len(map_df)} entries")
                except Exception as e:
                    logger.error(f"Failed to read mapping file: {e}")
            else:
                logger.warning(f"Mapping file not found at {PATHWAY_MAP}")

            # --- Pick best hits ---
            combined_results = []
            if not blast_df.empty:
                top_hits_idx = blast_df.groupby("qseqid")["bitscore"].idxmax()
                top_hits = blast_df.loc[top_hits_idx]
            else:
                top_hits = pd.DataFrame(columns=blast_df.columns)

            # --- Build results ---
            for i, orig_header in enumerate(headers):
                safe_h = orig_header.replace(" ", "_")
                top_hit_row = top_hits[top_hits["qseqid"] == safe_h]

                if not top_hit_row.empty:
                    row = top_hit_row.iloc[0]
                    top_hit = row["sseqid"]
                    bitscore = float(row["bitscore"])
                    pathways = "No Pathways Found"
                    
                    if not map_df.empty and "Entry" in map_df.columns and "Pathways" in map_df.columns:
                        paths = map_df.loc[map_df["Entry"] == top_hit, "Pathways"]
                        if (not paths.empty) and bitscore >= BITSCORE_THRESHOLD:
                            pathways = paths.values[0]
                else:
                    top_hit = "No Hit"
                    bitscore = 0.0
                    pathways = "No Pathways Found"

                gene_family = gene_preds[i] if i < len(gene_preds) else "Unknown"
                combined_results.append([orig_header, top_hit, bitscore, pathways, gene_family])

            # --- Save results ---
            result_df = pd.DataFrame(
                combined_results,
                columns=["Query/Header", "Top Hit", "Bitscore", "Pathways", "Predicted Gene Family"]
            )
            result_df.to_csv(RESULT_CSV, index=False)
            latest_predictions = combined_results

            # Store results in session and redirect to results page
            session['predictions'] = combined_results
            session['result_summary'] = {
                'total_sequences': len(sequences),
                'successful_predictions': len([p for p in gene_preds if p != "Unknown"]),
                'total_pathways_found': len([r for r in combined_results if r[3] != "No Pathways Found"])
            }

            logger.info(f"Processing complete. Results saved to {RESULT_CSV}")
            
            # Clean up temp files
            try:
                os.unlink(temp_fasta)
                os.unlink(blast_output)
            except:
                pass
            
            return redirect(url_for('results'))

        except Exception as e:
            logger.exception("Error during sequence processing")
            return f"⚠️ An error occurred: {str(e)}", 500

    return render_template("index.html", predictions=[])

@app.route("/download_csv")
def download_csv():
    global latest_predictions
    if not latest_predictions:
        return "⚠️ No predictions to download.", 400
    
    try:
        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerow(["Query/Header", "Top Hit", "Bitscore", "Pathways", "Predicted Gene Family"])
        writer.writerows(latest_predictions)
        output.seek(0)
        logger.info("CSV file generated for download")
        return send_file(
            io.BytesIO(output.getvalue().encode()),
            mimetype="text/csv",
            as_attachment=True,
            download_name="combined_predictions.csv"
        )
    except Exception as e:
        logger.exception("Error generating CSV")
        return f"⚠️ Error generating file: {str(e)}", 500

@app.route("/results")
def results():
    predictions = session.get('predictions', [])
    summary = session.get('result_summary', {})

    if not predictions:
        return redirect(url_for('index'))

    return render_template("results.html", predictions=predictions, summary=summary)

@app.route("/health")
def health():
    """Health check endpoint for Hugging Face Spaces"""
    return {
        "status": "healthy",
        "blast_available": BLAST_PATH is not None and os.path.exists(BLAST_PATH),
        "blast_db_exists": check_blast_db(),
        "model_loaded": model is not None,
        "esm_model_loaded": esm_model is not None
    }, 200

@app.route("/diagnostics")
def diagnostics():
    """Show system diagnostics (for debugging)"""
    return {
        "blast_path": BLAST_PATH,
        "blast_working": check_blast_availability(),
        "blast_db_path": BLAST_DB,
        "blast_db_exists": check_blast_db(),
        "pathway_map": PATHWAY_MAP,
        "pathway_map_exists": os.path.exists(PATHWAY_MAP),
        "model_loaded": model is not None,
        "encoder_loaded": encoder is not None,
        "esm_model_loaded": esm_model is not None,
        "base_dir": BASE_DIR,
        "current_dir": os.getcwd()
    }, 200

@app.route("/debug_blast")
def debug_blast():
    """Detailed BLAST debugging information"""
    # Check database files
    db_files = {}
    if BLAST_DB:
        db_files = {
            f"{BLAST_DB}.phr": os.path.exists(f"{BLAST_DB}.phr"),
            f"{BLAST_DB}.pin": os.path.exists(f"{BLAST_DB}.pin"),
            f"{BLAST_DB}.psq": os.path.exists(f"{BLAST_DB}.psq")
        }
    
    return {
        "blast_executable": {
            "path": BLAST_PATH,
            "exists": os.path.exists(BLAST_PATH) if BLAST_PATH else False,
            "accessible": check_blast_availability()
        },
        "blast_database": {
            "configured_path": BLAST_DB,
            "db_check_passed": check_blast_db(),
            "database_files": db_files
        },
        "models": {
            "lr_model": model is not None,
            "encoder": encoder is not None,
            "esm_model": esm_model is not None
        },
        "files": {
            "model_file": os.path.exists(MODEL_PATH),
            "encoder_file": os.path.exists(ENCODER_PATH),
            "pathway_map_file": os.path.exists(PATHWAY_MAP)
        }
    }, 200

if __name__ == "__main__":
    host = os.environ.get("HOST", "0.0.0.0")
    port = int(os.environ.get("PORT", 10000))
    debug = os.environ.get("DEBUG", "False").lower() == "true"
    
    logger.info(f"Starting Flask app on {host}:{port} (debug={debug})")
    logger.info(f"Model path: {MODEL_PATH}")
    logger.info(f"BLAST path: {BLAST_PATH}")
    
    app.run(host=host, port=port, debug=debug)
