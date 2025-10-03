import os
import io
import csv
import subprocess
import pandas as pd
import joblib
import numpy as np
import torch
import logging
from flask import Flask, request, render_template, send_file, session, redirect, url_for
from Bio.Seq import Seq
from transformers import EsmModel, EsmTokenizer

# logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[logging.FileHandler("app.log"), logging.StreamHandler()]
)

# ---------- Paths ----------
BASE_DIR = os.getcwd()  # Use current working directory for relative paths
MODEL_PATH = os.path.join(BASE_DIR, "best_model_LR.sav")
ENCODER_PATH = os.path.join(BASE_DIR, "label_encoder.pkl")

# Updated BLAST path for Linux
BLAST_PATH = os.path.join(BASE_DIR, "ncbi-blast-2.17.0+-aarch64-linux", "ncbi-blast-2.17.0+", "bin", "blastp")
BLAST_DB = os.path.join(BASE_DIR, "pathway_db")
MAPPING_FILE = os.path.join(BASE_DIR, "pathway_map.csv")
BITSCORE_THRESHOLD = 80.0
CONFIDENCE_THRESHOLD = 0.7  

RESULT_CSV = os.path.join(BASE_DIR, "combined_results.csv")

# Minimum length thresholds
MIN_PROTEIN_LENGTH = 50
MIN_DNA_LENGTH = 150

# ---------- Load Models ----------
model = None
encoder = None
esm_model = None
batch_converter = None

try:
    if os.path.exists(MODEL_PATH):
        model = joblib.load(MODEL_PATH)
        logging.info("LR model loaded successfully.")
except Exception:
    logging.exception("Failed to load LR model")

try:
    if os.path.exists(ENCODER_PATH):
        encoder = joblib.load(ENCODER_PATH)
        logging.info("Label encoder loaded.")
except Exception:
    logging.exception("Failed to load label encoder")

try:
    # Load the ESM2 model using Hugging Face Transformers
    model_name = "facebook/esm2_t6_8M_UR50D"
    esm_model = EsmModel.from_pretrained(model_name)
    tokenizer = EsmTokenizer.from_pretrained(model_name)
    
    esm_model.eval()
    esm_model = esm_model.to("cpu")  # Move model to CPU
    batch_converter = tokenizer
    logging.info("ESM2 model loaded and set to CPU.")
except Exception:
    logging.exception("Failed to load ESM model")

# ---------- DNA to Protein Helper ----------
def get_best_frame_protein(nuc_seq):
    nuc_seq = "".join(nuc_seq.split()).upper()
    trim_len = len(nuc_seq) - (len(nuc_seq) % 3)
    nuc_seq = nuc_seq[:trim_len]
    seq = Seq(nuc_seq)

    frames = [
        seq[i:].translate(to_stop=False) for i in range(3)
    ] + [
        seq.reverse_complement()[i:].translate(to_stop=False) for i in range(3)
    ]

    best_protein = ""
    for prot in frames:
        longest_orf = max(prot.split("*"), key=len)
        if len(longest_orf) > len(best_protein):
            best_protein = longest_orf
    return str(best_protein)

def detect_and_translate(seq):
    seq = seq.upper()
    # DNA heuristic
    if all(base in "ATGCN" for base in seq):
        logging.info("Detected DNA sequence → translating to protein.")
        if len(seq) < MIN_DNA_LENGTH:
            return None
        prot_seq = get_best_frame_protein(seq)
        if len(prot_seq) < MIN_PROTEIN_LENGTH:
            return None
        return prot_seq
    else:
        logging.info("Detected protein sequence → using as is.")
        if len(seq) < MIN_PROTEIN_LENGTH:
            return None
        return seq

# ---------- Functions ----------
def esm2_320_embed(sequence):
    if esm_model is None or batch_converter is None:
        raise RuntimeError("ESM model not available")
    
    # Tokenize the input sequence
    inputs = batch_converter(["seq1", sequence])
    with torch.no_grad():
        results = esm_model(inputs['input_ids'], return_contacts=False)
    
    token_representations = results["representations"][0]  # Use the correct layer
    seq_repr = token_representations[0, 1:-1].detach().cpu().numpy()
    return seq_repr.mean(axis=0)

def run_blast_and_get_dataframe(temp_fasta, blast_output):
    cmd = [
        BLAST_PATH,
        "-query", temp_fasta,
        "-db", BLAST_DB,
        "-out", blast_output,
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore"
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        logging.error("BLAST error:\n%s", proc.stderr)
    if not os.path.exists(blast_output) or os.path.getsize(blast_output) == 0:
        return pd.DataFrame(columns=["qseqid","sseqid","pident","length","evalue","bitscore"])
    cols = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]
    try:
        df = pd.read_csv(blast_output, sep="\t", names=cols)
        df["bitscore"] = df["bitscore"].astype(float)
        return df
    except Exception:
        logging.exception("Failed to read BLAST output")
        return pd.DataFrame(columns=cols)

# ---------- Flask app ----------
app = Flask(__name__)
app.secret_key = 'your-secret-key-here'  
latest_predictions = []

@app.after_request
def add_cache_control(response):
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
            return "⚠️ No valid sequences provided."

        # --- Convert DNA → Protein if needed & check length ---
        translated_sequences = []
        for seq in sequences:
            prot_seq = detect_and_translate(seq)
            translated_sequences.append(prot_seq)

        # --- Check if all sequences are too short ---
        if all(s is None for s in translated_sequences):
            if len(sequences) == 1:
                return "⚠️ Sequence is too short to predict."
            else:
                return "⚠️ Sequences are too short to predict."

        # Filter out None sequences
        filtered_sequences = [s for s in translated_sequences if s is not None]
        filtered_headers = [h for s,h in zip(translated_sequences, headers) if s is not None]

        # --- Gene family prediction ---
        try:
            features = [esm2_320_embed(seq) for seq in filtered_sequences]
            features = np.array(features)
            gene_preds = []
            if model is None:
                gene_preds = ["Unknown"] * len(filtered_sequences)
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
        except Exception:
            logging.exception("Gene family prediction failed")
            gene_preds = ["Unknown"] * len(filtered_sequences)

        # --- Save fasta for BLAST ---
        temp_fasta = os.path.join(os.getcwd(), "temp_input.fasta")
        header_map = {}
        with open(temp_fasta, "w") as f:
            for h, s in zip(filtered_headers, filtered_sequences):
                safe_h = h.replace(" ", "_")
                header_map[safe_h] = h
                f.write(f">{safe_h}\n{s}\n")

        # --- Run BLAST ---
        blast_output = os.path.join(os.getcwd(), "blast_results.txt")
        blast_df = run_blast_and_get_dataframe(temp_fasta, blast_output)

        # --- Mapping file ---
        map_df = pd.DataFrame()
        if os.path.exists(MAPPING_FILE):
            try:
                map_df = pd.read_csv(MAPPING_FILE, dtype=str)
            except Exception:
                logging.exception("Failed to read mapping file")

        # --- Best hits ---
        combined_results = []
        if not blast_df.empty:
            top_hits_idx = blast_df.groupby("qseqid")["bitscore"].idxmax()
            top_hits = blast_df.loc[top_hits_idx]
        else:
            top_hits = pd.DataFrame(columns=blast_df.columns)

        # --- Build results ---
        for i, orig_header in enumerate(filtered_headers):
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

        # Store results
        session['predictions'] = combined_results
        session['result_summary'] = {
            'total_sequences': len(filtered_sequences),
            'successful_predictions': len([p for p in gene_preds if p != "Unknown"]),
            'total_pathways_found': len([r for r in combined_results if r[3] != "No Pathways Found"])
        }

        return redirect(url_for('results'))

    return render_template("index.html", predictions=[])

@app.route("/download_csv")
def download_csv():
    global latest_predictions
    if not latest_predictions:
        return "⚠️ No predictions to download."
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(["Query/Header", "Top Hit", "Bitscore", "Pathways", "Predicted Gene Family"])
    writer.writerows(latest_predictions)
    output.seek(0)
    return send_file(
        io.BytesIO(output.getvalue().encode()),
        mimetype="text/csv",
        as_attachment=True,
        download_name="combined_predictions.csv"
    )

@app.route("/results")
def results():
    predictions = session.get('predictions', [])
    summary = session.get('result_summary', {})

    if not predictions:
        return redirect(url_for('index'))

    return render_template("results.html", predictions=predictions, summary=summary)

if __name__ == "__main__":
    app.run(debug=True)
