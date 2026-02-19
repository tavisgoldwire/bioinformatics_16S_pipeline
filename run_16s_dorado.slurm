#!/usr/bin/env bash
# ==============================================================================
# run_16s_dorado.slurm
# SLURM submission wrapper for nanopore_16s_dorado_pipeline.sh
# HiPerGator / hpg-turin partition, account duttonc
#
# Submit with:
#   sbatch run_16s_dorado.slurm
# ==============================================================================

# ── SLURM directives ───────────────────────────────────────────────────────────
#SBATCH --job-name=dorado_16s_pipeline
#SBATCH --account=duttonc
#SBATCH --qos=duttonc
#SBATCH --partition=hpg-turin
#SBATCH --gpus=2
#SBATCH --cpus-per-task=24
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=duttonc@ufl.edu
# Log path is set dynamically below; SLURM requires a literal path at submit time,
# so we use a timestamped name using %j (job ID).
#SBATCH --output=/blue/duttonc/duttonc/logs/dorado_16s_%j.log
#SBATCH --error=/blue/duttonc/duttonc/logs/dorado_16s_%j.err

# ==============================================================================
# USER-CONFIGURABLE VARIABLES — edit these before submitting
# ==============================================================================

# Path to directory containing .pod5 files
POD5_DIR="/blue/duttonc/duttonc/giraffe"

# Desired output directory (will be created if absent; use --force to overwrite)
OUT_DIR="/blue/duttonc/duttonc/giraffe/16s_pipeline_out"

# Dorado SUP model: either a short string resolved by dorado, or a full path.
# The full-path model you already have on HiPerGator:
MODEL="/blue/duttonc/shared_resources/dorado/dorado_basecalling_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
# If you want dorado to auto-resolve the model string instead, comment out above
# and uncomment below (requires internet access or a local model cache):
# MODEL="dna_r10.4.1_e8.2_400bps_sup@v5.0.0"

# Number of CPU threads to expose to dorado / samtools / pigz
THREADS="${SLURM_CPUS_PER_TASK}"   # inherits from #SBATCH --cpus-per-task

# Number of GPUs (should match --gpus above)
GPUS=2

# Keep unclassified reads? (true/false)
KEEP_UNCLASSIFIED="true"

# Demux mode: "after_basecall" (safer, recommended) or "during_basecall"
DEMUX_MODE="after_basecall"

# Trim primers with cutadapt after demux? (true/false)
# If true, cutadapt must be available in PATH.
TRIM_PRIMERS="false"

# Pass --force to overwrite an existing OUT_DIR (remove or comment out to protect results)
FORCE_FLAG=""        # set to "--force" to allow overwrite

# Path to the pipeline script itself
PIPELINE_SCRIPT="$(dirname "$0")/nanopore_16s_dorado_pipeline.sh"

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================
set -euo pipefail

# Ensure log directory exists (SLURM won't create it automatically)
mkdir -p "$(dirname "/blue/duttonc/duttonc/logs/placeholder")"

echo "============================================================"
echo " SLURM Job ID  : ${SLURM_JOB_ID}"
echo " Node          : $(hostname)"
echo " Start time    : $(date '+%Y-%m-%d %H:%M:%S')"
echo " POD5 dir      : ${POD5_DIR}"
echo " Output dir    : ${OUT_DIR}"
echo " Model         : ${MODEL}"
echo "============================================================"

# ── Load modules ───────────────────────────────────────────────────────────────
# Option A (preferred): load via HiPerGator module system
module load cuda/12.9.1
module load dorado

# Option B: if dorado is not a module, set DORADO_BIN to the executable path
# and export it; the pipeline script will use whichever is in PATH.
# export DORADO_BIN="/blue/duttonc/shared_resources/dorado/bin/dorado"
# export PATH="${DORADO_BIN%/*}:$PATH"

# Samtools is needed for BAM->FASTQ conversion; load it too.
module load samtools

# pigz for fast parallel compression (optional but recommended)
module load pigz 2>/dev/null || echo "[WARN] pigz module not available; will use gzip."

# ── Validate the pipeline script ──────────────────────────────────────────────
[[ -f "$PIPELINE_SCRIPT" ]] \
    || { echo "[ERROR] Pipeline script not found: $PIPELINE_SCRIPT" >&2; exit 1; }
chmod +x "$PIPELINE_SCRIPT"

# ── Redirect all subsequent output into a per-job log inside OUT_DIR ──────────
# (SLURM already captures stdout/stderr above, but this gives a clean copy.)
mkdir -p "${OUT_DIR}/logs"
exec > >(tee -a "${OUT_DIR}/logs/pipeline_${SLURM_JOB_ID}.log") 2>&1

# ── Print GPU visibility ───────────────────────────────────────────────────────
echo "[INFO] Checking GPU visibility..."
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader 2>/dev/null \
    || echo "[WARN] nvidia-smi failed; GPUs may still be usable via CUDA."

# ── Launch the pipeline ────────────────────────────────────────────────────────
echo "[INFO] Launching pipeline at $(date '+%Y-%m-%d %H:%M:%S')"

bash "$PIPELINE_SCRIPT" \
    --pod5_dir        "$POD5_DIR"           \
    --out_dir         "$OUT_DIR"            \
    --model           "$MODEL"              \
    --threads         "$THREADS"            \
    --gpus            "$GPUS"              \
    --keep_unclassified "$KEEP_UNCLASSIFIED" \
    --demux_mode      "$DEMUX_MODE"         \
    --trim_primers    "$TRIM_PRIMERS"       \
    ${FORCE_FLAG}

echo "[INFO] Pipeline finished at $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================"
echo " Job ${SLURM_JOB_ID} complete."
echo "============================================================"
