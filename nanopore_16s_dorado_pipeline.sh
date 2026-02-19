#!/usr/bin/env bash
# ==============================================================================
# nanopore_16s_dorado_pipeline.sh
# Production Oxford Nanopore 16S pipeline: dorado basecall + custom demux
# Compatible with HiPerGator (SLURM/GPU). Designed for Quick-16S Full-Length
# + SQK-LSK114 libraries; custom Zymo barcodes via minibar_and_barcodes.zip.
#
# Usage:
#   bash nanopore_16s_dorado_pipeline.sh \
#       --pod5_dir /path/to/pod5 \
#       --out_dir  /path/to/output \
#       --model    dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
#      [--threads 16] [--gpus 2] \
#      [--keep_unclassified true] \
#      [--demux_mode after_basecall|during_basecall] \
#      [--trim_primers false] \
#      [--force]
# ==============================================================================
set -euo pipefail

# ── Colour helpers ─────────────────────────────────────────────────────────────
RED='\033[0;31m'; YELLOW='\033[1;33m'; GREEN='\033[0;32m'
CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'
info()  { echo -e "${CYAN}[INFO]${NC}  $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $*"; }
error() { echo -e "${RED}[ERROR]${NC} $*" >&2; exit 1; }
ok()    { echo -e "${GREEN}[OK]${NC}    $*"; }

# ── Defaults ───────────────────────────────────────────────────────────────────
THREADS=16
GPUS=2
KEEP_UNCLASSIFIED=true
DEMUX_MODE="after_basecall"
TRIM_PRIMERS=false
FORCE=false
POD5_DIR=""
OUT_DIR=""
MODEL=""

BARCODE_ZIP_URL="https://zymo-microbiomics-service.s3.amazonaws.com/epiquest/epiquest_in4521/VUCKWUBPTZJFQRXS/rawdata/240903/minibar_and_barcodes.zip"

# ── Argument parsing ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --pod5_dir)            POD5_DIR="$2";            shift 2 ;;
        --out_dir)             OUT_DIR="$2";             shift 2 ;;
        --model)               MODEL="$2";               shift 2 ;;
        --threads)             THREADS="$2";             shift 2 ;;
        --gpus)                GPUS="$2";                shift 2 ;;
        --keep_unclassified)   KEEP_UNCLASSIFIED="$2";   shift 2 ;;
        --demux_mode)          DEMUX_MODE="$2";          shift 2 ;;
        --trim_primers)        TRIM_PRIMERS="$2";        shift 2 ;;
        --force)               FORCE=true;               shift   ;;
        *) error "Unknown argument: $1. Run with --help to see options." ;;
    esac
done

# ── Required arg validation ────────────────────────────────────────────────────
[[ -z "$POD5_DIR" ]] && error "--pod5_dir is required."
[[ -z "$OUT_DIR"  ]] && error "--out_dir is required."
[[ -z "$MODEL"    ]] && error "--model is required (e.g. dna_r10.4.1_e8.2_400bps_sup@v5.0.0 or full path)."

[[ "$DEMUX_MODE" =~ ^(after_basecall|during_basecall)$ ]] \
    || error "--demux_mode must be 'after_basecall' or 'during_basecall'. Got: $DEMUX_MODE"

# ── Validate POD5 input ────────────────────────────────────────────────────────
[[ -d "$POD5_DIR" ]] || error "POD5 directory not found: $POD5_DIR"
POD5_COUNT=$(find "$POD5_DIR" -maxdepth 2 -name "*.pod5" | wc -l)
[[ "$POD5_COUNT" -gt 0 ]] || error "No .pod5 files found under $POD5_DIR"
info "Found $POD5_COUNT POD5 file(s) in $POD5_DIR"

# ── Output layout ──────────────────────────────────────────────────────────────
LOGS_DIR="${OUT_DIR}/logs"
RES_DIR="${OUT_DIR}/resources"
BC_DIR="${OUT_DIR}/basecalled"
DEMUX_DIR="${OUT_DIR}/demux"
REPORTS_DIR="${OUT_DIR}/reports"
TRIMMED_DIR="${OUT_DIR}/demux_trimmed"

# Guard against overwriting existing results
if [[ -d "$OUT_DIR" && "$FORCE" == "false" ]]; then
    warn "Output directory already exists: $OUT_DIR"
    warn "Pass --force to overwrite. Exiting to protect existing results."
    exit 1
fi

mkdir -p "$LOGS_DIR" "$RES_DIR" "$BC_DIR" "$DEMUX_DIR" "$REPORTS_DIR"
mkdir -p "${REPORTS_DIR}/dorado_summary_files"

# ── Dependency checks ──────────────────────────────────────────────────────────
check_cmd() {
    command -v "$1" &>/dev/null \
        || error "Dependency not found: $1. Please load the appropriate module or add to PATH."
}

check_cmd dorado
check_cmd unzip
check_cmd awk
check_cmd find

# Prefer pigz over gzip for compression speed
if command -v pigz &>/dev/null; then
    GZIP_CMD="pigz"
    GZIP_THREADS="-p ${THREADS}"
    info "Using pigz for compression."
else
    GZIP_CMD="gzip"
    GZIP_THREADS=""
    warn "pigz not found; falling back to gzip (slower)."
fi

SAMTOOLS_AVAILABLE=false
command -v samtools &>/dev/null && SAMTOOLS_AVAILABLE=true && info "samtools found."

CUTADAPT_AVAILABLE=false
if [[ "$TRIM_PRIMERS" == "true" ]]; then
    command -v cutadapt &>/dev/null \
        || error "--trim_primers true requested but cutadapt is not in PATH."
    CUTADAPT_AVAILABLE=true
fi

# ── Lightweight dorado flag check ──────────────────────────────────────────────
info "Verifying dorado demux flags..."
DEMUX_HELP=$(dorado demux --help 2>&1 || true)

check_dorado_flag() {
    local flag="$1"
    echo "$DEMUX_HELP" | grep -q -- "$flag" \
        || error "dorado demux flag '$flag' not found in --help output. \
Your dorado version may differ. Check: dorado demux --help"
}

check_dorado_flag "--output-dir"
check_dorado_flag "--barcode-arrangement"
# --barcode-sequences name varies; check for both common spellings
if ! echo "$DEMUX_HELP" | grep -qE -- "--barcode-sequences|--barcode-seqs"; then
    error "Cannot find --barcode-sequences or --barcode-seqs in dorado demux --help. \
Please inspect: dorado demux --help and adjust the BARCODE_SEQ_FLAG in this script."
fi
# Set the correct flag name
if echo "$DEMUX_HELP" | grep -q -- "--barcode-sequences"; then
    BARCODE_SEQ_FLAG="--barcode-sequences"
else
    BARCODE_SEQ_FLAG="--barcode-seqs"
fi
info "Using barcode sequences flag: $BARCODE_SEQ_FLAG"

# ── Download & unzip barcode bundle ───────────────────────────────────────────
ZIP_PATH="${RES_DIR}/minibar_and_barcodes.zip"
EXTRACT_DIR="${RES_DIR}/barcode_bundle"

if [[ ! -f "$ZIP_PATH" ]]; then
    info "Downloading barcode bundle..."
    if command -v wget &>/dev/null; then
        wget -q -O "$ZIP_PATH" "$BARCODE_ZIP_URL"
    elif command -v curl &>/dev/null; then
        curl -sSL -o "$ZIP_PATH" "$BARCODE_ZIP_URL"
    else
        error "Neither wget nor curl found. Cannot download barcode bundle."
    fi
    ok "Download complete: $ZIP_PATH"
else
    info "Barcode bundle already present; skipping download."
fi

mkdir -p "$EXTRACT_DIR"
if [[ ! -f "${EXTRACT_DIR}/.extracted" ]]; then
    info "Unzipping barcode bundle..."
    unzip -o -q "$ZIP_PATH" -d "$EXTRACT_DIR"
    touch "${EXTRACT_DIR}/.extracted"
fi

# ── Auto-detect TOML and FASTA ────────────────────────────────────────────────
ARR_TOML=$(find "$EXTRACT_DIR" -type f -name "*.toml" | head -n 1)
BAR_FASTA=$(find "$EXTRACT_DIR" -type f \( -name "*.fa" -o -name "*.fasta" \) | head -n 1)

[[ -n "$ARR_TOML"  ]] || error "No .toml file found under $EXTRACT_DIR. \
Inspect the zip contents: unzip -l $ZIP_PATH"
[[ -n "$BAR_FASTA" ]] || error "No .fa or .fasta file found under $EXTRACT_DIR. \
Inspect the zip contents: unzip -l $ZIP_PATH"

ok "Detected barcode arrangement TOML : $ARR_TOML"
ok "Detected barcode sequences FASTA  : $BAR_FASTA"

# ── Provenance / versions ─────────────────────────────────────────────────────
VERSIONS_FILE="${REPORTS_DIR}/versions.txt"
{
    echo "===== Pipeline Provenance ====="
    echo "Date           : $(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo "Hostname       : $(hostname -f)"
    echo "SLURM_JOB_ID   : ${SLURM_JOB_ID:-not_set}"
    echo "Model          : $MODEL"
    echo "TOML file      : $ARR_TOML"
    echo "FASTA file     : $BAR_FASTA"
    echo "demux_mode     : $DEMUX_MODE"
    echo "trim_primers   : $TRIM_PRIMERS"
    echo "keep_unclassified : $KEEP_UNCLASSIFIED"
    echo ""
    echo "===== Tool Versions ====="
    echo -n "dorado         : "; dorado --version 2>&1 || true
    echo -n "pigz/gzip      : "; $GZIP_CMD --version 2>&1 | head -1 || true
    $SAMTOOLS_AVAILABLE && { echo -n "samtools       : "; samtools --version | head -1; } || true
    $CUTADAPT_AVAILABLE && { echo -n "cutadapt       : "; cutadapt --version; } || true
    echo "unzip          : $(unzip -v 2>&1 | head -1)"
    echo "bash           : $BASH_VERSION"
} > "$VERSIONS_FILE"
info "Versions written to $VERSIONS_FILE"

# ── Helper: merge BAM sub-dirs into FASTQ.gz per barcode ─────────────────────
# dorado demux writes one BAM per barcode into the output dir.
# We convert each to FASTQ.gz.
bam_to_fastq_gz() {
    local bam_dir="$1"
    local fq_dir="$2"
    mkdir -p "$fq_dir"

    find "$bam_dir" -maxdepth 1 -name "*.bam" | while read -r bam; do
        local name
        name=$(basename "$bam" .bam)
        local out_fq="${fq_dir}/${name}.fastq.gz"
        if [[ "$SAMTOOLS_AVAILABLE" == "true" ]]; then
            samtools fastq -@ "$THREADS" "$bam" | $GZIP_CMD $GZIP_THREADS -c > "$out_fq"
        else
            dorado summary "$bam" 2>/dev/null | tail -n +2 | awk '{print "@"$1"\n"$10"\n+\n"$11}' \
                | $GZIP_CMD $GZIP_THREADS -c > "$out_fq" \
                || warn "samtools not available and dorado summary fallback failed for $bam. Install samtools."
        fi
        ok "Converted: $bam -> $out_fq"
    done
}

# ── Helper: check if temp/scratch space is available ─────────────────────────
if [[ -n "${TMPDIR:-}" && -d "${TMPDIR:-}" ]]; then
    SCRATCH="${TMPDIR}/dorado_scratch_$$"
    info "Using local scratch at: $SCRATCH"
else
    SCRATCH="${OUT_DIR}/.scratch_$$"
    warn "TMPDIR not set; using $SCRATCH (may be slower on shared FS)."
fi
mkdir -p "$SCRATCH"
trap 'rm -rf "$SCRATCH"' EXIT

# ── BASECALLING ────────────────────────────────────────────────────────────────
MERGED_BAM="${BC_DIR}/basecalled_merged.bam"
MERGED_BAM_SCRATCH="${SCRATCH}/basecalled_merged.bam"

if [[ "$DEMUX_MODE" == "after_basecall" ]]; then
    # ─── Mode A: Basecall, then demux separately ─────────────────────────────
    info "=== STEP 1: Basecalling (SUP, no trim) ==="
    dorado basecaller \
        "$MODEL" \
        "$POD5_DIR" \
        --device cuda:all \
        --no-trim \
        > "$MERGED_BAM_SCRATCH"

    info "Copying basecalled BAM from scratch to $MERGED_BAM ..."
    cp "$MERGED_BAM_SCRATCH" "$MERGED_BAM"
    ok "Basecalling complete. BAM: $MERGED_BAM"

    # Optional: dorado summary
    if dorado summary --help &>/dev/null 2>&1; then
        dorado summary "$MERGED_BAM" \
            > "${REPORTS_DIR}/dorado_summary_files/sequencing_summary.txt" 2>/dev/null \
            || warn "dorado summary failed (non-fatal)."
    fi

    # ─── Mode A: Demux after basecall ─────────────────────────────────────────
    info "=== STEP 2: Demultiplexing (after basecall) ==="
    DEMUX_BAM_SCRATCH="${SCRATCH}/demux_bams"
    mkdir -p "$DEMUX_BAM_SCRATCH"

    DEMUX_UNCLASS_FLAG=""
    [[ "$KEEP_UNCLASSIFIED" == "true" ]] && DEMUX_UNCLASS_FLAG="--emit-fastq" || true
    # Note: dorado demux natively emits per-barcode BAMs; we convert below.
    # --emit-fastq is not universally supported; we use BAM -> fastq conversion.

    dorado demux \
        "$MERGED_BAM" \
        --output-dir "$DEMUX_BAM_SCRATCH" \
        --barcode-arrangement "$ARR_TOML" \
        ${BARCODE_SEQ_FLAG} "$BAR_FASTA" \
        --no-trim

    ok "Demux BAMs written to $DEMUX_BAM_SCRATCH"

    # Convert BAMs to per-barcode FASTQ.gz in DEMUX_DIR
    info "=== STEP 3: Converting demux BAMs to FASTQ.gz ==="
    find "$DEMUX_BAM_SCRATCH" -maxdepth 1 -name "*.bam" | while read -r bam; do
        local_name=$(basename "$bam" .bam)
        # Skip unclassified if not requested
        if [[ "$local_name" == *"unclassified"* && "$KEEP_UNCLASSIFIED" != "true" ]]; then
            info "Skipping unclassified (keep_unclassified=false)."
            continue
        fi
        bc_out="${DEMUX_DIR}/${local_name}"
        mkdir -p "$bc_out"
        out_fq="${bc_out}/${local_name}.fastq.gz"
        if [[ "$SAMTOOLS_AVAILABLE" == "true" ]]; then
            samtools fastq -@ "$THREADS" "$bam" | $GZIP_CMD $GZIP_THREADS -c > "$out_fq"
        else
            error "samtools is required to convert demux BAMs to FASTQ. Please load samtools module."
        fi
        ok "  -> $out_fq"
    done

else
    # ─── Mode B: Basecall + demux during basecall ─────────────────────────────
    info "=== STEP 1+2: Basecalling WITH inline demultiplexing ==="
    DEMUX_DURING_SCRATCH="${SCRATCH}/demux_during"
    mkdir -p "$DEMUX_DURING_SCRATCH"

    UNCLASS_FLAG=""
    [[ "$KEEP_UNCLASSIFIED" == "true" ]] && UNCLASS_FLAG="--emit-fastq" || true

    dorado basecaller \
        "$MODEL" \
        "$POD5_DIR" \
        --device cuda:all \
        --no-trim \
        --barcode-arrangement "$ARR_TOML" \
        ${BARCODE_SEQ_FLAG} "$BAR_FASTA" \
        --output-dir "$DEMUX_DURING_SCRATCH"

    ok "Basecall+demux complete. Output: $DEMUX_DURING_SCRATCH"

    info "=== STEP 3: Converting to FASTQ.gz and organising ==="
    find "$DEMUX_DURING_SCRATCH" -maxdepth 1 -name "*.bam" | while read -r bam; do
        local_name=$(basename "$bam" .bam)
        if [[ "$local_name" == *"unclassified"* && "$KEEP_UNCLASSIFIED" != "true" ]]; then
            continue
        fi
        bc_out="${DEMUX_DIR}/${local_name}"
        mkdir -p "$bc_out"
        out_fq="${bc_out}/${local_name}.fastq.gz"
        if [[ "$SAMTOOLS_AVAILABLE" == "true" ]]; then
            samtools fastq -@ "$THREADS" "$bam" | $GZIP_CMD $GZIP_THREADS -c > "$out_fq"
        else
            error "samtools required for BAM->FASTQ conversion."
        fi
        ok "  -> $out_fq"
    done
fi

# ── PRIMER TRIMMING (optional, after demux) ────────────────────────────────────
if [[ "$TRIM_PRIMERS" == "true" ]]; then
    info "=== STEP 4: Trimming primers with cutadapt ==="
    mkdir -p "$TRIMMED_DIR"

    # ── Primer sequences: adjust these to match your library ──────────────────
    # Quick-16S Full-Length typically uses 27F / 1492R
    # Replace with your actual primer sequences if different.
    FWD_PRIMER="AGAGTTTGATCMTGGCTCAG"   # 27F
    REV_PRIMER="TACGGYTACCTTGTTACGACTT"  # 1492R

    warn "Using default 27F/1492R primer sequences. Verify these match your kit!"

    find "$DEMUX_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r bc_dir; do
        bc_name=$(basename "$bc_dir")
        trimmed_bc_dir="${TRIMMED_DIR}/${bc_name}"
        mkdir -p "$trimmed_bc_dir"

        find "$bc_dir" -name "*.fastq.gz" | while read -r fq; do
            fq_name=$(basename "$fq")
            out_fq="${trimmed_bc_dir}/${fq_name}"
            cutadapt \
                -g "$FWD_PRIMER" \
                -a "$REV_PRIMER" \
                --discard-untrimmed \
                --minimum-length 100 \
                -j "$THREADS" \
                -o "$out_fq" \
                "$fq" \
                >> "${LOGS_DIR}/cutadapt_${bc_name}.log" 2>&1
            ok "Trimmed: $fq -> $out_fq"
        done
    done
    ok "Primer trimming complete. Trimmed reads in: $TRIMMED_DIR"
fi

# ── REPORTING ─────────────────────────────────────────────────────────────────
info "=== STEP 5: Generating read count report ==="
COUNTS_TSV="${REPORTS_DIR}/demux_read_counts.tsv"
printf "barcode_folder\treads\tbases\n" > "$COUNTS_TSV"

TOTAL_READS=0
declare -A BC_READS

find "$DEMUX_DIR" -mindepth 1 -maxdepth 1 -type d | sort | while read -r bc_dir; do
    bc_name=$(basename "$bc_dir")
    reads=0
    bases=0

    # Count reads and bases from all FASTQ.gz in this barcode folder
    while IFS= read -r fq; do
        fq_reads=$(zcat "$fq" | awk 'NR%4==1 {r++} END {print r+0}')
        fq_bases=$(zcat "$fq" | awk 'NR%4==2 {bases+=length($0)} END {print bases+0}')
        reads=$(( reads + fq_reads ))
        bases=$(( bases + fq_bases ))
    done < <(find "$bc_dir" -name "*.fastq.gz")

    printf "%s\t%d\t%d\n" "$bc_name" "$reads" "$bases" >> "$COUNTS_TSV"
done

ok "Read counts written to $COUNTS_TSV"

# ── Terminal summary ───────────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}════════════════════════════════════════${NC}"
echo -e "${BOLD}  PIPELINE SUMMARY${NC}"
echo -e "${BOLD}════════════════════════════════════════${NC}"

TOTAL_READS=$(awk 'NR>1 {sum+=$2} END {print sum+0}' "$COUNTS_TSV")
UNCLASS_READS=$(awk -F'\t' '$1~/unclassified/ {print $2}' "$COUNTS_TSV" | head -1)
UNCLASS_READS=${UNCLASS_READS:-0}

echo -e "Total reads           : ${BOLD}${TOTAL_READS}${NC}"
if [[ "$TOTAL_READS" -gt 0 ]]; then
    UNCLASS_PCT=$(awk -v u="$UNCLASS_READS" -v t="$TOTAL_READS" 'BEGIN {printf "%.1f", (u/t)*100}')
    echo -e "Unclassified reads    : ${UNCLASS_READS} (${UNCLASS_PCT}%)"
fi

echo ""
echo -e "${BOLD}Top 5 barcodes by read count:${NC}"
awk -F'\t' 'NR>1 && $1!~/unclassified/ {print $2"\t"$1}' "$COUNTS_TSV" \
    | sort -rn | head -5 \
    | awk -F'\t' '{printf "  %-30s %s reads\n", $2, $1}'

echo ""
echo -e "${GREEN}Pipeline complete.${NC} Outputs in: $OUT_DIR"
echo -e "${BOLD}════════════════════════════════════════${NC}"
echo -e "${GREEN}Pipeline complete.${NC} Outputs in: $OUT_DIR"
echo -e "${BOLD}════════════════════════════════════════${NC}"
