#!/bin/bash

#SBATCH --job-name=emu_S289_S384
#SBATCH --account=duttonc
#SBATCH --qos=duttonc
#SBATCH --partition=hpg-turin
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/blue/duttonc/duttonc/giraffe/superaccuracy/logs/emu_S289_S384_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=duttonc@ufl.edu

module load conda
source $(conda info --base)/etc/profile.d/conda.sh

# Activate/create conda environment and dependencies
if ! conda info --envs | grep -q "^emu_env[[:space:]]"; then
    conda create -y -n emu_env python=3.9
fi
conda activate emu_env
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -y emu
pip install osfclient

# Setup Emu database with OSF if needed
EMU_DB_DIR="/blue/duttonc/duttonc/databases/emu_db"
if [ ! -d "$EMU_DB_DIR" ] || [ ! -f "$EMU_DB_DIR/species_taxid.fasta" ]; then
    echo "Downloading Emu database from OSF..." >&2
    mkdir -p "$EMU_DB_DIR"
    export EMU_DATABASE_DIR="$EMU_DB_DIR"
    cd "$EMU_DB_DIR"
    osf -p 56uf7 fetch osfstorage/emu-prebuilt/emu.tar
    tar -xvf emu.tar
    rm emu.tar
    cd -
fi

INPUT_ROOT="/blue/duttonc/duttonc/giraffe/superaccuracy"
OUTPUT_ROOT="$INPUT_ROOT/emu_tax_S289_S384"
mkdir -p "$OUTPUT_ROOT"

# Loop ONLY over S289–S384, checking and skipping processed samples.
for i in $(seq 289 384); do
    SAMPLE="S${i}"
    SAMPLE_DIR="$INPUT_ROOT/$SAMPLE"
    OUTDIR="$OUTPUT_ROOT/$SAMPLE"
    MERGED_FASTQ="$INPUT_ROOT/${SAMPLE}_merged_for_emu.fastq.gz"
    OUTFILE="$OUTDIR/${SAMPLE}_merged_for_emu.fastq_rel-abundance.tsv"

    echo "Checking $SAMPLE, OUTFILE=$OUTFILE..."
    if [ -f "$OUTFILE" ]; then
        echo "Sample $SAMPLE already processed at $OUTFILE. Skipping."
        continue
    fi

    # Only process if sample folder exists
    if [ ! -d "$SAMPLE_DIR" ]; then
        echo "Sample folder $SAMPLE_DIR not found, skipping."
        continue
    fi

    mkdir -p "$OUTDIR"

    # Merge all FASTQ(.gz) files in current sample folder
    if ls "$SAMPLE_DIR"/*.fastq.gz 1> /dev/null 2>&1; then
        cat "$SAMPLE_DIR"/*.fastq.gz > "$MERGED_FASTQ"
    elif ls "$SAMPLE_DIR"/*.fastq 1> /dev/null 2>&1; then
        cat "$SAMPLE_DIR"/*.fastq | gzip > "$MERGED_FASTQ"
    else
        echo "No FASTQ files found for $SAMPLE, skipping."
        continue
    fi

    emu abundance "$MERGED_FASTQ" --db "$EMU_DB_DIR" --threads $SLURM_CPUS_PER_TASK --output-dir "$OUTDIR" --keep-counts
done

# Combine all outputs (TSV, with counts)
emu combine-outputs "$OUTPUT_ROOT" species --counts

# Also generate CSV for spreadsheet use
TSV_COMBINED="$OUTPUT_ROOT/emu-combined-species.tsv"
CSV_COMBINED="$OUTPUT_ROOT/emu-combined-species.csv"
if [ -f "$TSV_COMBINED" ]; then
    awk 'BEGIN{FS=OFS="\t"}{print}' "$TSV_COMBINED" | sed 's/\t/,/g' > "$CSV_COMBINED"
    echo "Combined CSV file: $CSV_COMBINED"
fi

echo "SUCCESS: S289–S384 processed. Combined table (counts) in $CSV_COMBINED and $TSV_COMBINED"
