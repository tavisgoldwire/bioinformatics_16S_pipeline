#!/bin/bash

#SBATCH --job-name=dorado_basecall
#SBATCH --account=duttonc
#SBATCH --qos=duttonc
#SBATCH --partition=hpg-turin
#SBATCH --gpus=2
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=dorado_basecall_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=duttonc@ufl.edu
#SBATCH --array=0-68

cd /blue/duttonc/duttonc/giraffe

# Ensure output directory exists
OUTPUT_DIR="superaccuracy"
mkdir -p "$OUTPUT_DIR"

# Gather .pod5 files in /blue/duttonc/duttonc/giraffe/
pod5_files=( *.pod5 )

POD5_FILE="${pod5_files[$SLURM_ARRAY_TASK_ID]}"
OUT_NAME="$OUTPUT_DIR/$(basename "$POD5_FILE" .pod5)_sup.fastq"

module load cuda/12.9.1
module load dorado

dorado basecaller --device cuda:all /blue/duttonc/shared_resources/dorado/dorado_basecalling_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 "$POD5_FILE" --no-trim --emit-fastq > "$OUT_NAME"
