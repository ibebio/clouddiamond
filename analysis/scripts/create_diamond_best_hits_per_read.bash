#!/bin/bash
# Parameters
CONDA_ENV=$1
CREATE_DIAMOND_BEST_HITS_PER_READ=$2
DIAMONDN_FILE=$3
OUTPUT_DIR=$4
ACCESSION=$5

# "python scripts/create_diamond_best_hits_per_read.py {input.diamondn_file} {params.output_dir} {params.accession}"
# Activate conda environment
source activate $CONDA_ENV

set -euo pipefail
# Run tool
$CREATE_DIAMOND_BEST_HITS_PER_READ $DIAMONDN_FILE $OUTPUT_DIR $ACCESSION

