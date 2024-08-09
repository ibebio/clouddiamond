#!/bin/bash
# Parameters
CONDA_ENV=$1
ANALYZE_ACCESSION=$2
DIAMONDN_FILE=$3
OUTPUT_DIR=$4
ACCESSION=$5
ACCESSION_XML=$6


# "python scripts/analyze_accession.py {input.diamondn_file} {params.output_dir} {params.accession} {input.accession_xml}"

# Activate conda environment
source activate $CONDA_ENV

set -euo pipefail
# Run tool
$ANALYZE_ACCESSION $DIAMONDN_FILE $OUTPUT_DIR $ACCESSION $ACCESSION_XML

