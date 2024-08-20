#!/bin/bash
set -euo pipefail
# Parameters
DMND_COV_STATS=$1
LENGTHS_FILE=$2
DIAMOND_HITS_PER_READ_FILE=$3
OUTPUT_DIR=$4
BINS=$5

###############################################################################
# Run full analysis
###############################################################################
zstdcat $DIAMOND_HITS_PER_READ_FILE | $DMND_COV_STATS \
        $LENGTHS_FILE \
        - \
        $OUTPUT_DIR/hitstats_bestperread \
        $BINS
