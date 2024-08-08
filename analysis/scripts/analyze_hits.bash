#!/bin/bash

# Parameters
DMND_COV_STATS=$1
LENGTHS_FILE=$2
DIAMOND_HIT_FILE=$3
BEST_HITS_999_FILE=$4
BEST_HITS_90_FILE=$5
OUTPUT_DIR=$6
BINS=$7

###############################################################################
# Run full analysis
###############################################################################
zstdcat $DIAMOND_HIT_FILE | $DMND_COV_STATS \
        $LENGTHS_FILE \
        - \
        $OUTPUT_DIR/hitstats_all \
        $BINS

###############################################################################
# Run 99.9% analysis
###############################################################################
zstdcat $BEST_HITS_999_FILE | $DMND_COV_STATS \
        $LENGTHS_FILE \
        - \
        $OUTPUT_DIR/hitstats_999 \
        $BINS

###############################################################################
# Run 90% analysis
###############################################################################
zstdcat $BEST_HITS_90_FILE | $DMND_COV_STATS \
        $LENGTHS_FILE \
        - \
        $OUTPUT_DIR/hitstats_90 \
        $BINS
