#!/bin/bash
# bash strict mode
# set -euo pipefail
###############################################################################
# This script is used to merge the hit tables from the different SRA accessions
# produced by the pipeline.
# It expects the following arguments:
# 1. The results directory, with subirectories for each SRA accession
# 2. The type of hit table to merge( all, 999, 90)
# 3. The output file name
# 4. Flag if additional metadata from the stats.txt file should be included.
#    If given, it is expected to be the string include_metadata.
# The script creates a merged hit table with an additional column for the 
# SRA accession.
###############################################################################

# To create all hits tables, the following code was used:

# Merge all hits tables

# cd /ebio/abt6_projects/clouddiamond/data/analysis/Ath-Rep/results

# # all
# /ebio/abt6_projects7/small_projects/ibezrukov2/cloud-diamond/analysis/scripts/merge_hits_tables.bash Ath-Rep all hitstats_Ath-Rep_all.tsv
# /ebio/abt6_projects7/small_projects/ibezrukov2/cloud-diamond/analysis/scripts/merge_hits_tables.bash Ath-Rep all hitstats_Ath-Rep_all_with_metadata.tsv include_metadata

# # 90
# /ebio/abt6_projects7/small_projects/ibezrukov2/cloud-diamond/analysis/scripts/merge_hits_tables.bash Ath-Rep 90 hitstats_Ath-Rep_90.tsv
# /ebio/abt6_projects7/small_projects/ibezrukov2/cloud-diamond/analysis/scripts/merge_hits_tables.bash Ath-Rep 90 hitstats_Ath-Rep_90_with_metadata.tsv include_metadata

# # 999
# /ebio/abt6_projects7/small_projects/ibezrukov2/cloud-diamond/analysis/scripts/merge_hits_tables.bash Ath-Rep 999 hitstats_Ath-Rep_999.tsv
# /ebio/abt6_projects7/small_projects/ibezrukov2/cloud-diamond/analysis/scripts/merge_hits_tables.bash Ath-Rep 999 hitstats_Ath-Rep_999_with_metadata.tsv include_metadata





###############################################################################
# Parse parameters
###############################################################################
# Check if the correct arguments were provided, either 3 or 4 and check 
# for the flag include_metadata
if [ "$#" -ne 3 ] && [ "$#" -ne 4 ]; then
    echo "Usage: $0 <results_dir> <hit_table_type> <output_file> [include_metadata]"
    exit 1
fi


# Assign the arguments to variables
RESULTS_DIR=$1
HIT_TABLE_TYPE=$2
OUTPUT_FILE=$3
# Check if the flag include_metadata was provided
if [ "$#" -eq 4 ]; then
    INCLUDE_METADATA=1
else
    INCLUDE_METADATA=0
fi


###############################################################################
# Merge hit tables
###############################################################################
# Create header for the output file
if [ $INCLUDE_METADATA -eq 1 ]; then
    # Add the following columns from the stats file:
    # num_reads       num_bases       num_hits        sequencer       year    lab  library_strategy        library_source
    echo -e "accession\tnum_reads\tnum_bases\tnum_hits\tsequencer\tyear\tlab\tlibrary_strategy\tlibrary_source\tnum_best_hits90\tnum_best_hits95\tnum_best_hits99\tnum_best_hits999\t$(head -n 1 $RESULTS_DIR/$(ls $RESULTS_DIR | head -n 1)/hitstats_$HIT_TABLE_TYPE.tsv)" > $OUTPUT_FILE
else
    echo -e "accession\t$(head -n 1 $RESULTS_DIR/$(ls $RESULTS_DIR | head -n 1)/hitstats_$HIT_TABLE_TYPE.tsv)" > $OUTPUT_FILE
fi

# Add the hit tables for all accessions to the output file
for ACCESSION_DIR in $RESULTS_DIR/*; do
    ACCESSION=$(basename $ACCESSION_DIR)
    if [ $INCLUDE_METADATA -eq 1 ]; then
        # Get the metadata from the stats file
        # we need the columns 2-15: num_reads, num_bases, num_hits, sequencer, year, lab, ecotype, genotype, library_strategy, library_source, num_best_hits90, num_best_hits95, num_best_hits99, num_best_hits999
        # set them as individual variables for each column
        NUM_READS="$(awk -F"\t" 'NR==2 {print $2}' $ACCESSION_DIR/${ACCESSION}_stats.txt)"
        NUM_BASES="$(awk -F"\t" 'NR==2 {print $3}' $ACCESSION_DIR/${ACCESSION}_stats.txt)"
        NUM_HITS="$(awk -F"\t" 'NR==2 {print $4}' $ACCESSION_DIR/${ACCESSION}_stats.txt)"
        SEQUENCER="$(awk -F"\t" 'NR==2 {print $5}' $ACCESSION_DIR/${ACCESSION}_stats.txt)"
        YEAR="$(awk -F"\t" 'NR==2 {print $6}' $ACCESSION_DIR/${ACCESSION}_stats.txt)"
        LAB="$(awk -F"\t" 'NR==2 {print $7}' $ACCESSION_DIR/${ACCESSION}_stats.txt)"
        LIBRARY_STRATEGY="$(awk -F"\t" 'NR==2 {print $10}' $ACCESSION_DIR/${ACCESSION}_stats.txt)"
        LIBRARY_SOURCE="$(awk -F"\t" 'NR==2 {print $11}' $ACCESSION_DIR/${ACCESSION}_stats.txt)"
        NUM_BEST_HITS90="$(awk -F"\t" 'NR==2 {print $12}' $ACCESSION_DIR/${ACCESSION}_stats.txt)"
        NUM_BEST_HITS95="$(awk -F"\t" 'NR==2 {print $13}' $ACCESSION_DIR/${ACCESSION}_stats.txt)"
        NUM_BEST_HITS99="$(awk -F"\t" 'NR==2 {print $14}' $ACCESSION_DIR/${ACCESSION}_stats.txt)"
        NUM_BEST_HITS999="$(awk -F"\t" 'NR==2 {print $15}' $ACCESSION_DIR/${ACCESSION}_stats.txt)"
        # Add the metadata to the hit table
        tail -n +2 $ACCESSION_DIR/hitstats_$HIT_TABLE_TYPE.tsv |\
        awk -v ACCESSION=$ACCESSION \
            -v NUM_READS=$NUM_READS \
            -v NUM_BASES=$NUM_BASES \
            -v NUM_HITS=$NUM_HITS \
            -v SEQUENCER="$SEQUENCER" \
            -v YEAR="$YEAR" \
            -v LAB="$LAB" \
            -v LIBRARY_STRATEGY="$LIBRARY_STRATEGY" \
            -v LIBRARY_SOURCE="$LIBRARY_SOURCE" \
            -v NUM_BEST_HITS90=$NUM_BEST_HITS90 \
            -v NUM_BEST_HITS95=$NUM_BEST_HITS95 \
            -v NUM_BEST_HITS99=$NUM_BEST_HITS99 \
            -v NUM_BEST_HITS999=$NUM_BEST_HITS999 \
             \
            '{print ACCESSION"\t"NUM_READS"\t"NUM_BASES"\t"NUM_HITS"\t"SEQUENCER"\t"YEAR"\t"LAB"\t"LIBRARY_STRATEGY"\t"LIBRARY_SOURCE"\t"NUM_BEST_HITS90"\t"NUM_BEST_HITS95"\t"NUM_BEST_HITS99"\t"NUM_BEST_HITS999"\t"$0}' >> $OUTPUT_FILE
    else
        # Add the hit table to the output file
        tail -n +2 $ACCESSION_DIR/hitstats_$HIT_TABLE_TYPE.tsv | awk -v ACCESSION=$ACCESSION '{print ACCESSION"\t"$0}' >> $OUTPUT_FILE
    fi
done


# for ACCESSION_DIR in $RESULTS_DIR/*; do
#     ACCESSION=$(basename $ACCESSION_DIR)
#     tail -n +2 $ACCESSION_DIR/hitstats_$HIT_TABLE_TYPE.tsv | \
#     awk -v ACCESSION=$ACCESSION '{print ACCESSION"\t"$0}' >> $OUTPUT_FILE
# done