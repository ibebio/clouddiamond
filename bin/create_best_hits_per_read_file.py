#!/usr/bin/env python3

# This script creates a bests hits per read file from a diamond output file. It takes the following inputs:

# 1. Diamond output file, either plain text or zstd compressed
# 2. Output folder, it will be created if it does not exist
# 3. SRA accession id, it will be used to name the output files

# Pandas is used to read the diamond output file and to find the best hits per read.
# The best hits are based on the bitscore column, for each read id, the best hit is the one with the highest bitscore.

# The input file has the following columns:
# seqid sseqid pident qstart qend sstart send evalue bitscore qlen slen cigar


import sys, os, re, argparse
import pandas as pd
import numpy as np
import zstandard as zstd

def parse_args():
    parser = argparse.ArgumentParser(description='Analyze diamond output file')
    parser.add_argument('diamond_file', help='Diamond output file')
    parser.add_argument('output_folder', help='Output folder')
    parser.add_argument('sra_accession', help='SRA accession id')
    return parser.parse_args()

def main():
    ###########################################################################
    # Initial setup
    ###########################################################################
    args = parse_args()
    # Verify the arguments
    if not os.path.exists(args.diamond_file):
        print('Diamond output file does not exist')
        sys.exit(1)
    diamond_file = args.diamond_file
    output_folder = args.output_folder
    sra_accession = args.sra_accession

    # Create output folder if it does not exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    ###########################################################################
    # Process diamond output file
    ###########################################################################

    # Read the diamond output file, either as plain text, or zstd compressed if it ends with .zst
    df = pd.read_csv(diamond_file, sep='\t', header=None)
    
    df.columns = ['seqid', 'sseqid', 'pident', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'cigar']

    # Now, for each seqid, find the best hit based on the bitscore column
    # Sort the data frame by seqid and bitscore
    df = df.sort_values(['seqid', 'bitscore'], ascending=[True, False])
    # Find the best hit for each seqid
    best_hits_for_read = df.groupby('seqid').first().reset_index()
    # save to csv
    best_hits_for_read.to_csv(os.path.join(output_folder, f'{sra_accession}_best_hits_per_read.diamondn.zst'), sep='\t', index=False, header=False)
  
if __name__ == '__main__':
    main()

