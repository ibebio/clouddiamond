#!/usr/bin/env python3

# This scripts analyzes the diamond output file. It takes the following inputs:
# 1. Diamond output file, either plain text or zstd compressed
# 2. Output folder, it will be created if it does not exist
# 3. SRA accession id, it will be used to name the output files
# 4. SRA XML file, it is used to get the number of reads, type of sequencer, year lab, and if available, ecotype

# The statics are computed using the pandas library.
# The input file has the following columns:
# seqid sseqid pident qstart qend sstart send evalue bitscore qlen slen cigar

# The following statistics are computed:
# 1. Number of reads
# 2. Number of hits
# 3. Sequencer type
# 4. Year when sequenced
#.5. Lab where sequenced
# 6. Ecotype
# 7. Number of best hits


import sys, os, re, argparse
import pandas as pd
import numpy as np
import zstandard as zstd

def parse_args():
    parser = argparse.ArgumentParser(description='Analyze diamond output file')
    parser.add_argument('diamond_file', help='Diamond output file')
    parser.add_argument('output_folder', help='Output folder')
    parser.add_argument('sra_accession', help='SRA accession id')
    parser.add_argument('sra_xml', help='SRA XML file')
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
    if not os.path.exists(args.sra_xml):
        print('SRA XML file does not exist')
        sys.exit(1)
    diamond_file = args.diamond_file
    output_folder = args.output_folder
    sra_accession = args.sra_accession
    sra_xml = args.sra_xml

    # Create output folder if it does not exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    
    ###########################################################################
    # Process SRA metadata from XML file
    ###########################################################################

    # Read the SRA XML file
    with open(sra_xml, 'r') as f:
        xml = f.read()
    # Get the number of reads, which is the "spots" attribute in RUN_SET
    num_reads = re.search(r'<RUN_SET.*?spots="(.*?)"', xml).group(1)
    print(f'Number of reads: {num_reads}')    

    # Get the number of bases, which is the "bases" attribute in RUN_SET
    num_bases = re.search(r'<RUN_SET.*?bases="(.*?)"', xml).group(1)
    print(f'Number of bases: {num_bases}')

    # Get the sequencer type
    sequencer = re.search(r'<INSTRUMENT_MODEL>(.*?)</INSTRUMENT_MODEL>', xml).group(1)
    print(f'Sequencer type: {sequencer}')

    # Get the year when sequenced, which is in the "published" attribute in RUN. The published attribute is "2022-07-14 22:42:20"
    year = re.search(r'<RUN.*?published="(.*?)"', xml).group(1)
    year = year.split('-')[0]

    print(f'Year sequenced: {year}')

    # Get the lab where sequenced. It is the center_name attribute in STUDY, like <STUDY accession="ERP138859" alias="9ec0e668-4413-45f3-9c24-051782dcfe7b" center_name="max planck institute for biology tuebingen">
    lab = re.search(r'<STUDY.*?center_name="(.*?)"', xml).group(1)
    print(f'Lab: {lab}')

    # Get the ecotype. This often does not exist, then set it to Not available. If it exists, it is in SAMPLE_ATTRIBUTE, like <SAMPLE_ATTRIBUTE><TAG>ecotype</TAG><VALUE>Col-0</VALUE></SAMPLE_ATTRIBUTE>    
    try:
        ecotype = re.search(r'<SAMPLE_ATTRIBUTE><TAG>ecotype</TAG><VALUE>(.*?)</VALUE></SAMPLE_ATTRIBUTE>', xml).group(1)
    except:
        ecotype = 'Not available'
    print(f'Ecotype: {ecotype}')

    # Get the genotype. This often does not exist, then set it to Not available. If it exists, it is in SAMPLE_ATTRIBUTE, like <SAMPLE_ATTRIBUTE><TAG>genotype</TAG><VALUE>Col-0</VALUE></SAMPLE_ATTRIBUTE> 
    try:
        genotype = re.search(r'<SAMPLE_ATTRIBUTE><TAG>genotype</TAG><VALUE>(.*?)</VALUE></SAMPLE_ATTRIBUTE>', xml).group(1)
    except:
        genotype = 'Not available'
    print(f'Genotype: {genotype}')

    # Get the library strategy. This is like    <LIBRARY_STRATEGY>ATAC-seq</LIBRARY_STRATEGY>
    try:
        library_strategy = re.search(r'<LIBRARY_STRATEGY>(.*?)</LIBRARY_STRATEGY>', xml).group(1)
    except:
        library_strategy = 'Not available'
    print(f'Library strategy: {library_strategy}')

    # Get the library source. This is like    <LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>
    try:
        library_source = re.search(r'<LIBRARY_SOURCE>(.*?)</LIBRARY_SOURCE>', xml).group(1)
    except:
        library_source = 'Not available'
    print(f'Library source: {library_source}')


    ###########################################################################
    # Process diamond output file
    ###########################################################################

    # Read the diamond output file, either as plain text, or zstd compressed if it ends with .zst
    df = pd.read_csv(diamond_file, sep='\t', header=None)
    
    df.columns = ['seqid', 'sseqid', 'pident', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'cigar']

    # Number of hits
    num_hits = df.shape[0]
    print(f'Number of hits: {num_hits}')

    ###########################################################################
    # Number of best hits
    # We use the following thresholds, defined as percentiles of the bitscore column
    # - 90% percentile  -- this is the data frame that is going to be stored, thus reduces the dataset to 10%
    #   of the original size. Since it is sorted, the other thresholds are easy to compute. We also store the correspoding bitscore values for each threshold.
    # - 95% percentile
    # - 99% percentile
    # - 99.9% percentile
    #
    # Save these best hits in an extra data frame
    thresholds = [90, 95, 99, 99.9]
    best_hits = {}
    best_hits_thresholds = {} # Store the bitscore values for each threshold
    num_best_hits = {}
    for threshold in thresholds:
        best_hits_thresholds[threshold] = np.percentile(df['bitscore'], threshold)
        best_hits[threshold] = df[df['bitscore'] >= best_hits_thresholds[threshold]]
        # Sort the best hits by bitscore
        best_hits[threshold] = best_hits[threshold].sort_values('bitscore', ascending=False)
        num_best_hits[threshold] = best_hits[threshold].shape[0]
        print(f'Number of best hits at {threshold} percentile: {num_best_hits[threshold]}')
        # Save the 90% and 99.9% best hits in a file
        if threshold == 90 or threshold == 99.9:
            best_hits[threshold].to_csv(os.path.join(output_folder, f'{sra_accession}_best_hits{threshold}.diamondn.zst'), sep='\t', index=False, header=False)
    
    # Save the bitscore thresholds from the best_hits dict in a file
    with open(os.path.join(output_folder, f'{sra_accession}_best_hits_thresholds.txt'), 'w') as f:
        # Write the header
        f.write('threshold\tbitscore\n')
        for threshold in thresholds:
            f.write(f'{threshold}\t{best_hits_thresholds[threshold]}\n')
    #
    ###########################################################################        
    

    # Save the statistics in a file as a table, and as a JSON file
    stats = pd.DataFrame({'sra_accession': [sra_accession], 'num_reads': [num_reads], 'num_bases': [num_bases], 'num_hits': [num_hits], 
                          'sequencer': [sequencer], 
                          'year': [year], 'lab': [lab], 'ecotype': [ecotype], 
                          'genotype': [genotype], 
                          'library_strategy': [library_strategy], 'library_source': [library_source], 
                          'num_best_hits90': [num_best_hits[90]], 'num_best_hits95': [num_best_hits[95]], 
                          'num_best_hits99': [num_best_hits[99]], 'num_best_hits999': [num_best_hits[99.9]]})
    
    stats.to_csv(os.path.join(output_folder, f'{sra_accession}_stats.txt'), sep='\t', index=False)
    stats.to_json(os.path.join(output_folder, f'{sra_accession}_stats.json'))
    

    ###################################
    # Compute statistics based on the sseqid column
    ###################################
    # Number of hits per sseqid
    hits_per_sseqid = df['sseqid'].value_counts()
    
    # Number of hits per sseqid, normalized by the number of bases * 1000 * 1000
    hits_per_sseqid_normalized = hits_per_sseqid / int(num_bases) * 1000 * 1000
    
    # Flag the sseqids that are best hits
    # Add best_hit column to the data frame
    hits_per_sseqid = hits_per_sseqid.to_frame()
    for threshold in thresholds:
        hits_per_sseqid[f'best_hit{threshold}'] = hits_per_sseqid.index.isin(best_hits[threshold]['sseqid'])
        # Convert the best_hit column to int
        hits_per_sseqid[f'best_hit{threshold}'] = hits_per_sseqid[f'best_hit{threshold}'].astype(int)
    
    # Join the two data frames
    hits_per_sseqid = pd.concat([hits_per_sseqid, hits_per_sseqid_normalized], axis=1)
    
    # Rename the columns
    hits_per_sseqid.columns = ['num_hits', 'best_hit90', 'best_hit95', 'best_hit99', 'best_hit999', 'num_hits_normalized']
    

    ###################################
    # Compute coverage statistics
    ###################################



    ###################################
    # Save the main stats data frame
    ###################################
    hits_per_sseqid.to_csv(os.path.join(output_folder, f'{sra_accession}_hits_per_sseqid.txt.zst'), sep='\t')

    

if __name__ == '__main__':
    main()

