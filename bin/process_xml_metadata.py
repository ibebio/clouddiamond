#!/usr/bin/env python3

import sys
import os
import re
import argparse

# Script to load an SRA xml file and try to extract the potential arabidopsis accession from the text or other fields
# The script will output the SRA accession, the new ecotype, whether it is a mix, whether it is a wild type, the original ecotype and the genotype
# The script will output the results to stdout.
# Was run like this:
# parallel --jobs 40 python3 ~/sp/cloud-diamond/bin/process_xml_metadata.py {} ~/sp/cloud-diamond/data/ecotypes.txt ::: metadata/*.xml >> ecotypes_genotypes.tsv

def parse_args():
    parser = argparse.ArgumentParser(description='Process SRA XML metadata')
    parser.add_argument('sra_xml', help='SRA XML file')
    parser.add_argument('accession_replacement', help='Accession replacement file')
    # Add a flag to only print the header
    parser.add_argument('--header', action='store_true', help='Only print the header')
    return parser.parse_args()

def main():
    args = parse_args()
    # Verify the arguments

    if args.header:
        print('sra_accession\tecotype\tisMix\tisWild\toriginal_ecotype\toriginal_genotype')
        sys.exit(0)

    if not os.path.exists(args.sra_xml):
        print('SRA XML file does not exist')
        sys.exit(1)
    if not os.path.exists(args.accession_replacement):
        print('Accession replacement file does not exist')
        sys.exit(1)
    sra_xml = args.sra_xml
    accession_replacement = args.accession_replacement

    # Read the SRA XML file
    with open(sra_xml, 'r') as f:
        xml = f.read()

    # Read the accession replacement file, which has a header.
    # The file has 4 , separated columns: Count,Original,Ecotype,isMix. Count is ignored here.
    # Read it as a dictionary, where the key is the original accession and the value is again a dictionary, with 2 keys:
    # the new accession and whether it is a mix with another accession(isMix). isMix can be either 1,0 or NA
    # If only the second column is present, the new accession is the same as the original one and isMix is 0
    # If the third column is present, but the fourth is not, the new accession is the second column and isMix is 0
    
    with open(accession_replacement, 'r') as f:
        next(f)
        accession_dict = {}
        for line in f:
            line = line.strip().split(',')
            if len(line) == 2:
                accession_dict[line[1]] = {'new_accession': line[1], 'isMix': 0}
            elif len(line) == 3:
                accession_dict[line[1]] = {'new_accession': line[2], 'isMix': 0}
            else:
                accession_dict[line[1]] = {'new_accession': line[2], 'isMix': line[3]}


    # Set the SRA accession as the basename without the xml extension
    sra_accession = os.path.basename(sra_xml).replace('.xml', '')

    # Get the ecotype from the TAG field, if it exists. If not, set to NA
    # the corresponding grep line is
    # grep -E -A 1 -i '<TAG>.*ecotype.*</TAG>' metadata/*.xml | grep -v 'TAG' | grep 'VALUE' |cut -f 2 -d '>' | cut -f 1 -d '<'
    # The xml excerpt looks like this:
    # <SAMPLE_ATTRIBUTE>
    #   <TAG>ecotype</TAG>
    #   <VALUE>Col-0</VALUE>
    #</SAMPLE_ATTRIBUTE>
    try:
        ecotype = re.search(r'<SAMPLE_ATTRIBUTE>\s*<TAG>.*?ecotype.*?</TAG>\s*<VALUE>\s*(.*?)\s*</VALUE>\s*</SAMPLE_ATTRIBUTE>', xml).group(1)
        # ecotype = re.search(r'<TAG>.*ecotype.*</TAG>.*<VALUE>(.*?)</VALUE>', xml, re.IGNORECASE).group(1)
    except AttributeError:
        ecotype = 'NA'
    
    # Perform the replacement of the ecotype, according to the dictionary
    if ecotype in accession_dict:
        new_ecotype = accession_dict[ecotype]['new_accession']
        isMix = accession_dict[ecotype]['isMix']
    else:
        new_ecotype = ecotype
        isMix = 'NA'
    


    # Parse the genotype value to determine whether it is a wild type. Check whether the string 'wild' is present in the genotype value.
    # If the genotype tag is not present, set it to NA
    # the corresponding grep line is
    #  grep -E -A 1 -i '<TAG>.*genotype.*</TAG>' metadata/*.xml | grep -v 'TAG' | grep 'VALUE' |cut -f 2 -d '>' | cut -f 1 -d '<' |grep -i wild
    try:
        genotype = re.search(r'<SAMPLE_ATTRIBUTE>\s*<TAG>.*?genotype.*?</TAG>\s*<VALUE>\s*(.*?)\s*</VALUE>\s*</SAMPLE_ATTRIBUTE>', xml).group(1)
        # genotype = re.search(r'<TAG>.*genotype.*</TAG>.*<VALUE>(.*?)</VALUE>', xml, re.IGNORECASE).group(1)
        if 'wild' in genotype.lower():
            isWild = 1
        else:
            isWild = 0
    except AttributeError:
        genotype = 'NA'
        isWild = 'NA'
    

    print(f'{sra_accession}\t{new_ecotype}\t{isMix}\t{isWild}\t{ecotype}\t{genotype}')


if __name__ == '__main__':
    main()
