#!/bin/bash
set -xv
# This script gets the SRA accession and the aligner name as input parameters
# It downloads the SRA file, converts it to fastq, aligns/diamondn's it to the reference genome
# The output are stats files

# Usage: process_accession.bash <SRA accession>
# Example: process_accession.bash SRR1234567

# Are running everything, or just one step?

echo "$(date) Processing accession $1"

##########################################################################################
# 0. Global configuration
##########################################################################################
SRA_DIR=/s3mnt/scratch
SRA_BUCKET=s3://diamond2-bucket1/scratch
DIAMOND_DB=/data/ref/Col-CEN/ColCEN_ATHILA.dmnd

# Check if the number of arguments is correct
if [ $# -lt 2 ] ; then
    echo "Usage: process_accession.bash <SRA accession> [step]"
    exit 1
fi

# Assign the arguments to variables
SRA_ACCESSION=$1
STEP=$2

FASTQ_DIR=$SRA_DIR/$SRA_ACCESSION/2_fastq

##########################################################################################
# 1. Download the SRA file
##########################################################################################
if [[ "${STEP}" == "" || "${STEP}" == "1_download" ]] ; then 
echo "$(date) 1. Download the SRA file"
# Create an accession directory
DOWNLOAD_DIR=$SRA_DIR/$SRA_ACCESSION/1_download
mkdir -p ${DOWNLOAD_DIR}

# Download the SRA file

# Check if the download was already done and successful (the flag file exists)
if [ ! -f ${SRA_DIR}/${SRA_ACCESSION}/1_download.done ]; then
    # Download the SRA file via AWS CLI
    aws s3 cp s3://sra-pub-run-odp/sra/${SRA_ACCESSION}/${SRA_ACCESSION} ${SRA_BUCKET}/${SRA_ACCESSION}/1_download/

    # Check if the download was successful
    if [ $? -ne 0 ]; then
        echo "Error: prefetch failed for ${SRA_ACCESSION}"
        exit 1
    fi

    # Create the flag file
    touch -f ${SRA_DIR}/${SRA_ACCESSION}/1_download.done
fi
fi # STEP

##########################################################################################
# 2. Convert the SRA file to fastq(s)
##########################################################################################
if [[ "${STEP}" == "" || "${STEP}" == "2_fastq" ]] ; then 
echo "$(date) 2. Convert the SRA file to fastq(s)"
# Create a directory for the fastq files
mkdir -p ${FASTQ_DIR}

# Check if the conversion was already done and successful (the flag file exists)
if [ ! -f ${SRA_DIR}/${SRA_ACCESSION}/2_fastq.done ]; then
    CURR_DIR=$(pwd)
    cd ${FASTQ_DIR}
    # Convert the SRA file to fastq(s)
    fastq-dump --split-3 --skip-technical --readids --clip -O ${FASTQ_DIR} ${DOWNLOAD_DIR}/${SRA_ACCESSION}
    # fasterq-dump --split-3 --skip-technical -O ${FASTQ_DIR} ${DOWNLOAD_DIR}/${SRA_ACCESSION}

    # pigz? not yet


    # Check if the conversion was successful
    if [ $? -ne 0 ]; then
        echo "Error: fastq-dump failed for ${SRA_ACCESSION}"
        exit 1
    fi
    cd ${CURR_DIR}
    # Create the flag file
    touch -f ${SRA_DIR}/{$SRA_ACCESSION}/2_fastq.done
fi
fi # STEP

##########################################################################################
# 3.  Join paired-end reads -- fastq-join
##########################################################################################
if [[ "${STEP}" == "" || "${STEP}" == "3_join" ]] ; then 
# Do we actually have paired-end reads?
if [[ -f ${FASTQ_DIR}/${SRA_ACCESSION}_1.fastq ]] && [[ -f ${FASTQ_DIR}/${SRA_ACCESSION}_2.fastq ]] ; then
    echo "$(date) 3. Join paired-end reads"
    # Create a directory for the joined fastq files
    JOIN_DIR=$SRA_DIR/$SRA_ACCESSION/3_join
    mkdir -p ${JOIN_DIR}

    # Check if the joining was already done and successful (the flag file exists)
    if [ ! -f ${SRA_DIR}/${SRA_ACCESSION}/3_join.done ]; then
        # Join the paired-end reads
        fastq-join ${FASTQ_DIR}/${SRA_ACCESSION}_1.fastq ${FASTQ_DIR}/${SRA_ACCESSION}_2.fastq -o ${JOIN_DIR}/${SRA_ACCESSION}_%.fastq
        # Check if the joining was successful
        if [ $? -ne 0 ]; then
            echo "Error: fastq-join failed for ${SRA_ACCESSION}"
            exit 1
        fi
        # concatenate all reads into one file
        cat ${JOIN_DIR}/${SRA_ACCESSION}_*.fastq > ${JOIN_DIR}/${SRA_ACCESSION}.merged.fastq
        rm ${JOIN_DIR}/${SRA_ACCESSION}_*.fastq
        
        # Create the flag file
        touch -f ${SRA_DIR}/${SRA_ACCESSION}/3_join.done
    fi
else # Just copy the fastq file to the join directory
    echo "$(date) 3. Copying the fastq file to the join directory"
    # Create a directory for the joined fastq files
    JOIN_DIR=$SRA_DIR/$SRA_ACCESSION/3_join
    mkdir -p ${JOIN_DIR}

    # Check if the joining was already done and successful (the flag file exists)
    if [ ! -f ${SRA_DIR}/${SRA_ACCESSION}/3_join.done ]; then
        # Copy the fastq file to the join directory
        cp ${FASTQ_DIR}/${SRA_ACCESSION}.fastq ${JOIN_DIR}/${SRA_ACCESSION}.merged.fastq

        # Check if the copying was successful
        if [ $? -ne 0 ]; then
            echo "Error: copying failed for ${SRA_ACCESSION}"
            exit 1
        fi

        # Create the flag file
        touch -f ${SRA_DIR}/${SRA_ACCESSION}/3_join.done
    fi
fi
fi # STEP


##########################################################################################
# 4.  Blast the fastq(s) agains the reference genome -- diamond blastn
##########################################################################################
if [[ "${STEP}" == "" || "${STEP}" == "4_diamondn" ]] ; then 
echo "$(date) 4. Blast the fastq(s) against the reference genome"
# Create a directory for the blast results
BLAST_DIR=$SRA_DIR/$SRA_ACCESSION/4_diamondn
mkdir -p ${BLAST_DIR}

# Check if the blast was already done and successful (the flag file exists)
CURR_DIR=$(pwd)
if [ ! -f ${SRA_DIR}/${SRA_ACCESSION}/4_diamondn.done ]; then
    # Blast the fastq(s) against the reference genome
    cd ${BLAST_DIR}

    diamond blastn --db ${DIAMOND_DB} -q ${SRA_DIR}/${SRA_ACCESSION}/3_join/${SRA_ACCESSION}.merged.fastq -o ${BLAST_DIR}/${SRA_ACCESSION}.diamondn \
            -f 6 qseqid sseqid pident qstart qend sstart send evalue bitscore qlen slen \
            --shape-mask 11111111111 \
            --minimizer-window 1 \
            --block-size 1.0 \
            -k0 \
            --evalue 0.0001 \
            
        
    # Check if the blast was successful
    if [ $? -ne 0 ]; then
        echo "Error: diamondn failed for ${SRA_ACCESSION}"
        exit 1
    fi
    cd ${CURR_DIR}
    # Create the flag file
    touch -f ${SRA_DIR}/${SRA_ACCESSION}/4_diamondn.done
fi
echo "$(date) Done processing accession $1"
fi # 4 STEP