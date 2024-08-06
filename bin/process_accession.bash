#!/bin/bash
#set -xv
set -o pipefail
set -e
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
PATH=${PATH}:/home/ec2-user/micromamba/envs/base/bin/:/home/ec2-user/sratoolkit.3.0.7-centos_linux64/bin
DELETE_PREV=1 # Delete data from previous steps. 3rd step joined fastq is not deleted, compressed later!

# Check if the number of arguments is correct
if [ $# -lt 1 ] ; then
    echo "Usage: process_accession.bash <SRA accession> [step]"
    exit 1
fi

# Assign the arguments to variables
SRA_ACCESSION=$1
STEP=$2

# Bail out on undefined variables, set this after the positional parameters have been assigned
set -u

DOWNLOAD_DIR=$SRA_DIR/$SRA_ACCESSION/1_download
FASTQ_DIR=$SRA_DIR/$SRA_ACCESSION/2_fastq
JOIN_DIR=$SRA_DIR/$SRA_ACCESSION/3_join
BLAST_DIR=$SRA_DIR/$SRA_ACCESSION/4_diamondn

##########################################################################################
# 1. Download the SRA file
##########################################################################################
if [[ "${STEP}" == "" || "${STEP}" == "1_download" ]] ; then 
echo "$(date) 1. Download the SRA file"
# Create an accession directory
mkdir -p ${DOWNLOAD_DIR}

# Download the SRA file

# Check if the download was already done and successful (the flag file exists)
if [ ! -f ${SRA_DIR}/${SRA_ACCESSION}/1_download.done ]; then
    # Download the SRA file via AWS CLI
    /home/ec2-user/micromamba/envs/base/bin/aws s3 cp s3://sra-pub-run-odp/sra/${SRA_ACCESSION}/${SRA_ACCESSION} ${SRA_BUCKET}/${SRA_ACCESSION}/1_download/

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
if [ ! -f ${SRA_DIR}/${SRA_ACCESSION}/2_fastq.done ] && [ -f ${SRA_DIR}/${SRA_ACCESSION}/1_download.done ] ; then
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
    touch -f ${SRA_DIR}/${SRA_ACCESSION}/2_fastq.done

    # Delete previous step data
    if [[ "${DELETE_PREV}" == "1" ]] ; then
        rm ${DOWNLOAD_DIR}/${SRA_ACCESSION}
    fi


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
    if [ ! -f ${SRA_DIR}/${SRA_ACCESSION}/3_join.done ] && [ -f ${SRA_DIR}/${SRA_ACCESSION}/2_fastq.done ] ; then
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
    if [ ! -f ${SRA_DIR}/${SRA_ACCESSION}/3_join.done ] && [ -f ${SRA_DIR}/${SRA_ACCESSION}/2_fastq.done ] ; then
        # Copy the fastq file to the join directory
        cp ${FASTQ_DIR}/${SRA_ACCESSION}_1.fastq ${JOIN_DIR}/${SRA_ACCESSION}.merged.fastq

        # Check if the copying was successful
        if [ $? -ne 0 ]; then
            echo "Error: copying failed for ${SRA_ACCESSION}"
            exit 1
        fi

        # Create the flag file
        touch -f ${SRA_DIR}/${SRA_ACCESSION}/3_join.done
    fi

    # Delete data from previous step
    if [[ "${DELETE_PREV}" == "1" ]] ; then
        rm ${FASTQ_DIR}/${SRA_ACCESSION}*.fastq
    fi

fi
fi # STEP


##########################################################################################
# 4.  Blast the fastq(s) agains the reference genome -- diamond blastn
##########################################################################################
if [[ "${STEP}" == "" || "${STEP}" == "4_diamondn" ]] ; then 
echo "$(date) 4. Blast the fastq(s) against the reference genome"
# Create a directory for the blast results
mkdir -p ${BLAST_DIR}

# Check if the blast was already done and successful (the flag file exists)
CURR_DIR=$(pwd)
if [ ! -f ${SRA_DIR}/${SRA_ACCESSION}/4_diamondn.done ] && [ -f ${SRA_DIR}/${SRA_ACCESSION}/3_join.done ] ; then
    # Blast the fastq(s) against the reference genome
    cd ${BLAST_DIR}

    # 1. ColCEN
    OUT_DIR=${BLAST_DIR}/ColCEN_ATHILA
    mkdir -p ${OUT_DIR}
    if [[ ! -f ${OUT_DIR}/4a_ColCEN_ATHILA.DONE ]] ; then
        diamond blastn --db ${DIAMOND_DB} -q ${SRA_DIR}/${SRA_ACCESSION}/3_join/${SRA_ACCESSION}.merged.fastq -o ${OUT_DIR}/${SRA_ACCESSION}.diamondn \
                -f 6 qseqid sseqid pident qstart qend sstart send evalue bitscore qlen slen cigar \
                --shape-mask 11111111111 \
                --minimizer-window 1 \
                --block-size 1.0 \
                --masking 0 \
                -k0 \
                --evalue 0.0001 
        # Check if the blast was successful
        if [ $? -ne 0 ]; then
            echo "Error: diamondn failed for ${SRA_ACCESSION}"
            exit 1
        fi
        touch -f ${OUT_DIR}/4a_ColCEN_ATHILA.DONE
    fi

    # 2. Ath-Rep
    OUT_DIR=${BLAST_DIR}/Ath-Rep
    mkdir -p ${OUT_DIR}
    if [[ ! -f ${OUT_DIR}/4b_Ath-Rep.DONE ]] ; then
        diamond blastn --db /data/ref/Ath-Rep/Ath-Rep.dmnd -q ${SRA_DIR}/${SRA_ACCESSION}/3_join/${SRA_ACCESSION}.merged.fastq -o ${OUT_DIR}/${SRA_ACCESSION}.diamondn \
                -f 6 qseqid sseqid pident qstart qend sstart send evalue bitscore qlen slen cigar \
                --shape-mask 11111111111 \
                --minimizer-window 1 \
                --block-size 1.0 \
                --masking 0 \
                -k0 \
                --evalue 0.001 
        # Check if the blast was successful
        if [ $? -ne 0 ]; then
            echo "Error: diamondn failed for ${SRA_ACCESSION}"
            exit 1
        fi
        touch -f ${OUT_DIR}/4b_Ath-Rep.DONE
    fi
    # 3. Col0-Telib
    OUT_DIR=${BLAST_DIR}/Col0-Telib
    mkdir -p ${OUT_DIR}
    if [[ ! -f ${OUT_DIR}/4c_Col0-Telib.DONE ]] ; then
        diamond blastn --db /data/ref/Col0-Telib/Col0-Telib.dmnd -q ${SRA_DIR}/${SRA_ACCESSION}/3_join/${SRA_ACCESSION}.merged.fastq -o ${OUT_DIR}/${SRA_ACCESSION}.diamondn \
                -f 6 qseqid sseqid pident qstart qend sstart send evalue bitscore qlen slen cigar \
                --shape-mask 11111111111 \
                --minimizer-window 1 \
                --block-size 1.0 \
                -k0 \
                --evalue 0.001 \
                --masking 0 
        # Check if the blast was successful
        if [ $? -ne 0 ]; then
            echo "Error: diamondn failed for ${SRA_ACCESSION}"
            exit 1
        fi
        touch -f ${OUT_DIR}/4c_Col0-Telib.DONE
    fi

    cd ${CURR_DIR}

    # Delete data from previous step
    #if [[ "${DELETE_PREV}" == "1" ]] ; then
    #    rm ${JOIN_DIR}/${SRA_ACCESSION}*.fastq
    #fi

    # Create the flag file
    touch -f ${SRA_DIR}/${SRA_ACCESSION}/4_diamondn.done
fi
echo "$(date) Done processing accession $1"
fi # 4 STEP


##########################################################################################
# 5.  Compress outputs
##########################################################################################
if [[ "${STEP}" == "" || "${STEP}" == "5_compress" ]] ; then 
echo "$(date) 5. Compress output and fastq merged input"

# Check if the compression was already done and successful
CURR_DIR=$(pwd)
if [ ! -f ${SRA_DIR}/${SRA_ACCESSION}/5_compress.done ] && [ -f ${SRA_DIR}/${SRA_ACCESSION}/4_diamondn.done ] ; then
 
    # Compress merged fastq
    if [[ ! -f ${SRA_DIR}/${SRA_ACCESSION}/3_join/5a_compress.DONE ]] ; then
        cat ${SRA_DIR}/${SRA_ACCESSION}/3_join/${SRA_ACCESSION}.merged.fastq | zstd -T0 --adapt -c - > ${SRA_DIR}/${SRA_ACCESSION}/3_join/${SRA_ACCESSION}.merged.fastq.zst
        if [ $? -ne 0 ]; then
            echo "Error: compression failed for ${SRA_ACCESSION} fastq"
            exit 1
        fi
        touch -f ${SRA_DIR}/${SRA_ACCESSION}/3_join/5a_compress.DONE
        rm ${SRA_DIR}/${SRA_ACCESSION}/3_join/${SRA_ACCESSION}.merged.fastq
    fi

    # Compress diamond  fastq
    if [[ ! -f ${BLAST_DIR}/ColCEN_ATHILA/5b_compress.DONE ]] ; then
        cat ${BLAST_DIR}/ColCEN_ATHILA/${SRA_ACCESSION}.diamondn | zstd -T0 --adapt -c - > ${BLAST_DIR}/ColCEN_ATHILA/${SRA_ACCESSION}.diamondn.zst
        if [ $? -ne 0 ]; then
            echo "Error: compression failed for ${SRA_ACCESSION} ColCEN_ATHILA"
            exit 1
        fi
        touch -f ${BLAST_DIR}/ColCEN_ATHILA/5b_compress.DONE
        rm ${BLAST_DIR}/ColCEN_ATHILA/${SRA_ACCESSION}.diamondn
    fi

    # Compress diamond  fastq
    if [[ ! -f ${BLAST_DIR}/Ath-Rep/5c_compress.DONE ]] ; then
        cat ${BLAST_DIR}/Ath-Rep/${SRA_ACCESSION}.diamondn | zstd -T0 --adapt -c - > ${BLAST_DIR}/Ath-Rep/${SRA_ACCESSION}.diamondn.zst
        if [ $? -ne 0 ]; then
            echo "Error: compression failed for ${SRA_ACCESSION} Ath-Rep"
            exit 1
        fi
        touch -f ${BLAST_DIR}/Ath-Rep/5c_compress.DONE
        rm ${BLAST_DIR}/Ath-Rep/${SRA_ACCESSION}.diamondn
    fi

    # Compress diamond  fastq
    if [[ ! -f ${BLAST_DIR}/Col0-Telib/5d_compress.DONE ]] ; then
        cat ${BLAST_DIR}/Col0-Telib/${SRA_ACCESSION}.diamondn | zstd -T0 --adapt -c - > ${BLAST_DIR}/Col0-Telib/${SRA_ACCESSION}.diamondn.zst
        if [ $? -ne 0 ]; then
            echo "Error: compression failed for ${SRA_ACCESSION} Col0-Telib"
            exit 1
        fi
        touch -f ${BLAST_DIR}/Col0-Telib/5d_compress.DONE
        rm ${BLAST_DIR}/Col0-Telib/${SRA_ACCESSION}.diamondn
    fi
    
    
    # Create the flag file
    touch -f ${SRA_DIR}/${SRA_ACCESSION}/5_compress.done

fi # .done

fi # Step
