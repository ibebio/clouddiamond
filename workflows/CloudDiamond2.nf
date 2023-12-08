params.sra_accessions_list = "accessions.txt"



////////////////////////////////////////////////////////////////////////////////
// 0. Global configuration, can/should go into nextflow config
////////////////////////////////////////////////////////////////////////////////
params.SRA_DIR="/s3mnt/scratch"
params.SRA_BUCKET="s3://diamond2-bucket1/scratch"
params.DIAMOND_DB="/data/ref/Col-CEN/ColCEN_ATHILA.dmnd"



////////////////////////////////////////////////////////////////////////////////
// 1. Download the SRA file
////////////////////////////////////////////////////////////////////////////////

process download_sra {
    tag { SRA_ACCESSION }
    label 'process_low'
    publishDir params.outdir, pattern: "*.DONE", 
    saveAs: { "${SRA_ACCESSION}/${it}" }

    input:
    val SRA_ACCESSION

    output:
    file "1_download.DONE"
    val SRA_ACCESSION

    script:
    """
    echo 1_download_sra $SRA_ACCESSION
    process_accession.bash $SRA_ACCESSION 1_download
    if [ \$? -eq 0 ] ; then
        touch -f 1_download.DONE
    fi
    """

}


process convert_sra_to_fastq {
    tag { SRA_ACCESSION }
    label 'process_low'
    publishDir params.outdir, pattern: "*.DONE", 
    saveAs: { "${SRA_ACCESSION}/${it}" }

    input:
    file download_sra_done
    val SRA_ACCESSION

    output:
    file "2_fastq.DONE"
    val SRA_ACCESSION

    script:
    """
    echo 2_convert_sra_to_fastq $SRA_ACCESSION
    process_accession.bash $SRA_ACCESSION 2_fastq
    if [ \$? -eq 0 ] ; then
        touch -f 2_fastq.DONE
    fi
    """

}


process join_paired_end_reads {
    tag { SRA_ACCESSION }
    label 'process_medium'
    publishDir params.outdir, pattern: "*.DONE", 
    saveAs: { "${SRA_ACCESSION}/${it}" }

    input:
    file convert_sra_to_fastq_done
    val SRA_ACCESSION

    output:
    file "3_join.DONE"
    val SRA_ACCESSION

    script:
    """
    echo 3_join_paired_end_reads $SRA_ACCESSION
    process_accession.bash $SRA_ACCESSION 3_join
    if [ \$? -eq 0 ] ; then
        touch -f 3_join.DONE
    fi
    """

}


process diamond_fastqs {
    tag { SRA_ACCESSION }
    label 'process_high'
    publishDir params.outdir, pattern: "*.DONE", 
    saveAs: { "${SRA_ACCESSION}/${it}" }

    input:
    file join_paired_end_reads_done
    val SRA_ACCESSION

    output:
    file "4_diamond_fastqs.DONE"
    val SRA_ACCESSION

    script:
    """
    echo 4_dimond_fastqs $SRA_ACCESSION
    process_accession.bash $SRA_ACCESSION 4_diamondn
    if [ \$? -eq 0 ] ; then
        touch -f 4_diamondn.DONE
    fi
    """


}


workflow CloudDiamond {
    Channel.fromPath(params.sra_accessions_list) \
    | splitCsv(header:false) \
    | map {line -> line[0] } \
    | download_sra \
    | convert_sra_to_fastq \
    | join_paired_end_reads \
    | diamond_fastqs
    
}



// // Step 0: Get a list of SRA accessions from a text file
// process GetAccessions {
//     input:
//     file accessionFile

//     output:
//     file "accessions.txt"

//     script:
//     """
//     # Command to extract SRA accessions from the text file and save them to accessions.txt
//     # Replace this command with the actual command to extract accessions from the file
//     grep "^SRA" ${accessionFile} > accessions.txt
//     """
// }

// // Step 1: Download the SRA file
// process DownloadSRA {
//     input:
//     file accessionFile

//     output:
//     file "sra_file.DONE"

//     script:
//     """
//     # Command to download the SRA file using the accession from accessions.txt
//     # Replace this command with the actual command to download the SRA file
//     wget -O sra_file.sra $(cat accessions.txt)
//     touch sra_file.DONE
//     """
// }

// // Step 2: Convert the SRA file to fastq
// process ConvertToFastq {
//     input:
//     file "sra_file.DONE"

//     output:
//     file "fastq_file.DONE"

//     script:
//     """
//     # Command to convert the SRA file to fastq format
//     # Replace this command with the actual command to convert the SRA file to fastq
//     fastq-dump sra_file.sra
//     touch fastq_file.DONE
//     """
// }

// // Step 3: Join paired end reads
// process JoinReads {
//     input:
//     file "fastq_file.DONE"

//     output:
//     file "joined_reads.DONE"

//     script:
//     """
//     # Command to join paired end reads
//     # Replace this command with the actual command to join the fastq reads
//     join-reads fastq_file.fastq
//     touch joined_reads.DONE
//     """
// }

// // Step 4: Blast fastqs against the reference genome
// process BlastFastq {
//     input:
//     file "joined_reads.DONE"

//     output:
//     file "blast_results.DONE"

//     script:
//     """
//     # Command to blast fastqs against the reference genome
//     # Replace this command with the actual command to perform the blast
//     blast fastq_file.fastq reference_genome.fasta
//     touch blast_results.DONE
//     """
// }

// // Define the workflow
// workflow {
//     // Get the list of SRA accessions
//     GetAccessions(accessionFile)

//     // Perform the following steps for each accession
//     Channel.fromPath("accessions.txt").set { accessions }

//     accessions.into { accession ->
//         // Download the SRA file
//         DownloadSRA(accession)

//         // Convert the SRA file to fastq
//         ConvertToFastq()

//         // Join paired end reads
//         JoinReads()

//         // Blast fastqs against the reference genome
//         BlastFastq()
//     }
// }
