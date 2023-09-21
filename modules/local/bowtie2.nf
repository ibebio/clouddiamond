process bowtie2 {
    tag { sample_id }
    label 'process_high'
    publishDir params.outdir, pattern: "{*.log, *.flagstat.txt}", saveAs: { "Samples/${sample_id}/${it}" }
    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.4' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:a0ffedb52808e102887f6ce600d092675bf3528a-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:a0ffedb52808e102887f6ce600d092675bf3528a-0"
    }

    input:
    tuple val(sample_id), path(fastq_files)
    path db

    output:
    tuple path('*.flagstat.txt')  , emit: TXT
    tuple path("*.log")                  , emit: LOGS

    script:
    """
    echo "#TRACE sample_id=${sample_id}"
    echo "#TRACE n_cpus=${task.cpus}"

    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`    

    # convert the incoming FASTQ file list to an array
    fastq_files=(${fastq_files})
    if [ \${#fastq_files[@]} == 2 ]; then
        # paired-end
        ( bowtie2 -p $task.cpus -x \$INDEX -1 \${fastq_files[0]} -2 \${fastq_files[1]} --met-file ${sample_id}.metrics.txt | samtools view -b -G 12 - > ${sample_id}.bam ) 3>&1 1>&2 2>&3 | tee ${sample_id}.align.log

    else
        # single-end
        ( bowtie2 -p $task.cpus -x \$INDEX -U \${fastq_files[0]} --met-file ${sample_id}.metrics.txt | samtools view -b -F 4 - > ${sample_id}.bam ) 3>&1 1>&2 2>&3 | tee ${sample_id}.align.log
    fi

    # run samtoools flagstat
    samtools flagstat ${sample_id}.bam > ${sample_id}.flagstat.txt
    """


}
