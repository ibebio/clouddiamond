/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)



// Validate input parameters
// WorkflowCloudDiamond.initialise(workflow, params, log)

// Define the sentinel for "done" signals // not sure if needed, copied from GEMaker
// DONE_SENTINEL = 1


/*
========================================================================================
    CONFIG FILES
========================================================================================
*/


/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
// def modules = params.modules.clone()


//
// MODULE: Local to the pipeline
//
// include { GET_SOFTWARE_VERSIONS as get_software_versions } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

// These 2 might me  needed, if the currently implemented nextflow temp cleanup is not suitable
// Module: clean_work_dirs
// include { clean_work_dirs } from '../modules/local/clean_work_dirs' addParams()

// Module: clean_work_files
// include { clean_work_files } from '../modules/local/clean_work_files' addParams()


// THIS IS THE SRA STUFF FROM GEMMaker
// ---------------------------------------------------------------------------------------

// Module: retrieve_sra_metadata
// Note: we have to pass in the work directory in pieces with a 'path:'
// prefix because Nextflow will recoginze the path, see that the
// diretory attributes have changed and won't use the cache. So, this is
// a bit of a hack to keep Nextflow from rerunning the retreive on a resume
// if it completed successfully.
include { retrieve_sra_metadata } from '../modules/local/retrieve_sra_metadata' addParams(workDirParent: "path:" + workflow.workDir.getParent(), workDirName: workflow.workDir.getName())

// Module: download_runs
// Note: we have to pass in the work directory in pieces with a 'path:'
// prefix because Nextflow will recoginze the path, see that the
// diretory attributes have changed and won't use the cache. So, this is
// a bit of a hack to keep Nextflow from rerunning the retreive on a resume
// if it completed successfully.
include { download_runs } from '../modules/local/download_runs' addParams(workDirParent: "path:" + workflow.workDir.getParent(), workDirName: workflow.workDir.getName())

// Module: failed_run_report
include { failed_run_report } from '../modules/local/failed_run_report' addParams()

// Module: fastq_dump
publish_pattern_fastq_dump = params.keep_retrieved_fastq
    ? "{*.fastq}"
    : "{none}"
include { fastq_dump } from '../modules/local/fastq_dump' addParams(publish_pattern_fastq_dump: publish_pattern_fastq_dump)


// Module: fastq_merge NEEDED?
include { fastq_merge } from '../modules/local/fastq_merge' addParams()
// ---------------------------------------------------------------------------------------
// END SRA STUFF FROM GEMMaker


// THIS IS THE DIAMOND ALIGNMENT STUFF

// include { diamond_blastx } from '../modules/local/diamond_blastx' addParams()

include { bowtie2 } from '../modules/local/bowtie2' addParams()


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow CloudDiamond {
    
    retrieve_sra_metadata(file(params.sra_accessions_list))
    SAMPLES_LIST = retrieve_sra_metadata.out.SRR2SRX
                .splitCsv()
                .groupTuple(by: 1)
                .map { [it[1], it[0].join(" "), "remote"] }
     

    // SAMPLES_BATCH = SAMPLES_LIST.buffer(size: 10, remainder: true)
    // SAMPLES_BATCH = SAMPLES_LIST
    // SAMPLES_BATCH.view()

    print(SAMPLES_LIST)
    download_runs(SAMPLES_LIST)
    SRA_FILES = download_runs.out.SRA_FILES

    fastq_dump(SRA_FILES)
    DOWNLOADED_FASTQ_FILES = fastq_dump.out.FASTQ_FILES

    fastq_merge(DOWNLOADED_FASTQ_FILES)
    MERGED_FASTQ_FILES = fastq_merge.out.FASTQ_FILES

    FAILED_SAMPLES = Channel.empty()
        .mix(
            download_runs.out.FAILED_SAMPLES.map { it[0] },
            fastq_dump.out.FAILED_SAMPLES.map { it[0] })


    FASTQ_FILES = MERGED_FASTQ_FILES

    // Run alignment
    if (params.aligner == 'diamond_blastx') {
        diamond_blastx(FASTQ_FILES, params.reference_db)
        ALIGNER_STATS = diamond_blastx.out.TXT
    } else if (params.aligner == 'bowtie2') {
        bowtie2(FASTQ_FILES, params.reference_db)
        ALIGNER_STATS = bowtie2.out.TXT
    } else {
        error "Unknown aligner: ${params.aligner}"
    }





}