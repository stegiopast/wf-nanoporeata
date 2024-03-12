#!/usr/bin/env nextflow

//Clean main to start coding

import groovy.json.JsonBuilder
@Grab('com.xlson.groovycsv:groovycsv:1.1')
import static com.xlson.groovycsv.CsvParser.parseCsv
import java.io.File
nextflow.enable.dsl = 2



include { minimap_alignment;minimap_transcript_alignment } from './lib/alignment.nf'
include { find_data_folders;make_directories;convert_gtf_to_df;create_feature_percentiles;see_progress} from './lib/initialize.nf'
include { activate_channels;fromPath_data;data_alignment_prep;reawake_next_round } from './lib/realtime_loop.nf'
include { merge_fastq;move_transcript_files;merge_table_of_all_folders;merge_salmon_annotation } from './lib/merge_files.nf'
include { readlength_check;infer_experiment } from './lib/qc.nf'
include { count_features;salmon_annotation } from './lib/quantify.nf'




// workflow module
workflow pipeline {
    take:
        checkup

    main:
        find_data_folders()
        make_directories(find_data_folders.out)
        convert_gtf_to_df(make_directories.out)
        create_feature_percentiles(make_directories.out)
        see_progress(make_directories.out)
        activate_channels(see_progress.out)
        fromPath_data(checkup)
        //string,string_array,data_string,data_seen_string
        string_collection = data_alignment_prep(fromPath_data.out)
        merge_fastq(string_collection.data_string.out, string_collection.data_seen_string.out, string_collection.string.out)
        readlength_check(string_collection.string, string_collection.string_array)

        //Genome alignment
        //bam_string,string_array
        bam_genome_collection = minimap_alignment(string_collection.string,string_collection.string_array)
        fc_collection = count_features(bam_genome_collection.bam_string,bam_genome_collection.string_array)
        merged_tables_collection = merge_table_of_all_folders(fc_collection.fc_string,fc_collection.string_array)
        infer_experiment(merged_tables_collection.fc_string)

        //Transcriptome alignment
        //bam_string,string_array
        bam_transcriptome_collection = minimap_transcript_alignment(string_collection.string,string_collection.string_array)
        moved_bams_collection = move_transcript_files(bam_transcriptome_collection.bam_string,bam_transcriptome_collection.string_array)
        salmon_collection = salmon_annotation(moved_bams_collection.string,moved_bams_collection.string_array)
        merged_salmon_annotation_collection = merge_salmon_annotation(salmon_collection.string,salmon_collection.string_array)
        
        //Restart new iteration
        reawake_next_round(fc_collection.fc_string,merged_salmon_annotation_collection.string,merged_tables_collection.fc_string,merged_tables_collection.string_array)

}

// entrypoint workflow
//WorkflowMain.initialise(workflow, params, log)
workflow {
    //Pinguscript.ping_start(nextflow, workflow, params)
    //Variables will be needed for configuration of the programm, maybe a switch to configuration file coming soon  
    //Folders of interest, annotation and alignment files
    log.info "General Threads: ${params.threads}"
    log.info "General folder: ${params.seq_data_folder}"
    
    params.data_folder = []
    params.sample_names = []
    //Define a path that contains the file you want to align your data to. Usually this file is in .fasta format
    log.info "Genome Fasta chosen: ${params.genome_fasta}"

    //Define a path that contains the transcriptome alignment file
    log.info "Transcriptome Fasta chosen: ${params.transcriptome_fasta}"

    //Define a path that contains the file you want to annotate your data to. Usually this file is in .gtf format 
    log.info "Gtf File chosen: ${params.genome_gtf}"

    //Define a directory in which all the scripts belonging to the application are stored. Especially fc.py and merge_fc.py are important for this script.
    log.info "Script directory: ${params.script_dir}"

    //Define a directory in which you want to memorize all your run specific data and the progress status of your data. This directory will be automatically created if it does not exist so far
    log.info "Run directory: ${params.output_dir}"

    //Conditions that will be compared in DESeq
    log.info "Conditions compared: ${params.conditions}"

    //Are samples barcoded
    log.info "Barcoded: ${params.barcoded}"

    //Metadata file
    log.info "Metadata: ${params.metadata}"
    

    //Reference BED file for RSeQC functions
    log.info "BED File chosen: ${params.bed_file}"

    //Maximal number of files per iteration
    log.info "Maximal batch size: ${params.batchsize}"
    //Channels are important for Activating and reactivating Processes at a given timepoint
    log.info "Type of Sample: ${params.drs}"

    checkup = Channel.empty()

    int continue_bool = 0

    int iteration = 0

    pipeline(checkup)
}

// workflow.onComplete {
//     Pinguscript.ping_complete(nextflow, workflow, params)
// }
// workflow.onError {
//     Pinguscript.ping_error(nextflow, workflow, params)
// }
