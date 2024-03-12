#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList
nextflow.enable.dsl = 2
nextflow.preview.recursion=true 

// seq_data_folder = "/home/charlotte/heackaton/HEK_HeLa_Experiment_folder/HEK_HeLa/20230323_1311_MN32609_FAQ51879_03abc106/fastq_pass/"
// output_dir = "~/out_dir"
// metadata = null
// barcoded = null
// genome_gtf = "/home/charlotte/heackaton/gencode.v43.primary_assembly.annotation.chr20.gtf"
// bed_file = null
// genome_fasta = null
// genome_index = "runMT-human_ont.mmi"
// transcriptome_fasta = null

include { ConvertGtfToDf; start_watch_path } from "./lib/initialize.nf"
include { MinimapIndex ; MinimapGenome ; MinimapTranscriptome } from "./lib/alignment.nf"
include { FeatureCount; CleanFeatureCountTable; UpdateFeatureCountTable; PublishFeatureCountTable; DESeq2Genome ; UpdateIterator } from "./lib/quantify.nf"
include { Salmon; CleanSalmonTable; UpdateSalmonTable; PublishSalmonTable } from "./lib/quantify.nf"

class IntObj{
    public int value = 0
}

iterator = new IntObj()
iterator.value = 1



workflow {
    index_output = MinimapIndex(file("${params.genome_fasta}"),file("${params.transcriptome_fasta}"))
    conversion_output = ConvertGtfToDf(file("${params.genome_gtf}"))
    conversion_output_channel = conversion_output.done
    index_output_channel = index_output.done
    ch_watched_final = start_watch_path()
    ch_watched_final.set{ samples }

    //Run gene alignment, annotation and quantification
    //minimap_index = MinimapIndex(samples,file("${params.genome_fasta}"),file("${params.transcriptome_fasta}"))
    minimap_genome_output= MinimapGenome(samples, file("${params.output_dir}/genome_index.mmi"), conversion_output_channel,index_output_channel)
    fc_output = FeatureCount(minimap_genome_output.aligned_bams, file("${params.genome_gtf}"))
    cleaned_fc = CleanFeatureCountTable(fc_output.single_fc_df, file("${params.metadata}"))
    UpdateFeatureCountTable.scan(cleaned_fc.clean_fc)
    PublishFeatureCountTable(UpdateFeatureCountTable.out)
    UpdateIterator(PublishFeatureCountTable.out)
    //DESeq2Genome(PublishFeatureCountTable.out.merged_csv,file("${params.metadata}"))


    //Run transcript alignment, annotation and quantification
    minimap_transcript_output = MinimapTranscriptome(samples, file("${params.output_dir}/transcriptome_index.mmi"),conversion_output_channel,index_output_channel)
    salmon_output = Salmon(minimap_transcript_output.aligned_bams, file("${params.genome_gtf}"), file("${params.transcriptome_fasta}"))
    cleaned_salmon = CleanSalmonTable(salmon_output.single_salmon_df, file("${params.metadata}"))
    UpdateSalmonTable.scan(cleaned_salmon.clean_salmon)
    PublishSalmonTable(UpdateSalmonTable.out)
    

}
//messages to display once the workflow has completed

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
}

//Rscript script.r "${new_input}" "${state}"