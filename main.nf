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

include { ConvertGtfToDf; CreateFeaturePercentiles; start_watch_path } from "./lib/initialize.nf"
include { MinimapIndex ; MinimapGenome ; MinimapTranscriptome } from "./lib/alignment.nf"
include { FeatureCount; CleanFeatureCountTable; UpdateFeatureCountTable; PublishFeatureCountTable; DESeq2Genome ; DESeq2Transcriptome; DTUanalysis; UpdateIterator ; UpdateIterator2 } from "./lib/quantify.nf"
include { Salmon; CleanSalmonTable; UpdateSalmonTable; PublishSalmonTable } from "./lib/quantify.nf"
include { RunDevelopmentEstimation; CountMappedReads; MergeMappedReadsTable; UpdateReadLengthDistribution } from "./lib/qc.nf"


workflow {
    feature_percentiles = CreateFeaturePercentiles(file("${params.bed_file}"))
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
    RunDevelopmentEstimation(UpdateIterator.out.run_statistics, PublishFeatureCountTable.out.merged_all, file("${params.output_dir}/inner_variability_plot.csv"),file("${params.output_dir}/inner_variability_per_sample.csv"),file("${params.output_dir}/exp_genes_counted_per_sample.csv"))
    DESeq2Genome(UpdateIterator.out.run_statistics,PublishFeatureCountTable.out.merged_all,file("${params.metadata}"))

    //Run read mapping quantification
    counted_reads = CountMappedReads(minimap_genome_output.aligned_bams, minimap_genome_output.aligned_bam_bais)
    MergeMappedReadsTable(counted_reads.read_counts,file("${params.metadata}"),file("${params.output_dir}/mapping_stats.txt"))

    //Read size distribution
    UpdateReadLengthDistribution(samples)
    
    //Process run time



    //Run transcript alignment, annotation and quantification
    minimap_transcript_output = MinimapTranscriptome(samples, file("${params.output_dir}/transcriptome_index.mmi"),conversion_output_channel,index_output_channel)
    salmon_output = Salmon(minimap_transcript_output.aligned_bams, file("${params.genome_gtf}"), file("${params.transcriptome_fasta}"))
    cleaned_salmon = CleanSalmonTable(salmon_output.single_salmon_df, file("${params.metadata}"))
    UpdateSalmonTable.scan(cleaned_salmon.clean_salmon)
    PublishSalmonTable(UpdateSalmonTable.out)
    UpdateIterator2(PublishSalmonTable.out)
    DESeq2Transcriptome(UpdateIterator2.out.run_statistics,PublishSalmonTable.out.merged_all,file("${params.metadata}"))
    DTUanalysis(DESeq2Transcriptome.out.dea_transcriptome_done,PublishSalmonTable.out.merged_all,file("${params.metadata}"), file("${params.output_dir}/converted_gtf.csv"))
    

}
//messages to display once the workflow has completed

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
}


//Rscript script.r "${new_input}" "${state}"