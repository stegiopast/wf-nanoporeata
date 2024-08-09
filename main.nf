#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList
import java.time.LocalDateTime
nextflow.enable.dsl = 2
nextflow.preview.recursion=true 

include { ConvertGtfToDf; CreateFeaturePercentiles; start_watch_path ; fetch_latest_bams ; CreateGenomeBamFilesForMerge; CreateTranscriptomeBamFilesForMerge; CopyBedfileAndMetadata } from "./lib/initialize.nf"
include { MinimapIndex ; MinimapGenome ; MinimapTranscriptome ; MinimapGenomeMergeBam ; MinimapTranscriptomeMergeBam; FindMergedGenomeBams; FindMergedTranscriptomeBams } from "./lib/alignment.nf"
include { FeatureCount; CleanFeatureCountTable; UpdateFeatureCountTable; PublishFeatureCountTable; DESeq2Genome ; DESeq2Transcriptome; DTUanalysis; UpdateIterator ; UpdateIterator2 } from "./lib/quantify.nf"
include { Salmon; CleanSalmonTable; UpdateSalmonTable; PublishSalmonTable } from "./lib/quantify.nf"
include { RunDevelopmentEstimation; CountMappedReads; MergeMappedReadsTable; DefineReadLengthDistribution; UpdateReadLengthDistribution; PublishReadLengthDistribution } from "./lib/qc.nf"


workflow {
    // Limit parralelsims globally
    CreateGenomeBamFilesForMerge(file("${params.metadata}"))
    CreateTranscriptomeBamFilesForMerge(file("${params.metadata}"))
    CopyBedfileAndMetadata(file("${params.bed_file}"), file("${params.metadata}"))
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

    //Merge single bams 
    
    //def latest_bams_genome = fetch_latest_bams("${params.output_dir}/bam_genome_merged/")
    //found_genome_bams = FindMergedGenomeBams(minimap_genome_output.aligned_bams)
    // def genome_bams = file("${params.output_dir}/bam_genome_merged/*.bam").collect()
    MinimapGenomeMergeBam(minimap_genome_output.aligned_bams, channel.fromPath("${params.output_dir}/bam_genome_merged/*.bam").collect())
    

    //Run read mapping quantification
    counted_reads = CountMappedReads(minimap_genome_output.aligned_bams, minimap_genome_output.aligned_bam_bais)
    MergeMappedReadsTable(counted_reads.read_counts,file("${params.metadata}"),file("${params.output_dir}/mapping_stats.txt"))

    //Read size distribution
    DefineReadLengthDistribution(samples,file("${params.metadata}"))
    UpdateReadLengthDistribution.scan(DefineReadLengthDistribution.out)
    PublishReadLengthDistribution(UpdateReadLengthDistribution.out)

    
    //Run transcript alignment, annotation and quantification
    minimap_transcript_output = MinimapTranscriptome(samples, file("${params.output_dir}/transcriptome_index.mmi"),conversion_output_channel,index_output_channel)
    salmon_output = Salmon(minimap_transcript_output.aligned_bams, file("${params.genome_gtf}"), file("${params.transcriptome_fasta}"))
    cleaned_salmon = CleanSalmonTable(salmon_output.single_salmon_df, file("${params.metadata}"))
    UpdateSalmonTable.scan(cleaned_salmon.clean_salmon)
    PublishSalmonTable(UpdateSalmonTable.out)
    UpdateIterator2(PublishSalmonTable.out)
    DESeq2Transcriptome(UpdateIterator2.out.run_statistics,PublishSalmonTable.out.merged_all,file("${params.metadata}"))
    DTUanalysis(DESeq2Transcriptome.out.dea_transcriptome_done,PublishSalmonTable.out.merged_all,file("${params.metadata}"), file("${params.output_dir}/converted_gtf.csv"))
    

    //def latest_bams_transcriptome = fetch_latest_bams("${params.output_dir}/bam_transcriptome_merged/")
    //found_transcriptome_bams = FindMergedTranscriptomeBams(minimap_transcript_output.aligned_bams)
    //def transcriptome_bams = fetch_latest_bams("${params.output_dir}/bam_transcriptome_merged/")
    MinimapTranscriptomeMergeBam(minimap_transcript_output.aligned_bams, channel.fromPath("${params.output_dir}/bam_transcriptome_merged/*.bam").collect())
    //Process run time
    // ProcessingTimeRegistration(file("${params.output_dir}/processing_time_table.csv"),minimap_genome_output.timestamp,UpdateIterator.out.timestamp,minimap_transcript_output.timestamp,UpdateIterator2.out.timestamp)

}
//messages to display once the workflow has completed

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
}


//Rscript script.r "${new_input}" "${state}"