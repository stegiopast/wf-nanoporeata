#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList
import java.time.LocalDateTime
nextflow.enable.dsl = 2
nextflow.preview.recursion=true 

include { ConvertGtfToDf; CreateFeaturePercentiles; start_watch_path ; fetch_latest_bams ; CreateGenomeBamFilesForMerge; CreateTranscriptomeBamFilesForMerge; CopyBedfileAndMetadata } from "./lib/initialize.nf"
include { MinimapIndex ; MinimapGenome ; MinimapTranscriptome ; MinimapGenomeMergeBam ; MinimapTranscriptomeMergeBam} from "./lib/alignment.nf"
include { FeatureCount; UpdateFeatureCountTable; PublishFeatureCountTable; DESeq2Genome ; DESeq2Transcriptome; DTUanalysis; UpdateIterator ; UpdateIterator2 } from "./lib/quantify.nf"
include { Salmon; UpdateSalmonTable; PublishSalmonTable } from "./lib/quantify.nf"
include { RunDevelopmentEstimation; CountMappedReads; MergeMappedReadsTable; DefineReadLengthDistribution; UpdateReadLengthDistribution; PublishReadLengthDistribution } from "./lib/qc.nf"


workflow {
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
    minimap_genome_output= MinimapGenome(samples, file("${params.out_dir}/genome_index.mmi"),file("${params.transcriptome_fasta}"), conversion_output_channel, index_output_channel)
    fc_output = FeatureCount(minimap_genome_output.aligned_bams, file("${params.genome_gtf}"), file("${params.metadata}"))
    UpdateFeatureCountTable.scan(fc_output.clean_fc)
    PublishFeatureCountTable(UpdateFeatureCountTable.out)
    UpdateIterator(PublishFeatureCountTable.out.merged_all)
    RunDevelopmentEstimation(UpdateIterator.out.run_statistics, PublishFeatureCountTable.out.merged_all, file("${params.out_dir}/inner_variability_plot.csv"),file("${params.out_dir}/inner_variability_per_sample.csv"),file("${params.out_dir}/exp_genes_counted_per_sample.csv"))
    DESeq2Genome(UpdateIterator.out.run_statistics,PublishFeatureCountTable.out.merged_all,file("${params.metadata}"))

    //Merge single genome aligned bams 
    MinimapGenomeMergeBam(minimap_genome_output.aligned_bams)//, Channel.fromPath("${params.out_dir}/bam_genome_merged/*"))

    //Run read mapping quantification
    counted_reads = CountMappedReads(minimap_genome_output.aligned_bams, minimap_genome_output.aligned_bam_bais)
    MergeMappedReadsTable(counted_reads.read_counts,file("${params.metadata}"),file("${params.out_dir}/mapping_stats.txt"))

    //Read size distribution
    DefineReadLengthDistribution(samples,file("${params.metadata}"))
    UpdateReadLengthDistribution.scan(DefineReadLengthDistribution.out)
    PublishReadLengthDistribution(UpdateReadLengthDistribution.out)

    
    //Run transcript alignment, annotation and quantification
    minimap_transcript_output = MinimapTranscriptome(samples, file("${params.out_dir}/transcriptome_index.mmi"),file("${params.transcriptome_fasta}"),conversion_output_channel,index_output_channel)
    salmon_output = Salmon(minimap_transcript_output.aligned_bams, file("${params.genome_gtf}"), file("${params.transcriptome_fasta}"), file("${params.metadata}"))

    UpdateSalmonTable.scan(salmon_output.clean_salmon)
    PublishSalmonTable(UpdateSalmonTable.out)
    UpdateIterator2(PublishSalmonTable.out.merged_all)
    DESeq2Transcriptome(UpdateIterator2.out.run_statistics,PublishSalmonTable.out.merged_all,file("${params.metadata}"))
    DTUanalysis(DESeq2Transcriptome.out.dea_transcriptome_done,PublishSalmonTable.out.merged_all,file("${params.metadata}"), file("${params.out_dir}/converted_gtf.csv"))
    

    //Merge single transcriptome aligned bams
    MinimapTranscriptomeMergeBam(minimap_transcript_output.aligned_bams)
    
    //Process run time
    // ProcessingTimeRegistration(file("${params.out_dir}/processing_time_table.csv"),minimap_genome_output.timestamp,UpdateIterator.out.timestamp,minimap_transcript_output.timestamp,UpdateIterator2.out.timestamp)

}

//messages to display once the workflow has completed
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
}