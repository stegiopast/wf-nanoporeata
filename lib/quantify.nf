import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList
import java.io.File




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                 FeatureCount Annotation                                                                 //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
process FeatureCount {
    label "nanoporeata"
    maxForks 1
    input:
        tuple val(ID), path(bam)
        path genome_gtf
    output:
        tuple val(ID),path("${ID}.${task.index}.csv"), emit: single_fc_df 
 
    script:
    """
    samtools index ${bam}
	featureCounts -a "${genome_gtf}" -F 'GTF' -L -T ${params.threads} -o ${ID}.${task.index}.csv ${bam} -T ${params.threads}
    """
}


process CleanFeatureCountTable {
    label "nanoporeata"
    publishDir "${params.output_dir}", mode: "copy"
    maxForks 1
    input:
        tuple val(ID),path(single_fc)
        path metadata
    output:
        path("feature_counts_${ID}_${task.index}.csv"), emit: clean_fc
        path("feature_counts_${ID}_${task.index}.csv")
    script:
    """
    python ${projectDir}/bin/initialize_fc_merge.py --input $single_fc --metadata $metadata
    mv merged_all_temp.csv feature_counts_${ID}_${task.index}.csv
    """

}


process UpdateFeatureCountTable {
    label "nanoporeata"
    maxForks 1
    input: 
        path output
    output:
        path "feature_counts_latest_${task.index}.csv"
    script:
       def new_table = output instanceof BlankSeparatedList ? output.first() : output
       def old_table = output instanceof BlankSeparatedList ? output.last() : "NOSTATE"
    """
    if [ $old_table == "NOSTATE" ]
    then 
        cat $new_table > feature_counts_latest_${task.index}.csv
    else
        python ${projectDir}/bin/continue_fc_merge.py --new_input $new_table --state $old_table
        mv merged_all_temp.csv feature_counts_latest_${task.index}.csv
    fi
    """
}

process PublishFeatureCountTable {
    label "nanoporeata"
    publishDir "${params.output_dir}", mode:"copy"
    maxForks 1
    input:
        path merged_csv
    output:
        path("merged_all.csv")
    script:
    """
    cat $merged_csv > merged_all.csv
    """

}
 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                 Salmon Annotation Transcriptome                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process Salmon{
    label "nanoporeata"
    maxForks 1
    input:
        tuple val(ID), path(bam)
        path genome_gtf
        path transcriptome_fasta

    output:
        tuple val(ID),path("salmon.${ID}.${task.index}"), emit: single_salmon_df 
    script:
    """
    samtools index ${bam}
    salmon quant -t ${params.transcriptome_fasta} -l A -a ${bam} -o salmon.${ID}.${task.index} -p ${params.threads}
    """
}


process CleanSalmonTable{
    label "nanoporeata"
    publishDir "${params.output_dir}", mode: "copy"
    maxForks 1
    input:
        tuple val(ID),path(single_salmon)
        path metadata
    output:
        path("salmon_${ID}_${task.index}.csv"), emit: clean_salmon
        path("salmon_${ID}_${task.index}.csv")
    script:
    """
    python ${projectDir}/bin/initialize_salmon_merge.py --input ${single_salmon}/quant.sf --metadata ${metadata}
    mv merged_all_temp.csv salmon_${ID}_${task.index}.csv
    """
}



process UpdateSalmonTable{
    label "nanoporeata"
    maxForks 1
    input: 
        path output
    output:
        path "salmon_latest_${task.index}.csv"
    script:
       def new_table = output instanceof BlankSeparatedList ? output.first() : output
       def old_table = output instanceof BlankSeparatedList ? output.last() : "NOSTATE"
    """
    if [ $old_table == "NOSTATE" ]
    then 
        cat $new_table > salmon_latest_${task.index}.csv
    else
        python ${projectDir}/bin/continue_salmon_merge.py --new_input $new_table --state $old_table
        mv merged_all_temp.csv salmon_latest_${task.index}.csv
    fi
    """

}

process PublishSalmonTable{
    label "nanoporeata"
    publishDir "${params.output_dir}", mode:"copy"
    maxForks 1
    input:
        path merged_csv
    output:
        path("salmon_merged_absolute.csv")
    script:
    """
    cat $merged_csv > salmon_merged_absolute.csv
    """
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                 DESeq2 Quantification                                                                   //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process DESeq2Genome {
    label "R"
    publishDir (path: "${params.output_dir}", mode: 'copy')
    maxForks 1
    cpus 3
    input:
        path count_table
        path metadata
    output:
        path("DDS_genes.R")
    script:
    """
    Rscript ${projectDir}/bin/nf_dea_function.R $metadata "gene" $count_table $task.cpus DDS_genes.R 
    """
}