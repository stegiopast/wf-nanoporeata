import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList
import java.io.File
import java.time.LocalDateTime




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
    cpus 4
    input:
        tuple val(ID), path(bam)
        path genome_gtf
    output:
        tuple val(ID),path("${ID}.${task.index}.csv"), emit: single_fc_df 
 
    script:
    """
    samtools index ${bam}
	featureCounts -a "${genome_gtf}" -F 'GTF' -L -T ${params.threads} -o ${ID}.${task.index}.csv ${bam} -T ${task.cpus}
    """
}


process CleanFeatureCountTable {
    label "nanoporeata"
    //publishDir "${params.output_dir}", mode: "copy"
    maxForks 1
    input:
        tuple val(ID),path(single_fc)
        path metadata
    output:
        path("feature_counts_${ID}_${task.index}.csv"), emit: clean_fc
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
        path("merged_all.csv"), emit: merged_all
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
    cpus 4
    input:
        tuple val(ID), path(bam)
        path genome_gtf
        path transcriptome_fasta

    output:
        tuple val(ID),path("salmon.${ID}.${task.index}"), emit: single_salmon_df 
    script:
    """
    samtools index ${bam}
    salmon quant -t ${params.transcriptome_fasta} -l A -a ${bam} -o salmon.${ID}.${task.index} -p ${task.cpus} --minAssignedFrags 0
    """
}


process CleanSalmonTable{
    label "nanoporeata"
    //publishDir "${params.output_dir}", mode: "copy"
    maxForks 1
    input:
        tuple val(ID),path(single_salmon)
        path metadata
    output:
        path("salmon_${ID}_${task.index}.csv"), emit: clean_salmon
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
        path "salmon_latest_${task.index}.csv", emit: merged_all
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
        path("salmon_merged_absolute.csv"), emit: merged_all
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
    cpus 4
    maxRetries 10
    input:
        val run_statistics
        path count_table
        path metadata

    when:
        run_statistics == 1 
    output:
        path '*gene.RData'
        val 1, emit: dea_genome_done
    script:
    """
    Rscript ${projectDir}/bin/dea_nextflow.R $metadata "gene" $count_table $task.cpus gene.RData 
    """
}


process DESeq2Transcriptome{
    label "R"
    publishDir (path: "${params.output_dir}", mode: 'copy')
    maxForks 1
    cpus 4
    maxRetries 10
    input:
        val run_statistics
        path count_table
        path metadata
    when:
        run_statistics == 1
    output:
        path('*transcript.RData')
        val 1, emit: dea_transcriptome_done

    script:
    """
    Rscript ${projectDir}/bin/dea_nextflow.R $metadata "transcript" $count_table $task.cpus transcript.RData 
    """
}


process DTUanalysis{
    label "R"
    publishDir (path: "${params.output_dir}", mode: 'copy')
    maxForks 1
    maxRetries 10 
    input:
        val transcriptome_done
        path count_table
        path metadata
        path gtf_file

    output:
        path('*transcript_usage.RData')
        

    script:
    """
    Rscript ${projectDir}/bin/dtu_nextflow.R $metadata $count_table $gtf_file transcript_usage.RData
    """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                 Update iterator                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


process UpdateIterator{
    publishDir "${params.output_dir}"
    input:
    path feature_count_merged

    output:
    val(run_statistics), emit: run_statistics
    path("*")

    script:
    if ((task.index % params.batchsize) == 0){
        run_statistics = 1
    }
    else{
        run_statistics = 0
    }
    """
    if [ ${task.index} -eq 1 ]
    then
        echo "" > inner_variability_plot.csv
        echo "" > inner_variability_per_sample.csv
        echo "" > exp_genes_counted_per_sample.csv
    else
        echo "Initialization has been performed successfully" > initialization_done.log
    fi
    """
}


process UpdateIterator2{
    publishDir "${params.output_dir}"
    input:
    path salmon_count_merged

    output:
    val(run_statistics), emit: run_statistics

    exec:
    if ((task.index % params.batchsize) == 0){
        run_statistics = 1
    }
    else{
        run_statistics = 0
    }
}