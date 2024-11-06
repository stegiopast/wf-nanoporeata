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
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}
    input:
        tuple val(ID), path(bam), val(full_bam)
        path genome_gtf
        path metadata
    output:
        //tuple val(ID),path("${ID}.${task.index}.csv"), emit: single_fc_df
        path("feature_counts_${ID}_${task.index}.csv"), emit: clean_fc 

    script:
    """
    samtools index ${bam}
    featureCounts -a "${genome_gtf}" -F 'GTF' -L -T ${params.threads} -o ${ID}.${task.index}.csv ${bam} -T ${task.cpus}
    python ${projectDir}/bin/initialize_fc_merge2.py --input ${ID}.${task.index}.csv --metadata $metadata
    mv merged_all_temp.csv feature_counts_${ID}_${task.index}.csv
    rm ${ID}.${task.index}.csv
    """
}




process UpdateFeatureCountTable {
    label "nanoporeata"
    maxForks 1
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}
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
        python ${projectDir}/bin/continue_fc_merge2.py --new_input $new_table --state $old_table
        mv merged_all_temp.csv feature_counts_latest_${task.index}.csv
        find . -type l -delete
    else
        python ${projectDir}/bin/continue_fc_merge2.py --new_input $new_table --state $old_table
        mv merged_all_temp.csv feature_counts_latest_${task.index}.csv
        find . -type l -delete
    fi
    """
}

process PublishFeatureCountTable {
    label "nanoporeata"
    publishDir "${params.out_dir}", mode:"move"
    maxForks 1
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}
    input:
        path merged_csv
    output:
        val("${params.out_dir}/merged_all.csv"), emit: merged_all
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
    cpus 4
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}
    input:
        tuple val(ID), path(bam), val(full_bam)
        path genome_gtf
        path transcriptome_fasta
        path metadata
    output:
        //tuple val(ID),path("salmon.${ID}.${task.index}"), emit: single_salmon_df 
        path("salmon_${ID}_${task.index}.csv"), emit: clean_salmon
    script:
    """
    samtools index ${bam}
    salmon quant -t ${params.transcriptome_fasta} -l A -a ${bam} -o salmon.${ID}.${task.index} -p ${task.cpus} --minAssignedFrags 0
    python ${projectDir}/bin/initialize_salmon_merge2.py --input salmon.${ID}.${task.index}/quant.sf --metadata ${metadata}
    mv merged_all_temp.csv salmon_${ID}_${task.index}.csv
    rm salmon.${ID}.${task.index} -r
    """
}


process UpdateSalmonTable{
    label "nanoporeata"
    maxForks 1
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}
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
        python ${projectDir}/bin/continue_salmon_merge2.py --new_input $new_table --state $old_table
        mv merged_all_temp.csv salmon_latest_${task.index}.csv
        find . -type l -delete
    else
        python ${projectDir}/bin/continue_salmon_merge2.py --new_input $new_table --state $old_table
        mv merged_all_temp.csv salmon_latest_${task.index}.csv
        find . -type l -delete
    fi
    """
}

process PublishSalmonTable{
    label "nanoporeata"
    publishDir "${params.out_dir}", mode:"move"
    maxForks 1
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}
    input:
        path merged_csv
    output:
        val("${params.out_dir}/salmon_merged_absolute.csv"), emit: merged_all
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
    publishDir (path: "${params.out_dir}", mode: 'copy')
    maxForks 1
    cpus 8
    maxRetries 10
    memory "10GB"
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}
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
    publishDir (path: "${params.out_dir}", mode: 'copy')
    maxForks 1
    cpus 4
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}
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
    publishDir (path: "${params.out_dir}", mode: 'copy')
    maxForks 1
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}
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
    publishDir(path: "${params.out_dir}", mode: "copy")
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}
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
    publishDir(path: "${params.out_dir}", mode: "copy")
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}
    input:
    path salmon_count_merged

    output:
    val(run_statistics), emit: run_statistics

    script:
    if ((task.index % params.batchsize) == 0){
        run_statistics = 1
    }
    else{
        run_statistics = 0
    }
    """
    echo "Done"
    """
}
