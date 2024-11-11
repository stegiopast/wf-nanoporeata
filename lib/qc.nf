import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList
import java.io.File
import java.time.LocalDateTime


process RunDevelopmentEstimation{
    publishDir(
        path: "${params.out_dir}",
        mode: 'move'
    )
    stageInMode "copy"
    label "nanoporeata"
    maxForks 1
    cpus 4
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}

    input:
        val run_statistics
        path merged_csv, stageAs: 'merged_all.csv'
        path inner_variability_plot, stageAs: 'inner_var.csv'
        path inner_variability_per_sample, stageAs: 'inner_var_p_sample.csv'
        path exp_genes_counted_per_sample, stageAs: 'exp_counts_p_sample.csv'
        
    when:
        run_statistics == 1
    
    output:
        path("exp_genes_counted_per_sample.csv")
        path("inner_variability_plot.csv")
        path("inner_variability_per_sample.csv")
    
    script:
    """
    python ${projectDir}/bin/infer_experiment_absolute_gene_amount.py -s merged_all.csv -o exp_counts_p_sample.csv
    python ${projectDir}/bin/infer_experiment_inner_variability.py -s merged_all.csv -d inner_var.csv -o inner_var_p_sample.csv
    mv exp_counts_p_sample.csv exp_genes_counted_per_sample.csv
    mv inner_var.csv inner_variability_plot.csv
    mv inner_var_p_sample.csv inner_variability_per_sample.csv
    """

}


process CountMappedReads{ 
    publishDir(
        path: "${params.out_dir}/mapped_reads/",
        mode: 'copy'
    )
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}
    input:
    tuple val(ID), path(bam)
    path(bai)

    output:
    tuple val(ID), path("mapped_read_counts_${ID}_${task.index}.csv"), emit: read_counts

    script:
    """
    bash ${projectDir}/bin/get_number_of_mapped_reads.sh $bam $ID mapped_read_counts_${ID}_${task.index}.csv 
    """
}

process MergeMappedReadsTable{
    publishDir(
        path: "${params.out_dir}",
        mode: 'copy'
    )
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}

    input:
    tuple val(ID), path(count_csv)
    path(metadata)
    path(construct_dir), stageAs: "construct_dir.csv"
    output:
    path("mapping_stats.txt")

    script:
    if (task.index == 1)
    """
    python ${projectDir}/bin/initialize_read_counts.py -i $count_csv -m $metadata -o mapping_stats.txt
    """
    else
    """
    python ${projectDir}/bin/continue_read_counts.py -i $count_csv -m $metadata -c construct_dir.csv -o mapping_stats.txt
    """
}

process DefineReadLengthDistribution{
    maxForks 1
    input:
    tuple val(ID), path(fastq_file)
    path(metadata)
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}
    output:
    path("${ID}_read_lengths_${task.index}.csv")

    script:
    """
    if [[ "$fastq_file" == *".gz" ]]
    then 
        zcat "$fastq_file" | awk 'NR%4==2' | awk '{ print length }' > ${ID}_read_length_file${task.index}.txt
    else
        cat "$fastq_file" | awk 'NR%4==2' | awk '{ print length }' > ${ID}_read_length_file${task.index}.txt
    fi
    cat ${ID}_read_length_file${task.index}.txt
    python ${projectDir}/bin/initialize_read_length2.py -s ${ID}_read_length_file${task.index}.txt -n ${ID} -m ${metadata} -o ${ID}_read_lengths_${task.index}.csv
    """
}

process UpdateReadLengthDistribution{
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}

    input:
    path(output)
    
    output:
    path("final_read_lengths_${task.index}.csv")

    script:
    def new_table = output instanceof BlankSeparatedList ? output.first() : output
    def old_table = output instanceof BlankSeparatedList ? output.last() : "NOSTATE"
    """
    if [ $old_table == "NOSTATE" ]
    then 
        python ${projectDir}/bin/continue_read_length2.py -n $new_table -s $old_table
        mv merged_all_readlengths_temp.csv final_read_lengths_${task.index}.csv
        find . -type l -delete
    else
        python ${projectDir}/bin/continue_read_length2.py -n $new_table -s $old_table
        mv merged_all_readlengths_temp.csv final_read_lengths_${task.index}.csv
        find . -type l -delete
    fi
    """
}

process PublishReadLengthDistribution{
    publishDir(path: "${params.out_dir}/ReadLengthFolder",
               mode: "move"
    )
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}
    input:
    path(final_read_lengths)

    output:
    path("*_read_length_pass.txt", optional:true)

    script:
    """
    python ${projectDir}/bin/publish_read_length.py -i $final_read_lengths
    """

}



// process ProcessingTimeRegistration{
//     publishDir "${params.out_dir}"
//     input:
//     path(time_table), stageAs: "previous_table.csv"
//     val(genome_process_start)
//     val(genome_process_end)
//     val(transcriptome_process_start)
//     val(transcriptome_process_end)

//     output:
//     path("processing_time_table.csv")

//     script:
//     long elapsed_time_genome_pipeline = genome_process_start.until(genome_process_end,ChronoUnit.SECONDS)
//     long elapsed_time_transcriptome_pipeline = transcriptome_process_start.until(transcriptome_process_end,ChronoUnit.SECONDS)
//     if (task.index == 1)
//     """
//     echo Tool,Iteration,Time > processing_time_table.csv
//     echo genome_pipeline,$task.index,$elapsed_time_genome_pipeline >> processing_time_table.csv
//     echo transcriptome_pipeline,$task.index,$elapsed_time_transcriptome_pipeline >> processing_time_table.csv
//     """
//     else
//     """
//     echo genome_pipeline,$task.index,$elapsed_time_genome_pipeline >> previous_table.csv
//     echo transcriptome_pipeline,$task.index,$elapsed_time_transcriptome_pipeline >> previous_table.csv
//     cp previous_table > processing_time_table.csv
//     """
// }





