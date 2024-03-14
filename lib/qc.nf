import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList
import java.io.File


process RunDevelopmentEstimation{
    label "nanoporeata"
    publishDir "${params.output_dir}"
    maxForks 1
    cpus 4

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
    publishDir "${params.output_dir}/mapped_reads/"
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
    publishDir "${params.output_dir}"
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

process UpdateReadLengthDistribution{
    publishDir "${params.output_dir}/ReadLengthFolder/"
    input:
    tuple val(ID), path(fastq_file)

    output:
    path("${ID}_read_lengths_pass.txt")

    script:
    def new_fastq = fastq_file instanceof BlankSeparatedList ? fastq_file.first() : fastq_file
    def old_fastq = fastq_file instanceof BlankSeparatedList ? fastq_file.last() : "NOSTATE"
    """
    if [ $old_fastq == "NOSTATE" ]
    then 
        echo Length > ${ID}_read_lengths_pass.txt
        bash ${projectDir}/bin/determine_read_length.sh $new_fastq ${ID}_read_lengths_pass.txt
    else
        cat $old_fastq > ${ID}_read_lengths_pass.txt
        bash ${projectDir}/bin/determine_read_length.sh $new_fastq ${ID}_read_lengths_pass.txt
    fi
    """
}


