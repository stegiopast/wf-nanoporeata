/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                          Size checkup                                                                   //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process readlength_check{
    input:
    val string
    val string_array

    script:
    println "Running readlength analysis"
    """
    for i in ${string}
    do
        if [ ${params.barcoded} -eq 1 ]
        then
            basis=\$(basename \$(dirname \${i}))
        else
            basis=\$(basename \$(dirname \$(dirname \$(dirname \${i}))))
            echo \$basis >> ${params.output_dir}error_logs/readlength_counting_error.log
        fi 
        bash ${params.script_dir}get_read_length_from_fastq.sh \$i \$basis ${params.output_dir}ReadLengthFolder || echo ${params.output_dir}error_logs/readlength_counting_error.log
    done 
    """
}

process infer_experiment{
    memory "4 GB"
    input:
    val fc_string

    script:
    if (fc_string != "")
    """
    workflow-glue infer_experiment_absolute_gene_amount -s ${params.output_dir}merged_all.csv -m "" -o ${params.output_dir}exp_genes_counted_per_sample.csv || echo "${fc_string}" >> ${params.output_dir}error_logs/quantify_genes_detected.log
    workflow-glue infer_experiment_inner_variability -s ${params.output_dir}merged_all.csv -m "" -d ${params.output_dir}inner_variability_plot.csv  -o ${params.output_dir}inner_variability_per_sample.csv || echo "${fc_string}" >> ${params.output_dir}error_logs/quantify_inner_variablity_detected.log
    """
    else
    """
    echo empty
    """


}