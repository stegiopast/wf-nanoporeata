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




// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// //                                                                                                                                                         // 
// //                                                                                                                                                         //
// //                                                                                                                                                         //
// //                                                                          Size checkup                                                                   //
// //                                                                                                                                                         //
// //                                                                                                                                                         //
// //                                                                                                                                                         //
// //                                                                                                                                                         //
// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// process readlength_check{
//     input:
//     val string
//     val string_array

//     script:
//     println "Running readlength analysis"
//     """
//     for i in ${string}
//     do
//         if [ ${params.barcoded} -eq 1 ]
//         then
//             basis=\$(basename \$(dirname \${i}))
//         else
//             basis=\$(basename \$(dirname \$(dirname \$(dirname \${i}))))
//             echo \$basis >> ${params.output_dir}error_logs/readlength_counting_error.log
//         fi 
//         bash ${params.script_dir}get_read_length_from_fastq.sh \$i \$basis ${params.output_dir}ReadLengthFolder || echo ${params.output_dir}error_logs/readlength_counting_error.log
//     done 
//     """
// }

// process infer_experiment{
//     memory "4 GB"
//     input:
//     val fc_string

//     script:
//     if (fc_string != "")
//     """
//     workflow-glue infer_experiment_absolute_gene_amount -s ${params.output_dir}merged_all.csv -m "" -o ${params.output_dir}exp_genes_counted_per_sample.csv || echo "${fc_string}" >> ${params.output_dir}error_logs/quantify_genes_detected.log
//     workflow-glue infer_experiment_inner_variability -s ${params.output_dir}merged_all.csv -m "" -d ${params.output_dir}inner_variability_plot.csv  -o ${params.output_dir}inner_variability_per_sample.csv || echo "${fc_string}" >> ${params.output_dir}error_logs/quantify_inner_variablity_detected.log
//     """
//     else
//     """
//     echo empty
//     """


// }