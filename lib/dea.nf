process DEA {
    input:
        file metadata
        file gene_or_trans
        file count_table
        file output_file
    output:
        file output_path, emit: output_path
    script:
        """
        Rscript ${metadata}/workflow_glue/R_scripts/nf_dea_function.R \
            ${gene_or_trans} \
            ${count_table} \
            ${params.threads} \
            ${output_file}
        """
}