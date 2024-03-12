process DTE {
    input:
        file metadata_path
        file counts_path
        file gtf_path
        file output_path
    output:
        file output_path, emit: output_path
    script:
        """
        Rscript ${projectDir}/workflow_glue/R_scripts/dte_nextflow.R \
            ${metadata_path} \
            ${counts_path} \
            ${gft_path} \
            ${output_path}
        """
}