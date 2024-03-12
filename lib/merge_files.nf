import groovy.json.JsonBuilder
@Grab('com.xlson.groovycsv:groovycsv:1.1')
import static com.xlson.groovycsv.CsvParser.parseCsv
import java.io.File;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                          Merge fastq                                                                    //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


process merge_fastq{
    input:
    val data_string
    val data_seen_string
    val string
    output:
    val data_string, emit: data_string

    script:
    if (string != "")
    """
    function fastqmergeParallel {
        par_basis=\$1
        run_dir=\$2
        data_to_process=""
        for i in ${string}
        do
            I=\$i
            if [ ${params.barcoded} -eq 1 ]
            then
            basis=\$(basename \$(dirname \$I))
            else
            basis=\$(basename \$(dirname \$(dirname \$(dirname \$I))))
            fi
            if [ \${basis} == \${par_basis} ]
            then
                data_to_process=\$data_to_process\$I" "
            fi
        done
        if [ ${params.barcoded} -eq 0 ]
        then
            if [ ! \$(ls ${params.seq_data_folder}\${par_basis}/*/fastq_pass/*.fastq.gz | wc -l) -eq 0 ]
                then
                zcat \${data_to_process} | gzip -c >> \${run_dir}\${par_basis}/merged_fastq/\${par_basis}.fastq.gz
                fi
            if [ ! \$(ls ${params.seq_data_folder}\${par_basis}/*/fastq_pass/*.fastq | wc -l) -eq 0 ]
                then
                cat \${data_to_process} >> \${run_dir}\${par_basis}/merged_fastq/\${par_basis}.fastq
                fi
        else
            if [ ! \$(ls ${params.seq_data_folder}*/*/fastq_pass/\${par_basis}/*.fastq.gz | wc -l) -eq 0 ]
            then
            zcat \${data_to_process} | gzip -c >> \${run_dir}\${par_basis}/merged_fastq/\${par_basis}.fastq.gz
            fi
            if [ ! \$(ls ${params.seq_data_folder}*/*/fastq_pass/\${par_basis}/*.fastq | wc -l) -eq 0 ]
            then
            cat \${data_to_process} >> \${run_dir}\${par_basis}/merged_fastq/\${par_basis}.fastq
            fi
        fi      
    }
    export -f fastqmergeParallel 
    start=\$(date +%s)
    if [ -f ${params.output_dir}corruption_cleanup.txt ]
    then
        cleanup=\$(cat ${params.output_dir}corruption_cleanup.txt)
        if [ \$cleanup -eq 1 ]
        then
            echo "0" > ${params.output_dir}corruption_cleanup.txt
        else
            parallel -v -u --env fastqmergeParallel --no-notice -j ${params.threads} fastqmergeParallel ::: ${data_string} ::: ${params.output_dir}
        fi
    else
        parallel -v -u --env fastqmergeParallel --no-notice -j ${params.threads} fastqmergeParallel ::: ${data_string} ::: ${params.output_dir}
    fi
    end=\$(date +%s)
    time=\$(echo "\$((\$end-\$start))")
    echo "fastq_merge,${iteration.value},\$time" >> ${params.output_dir}processing_time_table.csv
    """
    else
    """
    echo "No new data"
    """
}


process move_transcript_files{
    memory "4 GB"
    input: 
    val bam_string 
	val string_array 

    output:
    val bam_full_string, emit: bam_full_string
    val string_array, emit: string_array

    script:
    bam_full_string = ""
    for(i in 0..string_array.size()-1){
       if (string_array.size() > 0){
           element = string_array.get(i)
           if (params.barcoded == 1){
           basis = element.getParent().getName()
           }
           else {
               basis = element.getParent().getParent().getParent().getName()
           }
           filename= element.getName()

           converted_filename =  filename.minus(".${params.suffix}")
           outname= params.output_dir + basis + "/bam_files_transcripts/" + "full/" + converted_filename
           bam_full_string=bam_full_string + outname + ".bam" + " "
       }
    }
    move_transcripts_minimap2_running.value = 1
    if (bam_string != "")
        """
        for i in ${bam_string}
        do
            basis=\$(dirname \${i})
            filename=\$(basename \${i})
            cp \${i} \$basis/full/\$filename || echo "File does not exist: \$basis/full/\$filename" > ${params.output_dir}error_logs/copy_full_bam_failed.log
        done
        """
    else
        """
        echo empty
        """
}

process merge_table_of_all_folders{
    memory "4 GB"
    input:
    val fc_string
    val string_array
    
    output:
    val fc_string, emit: fc_string
    val string_array, emit: string_array 
 

    script:
    println "Feature count merging all files"
    
    //While loop is representing a threadlock situation
    if (fc_string != ""){
        data_string = ""
        for (i in 0..params.sample_names.size()-1){
	    data_string = data_string + params.output_dir + params.sample_names.get(i) + "/" + " "
        }
    }
    if (fc_string != "") 
        """
        workflow-glue merge_all_fc ${data_string} ${params.output_dir}merged_all.csv || echo "${fc_string}" >> ${params.output_dir}error_logs/feature_counts_merge_all_error.log
        """
    else
        """
        echo ""
        """

}


process merge_salmon_annotation{
    memory "4 GB"
    input:
    val string
    val string_array 
    
    output:
    val done, emit: done

    script:
    println "Salmon merging all files"  
    data_string = ""
    for (i in 0..params.sample_names.size()-1){
	data_string = data_string + params.output_dir + params.sample_names.get(i) + "/" + " "
    }
    if (string != "") 
        """
        workflow-glue merge_all_salmon ${data_string} ${params.output_dir}salmon_merged_tpm.csv ${params.output_dir}salmon_merged_absolute.csv || echo "${string}" >> ${params.output_dir}error_logs/salmon_merge_error.log
        """
    else
        """
        echo ""
        """
}


