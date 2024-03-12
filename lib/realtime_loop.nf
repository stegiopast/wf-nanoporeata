import groovy.json.JsonBuilder
@Grab('com.xlson.groovycsv:groovycsv:1.1')
import static com.xlson.groovycsv.CsvParser.parseCsv
import java.io.File;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                 Starting the analysis loop                                                              //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

process activate_channels{
     input:
     val done
     
     exec:
     //Activation of first processes by adding a value into a channel
     checkup << Channel.from(1)
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                        Detection of data in paths of interest                                                           //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//Process searches for all files in each folder of interest respectively 

process fromPath_data{
    memory "4 GB"
    input:
    val check

    output:
    val files, emit:files


    exec:
    println "Extracting files from Path"
    files = []
    if (params.barcoded == 1){    
        try{
            files << Channel.fromPath("${params.seq_data_folder}*/*/fastq_pass/*/*", type: "file").collect()
        }
        catch(Exception e){}
    }
    else{
        try{ 
            files << Channel.fromPath("${params.seq_data_folder}*/*/fastq_pass/*", type: "file").collect()
        }
        catch(Exception e){}  
    }
    done << Channel.from("yes")
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                       Preprocessing                                                                     //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Process compares value of data_seen_list with the list of all files that are present in the folders of 
//interest at each iteration and extracts new files. An array of new files and a concatenated string 
//with format "file1 data2 data3" is transported to the next processes.
//In case of an interruption of the pipeline and a subsequent corruption of the files 
//this process is also running a data cleanup

process data_alignment_prep{
        memory "4 GB"
        input:
        val files
        
        output:
        val string, emit: string 
        val string_array, emit:string_array
        val data_string, emit: data_string
        val data_seen_string, emit: data_seen_string
    
        script:
        string = ""
        string_array = []
        println "Preprocessing \n"
        println "Data seen before update: ${params.data_seen_list.size()} \n"
        println ""
        data_to_process_list = []
        if (files.val.get(0).getClass() == nextflow.util.ArrayBag){
        for(j in 0..files.val.size()-1){
            for(k in 0..files.val.get(j).size()-1){
                if(files.val.get(j).get(k).toString() in params.data_seen_list){
                }
                else if(files.val.get(j).get(k) in params.data_seen_list){
                }
                else{
                    data_to_process_list.add(files.val.get(j).get(k))
                }
            }
        }
        
        temp_list = []
        if (data_to_process_list.size > 1){
        for (q in 0..data_to_process_list.size()-1){
            if (params.barcoded == 1){    
               if (data_to_process_list.get(q).getParent().name in params.sample_names){
                   temp_list.add(data_to_process_list.get(q))
               }
            }
            else{
               if (data_to_process_list.get(q).getParent().getParent().getParent().name in params.sample_names){
                   temp_list.add(data_to_process_list.get(q))
               }
            }
        }
        data_to_process_list = temp_list
        }
        }
        if (suffix_determined.value == 0){
            if(data_to_process_list.size > 0){
                first_file = data_to_process_list.get(0).toString()
                first_file_parts = first_file.split('\\.').collect()
                last_part = first_file_parts.get(first_file_parts.size()-1)
                println "Last part: "
                println last_part
                println ""
                if (last_part == "gz"){
                    params.suffix = "fastq.gz"
                }
                else {
                    params.suffix = "fastq"
                }
                suffix_determined.value = 1
                println "Params suffix: "
                println params.suffix
                println ""
            }
        }
        println "Data to process: ${data_to_process_list.size()}"
        x_seen = []
        for(k in 0..Math.min(params.batchsize,data_to_process_list.size()-1)){
            if (data_to_process_list.size() > 1){
                x = Math.abs(new Random().nextInt() % (data_to_process_list.size()))
                //While loop is integrated to avoid double insertion of one file 
                while(x in x_seen){
                    x = Math.abs(new Random().nextInt() % (data_to_process_list.size()))
                }
                x_seen.add(x)
                params.data_seen_list.add(data_to_process_list.get(x))
                string = string + data_to_process_list.get(x) + " "
                string_array.add(data_to_process_list.get(x))
            }
            else if (data_to_process_list.size() == 1){
                params.data_seen_list.add(data_to_process_list.get(0))
                string = string + data_to_process_list.get(0) + " "
                string_array.add(data_to_process_list.get(0))   
            }
            else {
                string = ""
                string_array = []
            }
        }
        if (string != ""){
            iteration.value = iteration.value + 1
        }
        println ""
        println "New data: ${string_array.size()} \n"
        println "Data seen: ${params.data_seen_list.size()}"
        data_seen_string = ""
        data_string = ""  
        File file_tmp = new File("${params.output_dir}data_seen_tmp_processed.txt")
        file_tmp.write("")
        for(i in 0..params.data_seen_list.size()-1){
            if (params.data_seen_list.size() > 0){
                if (params.data_seen_list.get(i).toString() != ""){
                    data_seen_string = data_seen_string + params.data_seen_list.get(i).toString() + " "
                    file_tmp.append(params.data_seen_list.get(i).toString() + "\n")
                }
            }
        }
        for(i in 0..params.sample_names.size()-1){
        data_string = data_string + params.sample_names.get(i) + " "
        }
        """
        corruption_gene_bam=0
        corruption_transcript_bam=0
        bam_files_seen=" "
        bam_files_seen_transcripts=" "
        for ((i=1; i<=\$(cat ${params.output_dir}data_seen_tmp_processed.txt | wc -l); i++))
        do
            I=\$(sed -n \$((\$i))p ${params.output_dir}data_seen_tmp_processed.txt)
            if [ ${params.barcoded} -eq 1 ]
            then
                basis=\$(basename \$(dirname \${I}))
                filename=\$(basename \${I})
                converted_filename=\${filename/${params.suffix}/bam}
                outname=${params.output_dir}\${basis}/bam_files/\${converted_filename}
                outname_transcripts=${params.output_dir}\${basis}/bam_files_transcripts/\${converted_filename}
            else
                basis=\$(basename \$(dirname \$(dirname \$(dirname \${I}))))
                filename=\$(basename \${I})
                converted_filename=\${filename/${params.suffix}/bam}
                outname=${params.output_dir}\${basis}/bam_files/\${converted_filename}
                outname_transcripts=${params.output_dir}\${basis}/bam_files_transcripts/\${converted_filename}
            fi
            bam_files_seen=\${bam_files_seen}\${outname}" "
            bam_files_seen_transcripts=\${bam_files_seen_transcripts}\${outname_transcripts}" "
        done
        for i in ${data_string}
        do
            files_in_i=\$(ls ${params.output_dir}\${i}/bam_files/)
            files_in_i_transcripts=\$(ls ${params.output_dir}\${i}/bam_files_transcripts/)
            for j in \$files_in_i
                do
                    if [[ ! "\$bam_files_seen" =~ "${params.output_dir}\${i}/bam_files/\${j}" ]]
                    then
                        if [[ ! -d ${params.output_dir}\${i}/bam_files/\${j} ]]
                        then
                            rm ${params.output_dir}\${i}/bam_files/\${j} && echo ${params.output_dir}\${i}/bam_files/\${j} >> ${params.output_dir}error_logs/removed_bam_files.txt
                            corruption_gene_bam=1
                            echo "1" > ${params.output_dir}corruption_cleanup.txt
                        fi
                    fi
                done
            for j in \$files_in_i_transcripts
                do
                    if [[ ! "\$bam_files_seen_transcripts" =~ "${params.output_dir}\${i}/bam_files_transcripts/\${j}" ]]
                    then
                        if [[ ! -d ${params.output_dir}\${i}/bam_files_transcripts/\${j} ]]
                        then
                            rm ${params.output_dir}\${i}/bam_files_transcripts/\${j} && echo ${params.output_dir}\${i}/bam_files_transcripts/\${j} >> ${params.output_dir}error_logs/removed_bam_files_transcripts.txt
                            corruption_transcript_bam=1
                        fi
                    fi
                done
        done
        echo \$corruption_gene_bam >> ${params.output_dir}error_logs/corruption_gene.log
        function samtoolsParallelAll  {
            if [ ! \$(ls ${params.output_dir}\$1/bam_files/*.bam | wc -l) -eq 0 ]
            then
                par_basis=\$1
                run_dir=\$2
                sample=\${run_dir}\${par_basis}/bam_files/
                echo "Running samtoolsParallelAll" >> ${params.output_dir}error_logs/samtoolsParallelall_running.log
                sample_transcripts=\${run_dir}\${par_basis}/bam_files_transcripts/
                samtools merge \${run_dir}\${par_basis}/all_gene.bam \${sample}*.bam -f --threads ${params.threads} -c -p || echo "\${sample}" >> \${run_dir}error_logs/bam_file_genome_merge_after_corruption_error.log
                cp \${run_dir}\${par_basis}/all_gene.bam \${run_dir}bam_genome_merged/\${par_basis}.bam
                samtools index \${run_dir}bam_genome_merged/\${par_basis}.bam
                samtools merge \${run_dir}\${par_basis}/salmon/all.bam \${sample_transcripts}*.bam -f --threads ${params.threads} -c -p || echo "\${sample_transcripts}" >> \${run_dir}error_logs/bam_file_transcriptome_merge_after_corruption_error.log
                cp \${run_dir}\${par_basis}/salmon/all.bam \${run_dir}bam_transcriptome_merged/\${par_basis}.bam
                samtools index \${run_dir}bam_transcriptome_merged/\${par_basis}.bam
            fi
        }
        export -f samtoolsParallelAll
        function fastqmergeParallel {
        par_basis=\$1
        run_dir=\$2
        data_to_process=""
        for ((i=1; i<=\$(cat ${params.output_dir}data_seen_tmp_processed.txt | wc -l); i++))
        do
            I=\$(sed -n \$((\$i))p ${params.output_dir}data_seen_tmp_processed.txt)
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
                zcat \${data_to_process} | gzip > \${run_dir}\${par_basis}/merged_fastq/\${par_basis}.fastq.gz
                fi
            if [ ! \$(ls ${params.seq_data_folder}\${par_basis}/*/fastq_pass/*.fastq | wc -l) -eq 0 ]
                then
                cat \${data_to_process} > \${run_dir}\${par_basis}/merged_fastq/\${par_basis}.fastq
                fi
        else
            if [ ! \$(ls ${params.seq_data_folder}*/*/fastq_pass/\${par_basis}/*.fastq.gz | wc -l) -eq 0 ]
            then
            zcat \${data_to_process} | gzip > \${run_dir}\${par_basis}/merged_fastq/\${par_basis}.fastq.gz
            fi
            if [ ! \$(ls ${params.seq_data_folder}*/*/fastq_pass/\${par_basis}/*.fastq | wc -l) -eq 0 ]
            then
            cat \${data_to_process} > \${run_dir}\${par_basis}/merged_fastq/\${par_basis}.fastq
            fi
        fi      
        }
        export -f fastqmergeParallel 
       

        if [ -f ${params.output_dir}corruption_cleanup.txt ]
        then
            cleanup=\$(cat ${params.output_dir}corruption_cleanup.txt)
            if [ \$cleanup -eq 1 ]
            then 
            parallel -v -u --env samtoolsParallelAll --no-notice -j ${params.threads} samtoolsParallelAll ::: ${data_string} ::: ${params.output_dir}
            parallel -v -u --env fastqmergeParallel --no-notice -j ${params.threads} fastqmergeParallel ::: ${data_string} ::: ${params.output_dir}
            fi
        fi
        if [ ! -f ${params.output_dir}processing_time_table.csv ]
            then
                echo "Tool,Iteration,Time" >> ${params.output_dir}processing_time_table.csv
            fi
        iteration=${iteration.value}
        echo \$((\${iteration})) > ${params.output_dir}/current_iteration.txt
        """
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                      Reawake next iteration of data analysis                                                            //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


process reawake_next_round{
    input:
    val fc, emit: fc
    val salmon, emit: salmon
    val fc_string, emit:fc_string
    val string_array, emit: string_array 
    val data_string, emit: data_string

    exec:
    if (fc_string != ""){
    println "Feature Count merge done for current iteration"
    File file = new File("${params.output_dir}data_seen.txt")
    if (string_array.size() > 0){
        for (i in 0..string_array.size()-1){
            lines = file.readLines()
            file.write("${lines.get(0)}${string_array.get(i)} ")
        }
      }
    }
    if (fc_string == ""){
    sleep(30000)
    }
    def continue_text = new File("${params.output_dir}process_running.txt").getText()
    println "Continue ?"
    println continue_text
    println "____________"
    continue_bool.value = continue_text.toInteger()
    if (continue_bool.value == 1){
        checkup << Channel.from(1)
    }
    else {
        def continue_text2 = new File("${params.output_dir}process_running.txt")
        continue_text2.write "2"
        while(continue_bool.value != 1){
            def continue_text3 = new File("${params.output_dir}process_running.txt").getText()
            continue_bool.value = continue_text3.toInteger()
            if ( continue_bool.value == 1){
                println "Going to continue"
            }
            else {
                sleep(2000)
            }
        }
        checkup << Channel.from(1)
    }
}