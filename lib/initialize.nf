import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList
import java.io.File

process ConvertGtfToDf{ 
    publishDir "${params.output_dir}/"
    maxForks 3
    memory "8GB"
    maxRetries 10
    input:
    path(genome_gtf)

    output:
    path("./converted_gtf.csv")
    val 1, emit: done

    script:
    """
    python ${projectDir}/bin/convert_gtf_to_df.py -i ${genome_gtf} -o converted_gtf.csv
    """
}


def start_watch_path(){
    if (params.barcoded == 1){
        def ch_existing_input = Channel.fromPath("${params.seq_data_folder}/*/*/fastq_pass/*/*.fast*", type: "file") 
        | map { tuple(it.parent.name, it ) }
        | randomSample( 100000000 , 100 )

        def ch_watched = Channel.watchPath("${params.seq_data_folder}/*/*/fastq_pass/*/*.fast*", 'create,modify')      
        | until { file->file.name == 'STOP.fastq.gz' }
        | map { tuple(it.parent.name, it ) }
        
        ch_watched_final = ch_existing_input | concat(ch_watched)  
    }
    else{
        def ch_existing_input = Channel.fromPath("${params.seq_data_folder}/*/*/fastq_pass/*.fast*", type: "file") 
        //| randomSample( fromPath_counts )
        | map { tuple(it.parent.parent.parent.name, it ) }
        | randomSample( 100000000 , 100 )

        def ch_watched = Channel.watchPath("${params.seq_data_folder}/*/*/fastq_pass/*.fast*", 'create,modify')      
        | until { file->file.name == 'STOP.fastq.gz' }
        | map { tuple(it.parent.parent.parent.name, it ) }

        ch_watched_final = ch_existing_input | concat(ch_watched)     
    }
    return ch_watched_final
}

process CreateFeaturePercentiles{
    publishDir "${params.output_dir}/"
    maxForks 3
    memory "8GB"
    input:
        path genome_bed
    output:
        path "./g_percentiles.json", emit: feature_percentile_file
    script:
    """
    python ${projectDir}/bin/create_feature_percentiles.py --bed ${genome_bed} --output_dir ./
    """
}


// process find_data_folders(){
//     input: 

//     output:
//         val done, emit: done
//     exec:
//     if (params.barcoded == 1){
//         def text = new File(params.metadata).getText()
//         def data = parseCsv(text, separator: '\t')
//         data.each{line -> params.sample_names.add(line["Samples"])}
//     }
//     else {
//         def text = new File(params.metadata).getText()
//         def data = parseCsv(text, separator: '\t')
//         data.each{line -> params.sample_names.add(line["Samples"])}
//     }
//     done = 1
// }


// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// //                                                                                                                                                         // 
// //                                                                                                                                                         //
// //                                                                                                                                                         //
// //                                                            Make crucial Directories                                                                     //
// //                                                                                                                                                         //
// //                                                                                                                                                         //
// //                                                                                                                                                         //
// //                                                                                                                                                         //
// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// //Process creates necessary directories in each folder of interest respectively, performs minimap indexing, and performs a structure conversion of the gtf to a csv file
// process make_directories{
//     input:
//     val done
//     output:
// 	val finished, emit: finished
	
// 	script:
// 	string = ""
// 	for(i in 0..params.sample_names.size()-1){
// 	    string = string + params.output_dir + params.sample_names.get(i) + "/ "
// 	}
// 	println "Make directories in ${string}"
// 	finished = 1
// 	"""
//     if [ ! -d ${params.output_dir} ]
// 	then
// 	    mkdir ${params.output_dir}
// 	fi

//     if [ ! -d ${params.output_dir}r_objects ]
//     then   
//         mkdir ${params.output_dir}r_objects
//     fi

//     if [ ! -d ${params.output_dir}error_logs ]
//     then
//         mkdir ${params.output_dir}error_logs
//     fi

//     if [ ! -d ${params.output_dir}bam_genome_merged ]
//     then
//         mkdir ${params.output_dir}bam_genome_merged
//     fi
//     if [ ! -d ${params.output_dir}bam_transcriptome_merged ]
//     then
//         mkdir ${params.output_dir}bam_transcriptome_merged
//     fi
//     if [ ! -d ${params.output_dir}ReadLengthFolder ]
//     then
//         mkdir ${params.output_dir}ReadLengthFolder
//     fi
// 	for i in ${string}
// 	do
//       if [ ! -d \${i} ]
//       then
//          mkdir \${i}
//       fi

//       if [ ! -d \${i}bam_files ]
//       then
//         mkdir \${i}bam_files
//       fi

//       if [ ! -d \${i}bam_files_transcripts ]
//       then
//         mkdir \${i}bam_files_transcripts 
//       fi

//       if [ ! -d \${i}bam_files_transcripts/full ]
//       then
//             mkdir \${i}bam_files_transcripts/full
//       fi

// 	  if [ ! -d \${i}merged_fastq ]
// 	  then
// 	     mkdir \${i}merged_fastq
// 	  fi 

// 	  if [ ! -d \${i}merged_fc ]
// 	  then
// 	     mkdir \${i}merged_fc
// 	  fi 

//       if [ ! -d \${i}merged_fc_splice ]
// 	  then
// 	     mkdir \${i}merged_fc_splice
// 	  fi 

//       if [ ! -d \${i}single_fc ]
//       then 
//          mkdir \${i}single_fc
//       fi

//       if [ ! -d \${i}single_fc_splice ]
//       then 
//          mkdir \${i}single_fc_splice
//       fi
      
//       if [ ! -d \${i}fastqc ]
//       then
//          mkdir \${i}fastqc
//       fi
      
//       if [ ! -d \${i}salmon ]
//       then
//          mkdir \${i}salmon
//       fi

// 	done
//     echo 1 > ${params.output_dir}process_running.txt

//     if [ ${params.drs} -eq 1 ]
//     then
//         minimap2 --MD -ax splice -uf -k14 -d ${params.output_dir}MT-human_ont.mmi $params.genome_fasta || echo "Indiexing genome failed" >> ${params.output_dir}error_logs/index_genome.log
//         minimap2 --MD -ax map-ont -uf -k14 -d ${params.output_dir}MT-human_transcript_ont.mmi $params.transcriptome_fasta || echo "Indexing transcriptome failed" >> ${params.output_dir}error_logs/index_transcriptome.log
//     else
//         minimap2 --MD -ax splice -d ${params.output_dir}MT-human_ont.mmi $params.genome_fasta || echo "Indiexing genome failed" >> ${params.output_dir}error_logs/index_genome.log
//         minimap2 --MD -ax map-ont -d ${params.output_dir}MT-human_transcript_ont.mmi $params.transcriptome_fasta || echo "Indexing transcriptome failed" >> ${params.output_dir}error_logs/index_transcriptome.log
//     fi
//     """
// }


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                         // 
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                Feature percentiles and gtf conversion                                                   //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
//                                                                                                                                                         //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// process create_feature_percentiles{
//     input:
//     val finished

//     script:
//     """
//     workflow-glue createFeaturePercentiles --bed ${params.bed_file} --output_dir ${params.output_dir} || echo "Feature percentile preprocess not working" >> ${params.output_dir}error_logs/feature_percentile_preprocess_failed.log
//     """
// }


// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// //                                                                                                                                                         // 
// //                                                                                                                                                         //
// //                                                                                                                                                         //
// //                                                                Progress lookup in run folder                                                            //
// //                                                                                                                                                         //
// //                                                                                                                                                         //
// //                                                                                                                                                         //
// //                                                                                                                                                         //
// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// //Process reads file with read ids that have already been processed 
// //and checks for the current iteration of the nextflow pipeline
// process see_progress{
// 	input:
// 	val finished
	
// 	output:
// 	val done

// 	script:
// 	File file = new File("${params.output_dir}data_seen.txt")
// 	if(file.exists()){
//     		    println "Data seen exists"
//     		    def line = 0;
//     		    params.data_seen_list = []
//                 file.withReader{reader -> while((line = reader.readLine()) != null){
//      		    data = line.split(" ").collect{it as String}
//                 if (data.size() > 0){
//                     for (i in 0..data.size()-1){
//                         params.data_seen_list.add(data.get(i))
//                         }
//                     }
//                 }
//                 }
//                 if (params.data_seen_list.size() == 0){
//                     file.write(" ")
//                     params.data_seen_list = []
//                 }
	    
//     }
// 	else{
//    	    println "Data seen does not exist"
//    	    params.data_seen_list = []
//    	    file.write(" ")
// 	}
//     File file2 = new File("${params.output_dir}current_iteration.txt")
//     if(file2.exists()){
//     		    def line = 0;
//                 file2.withReader{reader -> while((line = reader.readLine()) != null){
//      		    data = line.split(" ").collect{it as int}
//                 if (data.size() > 0){
//                     iteration.value = data.get(0)
//                 }
//                 }
//         }
//     }
//     else
//     {
//         iteration.value = 0
//     }
//     println "Current iteration $iteration.value"
//     """
//     iteration=${iteration.value}
//     if [ \$iteration -gt 0 ]
//     then
//     iteration=\$((\${iteration}+1))
//     awk -F "," -v iter=\$iteration '{ if (\$2 < iter) print \$1,\$2,\$3 }' ${params.output_dir}processing_time_table.csv | sed 's/ /,/g;w ${params.output_dir}processing_time_table2.csv'
//     echo Tool,Iteration,Time > ${params.output_dir}processing_time_table3.csv & cat ${params.output_dir}processing_time_table2.csv >> ${params.output_dir}processing_time_table3.csv
//     cp ${params.output_dir}processing_time_table3.csv ${params.output_dir}processing_time_table.csv
//     rm ${params.output_dir}processing_time_table2.csv
//     rm ${params.output_dir}processing_time_table3.csv
//     fi
//     """


// }
