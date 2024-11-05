import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList
import java.io.File

process ConvertGtfToDf{ 
    publishDir(
        path: "${params.out_dir}/",
        mode: 'move'
    )
    maxForks 3
    memory "8GB"
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}
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

process CopyBedfileAndMetadata{ 
    publishDir(
        path: "${params.out_dir}/",
        mode: 'copy'
    )
    maxForks 3
    memory "8GB"
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 20 as long); return 'retry'}
    stageInMode "copy"
    input:
    path(bedfile)
    path(metadata)

    output:
    path("./basic_bed_file.bed")
    path("./metadata.tsv")
    val 1, emit: done

    script:
    """
    cp ${bedfile} basic_bed_file.bed
    cp ${metadata} metadata.tsv
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
    stageInMode "copy"
    publishDir(
        path: "${params.out_dir}/", 
        mode: 'move'
    )
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

def fetch_latest_bams(dir){
    files = file("${dir}/*.bam").collect()
    println files
    return files
}

process CreateGenomeBamFilesForMerge{
    publishDir(path: "${params.out_dir}/bam_genome_merged/", mode: "copy", overwrite: false)
    input:
    path(metadata)
    output:
    path("*.bam")
    """
    python ${projectDir}/bin/init_bam_merged.py -m ${metadata}
    """
}

process CreateTranscriptomeBamFilesForMerge{
    publishDir(path: "${params.out_dir}/bam_transcriptome_merged/", mode: "copy", overwrite: false)
    input:
    path(metadata)
    output:
    path("*.bam")
    """
    python ${projectDir}/bin/init_bam_merged.py -m ${metadata}
    """
}