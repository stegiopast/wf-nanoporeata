import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList
import java.io.File




process MinimapIndex {
    label "nanoporeata"
    publishDir "${params.output_dir}"
    maxForks 1
    memory "25GB"
    input:
    path genome_fasta
    path transcriptome_fasta

    output:
    path("*")
    val 1,emit: done

    script:
    """
    if [ ${params.drs} -eq 1 ]
    then
        minimap2 --MD -ax splice -uf -k14 -d genome_index.mmi $genome_fasta 
        minimap2 --MD -ax map-ont -uf -k14 -d transcriptome_index.mmi $transcriptome_fasta 
    else
        minimap2 --MD -ax splice -d genome_index.mmi $genome_fasta 
        minimap2 --MD -ax map-ont -d transcriptome_index.mmi $transcriptome_fasta 
    fi
    """
}





process MinimapGenome {
    label "nanoporeata"
       publishDir(
        path: "${params.output_dir}/${ID}/bam_files/",
        mode: 'copy',
    )
    maxForks 1
    memory "12 GB"
    maxRetries 10
    cpus 1
    input:
    tuple val(ID), path(fastq)
    path fasta
    val done
    val done2


    output: 
    tuple val(ID), path("genes_${ID}.out${task.index}.bam"), emit: aligned_bams 
  
    script:
    """
    if [ ${params.drs} -eq 1 ]
    then
        minimap2 --MD -ax splice -uf -k14 ${fasta} ${fastq} | samtools view -hbS -F 3844 | samtools sort > genes_${ID}.out${task.index}.bam
    else
        minimap2 --MD -ax splice ${fasta} ${fastq} | samtools view -hbS -F 3844 | samtools sort > genes_${ID}.out${task.index}.bam 
    fi
    """
}

process MinimapTranscriptome {
    label "nanoporeata"
       publishDir(
        path: "${params.output_dir}/${ID}/bam_files_transcripts/",
        mode: 'copy',
    )
    maxForks 1
    memory "12 GB"
    maxRetries 10
    cpus 1
    input:
    tuple val(ID), path(fastq)
    path fasta
    val done
    val done2


    output: 
    tuple val(ID), path("transcripts_${ID}.out${task.index}.bam"), emit: aligned_bams 

    script:
    """
    if [ ${params.drs} -eq 1 ]
    then
        minimap2 --MD -ax map-ont -uf -k14 ${fasta} ${fastq} | samtools view -hbS -F 3844 | samtools sort > transcripts_${ID}.out${task.index}.bam
    else
        minimap2 --MD -ax map-ont ${fasta} ${fastq} | samtools view -hbS -F 3844 | samtools sort > transcripts_${ID}.out${task.index}.bam
    fi
    """
}