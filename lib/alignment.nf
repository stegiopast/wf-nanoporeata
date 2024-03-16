import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList
import java.io.File
import java.time.LocalDateTime




process MinimapIndex {
    label "nanoporeata"
    publishDir "${params.output_dir}"
    maxForks 1
    maxRetries 10 
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
    memory "14GB"
    maxRetries 10
    cpus 4
    input:
    tuple val(ID), path(fastq)
    path fasta
    val done
    val done2


    output: 
    tuple val(ID), path("genes_${ID}.out${task.index}.bam"), emit: aligned_bams 
    path("genes_${ID}.out${task.index}.bam.bai"), emit:aligned_bam_bais

    script:
    """
    if [ ${params.drs} -eq 1 ]
    then
        minimap2 --MD -ax splice -uf -k14 -t ${task.cpus} ${fasta} ${fastq} | samtools view -hbS -F 3844 | samtools sort > genes_${ID}.out${task.index}.bam
        samtools index genes_${ID}.out${task.index}.bam
    else
        minimap2 --MD -ax splice -t ${task.cpus} ${fasta} ${fastq} | samtools view -hbS -F 3844 | samtools sort > genes_${ID}.out${task.index}.bam 
        samtools index genes_${ID}.out${task.index}.bam
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
    memory "14GB"
    maxRetries 10
    cpus 4
    input:
    tuple val(ID), path(fastq)
    path fasta
    val done
    val done2


    output: 
    tuple val(ID), path("transcripts_${ID}.out${task.index}.bam"), emit: aligned_bams 
    path("transcripts_${ID}.out${task.index}.bam.bai")

    script:
    """
    if [ ${params.drs} -eq 1 ]
    then
        minimap2 --MD -ax map-ont -uf -k14 -t ${task.cpus} ${fasta} ${fastq} | samtools view -hbS -F 3844 | samtools sort > transcripts_${ID}.out${task.index}.bam
        samtools index transcripts_${ID}.out${task.index}.bam
    else
        minimap2 --MD -ax map-ont -t ${task.cpus} ${fasta} ${fastq} | samtools view -hbS -F 3844 | samtools sort > transcripts_${ID}.out${task.index}.bam
        samtools index transcripts_${ID}.out${task.index}.bam
    fi
    """
}