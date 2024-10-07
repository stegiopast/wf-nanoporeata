import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList
import java.io.File
import java.time.LocalDateTime




process MinimapIndex {
    label "nanoporeata"
    publishDir(
        path: "${params.out_dir}"
    )
    maxForks 1
    maxRetries 10 
    errorStrategy 'retry'
    memory "25GB"
    input:
    path genome_fasta
    path transcriptome_fasta 

    output:
    path("genome_index.mmi")
    path("transcriptome_index.mmi")
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
        path: "${params.out_dir}/${ID}/bam_files/",
        mode: 'copy'
    )
    maxForks 1
    memory "14GB"
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry'}
    cpus 4
    input:
    tuple val(ID), path(fastq)
    path fasta, stageAs: 'genome_index.mmi'
    val done
    val done2


    output: 
    tuple val(ID), path("genes_${ID}.out${task.index}.bam"), emit: aligned_bams 
    path("genes_${ID}.out${task.index}.bam.bai"), emit:aligned_bam_bais

    script:
    """
    if [ ${params.drs} -eq 1 ]
    then
        minimap2 --MD -ax splice -uf -k14 -t ${task.cpus} genome_index.mmi ${fastq} > output.sam
        samtools view -hbS output.sam -F 3844 | samtools sort -o genes_${ID}.out${task.index}.bam
        rm output.sam
        #minimap2 --MD -ax splice -uf -k14 -t ${task.cpus} genome_index.mmi ${fastq} | samtools view -hb -F 3844 | samtools sort -o genes_${ID}.out${task.index}.bam
        samtools index genes_${ID}.out${task.index}.bam
    else
        #minimap2 --MD -ax splice -t ${task.cpus} genome_index.mmi ${fastq} | samtools view -hb -F 3844 | samtools sort -o genes_${ID}.out${task.index}.bam 
        minimap2 --MD -ax splice -t ${task.cpus} genome_index.mmi ${fastq} > output.sam
        samtools view -hbS output.sam -F 3844 | samtools sort -o genes_${ID}.out${task.index}.bam
        rm output.sam
        samtools index genes_${ID}.out${task.index}.bam
    fi
    """
}



process MinimapGenomeMergeBam{
    label "nanoporeata"
    cpus 2
    maxRetries 10
    maxForks 1
    memory "12GB"
    publishDir(
        path: "${params.out_dir}/bam_genome_merged/",
        mode: 'move'
    )
    input:
        tuple val(ID), path(bam)
        //file("${params.out_dir}/bam_genome_merged/*.bam")
        path bamFiles
    output:
        path("${ID}.bam")
        path("${ID}.bam.bai")
    script:
    """
    if [ -f ${ID}.bam ]
    then
        if [ -s ${ID}.bam ]
        then
            samtools merge -f --threads ${task.cpus} -o merged.bam ${ID}.bam ${bam} 
            samtools sort --threads ${task.cpus} merged.bam > merged_sorted.bam
            samtools view -h merged_sorted.bam > merged_sorted.sam
            grep '^@' merged_sorted.sam > merged_sorted_filtered.sam
            grep -v '^@' merged_sorted.sam | cut -f 1 | sort | uniq -u > unique_ids.txt
            grep -F -f unique_ids.txt merged_sorted.sam >> merged_sorted_filtered.sam
            samtools view -b merged_sorted_filtered.sam | samtools sort --threads ${task.cpus} - > final.bam
            rm -f merged.bam merged_sorted.bam merged_sorted.sam merged_sorted_filtered.sam unique_ids.txt barcode*.bam unclassified.bam
            mv final.bam ${ID}.bam
            samtools index ${ID}.bam
        else
            cp ${bam} final.bam
            rm -f barcode*.bam unclassified.bam
            mv final.bam ${ID}.bam
            samtools index ${ID}.bam
        fi
        
    else
        cp ${bam} final.bam
        samtools sort final.bam > final_sorted.bam
        rm -f final.bam barcode*.bam unclassified.bam
        mv final_sorted.bam ${ID}.bam
        samtools index ${ID}.bam
    fi
    """
}

process FindMergedGenomeBams{
    input:
        tuple val(ID), path(bam)
    output:
        tuple val(ID), path(bam), emit: aligned_bams
        path(files), emit: merged_bams
    script:
        files="${params.out_dir}/bam_genome_merged/${ID}.bam"
    """
    """
}

process MinimapTranscriptome {
    label "nanoporeata"
    publishDir(
        path: "${params.out_dir}/${ID}/bam_files_transcripts/",
        mode: 'copy',
    )
    maxForks 1
    memory "14GB"
    maxRetries 10
    errorStrategy {sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry'}
    cpus 4
    input:
    tuple val(ID), path(fastq)
    path fasta, stageAs: 'transcriptome_index.mmi'
    val done
    val done2

    output: 
    tuple val(ID), path("transcripts_${ID}.out${task.index}.bam"), emit: aligned_bams 
    path("transcripts_${ID}.out${task.index}.bam.bai")

    script:
    """
    if [ ${params.drs} -eq 1 ]
    then
        #minimap2 --MD -ax map-ont -uf -k14 -t ${task.cpus} transcriptome_index.mmi ${fastq} | samtools view -hbS -F 3844 | samtools sort -o transcripts_${ID}.out${task.index}.bam
        minimap2 --MD -ax map-ont -uf -k14 -t ${task.cpus} transcriptome_index.mmi ${fastq} > output.sam
        samtools view -hbS -F 3844 output.sam | samtools sort -o transcripts_${ID}.out${task.index}.bam
        rm output.sam
        samtools index transcripts_${ID}.out${task.index}.bam
    else
        #minimap2 --MD -ax map-ont -t ${task.cpus} transcriptome_index.mmi ${fastq} | samtools view -hbS -F 3844 | samtools sort -o transcripts_${ID}.out${task.index}.bam
        minimap2 --MD -ax map-ont -t ${task.cpus} transcriptome_index.mmi ${fastq} > output.sam
        samtools view -hbS -F 3844 output.sam | samtools sort -o transcripts_${ID}.out${task.index}.bam
        rm output.sam
        samtools index transcripts_${ID}.out${task.index}.bam
    fi
    """
}

process FindMergedTranscriptomeBams{
    input:
        tuple val(ID), path(bam)
    output:
        tuple val(ID), path(bam), emit: aligned_bams
        path(files), emit: merged_bams
    script:
        files="${params.out_dir}/bam_transcriptome_merged/${ID}.bam"
    """
    """
}


process MinimapTranscriptomeMergeBam{
    label "nanoporeata"
    cpus 2
    maxRetries 10
    maxForks 1
    memory "12GB"
    publishDir(
        path: "${params.out_dir}/bam_transcriptome_merged/",
        mode: 'move'
    )
    input:
        tuple val(ID), path(bam)
        //file("${params.out_dir}/bam_genome_merged/*.bam")
        path bamFiles
    output:
        path("${ID}.bam"), emit: bam_merged
        path("${ID}.bam.bai"), emit: bai_merged
    script:
    """
    if [ -f ${ID}.bam ]
    then
        if [ -s ${ID}.bam ]
        then
            samtools merge -f --threads ${task.cpus} -o merged.bam ${ID}.bam ${bam} 
            samtools sort --threads ${task.cpus} merged.bam > merged_sorted.bam
            samtools view -h merged_sorted.bam > merged_sorted.sam
            grep '^@' merged_sorted.sam > merged_sorted_filtered.sam
            grep -v '^@' merged_sorted.sam | cut -f 1 | sort | uniq -u > unique_ids.txt
            grep -F -f unique_ids.txt merged_sorted.sam >> merged_sorted_filtered.sam
            samtools view -b merged_sorted_filtered.sam | samtools sort --threads ${task.cpus} - > final.bam
            rm -f merged.bam merged_sorted.bam merged_sorted.sam merged_sorted_filtered.sam unique_ids.txt barcode*.bam unclassified.bam
            mv final.bam ${ID}.bam
            samtools index ${ID}.bam
        else
            cp ${bam} final.bam
            samtools sort final.bam > final_sorted.bam
            rm -f final.bam barcode*.bam unclassified.bam
            mv final_sorted.bam ${ID}.bam
            samtools index ${ID}.bam
        fi
    else
        cp ${bam} final.bam
        samtools sort final.bam > final_sorted.bam
        rm -f final.bam barcode*.bam unclassified.bam
        mv final_sorted.bam ${ID}.bam
        samtools index ${ID}.bam
    fi
    """
}
