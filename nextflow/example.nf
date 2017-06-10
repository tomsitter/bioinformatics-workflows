// Example sequence alignment workflow implemented in nextflow
 
//   © Copyright Government of Canada 2009-2016
//   Written by: Tom Sitter, for Agriculture and Agri-Food Canada


params.in = "$baseDir/data/samples/*.fastq"
params.genome = "$baseDir/data/genome.fa"

fastqs = Channel.fromPath { params.in }
                .ifEmpty { error "Cannot find any reads matching: ${params.in}" }

sequences = file(params.in)
genome_file = file(params.genome)

process bwa_map {
    publishDir "$baseDir/mapped_reads/"
    
    input:
    file genome_file
    file fastq from sequences

    output:
    file "${fastq.baseName}.bam" into mapped_reads

    """
    bwa mem $genome_file $fastq | samtools view -Sb -> ${fastq.baseName}.bam
    """
}

process samtools_sort {
    publishDir "$baseDir/sorted_reads" 
    
    input:
        file bamfile from mapped_reads
    output:
        file "sorted/$bamfile" into sorted_reads

    """
    mkdir sorted
    samtools sort -T ${bamfile.baseName} -O bam $bamfile > sorted/$bamfile
    """
}


process samtools_index {
    publishDir "$baseDir/sorted_reads/" 
    
    input:
        file sorted from sorted_reads
    output: 
        file "${sorted}.bai" 
    """
    samtools index $sorted
    """
}
