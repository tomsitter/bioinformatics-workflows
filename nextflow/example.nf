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
    input:
    file genome_file
    file fastq from sequences

    output:
    file "${fastq.baseName}.bam" into mapped_reads

    """
    bwa mem $baseDir/data/$genome_file $baseDir/data/samples/$fastq | samtools view -Sb -> ${fastq.baseName}.bam
    cp ${fastq.baseName}.bam $baseDir/mapped_reads/
    """
}

process samtools_sort {
    input:
        file bamfile from mapped_reads
    output:
        file "sorted/$bamfile" into sorted_reads

    """
    mkdir sorted
    samtools sort -T ${bamfile.baseName} -O bam $bamfile > sorted/$bamfile
    cp sorted/$bamfile $baseDir/sorted_reads/$bamfile
    """
}


process samtools_index {
    input:
        file sorted from sorted_reads

    """
    samtools index $sorted
    cp ${sorted}.bai $baseDir/sorted_reads/${sorted}.bai
    """
}

