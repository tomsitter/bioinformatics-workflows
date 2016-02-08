#!/bin/bash -ue
bwa mem /vagrant/bioinformatics-workflows/nextflow/data/genome.fa /vagrant/bioinformatics-workflows/nextflow/data/samples/C.fastq | samtools view -Sb -> C.bam
cp C.bam /vagrant/bioinformatics-workflows/nextflow/mapped_reads/
