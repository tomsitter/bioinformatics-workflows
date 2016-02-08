#!/bin/bash -ue
bwa mem /vagrant/bioinformatics-workflows/nextflow/data/genome.fa /vagrant/bioinformatics-workflows/nextflow/data/samples/A.fastq | samtools view -Sb -> A.bam
cp A.bam /vagrant/bioinformatics-workflows/nextflow/mapped_reads/
