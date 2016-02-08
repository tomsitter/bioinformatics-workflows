#!/bin/bash -ue
bwa mem /vagrant/bioinformatics-workflows/nextflow/data/genome.fa /vagrant/bioinformatics-workflows/nextflow/data/samples/B.fastq | samtools view -Sb -> B.bam
cp B.bam /vagrant/bioinformatics-workflows/nextflow/mapped_reads/
