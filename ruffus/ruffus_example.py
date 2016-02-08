''' This script is a simple sequence alignment workflow using ruffus

    © Copyright Government of Canada 2009-2016
    Written by: Tom Sitter, for Agriculture and Agri-Food Canada
'''

# All you need to get started
from ruffus import *
import os

# ruffus.cmdline can be used to pass variables in as arguments
starting_files = 'data/samples/*.fastq'
faFile = 'data/genome.fa'

# Directories must be created manually
@follows(mkdir('mapped_reads'))
@transform(starting_files, 
           suffix(".fastq"), 
           ".bam", 
           output_dir='mapped_reads')
def bwa_map(fqFile, bamFile):
    #python's subprocess.Popen() can be used here as well
    os.system(
        "bwa mem {original} {input} | samtools view -Sb - > {output}".format(
             original=faFile, input=fqFile, output=bamFile
        )
    )

    ''' The above syntax can also be shorted as follows
      os.system("bwa mem %s %s | samtools view -Sb - > %s". % faFile, fqFile, bamFile
    '''

@follows(mkdir('sorted_reads'))
@transform(bwa_map, 
           suffix(".bam"), 
           ".sorted.bam", 
           output_dir='sorted_reads')
def samtools_sort(bamFile, sortedBamFile):
    os.system(
        "samtools sort -T sorted_reads/{input} " \
        "-O bam {input} > {output}".format(
            input=bamFile, output=sortedBamFile
        )
    )

@transform(samtools_sort, 
           suffix(".sorted.bam"), 
           ".bai")
def samtools_index(bamFile, baiFile):
    os.system(
        "samtools index {input}".format( input=bamFile )
    )


# verbosity goes from 1(least) to 10(most)
# An example output from this command is below
pipeline_printout(verbose=10)

# The number of processes can be specified (default is 1)
pipeline_run(multiprocess=3)

'''
________________________________________
Tasks which will be run:

Task = "mkdir('sorted_reads')   before samtools_sort "
    "Make missing directories [/vagrant/ruffus/sorted_reads]"
    Multiple jobs Multiple outputs
       Make directories [sorted_reads]
         Job needs update: : Directories 'sorted_reads' is missing

Task = "mkdir('mapped_reads')   before bwa_map "
    "Make missing directories [/vagrant/ruffus/mapped_reads]"
    Multiple jobs Multiple outputs
       Make directories [mapped_reads]
         Job needs update: : Directories 'mapped_reads' is missing

Task = 'bwa_map'
    Multiple jobs Multiple outputs
       Job  = [data/samples/A.fastq
             -> mapped_reads/A.bam]
         Job needs update: ...
               Missing file [mapped_reads/A.bam]
       Job  = [data/samples/B.fastq
             -> mapped_reads/B.bam]
         Job needs update: ...
               Missing file [mapped_reads/B.bam]
       Job  = [data/samples/C.fastq
             -> mapped_reads/C.bam]
         Job needs update: ...
               Missing file [mapped_reads/C.bam]

Task = 'samtools_sort'
    Multiple jobs Multiple outputs
       Job  = [mapped_reads/A.bam
             -> sorted_reads/A.sorted.bam]
         Job needs update: ...
               Missing files [mapped_reads/A.bam, sorted_reads/A.sorted.bam]
       Job  = [mapped_reads/B.bam
             -> sorted_reads/B.sorted.bam]
         Job needs update: ...
               Missing files [mapped_reads/B.bam, sorted_reads/B.sorted.bam]
       Job  = [mapped_reads/C.bam
             -> sorted_reads/C.sorted.bam]
         Job needs update: ...
               Missing files [mapped_reads/C.bam, sorted_reads/C.sorted.bam]

Task = 'samtools_index'
    Multiple jobs Multiple outputs
       Job  = [sorted_reads/A.sorted.bam
             -> sorted_reads/A.bai]
         Job needs update: ...
               Missing files [sorted_reads/A.sorted.bam, sorted_reads/A.bai]
       Job  = [sorted_reads/B.sorted.bam
             -> sorted_reads/B.bai]
         Job needs update: ...
               Missing files [sorted_reads/B.sorted.bam, sorted_reads/B.bai]
       Job  = [sorted_reads/C.sorted.bam
             -> sorted_reads/C.bai]
         Job needs update: ...
               Missing files [sorted_reads/C.sorted.bam, sorted_reads/C.bai]

________________________________________

________________________________________
Tasks which will be run:


Task enters queue = "mkdir('sorted_reads')   before samtools_sort "
Task enters queue = "mkdir('mapped_reads')   before bwa_map "
Completed Task = "mkdir('sorted_reads')   before samtools_sort "
Completed Task = "mkdir('mapped_reads')   before bwa_map "
Task enters queue = 'bwa_map'
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::bwa_idx_load_from_disk] read 0 ALT contigs
[M::process] read 25000 sequences (2525000 bp)...
[M::process] read 25000 sequences (2525000 bp)...
[M::process] read 25000 sequences (2525000 bp)...
[M::mem_process_seqs] Processed 25000 reads in 2.777 CPU sec, 8.350 real sec
[M::mem_process_seqs] Processed 25000 reads in 2.764 CPU sec, 8.317 real sec
[M::mem_process_seqs] Processed 25000 reads in 2.809 CPU sec, 8.455 real sec
[main] Version: 0.7.12-r1039
[main] CMD: bwa mem data/genome.fa data/samples/B.fastq
[main] Real time: 12.544 sec; CPU: 2.914 sec
[main] Version: 0.7.12-r1039
[main] CMD: bwa mem data/genome.fa data/samples/A.fastq
[main] Real time: 12.767 sec; CPU: 2.919 sec
[main] Version: 0.7.12-r1039
[main] CMD: bwa mem data/genome.fa data/samples/C.fastq
[main] Real time: 12.779 sec; CPU: 2.962 sec
Completed Task = 'bwa_map'
Task enters queue = 'samtools_sort'
Completed Task = 'samtools_sort'
Task enters queue = 'samtools_index'
Completed Task = 'samtools_index'
'''