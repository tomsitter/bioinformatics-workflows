# nextflow installation

Homepage: http://www.nextflow.io/  
Github Repo: https://github.com/nextflow-io/nextflow

## Installation environment

Installation environment is a vagrant machine  [puphpet/centos65-x64](https://atlas.hashicorp.com/puphpet/boxes/centos65-x64)

CentOS release 6.7 (Final)  2.6.32-573.12.1.el6.x86_64

## Install java 

Nextflow supports Java 7 and 8, but if compiled for 8 cannot be run on Java 7

```
> yum install java
> java -version
openjdk version "1.8.0_71"
OpenJDK Runtime Environment (build 1.8.0_71-b15)
OpenJDK 64-Bit Server VM (build 25.71-b15, mixed mode)
```

## Download and Install nextflow
```
> wget -qO- get.nextflow.io | bash
```
This creates a nextflow executable in current directory.
We should move it somewhere to make it callable
```
> mv nextflow /usr/local/bin/nextflow
> nextflow -v
nextflow version 0.17.2.3480
```

## Usage

In order to run the example workflow you will need the bioinformatics tools installed.  
See the snakemake README.md for instructions on installing the tools using conda

Download the example.nf and the data directory to a workspace.  
For our example we will also need to make two directories to store results  
We can now run our workflow using the nextflow command  
```
> cd workspace
> mkdir mapped_reads sorted_reads
> nextflow example.nf
```

## Passing in parameters
```
> nextflow run pipeline.nf --db=/path/to/blast/db --query=/path/to/query.fasta
```
## Running on SGE

See http://www.nextflow.io/docs/latest/executor.html#sge-executor for more information

Create a file named nextflow.config in your work directory
Include the following contents

```
process {
  executor='sge'
  queue='<your execution queue>'
}
```

### Limiting resources for a task

There are a number of directives you can use to limit resource usage. The main ones are
- cpus
- queue
- memory
- penv
- time
- clusterOptions

Here is an example for limiting the cpus of a task

```
process big_job {

  cpus 8
  executor 'sge'

  """
  blastp -query input_sequence -num_threads ${task.cpus}
  """
}
```

### Specify working directory

The NXT_WORK environmental variable can be used to set the working directory
```
> export NXT_WORK='/scratch/nextflow'
```

### Logging all commands run

The ability to log each command run is not yet available. I have put in a feature request

A current workaround is to traverse the work/ directory and extract the .command.sh contents  
HOWEVER, this directory can contain the work from multiple runs so it is important to only navigate to directories for the run of interest.  
You can generate a trace file using the -with-trace flag. This file will list the working directory of each task in the current workflow.
```
> nextflow example.nf -with-trace
> cat trace.txt
task_id hash    native_id       name    status  exit    submit  duration        realtime        %cpu    rss     vmem   wcharr
2       6b/cf40ba       2823    bwa_map (2)     COMPLETED       0       2016-02-10 17:12:21.436 6.1s    2.9s    33.0%  97 BMBMB
1       24/0853ba       2824    bwa_map (1)     COMPLETED       0       2016-02-10 17:12:21.428 6.4s    3s      21.0%  97 BMBMB
3       e1/118e5d       3141    bwa_map (3)     COMPLETED       0       2016-02-10 17:12:27.532 5.8s    2.1s    18.0%  97 BMBMB
4       ed/cc5cb9       3208    samtools_sort (1)       COMPLETED       0       2016-02-10 17:12:27.832 5.7s    966ms  1.3 MBMB
5       38/900a7d       3415    samtools_sort (2)       COMPLETED       0       2016-02-10 17:12:33.320 5.6s    767ms  1.3 MBMB
7       3b/901f0a       3648    samtools_index (1)      COMPLETED       0       2016-02-10 17:12:38.898 187ms   80ms   0.7 KBMB
6       a4/116f00       3468    samtools_sort (3)       COMPLETED       0       2016-02-10 17:12:33.500 5.6s    830ms  1.2 MBMB
8       3c/40ad43       3701    samtools_index (2)      COMPLETED       0       2016-02-10 17:12:39.087 125ms   97ms   --
9       a3/50699b       3722    samtools_index (3)      COMPLETED       0       2016-02-10 17:12:39.121 120ms   34ms   --
```

The hash column has the working directory. A simple python script can parse this file and compile the commands.
