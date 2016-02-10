# nextflow installation
See https://github.com/nextflow-io/nextflow for more information

## Installation environment

Installation environment is a vagrant machine puphpet/centos65-x64

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

