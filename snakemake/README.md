# Snakemake

Johannes Köster and Sven Rahmann. *Snakemake—a scalable bioinformatics workflow engine*. Bioinformatics (2012) 28 (19): 2520-2522.
http://bioinformatics.oxfordjournals.org/content/28/19/2520.abstract

## Installation

Documentation: https://bitbucket.org/snakemake/snakemake/wiki/Documentation

Github repository: https://bitbucket.org/snakemake/snakemake/wiki/Home

For this tutorial it is recommended that you install everything through miniconda3
```
> wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
> bash Miniconda3-latest-Linux-x86_64.sh
> conda install -c bioconda snakemake
```

Next we will create a workspace and create and activate a virtual environment to work in
```
> mkdir workspace
> cd workspace
> conda create -n ruffusenv python
Fetching package metadata: ....
Solving package specifications: ............
Package plan for installation in environment /home/vagrant/miniconda3/envs/myenv:

The following NEW packages will be INSTALLED:

    openssl:    1.0.2f-0
    pip:        8.0.2-py35_0
    python:     3.5.1-0
    readline:   6.2-2
    setuptools: 19.6.2-py35_0
    sqlite:     3.9.2-0
    tk:         8.5.18-0
    wheel:      0.29.0-py35_0
    xz:         5.0.5-0
    zlib:       1.2.8-0
> source activate ruffusenv
prepending /home/vagrant/miniconda3/envs/myenv/bin to PATH
```

Now we can execute the pipeline. Download the data directory and Snakefile to your workspace and run
```
> snakemake
```
