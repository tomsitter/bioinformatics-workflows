# Ruffus

Goodstadt L. *Ruffus: a lightweight Python library for computational pipelines*. Bioinformatics. 2010 Nov 1;26(21):2778-9.
http://www.ncbi.nlm.nih.gov/pubmed/20847218

Online documentation available at: http://www.ruffus.org.uk/contents.html

Github repository available at: https://github.com/bunbun/ruffus

## Installation environment

Installation environment is a vagrant machine  [puphpet/centos65-x64](https://atlas.hashicorp.com/puphpet/boxes/centos65-x64)

CentOS release 6.7 (Final)  2.6.32-573.12.1.el6.x86_64

## Installation

ruffus is available as a Python 2 or 3 pip package.

If you wish to visualize your workflow you will also need Graphviz

```
> sudo pip install ruffus --upgrade
ruffus 2.6.3
> sudo yum install graphviz-python.x86_64
```

To run the ruffus workflow 
```
> python ruffus_example.py
```

_Note: This assumes you have the bioinformatics tools installed. See snakemake directory for those instructions_
