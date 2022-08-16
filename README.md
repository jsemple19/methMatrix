# methMatrix
Operations on methylation matrices from dSMF pipeline

## Installation
This package requires samtools to extract the methylation data for each read. I recommend installing [miniconda](https://docs.conda.io/en/latest/miniconda.html). And then install samtools with conda:
```
conda install -c bioconda samtools
```
Make sure the devtools package is installed and then install these packages from github:
```
#install.packages("devtools")
devtools::install_github("jsemple19/grangesutils")
devtools::install_github("jsemple19/methMatrix")
```




