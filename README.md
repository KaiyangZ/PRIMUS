# PRIMUS
Poisson scRNA Integration of Mixed Unknown Signals 

## Overview 
PRIMUS is a holistic clustering approach that identifies phenotypic cell groups from 
the scRNA-seq data while accounting for date source (e.g. patient, sample, dataset) -specific components as well as technical noises.

PRIMUS uses a bilinear Poisson regression model to simultaneously factorize the
expression data into the defined nuisance factors, undefined cellular phenotypes, and their
corresponding transcriptomic profiles.

## Installation
To run PRIMUS, open R and install PRIMUS from github:
```
devtools::install_github("KaiyangZ/PRIMUS")
```

## Usage
As input PRIMUS takes raw counts from scRNA-seq experiments,
a design matrix encoding the different nuisance factors, such as patient labels and
technical factors, and a vector of size factors. PRIMUS then uses a bilinear Poisson regression model to simultaneously factorize the raw data into the defined nuisance factors, undefined cellular phenotypes, and their corresponding transcriptomic profiles.
Check out this vignette for a quick start tutorial.
