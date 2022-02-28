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
technical factors, and a vector of size factors.
Check out [this vignette](https://htmlpreview.github.io/?https://github.com/KaiyangZ/PRIMUS/blob/master/vignettes/quickstart.html) for a quick start tutorial.

## Citation
Zhang K, Erkan EP, Jamalzadeh S, Dai J, Andersson N, Kaipio K, Lamminen T, Mansuri N, Huhtinen K, Carpén O, Hietanen S, Oikkonen J, Hynninen J, Virtanen A, Häkkinen A, Hautaniemi S, Vähärautio A. Longitudinal single-cell RNA-seq analysis reveals stress-promoted chemoresistance in metastatic ovarian cancer. Sci Adv. 2022 Feb 25;8(8):eabm1831. doi: 10.1126/sciadv.abm1831. Epub 2022 Feb 23. PMID: 35196078.
