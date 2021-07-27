[![Travis-CI Build Status](https://travis-ci.com/TomKellyGenetics/graphsim.svg?branch=master)](https://travis-ci.com/TomKellyGenetics/graphsim)
[![CircleCI build status](https://circleci.com/gh/TomKellyGenetics/graphsim.svg?style=svg)](https://circleci.com/gh/TomKellyGenetics/graphsim)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/TomKellyGenetics/graphsim?branch=master&svg=true)](https://ci.appveyor.com/project/TomKellyGenetics/graphsim)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/graphsim)](https://cran.r-project.org/package=graphsim)
[![Downloads](https://cranlogs.r-pkg.org/badges/graphsim)](https://CRAN.R-project.org/package=graphsim)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/graphsim?color=orange)](https://CRAN.R-project.org/package=graphsim)

[![DOI](https://zenodo.org/badge/67395230.svg)](https://zenodo.org/badge/latestdoi/67395230)
[![bioRxiv](https://img.shields.io/badge/biorXiv-10.1101%2F2020.03.02.972471-blue)](https://doi.org/10.1101/2020.03.02.972471)
[![status](https://joss.theoj.org/papers/96016c6a55d7f74bacebd187c6ededd6/status.svg)](https://joss.theoj.org/papers/10.21105/joss.02161)
[![](https://img.shields.io/badge/Altmetric-72-blue.svg)](https://www.altmetric.com/details/77053356)

--------------------------------------------------

# graphsim

## Version 1.0.2


###  Simulate Expression Data from 'igraph' Networks 

This package provides functions to develop simulated continuous data 
(e.g., gene expression) from a sigma covariance matrix derived from a 
graph structure in 'igraph' objects. Intended to extend 'mvtnorm' to 
take 'igraph' structures rather than sigma matrices as input. This 
allows the use of simulated data that correctly accounts for pathway
relationships and correlations. Here we present a versatile statistical 
framework to simulate correlated gene expression data from biological 
pathways, by sampling from a multivariate normal distribution derived 
from a graph structure. This package allows the simulation of biological
pathways from a graph structure based on a statistical model of 
gene expression, such as simulation of expression profiles that
of log-transformed and normalised data from microarray and RNA-Seq data.

#### Motivation

Network analysis of molecular biological pathways is important
for insights into biology and medical genetics. 
Gene expression profiles capture the regulatory state of a cell
and can be used to analyse complex molecular states with genome-scale data.
Biological pathways are more than simply sets of genes involved in functions,
they are rich in information of relationships defined by pathway structure.

Methods to infer biological pathways and gene regulatory networks from gene
expression data can be tested on  simulated datasets using this framework. This also allows for
pathway structures to be considered as a confounding variable when 
simulating gene expression data to test the performance of genomics analyses.

This package enables the generation of simulated gene expression datasets
containing pathway relationships from a known underlying network.
These simulated datasets can be used to evaluate various bioinformatics
methodologies, including statistical and network inference procedures.

Network analysis techniques have an important role in understanding
of biological pathways and interpretation of genomics studies.
Modelling biological pathways allows the evaluation of gene
regulatory network inference techniques (which so far rely on
experimental validation or resampling). This technique also
enables modelling datasets with correlated pathway-structures
to assess whether other genomics analysis techniques perform
as expected with the background of complex pathways.


## Installation

To install the latest release from CRAN:

```R
install.packages("graphsim")
```

To install the stable release of this package from github:

```R
# install.packages("devtools")
devtools::install_github("TomKellyGenetics/graphsim", ref = "master")
```

To get the development version of this package from github:

```R
# install.packages("devtools")
devtools::install_github("TomKellyGenetics/graphsim", ref = "dev")
```

## Usage

Please see the vignettes for demonstrations of this package on examples of simple simulated networks and the reactome pathway TGF-&beta; receptor signaling activates SMADs (R-HSA-2173789). An [article](https://doi.org/10.21105/joss.02161) with further details has been published in the 
_Journal of Open Source Software_.

A help menu can also be accessed within the R environment:

```
?graphsim
```

```
help("graphsim-package")
```

This will display a help page and link to documentation for each function.

--------------------------------------------------

## Citation

To cite package 'graphsim' in publications use:

>S. Thomas Kelly and Michael A. Black (2020). graphsim: Simulate Expression Data from
>'igraph' Networks. R package version 1.0.2.
>https://github.com/TomKellyGenetics/graphsim doi:10.5281/zenodo.3931288

A BibTeX entry for LaTeX users is:

```
  @Manual{,
    title = {{graphsim}: Simulate Expression Data from 'igraph' Networks },
    author = {S. Thomas Kelly and Michael A. Black},
    year = {2020},
    note = {R package version R package version 1.0.2.},
    url = {https://github.com/TomKellyGenetics/graphsim},
    doi = {10.5281/zenodo.3931288},
  }
```

Please also cite the publication describing use of this package where appropriate.

>Kelly, S.T. and Black, M.A. (2020). graphsim: An R package for simulating gene
>expression data from graph structures of biological pathways.
>_Journal of Open Source Software_, **5**(51), 2161, https://doi.org/10.21105/joss.02161


```
  @article{Kelly2020joss02161,
    doi = {10.21105/joss.02161},
    url ={https://doi.org/10.21105/joss.02161},
    year = {2020},
    publisher = {The Open Journal},
    volume = {5},
    number = {51},
    pages = {2161},
    author = {S. Thomas Kelly and Michael A. Black},
    title = {graphsim: An R package for simulating gene expression data from graph structures of biological pathways},
    journal = {Journal of Open Source Software} }
```

This article is also avaliable as a preprint.

>S. Thomas Kelly, Michael A. Black (2020)
> graphsim: An R package for simulating gene expression data from graph structures of biological pathways
> bioRxiv 2020.03.02.972471; doi:https://doi.org/10.1101/2020.03.02.972471

```
@article {Kelly2020.03.02.972471,
	author = {Kelly, S Thomas and Black, Michael A},
	title = {graphsim: An R package for simulating gene expression data from graph structures of biological pathways},
	elocation-id = {2020.03.02.972471},
	year = {2020},
	doi = {10.1101/2020.03.02.972471},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2020/03/04/2020.03.02.972471},
	eprint = {https://www.biorxiv.org/content/early/2020/03/04/2020.03.02.972471.full.pdf},
	journal = {bioRxiv}
}
```

--------------------------------------------------

## Contributions and Bug Reports

Please submit [issues](https://github.com/TomKellyGenetics/graphsim/issues) on GitHub to report
problems or suggest features. [Pull requests](https://github.com/TomKellyGenetics/graphsim/pulls)
to the `dev` branch on GitHub are also welcome to add features or correct problems. Please see
the [contributor guide](https://github.com/TomKellyGenetics/graphsim/blob/master/CONTRIBUTING.md) for more details.


