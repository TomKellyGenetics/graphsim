[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/graphsim)](https://cran.r-project.org/package=graphsim)
[![Travis-CI Build Status](https://travis-ci.org/TomKellyGenetics/graphsim.svg?branch=master)](https://travis-ci.org/TomKellyGenetics/graphsim)
[![CircleCI build status](https://circleci.com/gh/TomKellyGenetics/graphsim.svg?style=svg)](https://circleci.com/gh/TomKellyGenetics/graphsim)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/TomKellyGenetics/graphsim?branch=master&svg=true)](https://ci.appveyor.com/project/TomKellyGenetics/graphsim)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![codecov](https://codecov.io/gh/TomKellyGenetics/graphsim/branch/master/graph/badge.svg)](https://codecov.io/gh/TomKellyGenetics/graphsim)
[![DOI](https://zenodo.org/badge/67395230.svg)](https://zenodo.org/badge/latestdoi/67395230)

# graphsim

### Simulate Expression data from iGraph networks

Functions to develop simulated continuous data (e.g., gene expression) from a sigma covariance matrix derived from a graph structure in iGraph. Intended to extend mvtnorm to take iGraph structures rather than sigma matrices as input.

## Installation

To get the dependency package from CRAN:

```R
install.packages(c("igraph", "Matrix", "mvtnorm", "matrixcalc"))
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

Please see the vignettes for demonstrations of this package on examples of simple simulated networks and the reactome pathway TGF-&beta; receptor signaling activates SMADs (R-HSA-2173789). A manuscript with further details has been submitted for peer-review.

## Citation

To cite the graphsim package in publications use:

  S. Thomas Kelly and Michael A. Black (2020). graphsim: Simulate Expression data from iGraph networks 
  R package version 0.1.1. https://github.com/TomKellyGenetics/vioplot doi:10.5281/zenodo.1313986

A BibTeX entry for LaTeX users is

```
  @Manual{,
    title = {{graphsim}: Simulate Expression data from {iGraph} networks},
    author = {S. Thomas Kelly and Michael A. Black},
    year = {2020},
    note = {R package version 0.1.1},
    url = {https://github.com/TomKellyGenetics/graphsim},
    doi = {10.5281/zenodo.1313986},
  }
```

Please also acknowledge the manuscript describing use of this package once it is published. It is current avaliable as a preprint.

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
