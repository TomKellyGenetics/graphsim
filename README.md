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

To get the development version of this package from github:

```R
# install.packages("devtools")
devtools::install_github("TomKellyGenetics/graphsim")
```

## Vignettes

Please see the following demonstrations of this (and other iGraph extensions) on examples of simple simulated networks.


* https://rawgit.com/TomKellyGenetics/graphsim/master/vignettes/test_graph_convergent.html
* https://rawgit.com/TomKellyGenetics/graphsim/master/vignettes/test_graph_divergent.html
* https://rawgit.com/TomKellyGenetics/graphsim/master/vignettes/test_graph_reconvergent.html
* https://rawgit.com/TomKellyGenetics/graphsim/master/vignettes/test_graph_network.html
* https://rawgit.com/TomKellyGenetics/graphsim/master/vignettes/test_graph_network_inhibiting.html

## Citation

To cite the graphsim package in publications use:

  S. Thomas Kelly and Michael A. Black (2018). graphsim: Simulate Expression data from iGraph networks 
  R package version 0.1.0. https://github.com/TomKellyGenetics/vioplot doi:10.5281/zenodo.1313986

A BibTeX entry for LaTeX users is

```
  @Manual{,
    title = {{graphsim}: Simulate Expression data from {iGraph} networks},
    author = {S. Thomas Kelly and Michael A. Black},
    year = {2018},
    note = {R package version 0.1.0},
    url = {https://github.com/TomKellyGenetics/graphsim},
    doi = {10.5281/zenodo.1313986},
  }
```

Please also acknowledge the manuscript describing use of this package (currently in final preparations for submission).


