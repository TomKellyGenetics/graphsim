[![Travis-CI Build Status](https://travis-ci.org/TomKellyGenetics/graphsim.svg?branch=master)](https://travis-ci.org/TomKellyGenetics/graphsim)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/TomKellyGenetics/graphsim?branch=master&svg=true)](https://ci.appveyor.com/project/TomKellyGenetics/graphsim)

# graphsim: Simulate Expression data from iGraph networks

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
