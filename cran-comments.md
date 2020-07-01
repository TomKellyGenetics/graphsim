## Test environments
* ubuntu 14.04 (on travis-ci), R 3.6.2, R 4.0.0, R 4.0.2
* ubuntu 14.04 (on circle-ci), R 3.6.2
* win-builder (devel and release) Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Windows Server 2012 R2 x64 on appveyor)  R 3.6.2 Patched 
* rhub (release) Ubuntu Linux 16.04 LTS, R-release, GCC
* Fedora (devel) Linux, R-devel, clang, gfortran
* MacOS 8.6.0 R R 3.6.2
* MacOS 10.14.6 R 3.6.1, R 3.6.2, R 4.0.2

## Local R CMD check results

── R CMD check results ───────────────────────────────────── graphsim 1.0.0 ────
Duration: 1m 16.5s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

R CMD check succeeded

## Possible issues

Reactome.org urls have been checked and direct to the correct database. These are only needed for documentation.

graphsim/tests/figs is used exclusively for testing plotting functions with "vdiffr".

Link to local file works on GitHub.

Vignettes build with knitr::rmarkdown as used for this CRAN package: https://github.com/yixuan/prettydoc/

Building vignettes that contain real biological datasets are important to
demonstrate functionality of the package (requested by reviewers for JOSS).
Unfortunately, this means checks take considerable time to build vignettes.

## Release

This package has also been submitted to the Journal of Open-Source Software for peer-review. This release accompanies a revised version of the manuscript (currently under review). This manuscript has also been subitted to the biorxiv preprint server.

A minor release will be needed in the future to address reviewer comments and update citation information to the publication should it be accepted. 

## Vignettes

Vignettes are pregenerated to preserve headers, table of contents, and HTML style.

## Major Release

I will maintain this package at my current address <tom.kelly[at]riken.jp>
