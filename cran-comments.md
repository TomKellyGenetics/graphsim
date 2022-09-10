## Test environments
* ubuntu 14.04 (on travis-ci), R 3.6.2, R 4.0.0, R 4.0.2
* ubuntu 14.04 (on circle-ci), R 3.6.2
* win-builder (devel and release) Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Windows Server 2012 R2 x64 on appveyor)  R 3.6.2 Patched 
* rhub (release) Ubuntu Linux 16.04 LTS, R-release, GCC
* Fedora (devel) Linux, R-devel, clang, gfortran
* MacOS 8.6.0 R R 3.6.2
* MacOS Mojave 10.14.6 R 3.6.1, R 3.6.2, R 4.0.2
* MacOS Catalina 10.15.7 R 4.2.0

## Local R CMD check results

── R CMD check results ───────────────────────────────────── graphsim 1.0.3 ────
Duration: 1m 55.4s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

R CMD check succeeded

## Possible issues

Reactome.org urls have been checked and direct to the correct database. These are only needed for documentation.

graphsim/tests/figs is used exclusively for testing plotting functions with "vdiffr".

Link to local file works on GitHub. Links changed to "\\href"" calls

## Release

This package has also been accepted at the Journal of Open-Source Software for peer-review. This release updates citation information to the publication. 

## Vignettes

Vignettes build with knitr::rmarkdown as used for this CRAN package: https://github.com/yixuan/prettydoc/

Building vignettes that contain real biological datasets are important to
demonstrate functionality of the package (requested by reviewers for JOSS).

This means checks take considerable time to build vignettes. Pre-generated HTMLs
are provided and "eval=FALSE" has been added to the Rmarkdown version of examples
that take a long time to run.

Vignettes are pregenerated to preserve headers, table of contents, and HTML style.

Vignettes in Rmarkdown have been marked as eval=FALSE in some cases (to avoid checks timing out on CRAN).

Pre-generated results so that code and results match will be considered in the future.
This is an urgent release to ensure that checks pass on CRAN (without timing out)
and the (published) package is not archived and remains accessible.

## Minor Release

I will maintain this package at my current address <tom.kelly[at]riken.jp>
