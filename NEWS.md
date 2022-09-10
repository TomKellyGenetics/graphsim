# graphsim 1.0.3

Minor update to pass CRAN checks

- removed imports for unused packages

- removed dependency on deprecated sparse matrices

- correct links in citations and vignettes

# graphsim 1.0.2

Updates maintainer contact details.

- resolves vignette formatting #11 

- passes updated CRAN checks (links updated)

# graphsim 1.0.1

* Update citation to reflect acceptance at JOSS

* Update documentation (package help page, links and cross-references)

* Critical changes to vignettes to reduce build time (required for regular CRAN checks)

# graphsim 1.0.0

* Major stable release: note changes to results are possible! (legacy code should run without breaking)

* Expanded documentation and examples (consolidate into fewer vignettes for clarity)

* Resolves errors handling inhibiting edges

* Efficiently compute a state matrix from a vector of edge properties from paths

* Enables passing "sd" (standard deviation) to alter covariance of Sigma matrix

* Adds methods for computing using Laplacian matrices

* Adds function to compute simulations directly from an adjacency matrix

* Migrates computing states to sigma (these matrices include inhibitions)

# graphsim 0.1.1

* Initial CRAN release

* Unit testing for all functions

* Checking for compatible inputs

* Passing layout parameters to plotting function

# graphsim 0.1.1

* Full vignette of biological pathway

* Merge plot_directed from function from plot.igraph to avoid github dependencies

* Added example data for reactome pathways

# graphsim 0.1.0

* Initial version of the package for generating simulated data
