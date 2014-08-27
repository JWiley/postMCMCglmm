postMCMCglmm
============

R package for post estimation from MCMCglmm. Predicted probabilities from MCMCglmm probit models

The PDF manual is available under /inst/doc/.

Installation
------------

```
#install.packages("devtools")
require(devtools)

install_github("postMCMCglmm", "JWiley")
```

Citation
--------

  Joshua Wiley (2013). postMCMCglmm: Average marginal predicted
  probabilities from Bayesian ordered probit models. R package version
  0.1-2. doi: 10.5281/zenodo.11461.

[![DOI](https://zenodo.org/badge/5776/JWiley/postMCMCglmm.png)](http://dx.doi.org/10.5281/zenodo.11461)

Examples
--------

To see some demonstrations of models with MCMCglmm as well as how to
use the postMCMCglmm package, run:

```
# load the package
require(postMCMCglmm)

# see available demos
demo(package="postMCMCglmm")

# run one showing a simple, fixed effects example
demo("simpleMCMCglmm", package="postMCMCglmm")

# for a demonstration of new features in this package
demo("postMCMCglmm", package="postMCMCglmm")
```
