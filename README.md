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
