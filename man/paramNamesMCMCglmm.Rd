\name{paramNamesMCMCglmm}
\alias{paramNamesMCMCglmm}
\title{Extract the parameter names from an \code{MCMCglmm} object}
\usage{
  paramNamesMCMCglmm(object, ...)
}
\arguments{
  \item{object}{An \code{MCMCglmm} object}

  \item{\dots}{not used}
}
\value{
  A list with two elements: \item{fixed}{A character vector
  of the fixed effects parameter names} \item{random}{A
  character vector of the random effects parameter names}
}
\description{
  Simple function to extract the fixed and random effects
  parameter names from an \code{MCMCglmm} object.
}
\examples{
\dontrun{
  # a simple MCMCglmm model
  data(PlodiaPO)
  m <- MCMCglmm(PO ~ 1, random = ~ FSfamily, data = PlodiaPO, verbose=FALSE)

  # extract the parameter names
  paramNamesMCMCglmm(m)
}
}
\seealso{
  \code{\link{fixef.MCMCglmm}},
  \code{\link{ranef.MCMCglmm}}
}

