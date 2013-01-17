\name{predict2}
\alias{predict2}
\alias{predict2.MCMCglmm}
\title{Define a generic prediction function}
\usage{
  predict2(object, ...)

  \method{predict2}{MCMCglmm} (object, X, Z,
    use = c("all", "mean"), type = c("lp", "response"),
    ...)
}
\arguments{
  \item{object}{A model object to predict from}

  \item{\dots}{Additional arguments passed to the methods}

  \item{X}{The fixed effects design matrix. Can be the
  original or new data.}

  \item{Z}{The random effects design matrix. Can be the
  original or new data.}

  \item{use}{A character string. Use just the posterior
  \dQuote{mean} or \dQuote{all} posterior samples (the
  default).}

  \item{type}{A cahracter string. Either \dQuote{lp} for
  the linear predictor (the default) or \dQuote{response}
  for the predicted values on the response scale.}
}
\value{
  Either a matrix of the linear predictor if \code{type =
  "lp"} or a list of class MCMCglmmPredictedProbs if
  \code{type = "response"}
}
\description{
  This defintes a generic \code{predict2} function which is
  similar to the usual \code{predict} but can use different
  methods. In particular, the \code{MCMCglmm} method has
  features not available in the regular \code{predict}
  method for \code{MCMCglmm} objects.

  This is the main workhorse of the package. It is a
  \code{predict2} method for MCMCglmm objects. There are a
  few core arguments. The model \code{object} and design
  matrices, \code{X} (fixed effects) and \code{Z} (random
  effects). If \code{X} and \code{Z} are missing, it will
  attempt to fill them in from the model object (which
  optionally saves them). If \code{X} and \code{Z} are
  specified or \code{NULL}, they are not used. This is
  useful either for out of sample predictions or to use
  just the fixed effects. Note that these must be full
  design matrices, not data matrices.  For example, they
  must dummy code factors and include the intercept (if
  there was an intercept in the model).
}
\details{
  You can also use all posterior samples or just the mean.
  All is nice because it lets you construct highest
  posterior density (HPD) intervals around the predicted
  values, rather than just get an estimate. The mean is
  nice because if that is all you care about, it is much
  much faster. You can get either the linear predictor
  values or the response scale. However, response is
  currently only implemented for ordinal (probit) models.
  Theoretically it could be extended but the code is a
  pain.
}
\examples{
# to see available methods
methods(predict2)
\dontrun{
  ## Make me!
}
}
\seealso{
  \code{\link{summary.MCMCglmmPredictedProbs}},
  \code{\link{recycler}}
}

