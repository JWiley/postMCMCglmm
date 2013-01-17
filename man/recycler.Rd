\name{recycler}
\alias{recycler}
\title{Calculate change in predicted probabilities}
\usage{
  recycler(object, index = 2L, twiddle, values, ...)
}
\arguments{
  \item{object}{A \code{MCMCglmm} model object to use for
  recycled predictions.}

  \item{index}{An integer indicating the column of the
  fixed effects design matrix, X, to vary.  Defaults to
  2L.}

  \item{twiddle}{A twiddle value for continuous variables.
  Needs to be small enough for the scale of the predictor
  that a twiddle change is a reasonable approximation of
  taking the first derivative at a point. That is, a very
  small change.  If missing, reverts to .01.}

  \item{values}{Specific values to use for the varying
  predictor. These are primarily meant for discrete
  predictors rather than continuous ones.}

  \item{\dots}{Passed on to \code{predict2}}
}
\value{
  A list of class \code{MCMCglmmPredictedProbs} that are
  the differences in predicted probabilities for a one unit
  change (calculated from the twiddle value or between the
  discrete values supplied in \code{values}).
}
\description{
  \code{recycler} wraps many of the functions in
  \pkg{postMCMCglmm} to calculate the change in predicted
  probabilities for a twiddle change in the predictor, or
  for discrete predictors, it can use values so it is the
  change from 0 to 1 (for example). The result is a
  MCMCglmmPredictedProbs (of course a difference but still)
  object, so it can be summarized using the
  \code{MCMCglmmPredictedProbs} \code{summary} method. This
  gives average marginal recycled predicted probabilities,
  as well as highest posterior density intervals.
}
\examples{
\dontrun{
  ## Make me!
}
}
\seealso{
  \code{\link{summary.MCMCglmmPredictedProbs}},
  \code{\link{predict2.MCMCglmm}}
}

