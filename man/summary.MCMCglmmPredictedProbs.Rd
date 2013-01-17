\name{summary.MCMCglmmPredictedProbs}
\alias{summary.MCMCglmmPredictedProbs}
\title{Summary method for MCMCglmm predicted probabilities}
\usage{
  \method{summary}{MCMCglmmPredictedProbs} (object,
    marginalize = FALSE, level = 0.95, ...)
}
\arguments{
  \item{object}{A \code{MCMCglmmPredictedProbs} object to
  summarize}

  \item{marginalize}{Logical whether or not to marginalize
  by taking the row means. Defaults to \code{FALSE}.}

  \item{level}{A numeric value, the value to use when
  calculating HPD intervals. Defaults to .95.}

  \item{\dots}{Not currently used.}
}
\value{
  If HPD intervals are calculated, returns a list of
  matrices with the means, lower limit, and upper limit. If
  no HPD intervals are calculated, returns a list of means.
}
\description{
  Optionally first marginalizes across all observations by
  taking the row means. If the predicted values only used
  the posterior means, highest posterior density (HPD)
  intervals cannot be generated, so only the means are
  returned. Otherwise, it calculates the mean predicted
  probability, as well as the HPD interval. This can either
  be per observation or marginalized.
}
\examples{
\dontrun{
  ## Make me!
}
}
\seealso{
  \code{\link{predict2.MCMCglmm}}, \code{\link{recycler}}
}

