\name{summary.MCMCglmmPredictedLP}
\alias{summary.MCMCglmmPredictedLP}
\title{Summary method for MCMCglmm predicted values}
\usage{
  \method{summary}{MCMCglmmPredictedLP} (object,
    level = 0.95, ...)
}
\arguments{
  \item{object}{A \code{MCMCglmmPredictedLP} object to
  summarize}

  \item{level}{A numeric value, the value to use when
  calculating HPD intervals. Defaults to .95.}

  \item{\dots}{Not currently used.}
}
\value{
  If HPD intervals are calculated, returns a matrix with
  the means, lower limit, and upper limit. If no HPD
  intervals are calculated, returns a column vector of
  means.
}
\description{
  If the predicted values only used the posterior means,
  highest posterior density (HPD) intervals cannot be
  generated, so only the means are returned. Otherwise, it
  calculates the mean predicted value, as well as the HPD
  interval.
}
\examples{
\dontrun{
  ## Make me!
}
}
\seealso{
  \code{\link{predict2.MCMCglmm}}, \code{\link{recycler}}
}

