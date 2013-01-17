\name{jointPosterior}
\alias{jointPosterior}
\title{Plot the joint posterior from an MCMC}
\usage{
  jointPosterior(dat, x, y, xlab = "", ylab = "",
    probs = 0.95, plot = TRUE, topleftmargin = 0.2)
}
\arguments{
  \item{dat}{A data frame, typically of MCMC posterior
  samples.}

  \item{x}{A character string of the variable in \code{dat}
  to plot on the x-axis.}

  \item{y}{A character string of the variable in \code{dat}
  to plot on the y-axis.}

  \item{xlab}{A character string or expression containing
  the x-axis label. Expressions may be used for Greek or
  mathematical symbols. Defaults to blank.}

  \item{ylab}{A character string or expression containing
  the x-axis label. Expressions may be used for Greek or
  mathematical symbols. Defaults to blank.}

  \item{probs}{A numeric value of the probability to use
  for calculating the HPD region for the two dimensional
  density plot. Passed to \code{HPDregionplot} from the
  \pkg{emdbook} package.}

  \item{plot}{A logical whether the graph should be plotted
  or not.  Defaults to \code{TRUE}.}

  \item{topleftmargin}{A numeric value of the left margin
  of the upper univariate density plot. Need to adjust this
  to get the upper univariate density to line up with the
  bivariate density plot.}
}
\value{
  An invisible list of the univariate and bivariate density
  plots and the viewport used to put them all together.
  Primarily called for the side effect of creating a graph.
}
\description{
  A function to calculate the HPD region for a joint
  posterior, and plot the bivariate density with HPD
  region, and univariate densities on the margins using the
  \pkg{ggplot2} package.
}
\examples{
\dontrun{
  # sample data
  set.seed(10)
  dens2dtestdat <- as.data.frame(MASS::mvrnorm(4500, c(b1 = -.1, b2 = .05),
    Sigma = c(.05, .02)*matrix(c(1, -.5, -.5, 1), 2)*rep(c(.05, .02), each = 2)))
  d <- as.data.frame(mar2c$Sol[, 10:11]); colnames(d) <- c("b1", "b2")
  tmp <- as.data.frame(HPDregionplot(as.mcmc(d), n = 200)[[1]])
  jointPosterior(d, x = "b1", y = "b2", tmp, xlab = "Reactivity x Support", ylab = "Recovery x Support")

  # make the plot
  jointPosterior(dens2dtestdat, x = "b1", y = "b2", xlab = "Time x Constraint",
    ylab = bquote(Time^2 ~ x ~ Constraint))
  # clean up
  rm(dens2dtestdat)
}
}

