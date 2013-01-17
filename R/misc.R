#' Summary method for MCMCglmm predicted probabilities
#'
#' Optionally first marginalizes across all observations
#' by taking the row means. If the predicted values only used the posterior means,
#' highest posterior density (HPD) intervals cannot be generated, so only
#' the means are returned. Otherwise, it calculates the mean predicted probability,
#' as well as the HPD interval. This can either be per observation or marginalized.
#'
#' @param object A \code{MCMCglmmPredictedProbs} object to summarize
#' @param marginalize Logical whether or not to marginalize by taking the row means. Defaults to \code{FALSE}.
#' @param level A numeric value, the value to use when calculating HPD intervals. Defaults to .95.
#' @param \dots Not currently used.
#' @return If HPD intervals are calculated, returns a list of matrices
#'   with the means, lower limit, and upper limit. If no HPD intervals
#'   are calculated, returns a list of means.
#' @method summary MCMCglmmPredictedProbs
#' @export
#' @seealso \code{\link{predict2.MCMCglmm}}, \code{\link{recycler}}
#' @examples
#' \dontrun{
#'   ## Make me!
#' }
summary.MCMCglmmPredictedProbs <- function(object,
  marginalize = FALSE, level = .95, ...) {

  if (nrow(object[[1]]) > 1) {
    intervals <- TRUE
  } else {
    intervals <- FALSE
  }
  if (marginalize) {
    object <- lapply(object, rowMeans)
  }

  if (intervals) {
    res <- lapply(object, function(x) {
      if (!is.null(ncol(x))) {
        res <- cbind(colMeans(x), HPDinterval(mcmc(x), prob = level))
      } else {
        res <- cbind(mean(x), HPDinterval(mcmc(x), prob = level))
      }
      colnames(res) <- c("M", "LL", "UL")
      return(res)
    })
  } else {
    res <- lapply(object, mean)
  }
  return(res)
}

# simple little function to generate confusion matrices
# I use it for one column predicted class, one column actual
# then look at percentage correctly classified
confusion <- function(formula, data) {
  res <- xtabs(formula = formula, data = data)
  print(res)
  res/sum(res)
}


#' Plot the joint posterior from an MCMC
#'
#' A function to calculate the HPD region for a joint posterior,
#' and plot the bivariate density with HPD region, and univariate
#' densities on the margins using the \pkg{ggplot2} package.
#'
#' @param dat A data frame, typically of MCMC posterior samples.
#' @param x A character string of the variable in \code{dat} to plot on
#'   the x-axis.
#' @param y A character string of the variable in \code{dat} to plot on
#'   the y-axis.
#' @param xlab A character string or expression containing the
#'   x-axis label. Expressions may be used for Greek or mathematical
#'   symbols. Defaults to blank.
#' @param ylab A character string or expression containing the
#'   x-axis label. Expressions may be used for Greek or mathematical
#'   symbols. Defaults to blank.
#' @param probs A numeric value of the probability to use for calculating
#'   the HPD region for the two dimensional density plot. Passed to
#'   \code{HPDregionplot} from the \pkg{emdbook} package.
#' @param plot A logical whether the graph should be plotted or not.
#'   Defaults to \code{TRUE}.
#' @param topleftmargin A numeric value of the left margin of the
#'   upper univariate density plot. Need to adjust this to get the
#'   upper univariate density to line up with the bivariate density
#'   plot.
#' @return An invisible list of the univariate and bivariate density plots
#'   and the viewport used to put them all together. Primarily called for the
#'   side effect of creating a graph.
#' @importFrom grid unit grid.newpage pushViewport viewport grid.layout
#' @importFrom emdbook HPDregionplot
#' @examples
#' \dontrun{
#'   # sample data
#'   set.seed(10)
#'   dens2dtestdat <- as.data.frame(MASS::mvrnorm(4500, c(b1 = -.1, b2 = .05),
#'     Sigma = c(.05, .02)*matrix(c(1, -.5, -.5, 1), 2)*rep(c(.05, .02), each = 2)))
#'   d <- as.data.frame(mar2c$Sol[, 10:11]); colnames(d) <- c("b1", "b2")
#'   tmp <- as.data.frame(HPDregionplot(as.mcmc(d), n = 200)[[1]])
#'   jointPosterior(d, x = "b1", y = "b2", tmp, xlab = "Reactivity x Support", ylab = "Recovery x Support")
#'
#'   # make the plot
#'   jointPosterior(dens2dtestdat, x = "b1", y = "b2", xlab = "Time x Constraint",
#'     ylab = bquote(Time^2 ~ x ~ Constraint))
#'   # clean up
#'   rm(dens2dtestdat)
#' }
jointPosterior <- function(dat, x, y, xlab = "", ylab = "",
  probs = .95, plot=TRUE, topleftmargin = .2) {

  dat <- dat[, c(x, y)]

  cdat <- HPDregionplot(x = as.mcmc(dat), vars = 1:2, n = 200, prob = probs)
  cdat <- as.data.frame(cdat[[1]])

  pb1 <- ggplot(dat, aes_string(x = x)) + geom_density() +
    geom_rug() + theme_bw() +
    theme(axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = unit(c(.5, 0, 0, topleftmargin), "lines"))

  pb2 <- ggplot(dat, aes_string(x = y)) + geom_density() +
    geom_rug() + coord_flip() + theme_bw() +
    theme(axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = unit(c(0, .5, 0, 0), "lines"))

  pbivar <- ggplot(dat, aes_string(x = x, y = y)) +
    geom_vline(aes(x = 0), size = 1) +
    geom_hline(aes(y = 0), size = 1) +
    geom_point(alpha = .05) +
    geom_path(aes(x = x, y = y), cdat, colour = "blue", size = 1.5) +
    labs(x = xlab, y = ylab) +
    theme_bw() +
    theme(plot.margin = unit(rep(0, 4), "lines"))

  vp <- viewport(layout = grid.layout(2, 2,
    widths = unit(c(1.6, .4), "null"),
    heights = unit(c(.4, 1.6), "null")))

  if (plot) {
    grid.newpage()
    pushViewport(vp)

    print(pb1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(pb2, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
    print(pbivar, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  }

  bigres <- list(hist1 = pb1, hist2 = pb2, bivar = pbivar, viewport = vp)
  return(invisible(bigres))
}
