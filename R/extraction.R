#' Extract the parameter names from an \code{MCMCglmm} object
#'
#' Simple function to extract the fixed and random effects
#' parameter names from an \code{MCMCglmm} object.
#'
#' @param object An \code{MCMCglmm} object
#' @param \dots not used
#' @return A list with two elements:
#'   \item{fixed}{A character vector of the fixed effects parameter names}
#'   \item{random}{A character vector of the random effects parameter names}
#' @export
#' @seealso \code{\link{fixef.MCMCglmm}} \code{\link{ranef.MCMCglmm}}
#' @examples
#' \dontrun{
#'   # a simple MCMCglmm model
#'   data(PlodiaPO)
#'   m <- MCMCglmm(PO ~ 1, random = ~ FSfamily, data = PlodiaPO, verbose=FALSE)
#'
#'   # extract the parameter names
#'   paramNamesMCMCglmm(m)
#' }
paramNamesMCMCglmm <- function(object, ...) {
  fNames <- as.character(attr(terms(object$Fixed$formula), "variables"))[-c(1, 2)]

  rNames <- NULL

  if (!is.null(object$Random$formula)) {
    rNames <- as.character(attr(terms(object$Random$formula), "variables"))[-c(1)]
  }

  if ("(Intercept)" %in% colnames(object[["Sol"]])) {
    fNames <- c("(Intercept)", fNames)
  }

  list(fixed = fNames, random = rNames)
}

#' Extract the levels of factors used for random effects in \code{MCMCglmm} objects
#'
#' @param object An \code{MCMCglmm} model object
#' @param data The dataset used for the model
#' @param \dots Not currently used
#' @export
#' @seealso \code{\link{paramNamesMCMCglmm}} \code{\link{ranef.MCMCglmm}}
#' @examples
#' \dontrun{
#'   # a simple MCMCglmm model
#'   data(PlodiaPO)
#'   m <- MCMCglmm(PO ~ 1, random = ~ FSfamily, data = PlodiaPO, verbose=FALSE)
#'
#'   # extract the random effects levels
#'   ranefLevels(m, PlodiaPO)
#' }
ranefLevels <- function(object, data, ...) {
  n <- paramNamesMCMCglmm(object)$random
  res <- lapply(n, function(n) {
    levels(data[, n])
  })
  names(res) <- n
  return(res)
}

#' Internal function to extract the fixed or random effects from an \code{MCMCglmm} object
#'
#' Extracts the fixed or random effects portions from an MCMCglmm object.
#' Note for the random, these are the estimates themselves,
#' not the variability in the estimates. The \code{use} options let you
#' get either just the posterior mean or all the posterior samples.
#'
#' @param object An \code{MCMCglmm} object
#' @param use A character string indicating whether to return \dQuote{all} the
#'   posterior samples (the default) or only the \dQuote{mean} of them.
#' @param which A character string indicating whether to return the
#'   \dQuote{fixed} (the default) or \dQuote{random} effects.
#' @param \dots Not currently used.
#' @return A matrix of the posterior samples or means for the fixed or random effects.
#' @keywords internal
#' @seealso \code{\link{fixef.MCMCglmm}} \code{\link{ranef.MCMCglmm}}
#' @rdname extractEffects
.extractEffects <- function(object, use = c("all", "mean"),
  which = c("fixed", "random"), ...) {

  use <- match.arg(use)
  which <- match.arg(which)

  b <- object[["Sol"]]

  eff <- switch(which,
    fixed = {
      # this does not work because the intercept term contains
      # special characters (), that screw with the regular expressions
      # cannot use fixed matching because factor variables can
      # expand with arbitrarily named levels
      #unlist(lapply(paramNamesMCMCglmm(object)$fixed, function(n) {
      #  grep(paste0("^", n, ".*$"), colnames(b), value = TRUE)
      #}))
      object$X@Dimnames[[2]]
    },
    random = {
      unlist(lapply(paramNamesMCMCglmm(object)$random, function(n) {
        grep(paste0("^", n, "\\..*$"), colnames(b), value = TRUE)
      }))
    }
  )

  b <- b[, eff, drop=FALSE]

  switch(use,
    all = t(b),
    mean = as.matrix(colMeans(b)))
}

#' Extract fixed effects from an \code{MCMCglmm} object
#'
#' Function designed to extract the fixed effects from an
#' \code{MCMCglmm} model object. Can either extract all samples from the
#' fixed effects posteriors or return the posterior means.
#'
#' @param object An \code{MCMCglmm} model object to extract the effects from
#' @param use A character string indicating whether to extract
#'   all posterior samples or the mean of the posteriors. Defaults to
#'   "all".
#' @param \dots Arguments passed on to the worker function.
#' @return A matrix of the fixed effects
#' @importFrom nlme fixef
#' @S3method fixef MCMCglmm
#' @seealso \code{\link{ranef.MCMCglmm}}
#' @examples
#' \dontrun{
#'   # a simple MCMCglmm model
#'   data(PlodiaPO)
#'   m <- MCMCglmm(PO ~ 1, random= ~ FSfamily, data=PlodiaPO, verbose=FALSE)
#'
#'   # only extract average fixed effects
#'   fixef(m, use = "mean")
#'
#'   # histogram of posterior samples of fixed effects
#'   hist(fixef(m))
#'   # matches the mean
#'   rowMeans(fixef(m))
#' }
fixef.MCMCglmm <- function(object, use = c("all", "mean"), ...) {
  .extractEffects(object = object, use = use, which = "fixed", ...)
}

#' Extract random effects from an \code{MCMCglmm} object
#'
#' Function designed to extract the random effects from an
#' \code{MCMCglmm} model object. Can either extract all samples from the
#' random effects posteriors or return the posterior means.
#'
#' @param object An \code{MCMCglmm} model object to extract the effects from
#' @param use A character string indicating whether to extract
#'   all posterior samples or the mean of the posteriors. Defaults to
#'   "all".
#' @param \dots Arguments passed on to the worker function.
#' @return A matrix of the fixed effects
#' @importFrom nlme ranef
#' @S3method ranef MCMCglmm
#' @seealso \code{\link{fixef.MCMCglmm}}
#' @examples
#' \dontrun{
#'   # a simple MCMCglmm model
#'   data(PlodiaPO)
#'   m <- MCMCglmm(PO ~ 1, random= ~ FSfamily, data=PlodiaPO, pr=TRUE, verbose=FALSE)
#'
#'   # only extract average fixed effects
#'   head(ranef(m, use = "mean"))
#'
#'   # histogram of posterior samples of fixed effects
#'   hist(ranef(m)[1, ])
#'   # matches the mean
#'   rowMeans(ranef(m)[1:6, ])
#' }
ranef.MCMCglmm <- function(object, use = c("all", "mean"), ...) {
  .extractEffects(object = object, use = use, which = "random", ...)
}
