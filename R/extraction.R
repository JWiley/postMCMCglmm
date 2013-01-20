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
#' @seealso \code{\link{fixef.MCMCglmm}}, \code{\link{ranef.MCMCglmm}}
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
#' @seealso \code{\link{paramNamesMCMCglmm}}, \code{\link{ranef.MCMCglmm}}
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
#' @seealso \code{\link{fixef.MCMCglmm}}, \code{\link{ranef.MCMCglmm}}
#' @rdname extractEffects
.extractEffects <- function(object, use = c("all", "mean"),
  which = c("fixed", "random"), ...) {

  use <- match.arg(use)
  which <- match.arg(which)

  b <- as.matrix(object[["Sol"]])

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
      regex <- paste(paramNamesMCMCglmm(object)$random, collapse = "|")
      regex <- paste0("^(", regex, ")\\..*$")
      grep(regex, colnames(b), value = TRUE)
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
#' @method fixef MCMCglmm
#' @export
#' @export fixef
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
#' @method ranef MCMCglmm
#' @export
#' @export ranef
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



## stdranef(m[[1]], which = list(1), type = "lp")
## stdranef(m[[1]], which = list(2), type = "lp")
## stdranef(m[[1]], which = list(1, 2, c(1, 2)), type = "lp")
## stdranef(m[[1]], type = "lp")

## ## error
## #junk <- stdranef(m[[1]], which = list(1, 2, 3), type = "lp")

## stdranef(m[[1]], which = list("mrn", "pager"), type = "lp")

## # this does not work, check zero setting
## junk <- stdranef(m[[1]], type = "response")

stdranef <- function(object, which, type = c("lp", "response")) {
  type <- match.arg(type)

  if (is.null(object$Z)) stop("Z matrix must be saved")
  z <- object$Z

  re <- paramNamesMCMCglmm(object)$random

  if (missing(which)) which <- c(re, list(re))

  stopifnot(is.list(which))

  if (is.numeric(unlist(which))) {
    stopifnot(all(unlist(which) %in% seq_along(re)))
    which <- lapply(which, function(i) re[i])
  } else {
    stopifnot(all(unlist(which) %in% re))
  }

  index <- lapply(which, function(n) {
    n <- paste(n, collapse = "|")
    regex <- paste0("^(", n, ")\\..*$")
    index <- grep(regex, colnames(z))
    return(index)
  })

  # tricky, note that coefficients and predicted probabilities
  # do not come out with the same dimensions, the are transposed
  # samples on columns for coefs, samples on rows for predicted probs
  res <- switch(type,
    lp = {
      yhat <- lapply(index, function(i) ranef(object, use = "all")[i, ])
      lapply(yhat, function(m) matrix(apply(m, 2, var)))
    }, response = {
      yhat <- lapply(index, function(i) {
        tmp <- z
        tmp[, -i] <- 0L # bug here??
        predict2(object, X = NULL, Z = tmp, use = "all", type = type)
      })
      lapply(yhat, function(m) sapply(m, function(n) apply(n, 1, var)))
    })

  names(res) <- which

  M <- do.call(rbind, lapply(res, colMeans))

  finalres <- list(M = M, Data = res)
  class(finalres) <- "postMCMCglmmRE"

  return(finalres)
}

print.postMCMCglmmRE <- function(x, ...) {
  x <- x$M
  NextMethod(print, object = x)
}

