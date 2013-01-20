#' Define a generic prediction function
#'
#' This defintes a generic \code{predict2} function which
#' is similar to the usual \code{predict} but can use different
#' methods. In particular, the \code{MCMCglmm} method has features
#' not available in the regular \code{predict} method for \code{MCMCglmm}
#' objects.
#'
#' @param object A model object to predict from
#' @param \dots Additional arguments passed to the methods
#' @export
#' @rdname predict2
#' @examples
#' # to see available methods
#' methods(predict2)
predict2 <- function (object, ...) {
  UseMethod("predict2", object)
}

#' Get predicted values from a \code{MCMCglmm} object
#'
#' This is the main workhorse of the package.
#' It is a \code{predict2} method for MCMCglmm objects.
#' There are a few core arguments. The model \code{object}
#' and design matrices, \code{X} (fixed effects) and \code{Z} (random effects).
#' If \code{X} and \code{Z} are missing, it will attempt to fill
#' them in from the model object (which optionally saves them).
#' If \code{X} and \code{Z} are specified or \code{NULL}, they are not used.
#' This is useful either for out of sample predictions
#' or to use just the fixed effects. Note that these must be
#' full design matrices, not data matrices.  For example, they
#' must dummy code factors and include the intercept (if there
#' was an intercept in the model).
#'
#' You can also use all posterior samples or just the mean.
#' All is nice because it lets you construct highest posterior
#' density (HPD) intervals around the predicted values, rather
#' than just get an estimate. The mean is nice because if
#' that is all you care about, it is much much faster.
#' You can get either the linear predictor values or the response scale.
#' However, response is currently only implemented for ordinal (probit) models.
#' Theoretically it could be extended but the code is a pain.
#'
# @param object A \code{MCMCglmm} model object to use for prediction.
#' @param X The fixed effects design matrix. Can be the original or new data.
#' @param Z The random effects design matrix. Can be the original or new data.
#' @param use A character string. Use just the posterior \dQuote{mean} or
#'   \dQuote{all} posterior samples (the default).
#' @param type A cahracter string. Either \dQuote{lp} for the linear predictor
#'   (the default) or \dQuote{response} for the predicted values on the response scale.
# @param \dots Not currently used.
#' @return Either a matrix of the linear predictor if \code{type = "lp"} or a
#'   list of class MCMCglmmPredictedProbs if \code{type = "response"}
#' @method predict2 MCMCglmm
#' @export
#' @rdname predict2
#' @seealso \code{\link{summary.MCMCglmmPredictedProbs}}, \code{\link{recycler}}
#' @examples
#' \dontrun{
#'   ## Make me!
#' }
predict2.MCMCglmm <- function(object, X, Z, use = c("all", "mean"),
  type = c("lp", "response"), ...) {

  use <- match.arg(use)
  type <- match.arg(type)

  if (!is.null(object$X) && missing(X)) {
    X <- object$X
  }

  if (!is.null(object$Z) && missing(Z)) {
    Z <- object$Z
  }

  if (missing(X)) X <- NULL
  if (missing(Z)) Z <- NULL

  Xb <- Zu <- 0L

  b <- as.matrix(fixef(object, use = use))
  u <- as.matrix(ranef(object, use = use))

  if (!is.null(X)) Xb <- X %*% b
  if (!is.null(Z)) Zu <-  Z %*% u

  res <- t(as.matrix(Xb + Zu))

  dimnames(res) <- list(switch(use,
    all = paste0("Rep.", 1:nrow(res)), mean = NULL), NULL)

  if (type == "response") {
    if (all(object$family %in% c("ordinal"))) {
      stopifnot(length(unique(object$error.term)) == 1)

      stddev <- sqrt(unique(object$error.term) + 1)

      # cut points
      CP <- object$CP
      if (use == "mean") CP <- colMeans(CP)

      CP <- as.list(as.data.frame(CP))
      CP <- lapply(CP, function(x) {do.call("cbind", rep(list(x), dim(res)[2L]))})
      CP <- c(-Inf, 0, CP, Inf)
      # difference between cuts and predicted
      CP <- lapply(CP, `-` res)

      q <- vector("list", length(CP) - 2)
      for (i in 2:(length(CP) - 1)) {
        q[[i - 1]] <- mcmc(pnorm(CP[[i + 1]], 0, stddev) - pnorm(CP[[i]], 0, stddev))
      }
      q <- c(list(mcmc(1 - Reduce(`+`, q[1:(i - 1)]))), q)
      class(q) <- c("list", "MCMCglmmPredictedProbs")
      res <- q
    } else {
      stop("Function does not support response type for families beside ordinal")
    }
  }
  return(res)
}

## as.list.matrix <- function(x, ...) {
##   d <- dim(x)
##   ncols <- d[2L]
##   ic <- seq_len(ncols)

##   value <- vector("list", ncols)
##   for (i in ic) value[[i]] <- as.vector(x[, i])
##   names(value) <- dimnames(x)[[2L]]
##   return(value)
## }

## y <- list(1:10000)
## system.time(test1 <- do.call("cbind", rep(y, 4000)))
## system.time(test2 <- matrix(y[[1]])[, rep(1, 4000)])
## system.time(test3 <- matrix(rep(y[[1]], 4000), ncol=4000))

#' Calculate change in predicted probabilities
#'
#' \code{recycler} wraps many of the functions in \pkg{postMCMCglmm} to calculate
#' the change in predicted probabilities for a twiddle change in the predictor,
#' or for discrete predictors, it can use values so it is the change from 0 to 1
#' (for example). The result is a MCMCglmmPredictedProbs
#' (of course a difference but still) object, so it can be summarized using
#' the \code{MCMCglmmPredictedProbs} \code{summary} method.
#' This gives average marginal recycled predicted probabilities,
#' as well as highest posterior density intervals.
#'
#' @param object A \code{MCMCglmm} model object to use for recycled predictions.
#' @param index An integer indicating the column of the fixed effects design matrix, X,
#'   to vary.  Defaults to 2L.
#' @param twiddle A twiddle value for continuous variables. Needs to be small enough
#'   for the scale of the predictor that a twiddle change is a reasonable approximation
#'   of taking the first derivative at a point. That is, a very small change.
#'   If missing, reverts to .01.
#' @param values Specific values to use for the varying predictor. These are primarily
#'   meant for discrete predictors rather than continuous ones.
#' @param \dots Passed on to \code{predict2}
#' @seealso \code{\link{summary.MCMCglmmPredictedProbs}}, \code{\link{predict2.MCMCglmm}}
#' @return A list of class \code{MCMCglmmPredictedProbs} that are the differences
#'   in predicted probabilities for a one unit change (calculated from the
#'   twiddle value or between the discrete values supplied in \code{values}).
#' @export
#' @examples
#' \dontrun{
#'   ## Make me!
#' }
recycler <- function(object, index = 2L, twiddle, values, ...) {
  L1 <- missing(twiddle)
  L2 <- missing(values)
  if (!L1 && !L2) stop("Please specify either a twiddle value or actual values, not both")
  X1 <- X2 <- object$X

  stopifnot(index %in% seq.int(ncol(X1)) || index %in% colnames(X1))

  if (L1 && L2) {
    twiddle <- .01
  }
  if (L2) {
    X2[, index] <- X2[, index] + twiddle
  } else {
    stopifnot(length(values) == 2L)
    X1[, index] <- values[1]
    X2[, index] <- values[2]
    twiddle <- diff(values)
  }

  p1 <- predict2(object, X = X1, use = "all", type = "response", ...)
  p2 <- predict2(object, X = X2, use = "all", type = "response", ...)


  pdelta <- lapply(seq_along(p2), function(i) {
    (p2[[i]] - p1[[i]])/twiddle
  })

  class(pdelta) <- c("list", "MCMCglmmPredictedProbs")

  return(pdelta)
}
