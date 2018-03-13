# Copyright (c) 2015-2018 Hadrien Chauvin
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#####################

#' Model selection based on an \code{\link{fd_comb}} object.
#'
#' Wrappers around the functions of package \pkg{AICcmodavg} and extensions
#' inspired by it specifically tailored to nonlinear optimization.
#'
#' @param second.ord Whether to use AIC (\code{FALSE}) or AICc (\code{TRUE}).
#'   We discourage the use of AICc as this has been insufficiently dealt with
#'   within the statistical literature on generalized estimating equations.  If
#'   \code{TRUE}, the number of observations used for the small-sample
#'   adjustment is the number of clusters.  This clearly makes the adjustment
#'   bigger than it should be.
#' @param data The \code{\link{fd_data}} dataset on which optimization was
#'   performed.
#' @param ... (Potentially named) \code{\link{fd_nls}}/\code{link[stats]{nls}}
#'   objects.
#' @param .list Alternatively or additionally, a (potentially named) list of
#'   \code{\link{fd_nls}}/\code{link[stats]{nls}} objects can be provided.
#' @param vcov. Variance-covariance function to use on the free parameters. For
#'   example, the \pkg{sandwich} package offers an alternative to
#'   \code{\link[stats]{vcov}} and could be used for a Huber-White estimation.
#'   Alternatively, variance-covariance matrix.
#' @param object An \code{fd_modavg} object.
#' @param tab An \code{fd_aictab} object.
#' @param subset Used to subset \code{tab} the same way the \code{subset}
#'   argument of \code{\link[base]{subset}} works.
#' @param pos If \code{1}, returns the best model.  If \code{2}, returns the
#'   second-best model, etc.
#' @param verbose whether the name of the chosen model should be explicitly
#'   printed.
#'
#' @return \code{fd_aictab} returns a \code{data.frame} (see \emph{details}),
#' \code{fd_modavg} an \code{fd_modavg} object, \code{best} a "fit" object
#' (e.g., an \code{\link{fd_nls}} object).
#'
#' @details \code{fd_aictab} returns a \code{data.frame} of various summary
#' statistics concerning the models, ranked by increasing AIC or AICc.  The
#' output builds on the output of \code{\link[AICcmodavg]{aictab}}, with two
#' additional columns: \code{Effect} and \code{Q} (see the section on "effect
#' size" below for interpretation).
#'
#' The \code{fd_modavg} function works similarly to
#' \code{\link[AICcmodavg]{modavg}} function. It is based on a conceptualization
#' of AIC values as weights (Burnham and Anderson 2004). The weighting is done
#' both on the estimate and on the variance-covariance matrix.  \code{coef},
#' \code{vcov} and \code{weight} can be used to extract, respectively, the
#' estimated model-averaged coefficients, the variance-covariance matrix and the
#' weights used for the averaging.  Using QIC for correcting for variance
#' structure misspecification (Pan 2001) would have been better, but
#' implementation in R is sketchy and this is a fundamental matter beyond the
#' scope of this work.
#'
#' \code{best} returns the best (or second-best, ...) model of an
#' \code{fd_aictab} object.
#'
#' @section Effect size: The additional columns \code{Effect} and \code{Q} are
#'   to be interpreted as follows.
#'
#'   We think of a parameter as practically identifiable when its 95\%
#'   confidence interval does not include its lower or upper bounds.  This
#'   translates into a normalized effect size of roughly \eqn{2} (\eqn{p}-value
#'   < \eqn{0.05}, \code{*}).  Effect sizes greater than \eqn{2.6} (\eqn{p <
#'   0.01}, \code{**}) or even \eqn{3.3} (\eqn{p < 0.001}, \code{***}) are even
#'   more desirable. We think of a model as practically identifiable when the
#'   effect size of each of its parameters that have not reached the lower
#'   bounds is at least 2, and in general we quantify the practical
#'   identifiability of a model with the minimum of the effect sizes (we simply
#'   call it the \code{Effect} of the model).  A qualitative appraisal of this
#'   effect is given by \code{Q}.  Its value is either \code{""} (\code{Effect <
#'   2}), \code{"*"}, \code{"**"} or \code{"***"} depending on the minimum
#'   effect size.
#'
#' @name model-selection
#' @examples
#' control <- list(maxit=1)
#'
#' \dontrun{
#' control <- NULL
#' }
#'
#' data(FdSTyphimuriumWTC57)
#' dat <- cutoff(FdSTyphimuriumWTC57)
#' source(system.file("contrib", "nlsSA.R", package="fluodilution"))
#' CC <<- FdCommonConstraints
#' comb <- fd_comb(catcstr(attr(dat, "bestcstr"),
#'                         list(CC$`#delta_1100`, CC$`#delta_1111`)),
#'                 data=dat,
#'                 nlsfun="nlsSA",
#'                 control=control)
#'
#' tab <- fd_aictab(data = dat, .list = comb)
#' avg <- fd_modavg(tab)
#' invisible(vcov(avg))
#' weights(avg)
#' best(tab)
#'
#' @references Burnham KP, Anderson DR (2004). Multimodel inference:
#' understanding AIC and BIC in model selection. \emph{Sociological Methods
#' Research} \strong{33}: 261-304.
#'
#' Pan W (2001). Akaike's Information Criterion in Generalized Estimating
#' Equations. \emph{Biometrics} \strong{57} (1): 120-125.
#'
#' @section Required packages: Please install \pkg{AICcmodavg} for
#'   \code{fd_aictab}.
NULL

#' @rdname model-selection
#' @export
fd_aictab <- function (second.ord=FALSE, data, ..., .list = NULL) {
  nobs <- nlevels(data$Individual)

  arglst <- append(.list, list(...))
  if (length(arglst) == 0) stop("no model to select from")
  if (inherits(arglst[[1L]], "nls"))
    mdlsel <- arglst
  else {
    mdlsel <- Reduce(append, lapply(names(arglst), function (nm) {
      if (is.null(arglst[[nm]]) || length(arglst[[nm]]) == 0) {
        warning(
          "fd_aictab: model group '", nm, "' is null or empty ",
          "and won't be taken into account")
        NULL
      } else {
        setNames(arglst[[nm]], paste0(nm, ".", names(arglst[[nm]])))
      }
    }))
    if (length(mdlsel) == 0) stop("no model to select from")
  }
  if (is.null(names(mdlsel)))
    names(mdlsel) <- paste0("mod_", seq_along(mdlsel))

  excluded <- sapply(mdlsel, function (cur) {
    if (inherits(cur, "try-error") || is.null(cur$.flat)) {
      TRUE
    } else cur$.flat
  })
  mdlsel <- mdlsel[!as.logical(excluded)]

  effects_comb <- sapply(
    mdlsel,
    function (x) min((coef(x) - lower(model(x))) /
               sqrt(diag(vcov(x)))))

  selobj <- lapply(mdlsel, function (x) {
    class(x) <- "nls"
    return (x)
  })

  if (!requireNamespace("AICcmodavg", quietly=TRUE)) {
    stop("'AICcmodavg' required")
  }
  ans <- AICcmodavg::aictab(selobj, second.ord=second.ord, nobs=nobs)

  ans <- transform(transform(ans,
       Effect = effects_comb[as.character(Modnames)]),
             Q = cut(Effect, breaks=c(0, 2, 2.6, 3.3, Inf),
                   labels=c("", "*", "**", "***"),
                 include.lowest=TRUE))
  if (!is.null(nobs) & second.ord) {
    selected <- ans$K >= nobs
    ans$AICc[selected] <- ans$Delta_AICc[selected] <- Inf
  }
  ans <- ans[, c(1, 2, 3, 7, 9, 10)]
  return (structure(invisible(ans),
            fits=mdlsel,
            excluded=names(which(excluded)),
            class=c("fd_aictab", "data.frame")))
}

#' @rdname model-selection
#' @export
fd_modavg <- function (tab, vcov. = vcov) {
  if (NROW(tab) == 0) {
    return (NULL)
  }

  weights <- fd_weights(tab)

  models <- attr(tab, "fits")
  mdlsel <- match(as.character(tab$Modnames), names(models))

  all_coefs <- sapply(mdlsel,
       function (i) unlist(relisted_coef(models[[i]])))
  if (is.list(all_coefs))
    stop("relisted coefficients do not share a common structure")
  all_coefs <- t(all_coefs)
  avg_coef <- colSums(all_coefs * weights)

  avg_vcov <- cov(t(t(all_coefs) - avg_coef) * sqrt(weights))+
    Reduce("+", lapply(
      seq_along(weights),
      function (i)
        weights[[i]] * relisted_vcov(models[[mdlsel[i]]], vcov.)))

  structure(list(coefficients = avg_coef, vcov = avg_vcov,
           stddev = sqrt(diag(avg_vcov)),
           weights = weights,
           all_coefs = all_coefs),
        class=c("fd_modavg"))
}

#' @export
#' @rdname model-selection
vcov.fd_modavg <- function (object, ...) {
  object$vcov
}

#' @rdname model-selection
#' @export
weights.fd_modavg <- function (object, ...) {
  object$weights
}

#' @rdname model-selection
#' @export
best <- function (tab, subset=TRUE, pos=1, verbose = TRUE) {
  modname <- as.character(tab[eval(substitute(subset), tab), ]$Modnames[pos])
  if (verbose) {
    message("Best model:", modname)
  }
  ans <- attr(tab, "fit")[[modname]]
  ans$.picked <- modname
  ans
}

#' @export
#' @keywords internal
`[.fd_aictab` <- function (x, i, j, ...) {
  ans <- NextMethod()
  attr(ans, "fits") <- attr(x, "fits")[as.character(ans$Modnames)]
  attr(ans, "excluded") <- attr(x, "excluded")
  ans
}

#' @export
#' @keywords internal
fd_weights <- function (tab) {
  weights <- tab$LL - tab$K
  weights <- exp(weights - max(weights))
  setNames(weights / sum(weights), as.character(tab$Modnames))
}
