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

#' Minus log-likelihood of a fluorescence dilution model.
#'
#' This function returns a minus log-likelihood function to be used with general
#' maximum likelihood estimation packages such as \pkg{stats4} or \pkg{bbmle}.
#'
#' @param model An \code{\link{fd_model}} object.
#' @param data An \code{\link{fd_data}} dataset.
#' @param start Optional starting parameters.  The default behaviour
#'   (\code{NULL}) is to use the starting parameter provided by
#'   \code{\link{start}(model)}.
#' @param mgen Maximum number of generations to fit.  If \code{NULL}, what the
#'   proliferation model provides by default is used.
#' @param verbose If \code{TRUE} (default), print additional diagnostic
#'   messages.
#' @param weights A numeric vector of length equal to the number of rows of
#'   \code{data}. Gives the relative weights to put on the residuals.  Default
#'   (\code{NULL}) to equal weight.
#' @param stop_boundary Whether \code{\link{fd_predict}} should raise an error
#'   if the parameters cross the boundary.  By default (\code{FALSE}), just
#'   returns a very large number when it happens.
#'
#' @return Minus log-likelihood (same as \code{\link[stats]{logLik}} applied to
#'   an \code{\link{nls}} object).
#' @family optimization-related entries
#' @export
#'
#' @examples
#' library(stats4)
#' library(bbmle)
#'
#' # Create minuslogl function
#' data(FdSTyphimuriumWTC57)
#' dat <- subset(cutoff(FdSTyphimuriumWTC57), Individual == "140508.WT.C57")
#' mdl <- fd_model(dat, fmm="gaussian", proliferation="branching",
#'                 constraints = attr(dat, "bestcstr"))
#' fun <- fd_minuslogl(mdl, dat, verbose=FALSE)
#' print(fun)
#'
#' # This this function on lower parameters
#' do.call(fun, as.list(start(mdl)))
#'
#' # Optimize this minuslogl.  Notice that the 'start' parameter is not
#' # explicitly given as the arguments to 'fun' have default values
#' fun <- fd_minuslogl(mdl, dat, verbose=TRUE)
#' control <- list(maxit=1)
#' \dontrun{
#' control <- NULL
#' }
#' stats4::mle(fun, method="Nelder-Mead", control=control)
#'
#' # The user is advised to use GenSA for global search:
#' bbmle::mle2(fun, optimfun=GenSA::GenSA, control=control)
fd_minuslogl <- function (model, data, start = NULL, mgen = "guess",
                          verbose = TRUE,
                          weights = NULL,
                          stop_boundary = FALSE) {
  if (verbose == TRUE) verbose <- 5
  if (mgen == "guess") {
    mgen <- guess_mGen(data)
    if (verbose > 0) cat("mgen guessed: ", mgen, "\n")
    if (mgen > 20) {
      warning("guessed 'mgen' greater than 16, capped at 16")
      mgen <- 16
    }
  }
  if (is.null(mgen)) mgen <- "NULL"
  if (is.null(start)) start <- match.fun("start")(model)
  if (verbose > 0) {
    cat("Starting values:\n")
    print(rbind(start = start,
                lower = lower(model),
                upper = upper(model)))
  }

  # One first run to make sure everything is alright (also useful
  # when partial != NULL)
  invisible(fd_predict(model, start(model), data, mgen = mgen))

  N <- NROW(data)
  if (is.null(weights)) weights <- rep_len(1, N)
  zw <- weights == 0

  # nolint start
  shift <- -N * (log(2 * pi) + 1 - log(N) - sum(log(weights + zw))) / 2.0
  # (used in an eval)
  # nolint end

  report <- local({
    lastReportTime <- proc.time()[3]
    bestAns <- NA
    function (ans) {
      if (verbose & proc.time()[3] > lastReportTime + verbose) {
        lastReportTime <<- proc.time()[3]
        if (!is.finite(bestAns) || (is.finite(bestAns) && ans < bestAns)) {
          bestAns <<- ans
        }
        cat(sep="", "*** Best: ", bestAns,
          " (Current: ", ans, ") ***\n")
      }
      ans
    }
  })

  # Avoid expensive calls to fd_data
  data <- as.data.frame(data)

  eval(parse(text=paste0("function (\n",
       paste(paste("  ", names(start), "=", start),
         collapse=",\n"), "\n) {\n",
       "  param <- c(\n",
       paste(paste("    ", names(start), "=", names(start)),
         collapse=",\n"),
       "\n  )\n",
       "  pred <- fd_predict(model, param, data, mgen = ", mgen, ", ",
       "stop_boundary = ", stop_boundary, ")\n",
       "  pred$fitted <- pred$y\n",
       "  pred$y <- NULL\n",
       "  pred <- merge(pred, data)\n",
       "  report(-(-N * log(sum(weights * (pred$y - pred$fitted)^2)) /\n",
       "    2 + shift))\n",
       "}")))
}

#' @title Ordinary Least Square (OLS) optimization of a fluorescence dilution
#'   model.
#'
#' @description \code{fd_nls} acts as a wrapper around an
#' \code{\link[stats]{nls}}-like optimizer with additional checks and
#' capabilities.  Other functions support the activity of \code{fd_nls}.
#'
#' @param globalname The name of an \code{\link{fd_model}} object sitting in the
#'   global environment (that is, a \code{\link{character}} object).  See
#'   \code{\link{fd_model-functions}} for the rationale.
#' @param data The (experimental) \code{\link{fd_data}} dataset.
#' @param newdata The new \code{\link{fd_data}} dataset for which to make
#'   predictions. Notice that, contrary to \code{\link[stats]{predict.nls}},
#'   here \code{newdata} cannot be missing.
#' @param mgen The maximum number of generations (integer).  If \code{NULL}, the
#'   default of the proliferation model is considered.  If \code{"guess"},
#'   guessed through \code{guess_mGen}.
#' @param loop If \code{mgen} is not enough (there is a lot of cells in the
#'   maximum number of generations), number of additional optimization rounds
#'   that \code{fd_nls} can undertake with ever larger \code{mgen}.  This
#'   parameter is important for automatic combing, e.g. with
#'   \code{\link{fd_comb}}.
#' @param trace Whether to output additional trace information.
#' @param control Control parameters to pass along to the optimizer.
#' @param start Alternative starting coefficients to \code{\link{start}}.
#' @param nlsfun A user-specified optimization function.  Must (at least) accept
#'   the same arguments as \code{\link[stats]{nls}} (default).  The
#'   "contributed" \code{nlsSA} function can also be used (see details).
#' @param object An \code{fd_nls}, \code{\link[stats]{nls}},
#'   \code{\link[nlme]{nlsList}} or \code{\link[nlme]{nlme}} object.
#' @param drop Whether to return the relisted parameters directly if there is
#'   only one set of coefficients or a one-element list containing the relisted
#'   parameters.
#' @param vcov. Variance-covariance function to use on the free parameters. For
#'   example, the \pkg{sandwich} package offers an alternative to
#'   \code{\link[stats]{vcov}} and could be used for a Huber-White estimation.
#'   Alternatively, variance-covariance matrix to use.
#' @param n Number of Monte-Carlo simulations to run.
#' @param model An \code{\link{fd_model}} object.  Must be specified if
#'   \code{object} is not an \code{fd_nls} object (i.e. \code{model(object) ==
#'   NULL}).
#' @param ... For \code{fd_nls}, additional arguments to pass along to the
#'   optimizer. For \code{fd_residuals}, the \code{fd_nls},
#'   \code{\link[stats]{nls}}, \code{\link[nlme]{nlsList}} or
#'   \code{\link[nlme]{nlme}} objects to find the residuals of (the arguments
#'   have to be named).  For \code{relisted_fit}, additional parameters to pass
#'   along to \code{relisted_vcov}. Not used by \code{predict} or \code{vcov}.
#' @param .list Named list of additional \code{fd_nls},
#'   \code{\link[stats]{nls}}, \code{\link[nlme]{nlsList}} or
#'   \code{\link[nlme]{nlme}} objects.
#'
#' @return \code{fd_nls} returns an \code{fd_nls} object, \code{model} an
#'   \code{fd_model} object, \code{fd_residuals} returns a \code{data.frame},
#'   \code{relisted_coef} returns a list of structured parameter values,
#'   \code{relisted_vcov} returns a variance-covariance matrix,
#'   \code{relisted_fit} an object with which \code{coef} and \code{vcov} can be
#'   used (see \emph{details}).
#'
#' @details \code{fd_nls} returns an \code{fd_nls} object, inheriting \code{nls}
#' (the return class of \code{\link[stats]{nls}}), augmented by model
#' information.  Therefore, the functions that can be used on an \code{nls}
#' object can be used on an \code{fd_nls} object, such as \code{coef},
#' \code{predict} or \code{summary}.
#'
#' \code{model} returns the \code{fd_model} object that was used to perform the
#' optimization. \code{fd_residuals} returns \code{data} with the content of the
#' \code{y} column replaced by the predicted values.  \code{predict} returns a
#' numeric vector of the predicted values for the \code{\link{fd_data}} dataset
#' \code{newdata}.  \code{relisted_coef} is a wrapper around
#' \code{\link{relist.fd_model}} and return the relisted optimal coefficients.
#' In some cases, e.g. \code{\link[nlme]{nlme}} and \code{\link[nlme]{nlsList}},
#' instead of a vector, \code{coef} returns a matrix with the coefficients in
#' column and the group levels in rows. In this case, \code{relisted_coef}
#' returns a list of relisted coefficients.  If \code{drop=FALSE}, in the case
#' only one set of coefficients is returned, a one-element list of relisted
#' coefficients is returned instead of the relisted coefficients directly.
#'
#' \code{relisted_vcov} uses a Monte Carlo simulation, dependent on package
#' \pkg{mvtnorm}, to get to the variance-covariance matrix of the relisted
#' parameters from the variance-covariance matrix of the free parameters
#' returned by \code{vcov} (using singular value decomposition).  The rows and
#' columns are ordered as per\verb{
#' }\code{names(unlist(relisted_coef(object)))}.
#'
#' \code{relisted_fit} calls both \code{relisted_coef} and \code{relisted_vcov}
#' and returns an object with which \code{coef} and \code{vcov} can be used.
#'
#' @section Remark: Instead of \code{nls}, \code{GenSA} can be used for
#'   optimization through the "contributed" wrapper \code{nlsSA}.  \code{nlsSA}
#'   can be sourced with \verb{
#' }\code{source(system.file("contrib", "nlsSA.R",
#'   package="fluodilution"))}.
#'
#' @export
#' @family optimization-related entries
#'
#' @examples
#' data(FdSTyphimuriumWTC57)
#' dat <- cutoff(FdSTyphimuriumWTC57)
#' mdl <<- fd_model(data = dat, constraints = attr(dat, "bestcstr"))
#' guess_mGen(dat)
#' control <- list(maxit = 1L)
#' \dontrun{
#' control <- NULL
#' }
#'
#' \dontrun{
#' fit <- fd_nls("mdl", dat,
#'               algorithm="port", control=control)
#' # Error in match.fun(nlsfun)(fd_formula(globalname, mgen = mgen), data,  :
#' #     Convergence failure: iteration limit reached without convergence (10)
#' }
#'
#' # Alternatively, GenSA can be used through nlsSA:
#' source(system.file("contrib", "nlsSA.R", package="fluodilution"))
#' fit <- fd_nls("mdl", dat, nlsfun=nlsSA, control=control)
#'
#' \dontrun{
#' # This gives a summary for free parameters
#' print(summary(fit))
#' }
#'
#' # To go back to the unconstrained parameters, one can use
#' # relisted_fit or relisted_coef/relisted_vcov separately
#' relfit <- relisted_fit(fit)
#' invisible(relisted_coef(fit))
#' invisible(relisted_vcov(fit))
#'
#' # It is also possible to use a sandwich estimator instead
#' relfit <- relisted_fit(fit, vcov. = sandwich::sandwich)
#' print(lmtest::coeftest(relfit))
#'
#' @section Required packages: Please install \pkg{mvtnorm} for
#'   \code{relisted_vcov} and \code{relisted_fit}.
fd_nls <- function (globalname, data, mgen = "guess",
                    start=NULL,
                    loop=3, trace=TRUE, control=NULL, ...,
                    nlsfun=nls) {
  mdl_obj <- get_global_model(globalname)
  ctl <- control
  if (mgen == "guess") {
    mgen <- guess_mGen(data)
    if (trace)
      cat("mgen guessed: ", mgen, "\n")
    if (mgen > 16) {
      warning("guessed 'mgen' greater than 16, capped at 16")
      mgen <- 16
      loop <- 0
    }
  }
  if (is.null(start)) start <- match.fun("start")(mdl_obj)
  if (trace) {
    cat("Starting values:\n")
    print(rbind(start = start,
                lower = lower(mdl_obj),
                upper = upper(mdl_obj)))
  }
  ans <- match.fun(nlsfun)(fd_formula(globalname, mgen = mgen),
                           data,
                           na.action = na.fail,
                           start = start,
                           lower = lower(mdl_obj),
                           upper = upper(mdl_obj),
                           control = ctl,
                           trace = trace,
                           ...)
  ans$.modelname <- globalname
  ans$.model <- fd_clone(mdl_obj)
  ans$.mgen <- mgen
  ans$.flat <- any(is.na(ans$m$incr()))
  rel <- relisted_coef(ans)$pro
  gen <- reshape2::melt(
    ans$.model$pro$model(
      rel,
      times=sort(unique(data$Time)),
      mgen=mgen
    )$live_pop
  )
  if (any(gen$value[gen$Var2 == mgen + 1] > 0.1)) {
    if (loop <= 0) {
      warning("'mgen' (", mgen, ") is probably not big enough")
      ans$.mGen_hit <- max(gen$value[gen$Var2 == mgen + 1])
    } else {
      nm <- get_name(mdl_obj$pro)
      if (nm == "branching")
        mGen_new <-
          ceiling(max(rel$p / rel$f$mm * max(data$Time), mgen) * 1.2)
      else if (nm == "cyton")
        mGen_new <- ceiling(max(max(data$Time) / rel$f$mm, mgen) * 1.2)
      else
        mGen_new <- ceiling(mgen * 1.2)
      message("'mgen' (", mgen, ") is probably not big enough, ",
              "fitting again with mgen=", mGen_new)

      nm <- names(match.call()[-1])
      for (x in names(list(...))) assign(x, list(...)[[x]])
      lst <- setNames(lapply(nm, function (x) get(x)), nm)
      lst$loop <- loop - 1
      lst$mgen <- mGen_new
      do.call(Recall, lst)
    }
  } else {
    ans$.mGen_hit <- 0
  }
  hit_boundary <- coef(ans) <= lower(mdl_obj) | coef(ans) >= upper(mdl_obj)
  if (any(hit_boundary)) {
    message(paste(names(coef(ans))[hit_boundary], collapse=", "),
            " hit the boundary: you probably want to use more constraints")
  }
  class(ans) <- c("fd_nls", class(ans))
  ans
}

#' @export
#' @rdname fd_nls
model <- function (object) {
  ans <- NULL
  if (is.list(object) || is.environment(object)) ans <- object$.model
  if (is.null(ans)) ans <- attr(object, ".model")
  ans
}

#' @export
#' @rdname fd_nls
fd_residuals <- function (data, ..., .list = NULL) {
  arglst <- append(list(...), .list)
  if (is.null(names(arglst)))
    names(arglst) <- paste0("mod_", seq_along(arglst))
  if (length(arglst) == 0) {
    stop("fd_residuals: no object was provided")
  }
  df_res <- plyr::rbind.fill(lapply(names(arglst), function (nm) {
    if (inherits(data, "fd_data")) dat_cur <- data
    else {
      dat_cur <- data[[nm]]
      if (is.null(dat_cur))
        stop("'", nm, "' not found in 'data' list")
    }
    if (inherits(arglst[[nm]], "lmList")) {
      plyr::rbind.fill(lapply(names(arglst[[nm]]), function (ind) {
        transform(subset.data.frame(dat_cur, Individual == ind),
                  resid = resid(arglst[[nm]][[ind]], mode="p"),
                  fitted = fitted(arglst[[nm]][[ind]]),
                  algo = nm)
      }))
    } else {
      cat("CLAZZ :", nm, "\n")
      print(class(arglst[[nm]]))
      print(class(data))
      print(class(dat_cur))
      print(c(
        dat_cur = NROW(dat_cur),
        resid = length(resid(arglst[[nm]], mode="p")),
        fitted = length(fitted(arglst[[nm]]))
      ))
      transform(dat_cur,
                resid = resid(arglst[[nm]], mode="p"),
                fitted = fitted(arglst[[nm]]),
                algo = nm)
    }
  }))
  if (!inherits(data, "fd_data"))
    data <- do.call(rbind, data)
  data$Id <- 1:NROW(data)
  df_res <- merge(data, df_res)
  df_res[order(df_res$Id), ]
}

#' @export
#' @rdname fd_nls
predict.fd_nls <- function (object,
              newdata = stop("'newdata' is missing"),
              model = NULL, ...) {
  newdata$Id <- 1:NROW(newdata)
  if (inherits(object, "lmList")) {
    ans <- plyr::rbind.fill(lapply(
      levels(newdata$Individual),
      function (nm) {
        fit <- object[[nm]]
        categories <- names(relisted_coef(fit, drop=FALSE)[[1L]]$pro)
        data <- subset(newdata, Category %in% categories)
        if (is.null(fit)) return (NULL)
        mdl <- model
        if (is.null(mdl)) mdl <- match.fun("model")(fit)
        if (is.null(mdl)) stop("no model specified")
        dat <- subset(data, Individual == nm)
        co <- relist(coef(fit), mdl)
        mdl <- update(mdl, data=dat)
        co <- fd_freepar(co, mdl)
        fd_predict(model = mdl, param = co,
                   data = dat,
                   mgen = fit$.mgen)
      }
    ))
  } else if (inherits(object, "lme")) {
    categories <- names(relisted_coef(object, drop=FALSE)[[1L]]$pro)
    data <- subset(newdata, Category %in% categories)
    ans <- plyr::rbind.fill(lapply(
      levels(newdata$Individual),
      function (nm) {
        if (!(nm %in% rownames(coef(object)))) return (NULL)
        co <- as.matrix(coef(object))[nm, ]
        mdl <- model
        if (is.null(mdl)) mdl <- match.fun("model")(object)
        if (is.null(mdl)) stop("no model specified")
        dat <- subset(data, Individual == nm)
        co <- relist(co, mdl)
        mdl <- update(mdl, data=dat)
        co <- fd_freepar(co, mdl)
        fd_predict(model = mdl,
                   param = co,
                   data = dat,
                   mgen = object$.mgen)
      }
    ))
  } else {
    categories <- names(relisted_coef(object, drop=FALSE)[[1L]]$pro)
    data <- fd_data(subset(newdata, Category %in% categories),
            categories = categories)
    mdl <- model
    if (is.null(mdl)) mdl <- match.fun("model")(object)
    if (is.null(mdl)) stop("no model specified")
    co <- relist(coef(object), mdl)
    mdl <- update(mdl, data=data)
    co <- fd_freepar(co, mdl)
    ans <- fd_predict(model = mdl,
                      param = co,
                      data = data,
                      mgen = object$.mgen)
  }
  newdata$y <- NA
  newdata$y[ans$Id] <- ans$y
  newdata$y
}

#' @export
#' @rdname fd_nls
relisted_coef <- function (object, drop=TRUE) {
  co <- rbind(coef(object))
  ans <- setNames(
    lapply(
      1:NROW(co),
      function (i) relist(co[i, ], object$.model)),
    rownames(co)
  )
  if (drop && length(ans) == 1) ans[[1L]]
  else ans
}

#' @export
#' @rdname fd_nls
relisted_vcov <- function (object, vcov. = vcov, n=1000) {
  if (is.function (vcov.)) vcov. <- vcov.(object)
  if (!requireNamespace("mvtnorm", quietly=TRUE))
    stop("package 'mvtnorm' required for 'relisted_vcov'")

  # Actually, bounding does not make much sense...  It will distorts the
  # covariance matrix, that's all.  Bounding would be interesting, however,
  # when you want to do some Monte Carlo analysis.

  # nolint start
#   rand <- tmvtnorm::rtmvnorm(n=n, lower = lower(object),
#           upper = upper(object),
#               mean = coef(object),
#            sigma = vcov.,
#            algorithm="rejection")
  # nolint end

  # Not possible to use method "eigen" in every case, but "svd" works.
  # "chol": when rank deficient, complain.
  rand <-
    mvtnorm::rmvnorm(n=n, mean=colMeans(rbind(coef(object))),
             sigma=vcov., method="svd")
  relisted <- t(sapply(1:NROW(rand), function (i) {
    unlist(relist(rand[i, ], model(object)))
  }))
  cov(relisted)
}

#' @export
#' @rdname fd_nls
relisted_fit <- function (object, ...) {
  structure(list(
    coefficients =
      rowMeans(sapply(relisted_coef(object, drop=FALSE), unlist)),
    vcov = relisted_vcov(object, ...)
  ), class="relisted_fit")
}

#' @export
#' @rdname fd_nls
vcov.relisted_fit <- function (object, ...) {
  object$vcov
}


# mgen --------------------------------------------------------------------

#' @export
#' @rdname fd_nls
guess_mGen <- function (data) {
  rates <- fd_rates(data, cut_first = TRUE)
  rate_hists <- quantile(
    -rates$Rate[rates$Weight == "hist" & rates$Rate < 0 &
                is.finite(rates$Rate)],
    0.9)
  rate_Ns <- quantile(
    rates$Rate[rates$Weight == "N" & rates$Rate > 0 &
               is.finite(rates$Rate)],
    0.9)
  r <- max(rate_hists, rate_Ns, na.rm=TRUE)

  ans <- ceiling(max(data$Time) * r / log(2))
  if (!is.finite(ans))
    stop("internal error")
  ans
}

#' @export
#' @keywords internal
fd_rates <- function (data, cut_first = FALSE, group_hists = TRUE) {
  if (cut_first) {
    ans_hists <- plyr::ddply(
      subset.data.frame(data, Weight == "hist" & a > 0),
      setdiff(colnames(data), c("a", "b", "y")),
      function (df) {
        df2 <- df[1, ]
        cutoff <- exp(attr(data, "fmm")$m0[df$Inoculum[1]]) / 2.0
        df2$y <- with(
          df,
          exp(sum(y[y < cutoff] * (log(a) + log(b)) / 2.0) /
            sum(y[y < cutoff])))
        df2[, setdiff(colnames(df2), c("a", "b"))]
      })
  } else {
    ans_hists <- plyr::ddply(
      subset.data.frame(data, Weight == "hist"),
      setdiff(colnames(data), c("a", "b", "y")),
      plyr::summarise,
      y = exp(sum(y * (log(a) + log(b)) / 2.0) / sum(y)))
  }
  ans_hists <- plyr::ddply(
    ans_hists, setdiff(colnames(ans_hists),
               if (group_hists)
                 c("Timepoint", "Time", "y",
                 "Individual", "Inoculum")
               else c("Timepoint", "Time", "y")),
    function (df) {
      df2 <- plyr::ddply(df, "Time", plyr::summarise,
                 y = exp(mean(log(y))))
      if (NROW(df2) == 1)
        warning("'group_hists == FALSE' but some histograms ",
            "stand alone so rates cannot be determined.")
      df2 <- df2[order(df2$Time), ]
      rate <- diff(log(df2$y)) / diff(df2$Time)
      rate <- c(rate, rate[length(rate)])
      transform(df, Rate = rate[findInterval(Time, df2$Time)])
    }
  )

  ans_Ns <- subset.data.frame(data, Weight == "N")
  ans_Ns <- ans_Ns[, setdiff(colnames(ans_Ns), c("a", "b"))]
  ans_Ns <- plyr::ddply(
    ans_Ns,
    setdiff(colnames(ans_Ns), c("y", "Timepoint", "Time")),
    function (df) {
      df2 <- plyr::ddply(df, "Time", plyr::summarise, y = 10^mean(y))
      df2 <- df2[order(df2$Time), ]
      rate <- diff(log(df2$y)) / diff(df2$Time)
      rate <- c(rate, rate[length(rate)])
      transform(df, Rate = rate[findInterval(Time, df2$Time)])
    }
  )

  rbind(ans_hists, ans_Ns)
}


# Comb --------------------------------------------------------------------



#' Comb through many constraints using \code{\link{fd_nls}}.
#'
#'
#' Comb through many constraints using \code{\link{fd_nls}}.  The result can
#' then be fed to \code{\link{fd_aictab}}.
#'
#' @param cstr A list of \code{\link{constraints}} to comb through
#' @param fmm A finite mixture model, same format as with
#'   \code{\link{fd_model}}.
#' @param proliferation A proliferation model, same format as with
#'   \code{\link{fd_model}}.
#' @param data The \code{\link{fd_data}} dataset on which to perform the search.
#' @param continue If \code{TRUE} (default), an error in an optimization does
#'   not stop the whole search but simply results in the corresponding entry of
#'   the return value to be \code{NULL}.  If \code{FALSE}, an error is raised
#'   instead.
#' @param mc.cores Number of cores to use for parallel search using
#'   \pkg{parallel}. \code{mc.cores > 1} fails on Windows machines.
#' @param file Output \code{file} (or file name) for parallel search (during
#'   which the console output is silenced to avoid conflicts).  The output file
#'   can be read in real time on POSIX machines using \code{tail -f <file>}.
#' @param globalname Global name to use for storing the \code{\link{fd_model}}
#'   object. See \code{\link{fd_model-functions}} for more details.
#' @param ... Additional parameters to pass along to \code{\link{fd_nls}}.
#'
#' @return A named list, each entry being either the result of the optimization
#'   (an \code{fd_nls} object) or \code{NULL} if the optimization failed, using
#'   the constraints \code{cstr} and in the same order as the constraints
#'   \code{cstr}.  The names are the names of \code{cstr}, if provided, or
#'   \code{constr_<i>}, with \code{<i>} the index of the constraint if no name
#'   was made available.
#' @export
#' @family optimization-related entries
#' @seealso \code{\link{fd_nls}}, \code{\link{fd_aictab}}
#'
#' @examples
#' # Shared constraints
#' CC <<- FdCommonConstraints
#' `#macr` <- ~CC$`#noss` + pro:all:{res <- 0} + fmm:{c0 <- 0} + CC$`#mm_xyxz`
#'
#' control <- list(maxit = 1L)
#' \dontrun{
#' control <- list(maxit = 100L, simple.function = TRUE)
#' }
#'
#' source(system.file("contrib", "nlsSA.R", package="fluodilution"))
#' data(FdSTyphimuriumWTC57)
#' ans <- fd_comb(lapply(list("1100" = CC$`#delta_1100`,
#'                            "1111" = CC$`#delta_1111`),
#'                       catcstr, `#macr`),
#'                data = cutoff(FdSTyphimuriumWTC57),
#'                control = control,
#'                nlsfun = nlsSA,
#'                trace = TRUE)
#'
#' summary(ans[[1L]])
#' summary(ans[[2L]])
#'
#' @section Required packages: Please install \code{parallel}.
fd_comb <- function (cstr, fmm="gaussian", proliferation="branching",
                     data, continue = TRUE, mc.cores=1, file=NULL,
                     globalname = "mdl", ...) {
  if (!is.null(file))
    file.create(file)
  if (length(cstr) == 0) return (NULL)
  if (is.null(names(cstr))) {
    names(cstr) <- paste0("constr_", seq_along(cstr))
  }
  if (mc.cores > 1) {
    if (!requireNamespace("parallel", quietly=TRUE))
      stop("'parallel' must be installed for 'mc.cores > 0'")
    message("Using parallel computing with ", mc.cores, " cores")
    fn <- function (...) parallel::mclapply(mc.cores = mc.cores, ...)
  } else {
    fn <- lapply
  }
  ans <- setNames(fn(seq_along(cstr),
          function (cur) {
    if (!is.null(file)) {
      if (is.character(file))
        fileConn <- file(file)
      else
        fileCon <- file
      sink(fileConn, append=TRUE)
      sink(fileConn, append=TRUE, type="message")
    }
    tryCatch({
      cat(sep="", get_name(proliferation), ": ", names(cstr)[cur], "\n")
      dat <- data
      if (continue) funmdl <- dplyr::failwith(default=NULL, fd_model)
      else funmdl <- fd_model
      assign(globalname, funmdl(fmm=fmm,
               proliferation=proliferation,
          data=dat,
          constraints = cstr[[cur]],
          boxed=TRUE), envir = .GlobalEnv)
      if (is.null(get(globalname, .GlobalEnv))) return (NULL)
      if (continue) fun <- dplyr::failwith(default=NULL, fd_nls)
      else fun <- fd_nls
      fit <- fun(globalname, dat, ...)
      return (fit)
    },
    finally = {
      if (!is.null(file)) {
        sink(NULL)
        sink(NULL, type="message")
      }
    })
  }),
  names(cstr))
  return (structure(ans, class="fd_comb"))
}