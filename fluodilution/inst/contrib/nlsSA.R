# Copyright (c) 2015-2018 Hadrien Chauvin
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#####################

library(stats)

# Borrowed from R Core ----------------------------------------------------

# Exactly the same as the non-exported function in stats (stats:::nlsModel).
# However, we reproduce it here for backward compatibility.
# The same way, nlsGeneral is a blatant rip-off of stats::nls.

# The license for nls in R core is the GNU General Public License version 2
# or later, the same that applies to this file (by virtue of the GPL license).

# The 'stats' package has been authored by the R Core Team and
# contributors worldwide.

nlsModelGeneral <- function (form, data, start, wts, upper = NULL)
{
  thisEnv <- environment()
  env <- new.env(hash = TRUE, parent = environment(form))
  for (i in names(data)) assign(i, data[[i]], envir = env)
  ind <- as.list(start)
  parLength <- 0
  for (i in names(ind)) {
    temp <- start[[i]]
    storage.mode(temp) <- "double"
    assign(i, temp, envir = env)
    ind[[i]] <- parLength + seq_along(start[[i]])
    parLength <- parLength + length(start[[i]])
  }
  getPars.noVarying <- function() unlist(setNames(lapply(names(ind),
    get, envir = env), names(ind)))
  getPars <- getPars.noVarying
  internalPars <- getPars()
  if (!is.null(upper)) upper <- rep_len(upper, parLength)
  useParams <- rep(TRUE, parLength)
  lhs <- eval(form[[2L]], envir = env)
  rhs <- eval(form[[3L]], envir = env)
  .swts <- if (!missing(wts) && length(wts)) sqrt(wts) else rep_len(1, length(rhs))
  assign(".swts", .swts, envir = env)
  resid <- .swts * (lhs - rhs)
  dev <- sum(resid^2)
  if (is.null(attr(rhs, "gradient"))) {
    getRHS.noVarying <- function() {
      if (is.null(upper))
        numericDeriv(form[[3L]], names(ind), env)
      else numericDeriv(form[[3L]], names(ind), env, ifelse(internalPars <
        upper, 1, -1))
    }
    getRHS <- getRHS.noVarying
    rhs <- getRHS()
  }
  else {
    getRHS.noVarying <- function() eval(form[[3L]], envir = env)
    getRHS <- getRHS.noVarying
  }
  dimGrad <- dim(attr(rhs, "gradient"))
  marg <- length(dimGrad)
  if (marg > 0L) {
    gradSetArgs <- vector("list", marg + 1L)
    for (i in 2L:marg) gradSetArgs[[i]] <- rep(TRUE, dimGrad[i -
      1])
    useParams <- rep(TRUE, dimGrad[marg])
  }
  else {
    gradSetArgs <- vector("list", 2L)
    useParams <- rep(TRUE, length(attr(rhs, "gradient")))
  }
  npar <- length(useParams)
  gradSetArgs[[1L]] <- (~attr(ans, "gradient"))[[2L]]
  gradCall <- switch(length(gradSetArgs) - 1L, call("[", gradSetArgs[[1L]],
    gradSetArgs[[2L]], drop = FALSE), call("[", gradSetArgs[[1L]],
    gradSetArgs[[2L]], gradSetArgs[[2L]], drop = FALSE),
    call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]],
      gradSetArgs[[3L]], drop = FALSE), call("[", gradSetArgs[[1L]],
      gradSetArgs[[2L]], gradSetArgs[[2L]], gradSetArgs[[3L]],
      gradSetArgs[[4L]], drop = FALSE))
  getRHS.varying <- function() {
    ans <- getRHS.noVarying()
    attr(ans, "gradient") <- eval(gradCall)
    ans
  }
  QR <- qr(.swts * attr(rhs, "gradient"))
  qrDim <- min(dim(QR$qr))
  if (QR$rank < qrDim) {
    message("singular gradient matrix at initial parameter estimates")
  }

  getPars.varying <- function() unlist(setNames(lapply(names(ind),
    get, envir = env), names(ind)))[useParams]
  setPars.noVarying <- function(newPars) {
    assign("internalPars", newPars, envir = thisEnv)
    for (i in names(ind)) assign(i, unname(newPars[ind[[i]]]),
      envir = env)
  }
  setPars.varying <- function(newPars) {
    internalPars[useParams] <- newPars
    for (i in names(ind)) assign(i, unname(internalPars[ind[[i]]]),
      envir = env)
  }
  setPars <- setPars.noVarying
  on.exit(remove(i, data, parLength, start, temp, m))
  m <- list(resid = function() resid, fitted = function() rhs,
    formula = function() form, deviance = function() dev,
    lhs = function() lhs, gradient = function() .swts * attr(rhs,
      "gradient"), conv = function() {
      if (npar == 0) return(0)
      rr <- qr.qty(QR, resid)
      sqrt(sum(rr[1:npar]^2) /sum(rr[-(1:npar)]^2))
    }, incr = function() qr.coef(QR, resid), setVarying = function(vary = rep(TRUE,
      length(useParams))) {
          assign("useParams", if (is.character(vary)) {
                temp <- logical(length(useParams))
                temp[unlist(ind[vary])] <- TRUE
                temp
            } else if (is.logical(vary) && length(vary) != length(useParams)) stop("setVarying : 'vary' length must match length of parameters") else {
                vary
            }, envir = thisEnv)
            gradCall[[length(gradCall) - 1L]] <<- useParams
            if (all(useParams)) {
                assign("setPars", setPars.noVarying, envir = thisEnv)
                assign("getPars", getPars.noVarying, envir = thisEnv)
                assign("getRHS", getRHS.noVarying, envir = thisEnv)
                assign("npar", length(useParams), envir = thisEnv)
            } else {
                assign("setPars", setPars.varying, envir = thisEnv)
                assign("getPars", getPars.varying, envir = thisEnv)
                assign("getRHS", getRHS.varying, envir = thisEnv)
                assign("npar", length(seq_along(useParams)[useParams]),
                  envir = thisEnv)
            }
        }, setPars = function(newPars) {
            setPars(newPars)
            assign("resid", .swts * (lhs - assign("rhs", getRHS(),
                envir = thisEnv)), envir = thisEnv)
            assign("dev", sum(resid^2), envir = thisEnv)
            assign("QR", qr(.swts * attr(rhs, "gradient")), envir = thisEnv)
            return(QR$rank < min(dim(QR$qr)))
        }, getPars = function() getPars(), getAllPars = function() getPars(),
        getEnv = function() env, trace = function() {
            cat(format(dev), ": ", format(getPars()))
            cat("\n")
        }, Rmat = function() qr.R(QR), predict = function(newdata = list(),
            qr = FALSE) eval(form[[3L]], as.list(newdata), env))
    class(m) <- "nlsModel"
    m
}

SA_nls <- function () {
    if (!requireNamespace("GenSA", quietly=TRUE))
        stop("package 'GenSA' has to be installed for nlsSA to work")
    structure(function (...) {
            ans <- GenSA::GenSA(...)
            ans$convInfo <- list(isConv = TRUE,
                                 finIter = NA,
                                 finTol = NA,
                                 stopCode = NA,
                                 stopMessage = NA)
            return (ans)
        }, general = TRUE, supportBoxing = TRUE, name = "SA")
}

LM_nls <- function () {
    if (!requireNamespace("minpack.lm", quietly=TRUE))
        stop("'minpack.lm' required for 'nlsLM'")
    structure(minpack.lm::nls.lm, general = FALSE, supportBoxing = TRUE, name = "LM")
}

optimx_nls <- function (method) {
    structure(function (...) {
        ans <- optimx(method = method, ...)
        list(ans = ans,
             par = coef(ans)[1, ],
             convInfo = list(isConv = TRUE,
                             finIter = NA,
                             finTol = NA,
                             stopCode = NA,
                             stopMessage = NA))
    }, general = TRUE, supportBoxing = TRUE, name = "optimx")
}

#' @inheritParams stats::nls
#' @export
nlsSA <- function (...) {
    nlsGeneral(algorithm="SA", ...)
}

nlsLM <- function (...) {
    nlsGeneral(algorithm="LM", ...)
}

nlsOptimx <- function (method, ...) {
    nlsGeneral(algorithm=optimx_nls(method), ...)
}

nlsCombine <- function (algorithms, controls, save = TRUE, ...) {
  if (!is.list(algorithms) || !is.list(controls) ||
    length(algorithms) != length(controls)) {
    stop("'algorithms' and 'controls' must both be lists, of the same length")
    }
  fits <- lapply(seq_along(algorithms), function (i) {
    nlsGeneral(algorithm = algorithms[[i]], control = controls[[i]],
               ...)
  })

  best <- which.max(sapply(fits, logLik))
  ans <- fits[[best]]
  if (save) {
    attr(ans, "fits") <- fits
  }
  attr(ans, "best") <- best
  attr(ans, "algorithms") <- sapply(fits, function (cur) cur$call$algorithm)
  ans
}

nlsGeneral <- function (
  formula, data = parent.frame(), start, jac = NULL,
  algorithm = "SA", control = list(),
  trace = FALSE, subset, weights, na.action,
  model = FALSE, lower = -Inf, upper = Inf, ...)
{
  formula <- as.formula(formula)
  if (!is.list(data) && !is.environment(data))
      stop("'data' must be a list or an environment")
  mf <- match.call()
  varNames <- all.vars(formula)

  if (is.null(start)) {
      if (is.null(lower) || is.null(upper) ||
            any(!is.finite(lower) || any(!is.finite(upper))))
        stop("start is not specified, but no finite boundary is given")
     start <- (lower + upper) / 2.0
  }

  # ALGORITHM
  if (is.function(algorithm)) {
      algorithmFunc <- algorithm
      algorithm <- attr(algorithmFunc, "name")
  }
  else if (is.character(algorithm)) {
      if (algorithm %in% c("default", "plinear", "port")) {
          # Go pack to the standard nls
          nm <- names(mf[-1])
          lst <- setNames(lapply(nm, function (x) get(x)), nm)
          ans <- do.call(nls, lst)
          ans$call <- mf
          ans$data <- substitute(data)
          return (ans)
      } else if (algorithm %in% c("brute-force", "grid-search", "random-search",
                                  "plinear-brute", "plinear-random")) {
          if (!requireNamespace("nls2", quietly=TRUE))
              stop("'nls2' required for NLS algorithm ", algorithm)
          nm <- names(mf[-1])
          lst <- setNames(lapply(nm, function (x) get(x)), nm)
          if (is.null(lst$lower) || is.null(lst$upper) ||
                any(!is.finite(lst$lower)) || any(!is.finite(lst$upper))) {
              stop("for an 'nls2' algorithm, both 'lower' and 'upper' must be finite")
          }
          if (length(lst$lower) != length(lst$upper))
              stop("'lower' and 'upper' must have the same length")
          lst$start <- as.data.frame(rbind(lst$lower, lst$upper))
          lst$lower <- lst$upper <- NULL
          ans <- do.call(nls2::nls2, lst)
          ans$call <- mf
          ans$data <- substitute(data)
          return (ans)
      } else {
          algorithmFunc <- match.fun(paste0(algorithm, "_nls"))()
      }
  } else
      stop("'algorithm' is neither a string or a function")

  if (length(formula) == 2L) {
    formula[[3L]] <- formula[[2L]]
    formula[[2L]] <- 0
  }

  form2 <- formula
  form2[[2L]] <- 0
  varNamesRHS <- all.vars(form2)
  mWeights <- missing(weights)

  ## if trace = TRUE, set nls.lm.control$nprint = 1
  #if (trace) control$nprint <- 1

  pnames <- if (missing(start)) {
    if (!is.null(attr(data, "parameters"))) {
      names(attr(data, "parameters"))
    }
    else {
      cll <- formula[[length(formula)]]
      func <- get(as.character(cll[[1L]]))
      if (!is.null(pn <- attr(func, "pnames")))
        as.character(as.list(match.call(func, call = cll))[-1L][pn])
    }
  } else names(start)

  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()
  if (length(pnames)) varNames <- varNames[is.na(match(varNames, pnames))]
  lenVar <- function(var) tryCatch(length(eval(as.name(var), data, env)), error = function(e) -1)

  if (length(varNames)) {
    n <- sapply(varNames, lenVar)
    if (any(not.there <- n == -1)) {
      nnn <- names(n[not.there])
      if (missing(start)) {
        warning("No starting values specified for some parameters.\n",
                "Initializing ", paste(sQuote(nnn), collapse = ", "),
                " to '1.'.\n", "Consider specifying 'start' or using a selfStart model")
        start <- as.list(rep(1, length(nnn)))
        names(start) <- nnn
        varNames <- varNames[i <- is.na(match(varNames, nnn))]
        n <- n[i]
      }
      else stop("parameters without starting value in 'data': ", paste(nnn, collapse = ", "))
    }
  } else {
    if (length(pnames) && any((np <- sapply(pnames, lenVar)) == -1)) {
      message("fitting parameters ", paste(sQuote(pnames[np == -1]), collapse = ", "), " without any variables")
      n <- integer()
    } else stop("no parameters to fit")
  }

  respLength <- length(eval(formula[[2L]], data, env))

  if (length(n) > 0L) {
    varIndex <- n%%respLength == 0
    if (is.list(data) && diff(range(n[names(n) %in% names(data)])) > 0) {
      mf <- data
      if (!missing(subset)) warning("argument 'subset' will be ignored")
      if (!missing(na.action)) warning("argument 'na.action' will be ignored")
      if (missing(start)) start <- getInitial(formula, mf)
      startEnv <- new.env(hash = FALSE, parent = environment(formula))
      for (i in names(start)) assign(i, start[[i]], envir = startEnv)
      rhs <- eval(formula[[3L]], data, startEnv)
      n <- NROW(rhs)
      wts <- if (mWeights) rep(1, n) else eval(substitute(weights), data, environment(formula))
    } else {
      mf$formula <- as.formula(paste("~", paste(varNames[varIndex], collapse = "+")), env = environment(formula))
      mf$start <- mf$control <- mf$algorithm <- mf$trace <- mf$model <- NULL
      mf$lower <- mf$upper <- NULL
      mf[[1L]] <- as.name("model.frame")
      mf <- eval.parent(mf)
      n <- nrow(mf)
      mf <- as.list(mf)
      wts <- if (!mWeights) model.weights(mf) else rep(1, n)
    }
    if (any(wts < 0 | is.na(wts))) stop("missing or negative weights not allowed")
  } else {
    varIndex <- logical()
    mf <- list(0)
    wts <- numeric()
  }

  if (missing(start)) start <- getInitial(formula, mf)
  for (var in varNames[!varIndex]) mf[[var]] <- eval(as.name(var), data, env)
  varNamesRHS <- varNamesRHS[varNamesRHS %in% varNames[varIndex]]

  mf <- c(mf, start)
  lhs <- eval(formula[[2L]], envir = mf)
  m <- match(names(start), names(mf))
  .swts <- if (!missing(wts) && length(wts)) sqrt(wts)

  ###
  nlsctl <- list(report.frq = 10)
  commonctlnames <- intersect(names(control), names(nlsctl))
  nlsctl[commonctlnames] <- control[commonctlnames]
  optimctl <- control[setdiff(names(control), names(nlsctl))]

  FCT_RES <- local({
      lastReportTime <- proc.time()[3]
      startTime <- proc.time()[3]
      nbCalls <- 0
      bestSumRes <- NA
      function (par) {
          nbCalls <<- nbCalls + 1
        mf[m] <- par
        rhs <- eval(formula[[3L]], envir = mf, environment(formula))
        res <- lhs - rhs
        res <- .swts * res
        if (trace && proc.time()[3] > lastReportTime + nlsctl$report.frq) {
              lastReportTime <<- proc.time()[3]
              curRes <- sum(res^2)
              if (!is.finite(bestSumRes) || (is.finite(curRes) && curRes < bestSumRes))
                  bestSumRes <<- curRes
              cat(sep="", "Best sum of square: ", bestSumRes,
                  " (current: ", curRes, "; avg exec time: ",
                  (proc.time()[3] - startTime) / nbCalls,
                  "s)\n")
              if (nbCalls > 30) {
                  nbCalls <- 0
                  startTime <- proc.time()[3]
              }
        }
        res
      }
  })
  if (attr(algorithmFunc, "general")) {
      FCT <- function(par) {
        sum(FCT_RES(par)^2)
      }
  } else {
      FCT <- FCT_RES
  }

  if (attr(algorithmFunc, "supportBoxing")) {
      if (length(lower) == 1 && lower == -Inf) lower <- rep(-Inf, length(start))
      if (length(upper) == 1 && upper == Inf) upper <- rep(Inf, length(start))
      NLS <- algorithmFunc(par = start, fn = FCT, control = optimctl,
                lower = lower, upper = upper, ...)
  } else {
      NLS <- algorithmFunc(par = start, fn = FCT, control = optimctl,
                       ...)
  }

  ## previous versions before bounds were included
  # NLS <- nls.lm(par = start, fn = FCT, jac = jac, control = control, ...)
  # That's because GenSA does not set names
  start <- setNames(NLS$par, names(start))

  ##pass optimized parameters to 'nlsModel'
  m <- nlsModelGeneral(formula, mf, start, wts)

  convInfo <- NLS$convInfo
  nls.out <- list(m = m,
                  convInfo = convInfo,
                  data = substitute(data),
                  call = match.call())

  nls.out$call$algorithm <- algorithm
  nls.out$opt <- NLS
  ## need to use '$tol' parameter from nls.control to make 'predict.nls' work
  #nls.out$call$control <- nls.control()
  nls.out$call$trace <- FALSE
  nls.out$na.action <- attr(mf, "na.action")
  nls.out$dataClasses <- attr(attr(mf, "terms"), "dataClasses")[varNamesRHS]
  if (model)
    nls.out$model <- mf
  if (!mWeights)
    nls.out$weights <- wts
  nls.out$control <- control
  class(nls.out) <- c("nls", "nlsSA")
  nls.out
}


# Borrowed from nlme ------------------------------------------------------

nlsListGeneral <-
  ## A list of nls objects
  function(model, data, start, control, level, subset, na.action = na.fail,
           pool = TRUE, ...) {
      UseMethod("nlsListGeneral")
  }

nlsListGeneral.selfStart <-
  function (model, data, start, control, level, subset, na.action = na.fail,
            pool = TRUE, ...)
{
  mCall <- as.list(match.call())[-1]
  if (!inherits(data, "groupedData")) {
    stop("second argument must be a groupedData object")
  }
  marg <- substitute(model)
  if (mode(marg) != "name") {
    stop("cannot use an anonymous function for the model")
  }
  # Build up a call to the model function
  m <- call(as.character(marg))
  args <- lapply(names(formals(eval(marg))), as.name)
  args[[1L]] <- getCovariateFormula(data)[[2L]]
  m[1 + seq_along(args)] <- args
  form <- formula(data)
  form[[3L]][[2L]] <- m
  mCall$model <- form
  do.call("nlsListGeneral.formula", mCall)
  }

#' @inheritParams stats:::nls
#' @export
nlsListSA <- function (model, ...) {
    nlsListGeneral(model, algorithm="SA", ...)
}

nlsListLM <- function (model, ...) {
    nlsListGeneral(model, algorithm="LM", ...)
}

nlsListOptimix <- function (model, method, ...) {
    nlsListGeneral(model, algorithm=optimx_nls(method), ...)
}

nlsListGeneral.formula <-
  function(model, data, start = NULL, control, level, subset,
           na.action = na.fail, pool = TRUE, algorithm = "SA", ...) {

  Call <- match.call()
  if (!missing(subset)) {
    data <-
      data[eval(asOneSidedFormula(Call[["subset"]])[[2L]], data),, drop = FALSE]
  }
  if (!inherits(data, "data.frame")) data <- as.data.frame(data)
  data <- na.action(data)
  if (is.null(grpForm <- getGroupsFormula(model))) {
    if (inherits(data, "groupedData")) {
      if (missing(level)) level <- length(getGroupsFormula(data, asList = TRUE))
      else if (length(level) > 1) {
        stop("multiple levels not allowed")
      }
      groups <- getGroups(data, level = level)[drop = TRUE]
      grpForm <- getGroupsFormula(data)
    } else {
      stop("'data' must be a \"groupedData\" object if 'formula' does not include groups")
    }
  } else {
    if (missing(level)) {
      level <- length(getGroupsFormula(model, asList = TRUE))
    } else if (length(level) > 1) {
      stop("multiple levels not allowed")
    }
    model <- eval(parse(text = paste(paste(deparse(model[[2L]]), collapse=" "),
                        paste(deparse(getCovariateFormula(model)[[2L]]), collapse=" "),
            sep = "~")))
    groups <- getGroups(data, form = grpForm, level = level)[drop = TRUE]
  }
  if (is.null(start) && is.null(attr(data, "parameters"))) {
    ## no starting values
    ## checking for old-style selfStart functions
    FUN <- eval(model[[3L]][[1L]])
    if (is.function(FUN) && class(FUN) != "selfStart" &&
        !is.null(attr(FUN, "initial"))) {
      stop("old-style self-starting model functions\nare no longer supported.\nNew selfStart functions are available.\nUse\n  SSfpl instead of fpl,\n  SSfol instead of first.order.log,\n  SSbiexp instead of biexp,\n  SSlogis instead of logistic.\nIf writing your own selfStart model, see\n  \"help(selfStart)\"\nfor the new form of the \"initial\" attribute.")
    }
  }

  controlvals <- control
  val <- lapply(split(data, groups),
        function(dat, formula, start, control, first = TRUE) {
                  ans <- try({
                    data <- as.data.frame(dat)
                      nlsGeneral(formula = formula, algorithm = algorithm,
                                 start = start,
                            data = data, control = control,
                            ...)
                  })
                  if (inherits(ans, "try-error"))
                    NULL
                  else ans
        }, formula = model, start = start, control = controlvals)
  if (inherits(data, "groupedData")) {
    ## saving labels and units for plots
    attr(val, "units") <- attr(data, "units")
    attr(val, "labels") <- attr(data, "labels")
    attr(val, "outer") <- attr(data, "outer")
  }
  attr(val, "dims") <- list(N = nrow(data), M = length(val))
  attr(val, "call") <- Call
  attr(val,"groups") <- ordered(groups, levels = names(val))
  attr(val, "origOrder") <- match(unique(as.character(groups)), names(val))
  attr(val, "pool") <- pool
  attr(val, "groupsForm") <- grpForm
  class(val) <- c("nlsListSA", "nlsList", "lmList")
  val
}







