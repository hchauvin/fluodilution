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

# fd_model and related ----------------------------------------------------

#' Create and maintain a fluorescence dilution (FD) model.
#'
#' These functions create and maintain a fluorescence dilution (FD) model,
#' linking together a Finite Mixture Model (FMM), a proliferation model and a
#' set of constraints.
#'
#' @param fmm A Finite Mixture Model (FMM), see \code{\link{finite-mixture}}.
#' @param proliferation A proliferation model, see \code{\link{proliferation}}.
#' @param data An optional \code{\link{fd_data}} used to construct the FMM and
#'   the proliferation model.  If not specified, default values are used.
#' @param constraints A \code{cstr} object, as produced by
#'   \code{\link{cstrlist}} or an \code{expression} if only the
#'   \code{constraint} component is set (see \code{\link{constraints}}).
#' @param partial A numeric vector giving the index (not the name) of the
#'   categories in the dataset to consider for fitting.  By default, all the
#'   categories are considered.  Currently, the categories can only be
#'   contiguous and include \code{1} (e.g., \code{1:3} but not \code{2L:3L}).
#'   \code{partial} helps to shave off some computation time.
#' @param boxed Whether the parameters of the model should be boxed by upper and
#'   lower bounds (\code{TRUE}) or "unboxed" by an \code{atanh} transformation
#'   so as to be used with optimization methods that do not allow the
#'   specification of a boundary (such as the default algorithm of
#'   \code{\link[stats]{nls}}).  However, if the likely optimal parameters are
#'   far from the boundary, \code{boxed=TRUE} can in general be used (for
#'   example with \code{\link[nlme]{nlme}}).
#' @param process A function that is used internally by
#'   \code{\link{fd_predict}}.  For advanced use, the default
#'   \code{fd_process_default} can be overridden to accommodate for, e.g.,
#'   "weights" other than fluorescence histograms and cell counts (see source
#'   code of \verb{
#'   }\code{fd_process_default}, exported but not documented, for
#'   format).
#' @param memoise Whether the Finite Mixture Model stage should be memoised with
#'   package \pkg{memoise}.  Memoisation drastically improves performance but
#'   should probably be turned off, out of memory concerns, when optimizing
#'   autofluorescence level and specific number of molecules with FMMs
#'   \code{\link{fd_fmm_af}} and \code{\link{fd_fmm_af_bp}}.
#' @param object An \code{fd_model} object.
#' @param clean Whether \code{fd_clean} should be called as well.
#' @param addcstr Additional constraints to be added (using
#'   \code{\link{catcstr}}) to the current set of constraints.  Either a
#'   \code{\link{cstrlist}} or an \code{expression} if only the
#'   \code{constraint} component is to be updated.
#' @param start Numeric vector.  Overrides the starting parameters available
#'   through\verb{
#' }\code{\link{start}(object)}.
#' @param ... Not used.
#'
#' @return \code{fd_model} returns an \code{fd_model} object (an
#'   \code{\link[base]{environment}}). \code{fd_clean} modifies the model passed
#'   as a parameter and does not have a meaningful return value.
#'   \code{fd_clone}, on the other hand, returns a model that is not identical
#'   to the model passed along as argument.  \code{update} also returns a new
#'   environment.
#' @seealso \code{\link{fd_model-functions}} for functions to be used with an
#'   \code{fd_model} object, \code{\link{fd_simulate}} for simulating a
#'   population, \code{\link{fd_nls}} for fitting an FD model.
#' @export
#'
#' @details \code{fd_model} creates a new fluorescence dilution, two-step
#' hierarchical model. This model can then be optimized with
#' \code{\link{fd_nls}}.  Predicting experimental values necessitate the
#' calculation of internal quantities that can be bulky in memory:
#' \code{fd_clean} provides for their removal.  Because models are environments,
#' copying them can only be done through \code{fd_clone}.  \code{update} updates
#' an existing model by changing some of its arguments.
#'
#' @section From full parameter list to vector of free parameters: The
#'   parametrization of an FD model can be conveniently made with nested lists
#'   (see \code{\link{constraints}}).  More specifically, the top-level list has
#'   two elements named \code{fmm} and \code{pro} (for the format of those two
#'   elements, see respectively \code{\link{finite-mixture}} and
#'   \code{\link{proliferation}}).  However, multidimensional optimization
#'   functions only take vectors as arguments (the "free" parameters).
#'   Therefore, for simple problems, the starting parameters and lower/upper
#'   bounds are \code{unlist}ed and passed as arguments and the model is
#'   evaluated using a \code{relist}ing.  Here, we needed to introduce further
#'   steps along the way.
#'
#'   To anchor the discussion, we call the parametrization with which the model
#'   is eventually evaluated the "natural parametrization".  A first step away
#'   from this parametrization and towards the vector of free parameters
#'   involves reparametrizing some parameters to allow the specification of
#'   bounds independent from each other and rescaling to speed up convergence.
#'   In the case of the branching processes, because the proportions must sum to
#'   \code{1}, i.e. \code{p + res + d = 1}, with \code{p} the proportion of
#'   cells that divide, \code{res} the proportion of non-growers and \code{d}
#'   the proportion of cells that die or leave the compartment, we reparametrize
#'   \code{p} and \code{res} to \code{`p'` = p + res} and \code{`res'` = res /
#'   (p + res)}, allowing us to give the lower bounds \code{c(`p'` = 0, `res'` =
#'   0)} and upper bounds \code{c(`p'` = 1, `res'` = 1)}. Alternatively,
#'   Lagrange multipliers could have been used, but this leads to instability in
#'   such large nonlinear problems.
#'
#'   From this reparametrized version, we go on to constraining with
#'   \code{\link{constraints}}.  Constraining allows to reduce the dimension of
#'   a problem that would be too complicated to solve, even using global search
#'   (the "curse" of dimensionality), and in any case likely to be
#'   overparametrized in practice (henceforth not identifiable in practice).
#'
#'   The free parameters that are not constrained are then \code{unlist}ed, and
#'   it is this named vector that is given to optimizers and returned by
#'   \code{\link{start}}, \code{\link{lower}} and \code{\link{upper}}.
#'
#' In short, the steps from the natural parametrization to free parameters are:
#' \preformatted{
#' Action                Object
#' --------------        --------------
#'                       (1) natural
#' transformation   ->   (2) reparametrized
#' constraining     ->   (3) constrained
#' unlisting        ->   (4) unlisted
#' }
#' and naturally when evaluating the model the reverse procedure is applied:
#' \preformatted{
#' Action                        Object
#' ----------------------        --------------
#'                               (4) unlisted
#' relisting                ->   (3) constrained
#' expansion                ->   (2) reparametrized
#' inverse transformation   ->   (1) natural
#' }
#'
#' @section Warning:
#' The user should check with \code{\link{lower}} and \code{\link{upper}}
#' that the boundary
#' suits their needs.  If they do not, the \code{constraints} argument
#' should be a
#' \code{\link{cstrlist}} object.
#'
#' @examples
#' # Create a new model with default values
#' mdl <- fd_model()
#' print(mdl)
#' print(summary(mdl))
#'
#' # Original size
#' sum(sapply(ls(mdl), function (nm) object.size(mdl[[nm]])))
#' sim <- fd_simulate(fd_draw_unif(mdl), c(2, 6, 12))
#'
#' # Do calculations
#' fd_minuslogl(mdl, sim, verbose=FALSE)()
#' sum(sapply(ls(mdl), function (nm) object.size(mdl[[nm]])))
#'
#' # Clean up a bit
#' fd_clean(mdl)
#' # Reduced size
#' sum(sapply(ls(mdl), function (nm) object.size(mdl[[nm]])))
#'
#' # Clone
#' mdl2 <- fd_clone(mdl)
#' identical(mdl, mdl2)  # FALSE
#'
#' # Update model (updating clones as well)
#' mdl3 <- update(mdl2, proliferation="cyton")
#' print(mdl2$pro)
#' print(mdl3$pro)
#' identical(mdl2, mdl3)  # FALSE
fd_model <- function (fmm = "gaussian",
                      proliferation = "branching", data = NULL,
                      constraints = NULL,
                      partial = NULL, boxed = TRUE,
                      process = "fd_process_default",
                      memoise = TRUE) {
  if (!is.null(data) && !inherits(data, "fd_data"))
    stop("'data' his not an 'fd_data'")
  if (!is.logical(boxed))
    stop("'boxed' should be either TRUE or FALSE")
  if (!is.cstr(constraints))
    stop("'constraints' should be a valid constraint")
  if (!is.null(partial) && !is.numeric(partial))
    stop("'partial' should be numeric")

  # Process fmm and proliferation to a common format
  ret <- new.env()
  ret$fmm <- get_fmm(fmm, data)
  ret$pro <- get_proliferation(proliferation, data)

  if (is.null(process)) {
    # Backward compatibility
    process <- "fd_process_default"
  }
  if (is.character(process)) ret$processname <- process
  else ret$processname <- "(function)"
  ret$process <- process

  # Manage transformations
  ret$boxed <- boxed
  ret$trans <- function (object) {
    list(fmm = ret$fmm$trans(object$fmm),
       pro = ret$pro$trans(object$pro))
  }
  ret$trans_inverse <- function (object) {
    list(fmm = ret$fmm$trans_inverse(object$fmm),
       pro = ret$pro$trans_inverse(object$pro))
  }

  # Manage constraints
  if (is.null(constraints)) constraints <- cstrlist()
  else if (!is.list(constraints))
    constraints <- cstrlist(constraints)
  ret$form <- constraints
  ret$partial <- partial
  for (form in c("constraints", "start", "lower", "upper")) {
    ret$form[[form]] <- ret$fmm$constrain(ret$form[[form]], form)
    ret$form[[form]] <- ret$pro$constrain(ret$form[[form]], form)
    if (form != "constraints") {
      ret[[form]] <- applycstr(applycstr(ret$start, ret$form[[form]]),
                   ret$form$constraints)
      if (length(ret[[form]]) != 2 ||
          !all(names(ret[[form]]) %in% c("fmm", "pro")))
        stop("the constraints must contain 'fmm' and 'pro' ",
           "and nothing else")
      if (!all(names(ret[[form]]) == c("fmm", "pro")))
        ret[[form]] <- ret[[form]][2L:1L]
      if (any(all.equal(names(ret[[form]]$pro)[names(ret[[form]]$pro) %in%
                     ret$pro$categories],
          ret$pro$categories[ret$pro$categories %in%
                      names(ret[[form]]$pro)]) != TRUE))
        stop("'constraints': categories in wrong order")
      if (!is.null(ret$partial)) {
        if (length(ret$partial) == 0)
          stop("'partial' is of length 0")
        if (!all(ret$partial %in% seq_along(ret[[form]]$pro)))
          stop("'partial': some elements are not present in ",
             "constraints")
        if (any(duplicated(ret$partial)))
          stop("'partial': some elements are duplicated")
        ret[[form]]$pro <- ret[[form]]$pro[ret$partial]
      }
    }
  }

  # Get free parameters
  ret$parsedcstr <- parsecstr(ret$form$constraints)
  ret$freepar <- setdiff(names(unlist(ret$start)),
               make.names(ret$parsedcstr$prefix))
  if (length(ret$freepar) == 0) stop("no free parameter")

  # Start, lower and upper with only free parameters
  ret$boxed_lower <- getfreepar(ret, ret$lower, boxed=TRUE)
  ret$boxed_upper <- getfreepar(ret, ret$upper, boxed=TRUE)
  ret$free_lower <- getfreepar(ret, ret$lower, boxed=ret$boxed)
  ret$free_upper <- getfreepar(ret, ret$upper, boxed=ret$boxed)
  ret$free_start <- getfreepar(ret, ret$start, boxed=ret$boxed)

  comp <- ret$free_start >= ret$free_lower & ret$free_start <= ret$free_upper
  if (any(is.na(comp)) || !all(comp))
    stop("'start' is not within model boundary")

  # Check that all the constraints make sense
  if (any(all.equal(names(unlist(ret$lower)),
            names(unlist(ret$start))) != TRUE) ||
      any(all.equal(names(unlist(ret$upper)),
              names(unlist(ret$start))) != TRUE)) {
    list_diff <- function (change, base) {
      ans <- c()

      a <- change[!(change %in% base)]
      if (length(a) > 0) ans <- c(ans, paste0("+", a))

      b <- base[!(base %in% change)]
      if (length(b) > 0) ans <- c(ans, paste0("-", b))

      return (ans)
    }
    stop(
      "lower, upper and start have different structure:\n",
      "start has names ", paste0(collapse=", ", names(unlist(ret$start))),
      "]; lower has names (+/- start) [",
      paste0(
        collapse=", ",
        list_diff(
          names(unlist(ret$lower)),
          base = names(unlist(ret$start)))),
      "]; upper has names (+/- start [",
      paste0(
        collapse=", ",
        list_diff(
          names(unlist(ret$upper)),
          base = names(unlist(ret$start)))),
      "]")
  }
  if (!all(names(unlist(ret$constraints)) %in% names(unlist(ret$start)))) {
    orphan <- setdiff(names(unlist(ret$constraints)),
              names(unlist(ret$start)))
    warning("some parameters in 'constraints' do not exist in 'start': ",
        paste(orphan, collapse = " "))
  }

  ret$call <- match.call()
  ret$memoise <- memoise

  return (structure(ret, class="fd_model"))
}

#' @export
#' @rdname fd_model
fd_clean <- function (object) {
  model <- object
  model$parsedcstr <- NULL
  model$assigncstr <- NULL
  model$mem_diff_cdf <- NULL
}

#' @export
#' @rdname fd_model
fd_clone <- function (object, clean=TRUE) {
  model <- object
  if (!inherits(model, "fd_model"))
    stop("'model' is not an 'fd_model'")
  if (!is.logical(clean))
    stop("'clean' must be logical")
  mdl <- new.env()
  for (n in ls(model, all.names=TRUE)) assign(n, get(n, model), mdl)
  if (clean) fd_clean(mdl)
  structure(mdl, class="fd_model")
}

#' @export
#' @rdname fd_model
update.fd_model <- function (object, data = NULL,
                             fmm, proliferation,
                             addcstr, start, partial, boxed,
                             process,
                             ...) {
  model <- object
  form <- model$form
  if (!missing(addcstr)) {
    if (!is.cstr(addcstr))
      stop("'addcstr' should be a valid constraint")
    if (is.language(addcstr))
      addcstr <- cstrlist(constraints = addcstr)
    form$constraints <- catcstr(form$constraints, addcstr$constraints)
    form$start <- catcstr(form$start, addcstr$start)
    form$lower <- catcstr(form$lower, addcstr$lower)
    form$upper <- catcstr(form$upper, addcstr$upper)
  }
  if (missing(fmm))
    fmm <- get_fmm(model$fmm$name, data)
  if (missing(proliferation))
    proliferation <- get_proliferation(model$pro$name, data)
  if (missing(partial))
    partial <- model$partial
  if (missing(boxed))
    boxed <- model$boxed
  if (missing(process))
    process <- model$process
  mdl <- fd_model(fmm=fmm,
          proliferation=proliferation,
          data=data,
          constraints = form,
          partial = partial,
          boxed=boxed,
          process = process)
  if (!missing(start)) {
    if (!is.null(coef(start)))
      st_ <- colMeans(rbind(coef(start)))
    else st_ <- start
    st <- start(mdl)
    if (!all(names(st_) %in% names(st)))
      message("some 'start' params are constrained")
    st[intersect(names(st_), names(st))] <-
      st_[intersect(names(st_), names(st))]
  } else {
    st <- match.fun("start")(mdl)
  }
  mdl$start <- mdl$trans(relist(st, mdl))
  mdl$free_start <- getfreepar(mdl, mdl$start)
  mdl
}

#' @keywords internal
#' @export
print.fd_model <- function (x, ...) {
  model <- x
  cat(sep="", "Call: ", paste(deparse(model$call), collapse="\n"), "\n\n")
  cat(sep="", "'fd_model' object\n")
}

#' @keywords internal
#' @export
summary.fd_model <- function (object, ...) {
  model <- object
  ans <- list()
  ans$call <- model$call
  cstr <- model$parsedcstr$value
  if (!is.null(cstr)) {
    names(cstr) <- make.names(model$parsedcstr$prefix)
  }
  start <- unlist(model$start)
  start_nm <- names(start)
  ans$boxed <- model$boxed
  ans$partial <- model$partial
  ans$processname <- model$processname
  ans$params <- data.frame(Start = start,
               Lower = unlist(model$lower)[start_nm],
               Upper = unlist(model$upper)[start_nm])
  if (!is.null(cstr))
    ans$params$Constraints <- cstr[start_nm]
  ans$orphan <- setdiff(names(cstr), start_nm)
  ans$fmm <- model$fmm
  ans$pro <- model$pro
  structure(ans, class = "summary.fd_model")
}

#' @keywords internal
#' @export
print.summary.fd_model <- function (x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"),
    "\n", sep = "")

  cat("\nBoxed:", x$boxed, "\n")
  cat("Partial:",
    if (length(x$partial) == 0) "none"
    else paste(x$partial, collapse=","),
    "\n")
  cat("Process:", x$processname, "\n")
  cat("\nFMM:\n")
  print(x$fmm)
  cat("\nProliferation:\n")
  print(x$pro)

  cat("\nParameters:\n")
  print(x$params)
  if (length(x$orphan) > 0L)
    cat("\n(Constraint orphans: ", x$orphan, ")\n")
}

# Misc functions ----------------------------------------------------------

#' @name fd_model-functions
#'
#' @title Miscellaneous functions for the optimization of an FD model.
#'
#' @description These functions can be used when calling nonlinear optimization
#'   functions such as \code{\link[stats]{nls}}, \code{\link[nlme]{nlsList}},
#'   \code{\link[nlme]{nlme}} or \code{\link[nlme]{gnls}} (alternatively, the
#'   wrapper \code{\link{fd_nls}} can be used).
#'
#' @param x An \code{\link{fd_model}} object.
#' @param skeleton An \code{\link{fd_model}} object.
#' @param mgen Maximum number of generations (if \code{NULL}, the default of the
#'   proliferation model is taken).
#' @param globalname The name of an \code{\link{fd_model}} object sitting in the
#'   global environment (that is, a \code{\link{character}} object).  See below
#'   for the rationale.
#' @param flesh A named vector of free parameters, such as returned by
#'   \code{start(skeleton)}.
#' @param free Should the starting, lower or upper values returned by
#'   \code{start}, \code{lower} and \code{upper} be given only for the free
#'   parameters (in this case a named vector is returned) or for all the
#'   parameters, including the constrained ones (in this case a list is
#'   returned)?
#' @param structured An unconstrained, transformed and structured list of
#'   parameters (e.g., as returned by \code{relist}).
#' @param ... Not used.
#'
#' @return \code{fd_formula} builds a formula suitable for using with
#' \code{\link[stats]{nls}}, \code{\link[nlme]{nlsList}},
#' \code{\link[nlme]{nlme}} or \code{\link[nlme]{gnls}}. \code{start} returns
#' the default starting parameters for the FD model, either as a named vector of
#' free parameters (the default) or an unconstrained, transformed and structured
#' list of parameters.  \code{lower} and \code{upper} are two new generic
#' functions, aligned with the way \code{start} is defined in the
#' \code{\link{stats}} package (originally with time-series in mind), and for an
#' \code{fd_model} return the lower and upper bounds along the same modalities.
#' \code{relist} takes a named vector of free parameters \code{flesh} such as
#' the ones returned by \code{start}, \code{lower} and \code{upper} (when
#' \code{free=TRUE}) and turns them into an unconstrained, transformed and
#' structured list of parameters. \code{fd_freepar} performs the inverse
#' operations and from a the result of \code{relist} returns a named vector of
#' free parameters.
#'
#' @section Why a "global name" instead of a regular object:
#'
#'   Because \code{\link[stats]{nls}}-like functions and
#'   \code{\link[nlme]{nlme}} use different scoping rules, great care must be
#'   exercised concerning the scoping of the functions called in a formula.
#'   Indeed, for \code{nls}, the functions are first looked at in the
#'   environment of the formula, then in the \code{data} parameter if it is a
#'   list or an environment, and finally in the global environment.  For
#'   \code{nlme}, however, the scope is \emph{always} the global environment
#'   (Lumley 2003) and \code{data} can only be a \code{data.frame} for obvious
#'   reasons (it is cut down in many pieces according to the various levels of
#'   grouping).  In the end, to stress this point, we decided that not an
#'   \code{fd_model} object should be passed along to \code{formula} but the
#'   \emph{name} of an \code{fd_model} (that is, a \code{\link{character}}
#'   object) sitting in the global environment.  This implementation, of course,
#'   has many drawbacks, including the mandatory use of global variables, but it
#'   seems to be the only sensible one when using such complex nonlinear models
#'   within the formula paradigm.  To avoid such a contortion, general
#'   likelihood optimizers such as \code{\link[stats4]{mle}} or
#'   \code{\link[bbmle]{mle2}} could well be used instead of \code{nls} (see
#'   \code{\link{fd_minuslogl}}).  Unfortunately, no such alternative exists for
#'   either \code{\link[nlme]{nlme}} or \code{\link[nlme]{gnls}}.
#'
#' @examples
#' data(FdSTyphimuriumWTC57)
#' dat <- cutoff(FdSTyphimuriumWTC57)
#' mdl <<- fd_model(data=dat, constraints=attr(dat, "bestcstr"))
#'
#' rbind(start = unlist(relist(start(mdl), mdl)),
#'       lower = unlist(relist(lower(mdl), mdl)),
#'       upper = unlist(relist(upper(mdl), mdl)))
#'
#' # fd_formula can be fed directly to, e.g., nls
#' # (alternatively, fd_nls can be used: it is essentially a wrapper
#' # around nls with additional checks and an automatic, albeit far from
#' # perfect, determination of mgen)
#' fd_formula("mdl")
#'
#' control <- list(maxit = 1L)
#' \dontrun{
#' control <- NULL
#' }
#'
#' # Let's use the "port" algorithm of 'nls'
#' \dontrun{
#' fit <- nls(fd_formula("mdl"), data = FdSTyphimuriumWTC57,
#'            algorithm="port", start=start(mdl),
#'            lower=lower(mdl), upper=upper(mdl),
#'            control=control)
#' # Error in nls(fd_formula("mdl"), data = FdSTyphimuriumWTC57,
#' # algorithm = "port",  :
#' #     Convergence failure: iteration limit reached without convergence (10)
#' }
#'
#' # As it is likely to fail, we can use the global search
#' # algorithm 'GenSA' as well, again in an 'nls'-like framework
#' source(system.file("contrib", "nlsSA.R", package="fluodilution"))
#' fit <- nlsSA(fd_formula("mdl"), data = FdSTyphimuriumWTC57,
#'              start=start(mdl),
#'              lower=lower(mdl), upper=upper(mdl),
#'              control=control)
#'
#' # A third option is to use fd_minuslogl, see the relevant
#' # page for an example.
#'
#' @references Lumley T (2003).  \emph{Standard nonstandard evaluation rules.}
#' \url{http://developer.r-project.org/nonstandard-eval.pdf}
NULL

#' @export
#' @rdname fd_model-functions
fd_formula <- function (globalname, mgen = NULL) {
  object <- get_global_model(globalname)
  colnames <- c("Category", "Time", "Timepoint", "a", "b", "y",
"Inoculum", "Weight", "Type")
  if (is.null(mgen)) mgen <- "NULL"
  rhs <- paste0("fd_predict_mat(", globalname, ")(",
    "data.frame(",
    paste(paste0(setdiff(colnames, "y"),
           " = ", setdiff(colnames, "y")),
        collapse=", "),
    "), cbind(", paste(names(start(object)), collapse=", "),
    "), mgen = ", mgen, ")")
  ans <- as.formula(paste0("y ~ ", rhs))
  environment(ans) <- globalenv()
  ans
}

#' @export
#' @rdname fd_model-functions
start.fd_model <- function (x, free = TRUE, ...) {
  object <- x
  if (!inherits(object, "fd_model"))
    stop("'object' is not an 'fd_model'")
  if (free) object$free_start else object$trans_inverse(object$start)
}

#' @export
#' @rdname fd_model-functions
lower <- function (x, ...) UseMethod("lower")

#' @export
#' @rdname fd_model-functions
lower.fd_model <- function (x, free = TRUE, ...) {
  object <- x
  if (!inherits(object, "fd_model"))
    stop("'object' is not an 'fd_model'")
  if (free) object$free_lower else object$trans_inverse(object$lower)
}

#' @export
#' @rdname fd_model-functions
upper <- function (x, ...) UseMethod("upper")

#' @export
#' @rdname fd_model-functions
upper.fd_model <- function (x, free = TRUE, ...) {
  object <- x
  if (free) object$free_upper else object$trans_inverse(object$upper)
}

#' @export
#' @rdname fd_model-functions
relist.fd_model <- function (flesh, skeleton) {
  param <- flesh
  model <- skeleton
  if (is.null(names(param))) {
    names(param) <- model$freepar
  } else if (any(all.equal(names(param), names(model$boxed_lower)) != TRUE)) {
    stop("'param' is not structured according to the model (wrong names)")
  }
  if (!model$boxed) {
    param <- box.param(param, model$boxed_lower, model$boxed_upper)
  }
  object <- model$start
  eval(fd_assigncstr(model))
  model$trans_inverse(object)
}

#' @export
#' @rdname fd_model-functions
fd_freepar <- function (structured, x) {
  object <- x
  structured <- object$trans(structured)
  structured <- applycstr(structured, object$form$constraints)
  fp <- unlist(structured)[object$freepar]
  if (!object$boxed) {
    fp <- unbox.param(fp, object$boxed_lower, object$boxed_upper)
  }
  return (fp)
}

# Fitting -----------------------------------------------------------------

#' @export
#' @keywords internal
fd_process_default <- function (env) {
  # Predict hists, if any
  for (i in 1:2) {
    type <- c("hists", "hists_lost")[i]
    if (NROW(env$splitted[[type]]) > 0) {
      # First, you need to get to the CDF
      cdf <- mem_diff_cdf(env$model, env$splitted[[type]],
                env$psi, env$mgen)

      # Ugly hack
      model_pop <- apply(
        simplify2array(env$prolif[[c("live_pop", "lost_pop")[i]]]),
        c(1, 3),
        function (x) list(x))
      matching <- match(levels(env$splitted[[type]]$Category),
                env$model$pro$categories)
      cur_pred <- t(sapply(
        model_pop[cbind(findInterval(env$splitted[[type]]$Time,
                       env$times),
                matching[env$splitted[[type]]$Category])],
        "[[", 1))

      # Then mix
      env$splitted[[type]]$y <-
        env$model$fmm$htrans$trans(rowSums(cur_pred * cdf))
    }
  }

  # Predict proportions, if any
  for (i in 1:2) {
    type <- c("props", "props_lost")[i]
    if (NROW(env$splitted[[type]]) > 0) {
      model_pop <- simplify2array(
        env$prolif[[c("live_pop", "lost_pop")[i]]])
      should_mGen <- max(env$splitted[[type]]$a)
      if (should_mGen > NCOL(model_pop) - 1) {
        stop("'mgen' not big enough for 'prop': should be at least ",
           should_mGen)
      }
      matching <- match(levels(env$splitted[[type]]$Category),
                        env$model$pro$categories)
      env$splitted[[type]]$y <-
        model_pop[cbind(
          findInterval(env$splitted[[type]]$Time, env$times),
          env$splitted[[type]]$a + 1,
          matching[env$splitted[[type]]$Category])]

    }
  }

  # Do the same with cell counts, if any
  for (i in 1:2) {
    type <- c("Ns", "Ns_lost")[i]
    if (NROW(env$splitted[[type]]) > 0) {
      matching <- match(levels(env$splitted[[type]]$Category),
                        env$model$pro$categories)
      Pred_Ns <- env$prolif[[type]][
        cbind(findInterval(env$splitted$Ns$Time, env$times),
              matching[env$splitted$Ns$Category])]
      LPred_Ns <- env$model$fmm$cctrans$trans(Pred_Ns * 10^env$psi$c0)
      if (length(LPred_Ns) != NROW(env$splitted[[type]])) {
        stop("internal error")
      }
      env$splitted[[type]]$y <- LPred_Ns
    }
  }
}

#' Predict the histograms and cell counts of an \code{fd_data} dataset from a
#' set of model parameters.
#'
#' Predict the histograms and cell counts of an \code{fd_data} dataset from a
#' set of model parameters. \code{fd_predict} is used internally when
#' \code{\link{fd_formula}}, \code{\link{fd_simulate}} or
#' \code{\link{fd_minuslogl}} are called, but it can also be called directly.
#'
#' @param model An \code{\link{fd_model}} object.
#' @param param A named vector of free parameters, such as returned by
#'   \code{\link{start}}.
#' @param data An \code{\link{fd_data}} dataset.
#' @param mgen Maximum number of generations.  If \code{NULL}, taken from the
#'   default found in the proliferation model.
#' @param stop_boundary Whether the function should stop when \code{param} is
#'   out of bounds or simply return very large predicted values.
#' @return The original \code{fd_data}, not necessarily in the same order, with
#'   the \code{"y"} column filled with the predicted values.
#' @family simulation-related entries
#' @export
#' @examples
#' library(dplyr)
#'
#' data(FdSTyphimuriumWTC57)
#' dat <- cutoff(FdSTyphimuriumWTC57)
#' dat <- dat[order(dat$Type, dat$Inoculum), ]
#' mdl <- fd_model(data=dat, constraints=attr(dat, "bestcstr"))
#' pred <- fd_predict(mdl, start(mdl), dat)
#' if (any(all.equal(pred %>% dplyr::select(-y),
#'                   as.data.frame(dat) %>% dplyr::select(-y)) != TRUE))
#'   stop("should be equal")
#' dat_residuals <- pred$y - dat$y
fd_predict <- function (model, param, data, mgen=NULL, stop_boundary = FALSE) {
  # Check parameters
  if (any(!is.finite(param))) {
    warning("'theta' is not finite")
    data$y <- Inf
    return (data)
  }
  if (any(param < lower.fd_model(model) | param > upper.fd_model(model))) {
    if (stop_boundary)
      stop("parameters out of bounds")
    return (transform(data, y = 1e10))
  }
  if (inherits(data, "fd_data")) {
    data <- as.data.frame(data)
  }
  if (!all(levels(data$Category) %in% model$pro$categories))
    stop("the data do not have the right categories")

  # Relist and apply model
  env <- new.env()

  env$model <- model

  env$relisted <- relist.fd_model(param, model)
  env$theta <- env$relisted$pro
  env$psi <- env$relisted$fmm

  env$times <- sort(unique(data$Time))
  env$mgen <- mgen
  env$prolif <- model$pro$model(env$theta, env$times, env$mgen)

  # Sorting: of course additional processing time, but we cannot take the risk
  data <- data[order(data$Type, data$Inoculum), ]
  env$splitted <- split(data, data$Type)

  # Process data
  match.fun(model$process)(env)

  # Slightly better, in terms of performance, not to unsplit (even if one
  # then has to sort the data)
  data$y <- unlist(lapply(env$splitted, function (cur) cur$y))

  return (data)
}


# Internal ----------------------------------------------------------------

getfreepar <- function (model, structured, boxed=model$boxed) {
  structured <- applycstr(structured, model$form$constraints)
  fp <- unlist(structured)[model$freepar]
  if (!boxed) {
    fp <- unbox.param(fp, model$boxed_lower, model$boxed_upper)
  }
  return (fp)
}

fd_assigncstr <- function (model) {
  if (is.null(model$assigncstr)) {
    if (is.null(model$parsedcstr)) {
      model$parsedcstr <- parsecstr(model$form$constraints)
    }
    model$assigncstr <- paste(
      paste0("object$",
        gsub("\\.", "\\$", model$freepar),
        " <- param[['", model$freepar, "']]", collapse="\n"),
      "\n", assigncstr("object", model$parsedcstr)
    )
    if (!is.null(model$partial)) {
      model$assigncstr <- paste0(
        model$assigncstr,
        "\nobject$pro <- object$pro[", deparse(model$partial), "]\n"
      )
    }
    model$assigncstr %<>% parse(text=.)
  }
  model$assigncstr
}

mem_diff_cdf <- function (model, ...) {
  if (is.null(model$mem_diff_cdf)) {
    model$mem_diff_cdf <- function (df, psi, mgen) {
      if (is.null(mgen)) mgen <- model$pro$mgen
      do.call(rbind, lapply(unique(df$Inoculum), function (inoc) {
        ss <- subset.data.frame(df, Inoculum == inoc)

        # as.character: additional safeguard...
        psi$i <- as.character(inoc)
        model$fmm$diff_cdf(
          psi = psi,
          a = ss$a, b = ss$b,
          gen = 0:mgen
        )
      }))
    }
    if (is.null(model$memoise)) {
      # Backward compatibility
      model$memoise <- TRUE
    }
    if (model$memoise)
      model$mem_diff_cdf <- memoise::memoise(model$mem_diff_cdf)
  }
  model$mem_diff_cdf(...)
}

#' @export
#' @keywords internal
fd_predict_mat <- function (model) {
  structure(function (data, params, mgen=NULL) {
    if (!is.matrix(params)) stop("'params' must be a matrix")

    # To preserve the original order
    data$Id <- 1:NROW(data)

    # First, group the thetas
    uparams <- unique(params)

    df_pred <- do.call(rbind, lapply(1:NROW(uparams), function (i) {
      param <- uparams[i, ]

      rows <- colSums(t(params) != param) == 0
      fd_predict(model, param, data=data[rows, ], mgen=mgen)
    }))

    df_pred[order(df_pred$Id), "y"]
  },
  class="fd_predict_mat",
  .model = model)
}

get_global_model <- function (globalname) {
  if (!is.character(globalname) || !exists(globalname, envir=globalenv()))
    stop("'globalname' must be a character string and sits in ",
       "the global environment")
  object <- get(globalname, envir = globalenv())
  if (!inherits(object, "fd_model"))
    stop("'", globalname, "' is not an 'fd_model'")
  object
}
