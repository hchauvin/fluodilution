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

#' Constrain some parameters to specific values.
#'
#' Fully-fledged fluorescence dilution models feature a large amount of
#' parameters, not all of which are practically identifiable to an acceptable
#' degree.  Therefore, it is important to constrain them: for instance, by
#' requiring that a gamma distribution should have a coefficient of variation of
#' 1, turning it into an exponential distribution, or that different
#' compartments share the same time to division.  We devised a set of tools to
#' help this constraining.
#'
#' \code{cstrlist} is used to bring together the specification of overall
#' constraints and the starting values and lower and upper bounds (they are also
#' given using constraints).  The result can be used with, e.g.,
#' \code{\link{fd_model}}.  \code{fixcstr} is an alternative way to constructing
#' constraints: it returns the original set of constraints augmented by the
#' fixing of some of the parameters to the values found in \code{fit}.
#' \code{catcstr} is used to concatenate together an arbitrary number of
#' constraints (it does not work on a \code{cstrlist} object).
#'
#' @param constraints A set of constraints (format detailed below) that should
#'   apply to the starting parameters in an optimization, the lower and upper
#'   bounds and the result of the optimization.
#' @param start A set of constraints that apply only to the starting parameters.
#' @param lower A set of constraints that apply only to the lower bounds.
#' @param upper A set of constraints that apply only to the upper bounds.
#' @param cstr A set of constraints (format detailed below) or a \code{cstrlist}
#'   object.
#' @param fit Either an object the \code{\link[stats]{coef}} method can be used
#'   with (e.g., the result of \code{nls}) or a named vector of coefficients.
#' @param before Whether the new constraining should apply before (default) or
#'   after the current constraining.
#' @param ... Constraints (or lists of constraints) to "concatenate" together.
#'   The leftmost is concatenated before the rightmost and thus would apply
#'   before.  If lists are given, they are "combinatorially" concatenated (as
#'   with a Cartesian product).
#' @param drop If the returned list of constraints have only one element, return
#'   this element instead of a one-element list.
#'
#' @return \code{makecstr} returns a \code{cstrlist} object suitable for use
#'   with, e.g., \code{\link{fd_model}}.  \code{fixcstr} returns either the new
#'   set of constraints or a \code{cstrlist} object, depending on the type of
#'   argument \code{cstr}.  \code{catcstr} returns the resulting concatenation
#'   (a new set of constraints).
#'
#' @section Constraint format:
#' In EBNF,
#' \preformatted{
#'    constraint ::= [ term [ '+' ] constraint ]
#'          term ::= [ prefix ':' ] '{' assignments '}'
#'                 | name
#'                 | 'NULL'
#'                 | [ prefix ':' ] '(' constraint ')'
#'        prefix ::= prefix_atom [ ':' prefix ]
#'   prefix_atom ::= language
#'                 | '(' language {'/' language} ')'
#'   assignments ::= language '<-' (language | '..free..')
#'                   [ ';' assignments ]}
#'
#' This format is specifically tailored to constraining the value of nested
#' lists with structures that are somewhat "parallel".  For instance, in
#' fluorescence dilution, the parameters for the two-tier model are structured
#' as a list with two elements, \code{fmm} (parameters for the finite mixture
#' model) and \code{pro} (parameters for the proliferation model). \code{pro} is
#' structured as a named list, with each element giving the parameters for a
#' specific compartment, in their flow order. In turn, each compartment is
#' parametrized by gamma distributions, which are represented as lists of three
#' elements (a detailed account of this structure is given in
#' \code{\link{proliferation}}, and for FMM in \code{\link{finite-mixture}}).
#' Therefore, the constraining can be "packed" more efficiently, benefiting from
#' this "parallel" organization.
#'
#' In the end, with nested lists,
#' constraining a set of parameters \code{params} amount to using such
#' expressions:
#' \verb{
#' }\code{params$pro$One$f0$delta <- 1}.
#'
#' Such a constraining can be written in our format as
#' \code{pro:One:f0:{delta <- 1}} or
#' alternatively \code{pro:One:{f0$delta <- 1}}: the left-hand side is
#' simply prefixed.
#'
#' To make expressions lighter, it is possible to combine
#' \verb{
#' }\code{params$pro$One$f0$delta <- 1; params$pro$One$g0$delta <- 1}
#' \verb{
#' } into \code{pro:One:(f0/g0):{delta <- 1}}.
#'
#' For the proliferation model, the \code{all} prefix is expanded:
#' for three compartments \code{One}, \code{Two} and \code{Three},
#' \code{pro:all:(f0/g0):{delta <- 1}} expands into
#' \code{pro:(One/Two/Three):(f0/g0):{delta <- 1}}.
#'
#' It is also possible to set more than one parameter inside the brackets:
#' \verb{
#' }\code{pro:all:f0:{delta <- 1; ss <- 0.5}}.
#'
#' Previously constrained parameters can be "unconstrained" and constraints
#' linked together
#' using \code{+}:
#' \code{pro:all:f0:{delta <- 1; ss <- 0.5} +
#' pro:One:f0:{delta <- ..free..}}.
#'
#' A parameter that has been previously freed can be constrained again:
#' \preformatted{pro:all:f0:{delta <- 1; ss <- 0.5} +
#' pro:One:f0:{delta <- ..free..} +
#' pro:one:f0:{delta <- 0.5}}
#' and in general the rightmost value prevails.  For instance,
#' in the following \code{delta} is constrained to \code{0.5} in
#' compartment \code{One}
#' but to \code{1} in every other compartment:
#' \verb{
#' }\code{pro:all:f0:{delta <- 1; ss <- 0.5} +
#' pro:One:f0:{delta <- 0.5}}.
#'
#' To bind parameters together, it is possible to use the special
#' placeholders \code{.L1}
#' (stands for the list containing the left-hand parameter in the assignment),
#' \code{.L2} (the list just above) and \code{.L3} (one additional level up).
#' For example,
#' to bind probabilities together one can write
#' \code{pro:all:{p0 <- .L1$p}},
#' to bind distributions inside the same compartment
#' \verb{
#' }\code{pro:all:{f0$mm <- .L2$f$mm}} or \code{pro:all:f0:{mm <- .L2$f$mm}}
#' \verb{
#' } and
#' to bind distributions across compartments
#' \verb{
#' }\code{pro:all:{f$mm <- .L3$One$f$mm} + pro:One:f:{mm <- ..free..}}.
#'
#' When a term is a name beginning with \code{#} or in the form
#' \code{listname$name}, it is expanded in the
#' global environment:
#' for instance, if \code{`#noss`}
#' is defined in the global environment as
#' \verb{
#' }\code{`#noss` <- ~ pro:all:(f0/g0/f/g):{ss <- 0.5}},
#' \verb{
#' } then
#' \code{`#noss` + pro:One:(f0/g0):{delta <- 1}} is expanded to
#' \verb{
#' }\code{pro:all:(f0/g0/f/g):{ss <- 0.5} + pro:One:(f0/g0):{delta <- 1}}
#' (see the \emph{examples} section for more).
#'
#' Notice also that
#' a tilde \code{~} is used to indicate \emph{R} not to evaluate the expression
#' on the right of it
#' (it is a valid formula but not a meaningful statement).
#'
#' We believe that this format offers a very powerful, flexible and terse
#' approach to constraining complex and large nonlinear models.  It has
#' applications well
#' beyond the realm of fluorescence dilution.
#'
#' @seealso \code{\link{finite-mixture}}, \code{\link{proliferation}},
#' \code{\link{fd_model}}.
#'
#' @examples
#' # Constrain all 'delta' to be 1
#' `#c1` <- ~ pro:all:(f/g/f0/g0):{delta <- 1}
#'
#' # It is also possible to use a set of often-used constraints
#' CC <<- FdCommonConstraints
#'
#' # The constrains can be combined together in two different ways
#' `#c2` <- ~ `#c1` + CC$`#delta_1111`
#' `#c3` <- catcstr(`#c2`, CC$`#noss`, ~ pro:all:{p0 <- .L1$p})
#'
#' # Starting parameters and lower/upper bounds are also
#' # specified using constraints
#' cstrlist(constraints = `#c3`,
#'          start = ~ pro:all:{p <- 0.5},
#'          lower = ~ pro:all:{p <- 0.1},
#'          upper = ~ pro:all:{p <- 0.95})
#'
#' # An alternative way to constructing constraints uses 'fixcstr'
#' `#c4` <- fixcstr(`#c3`, c(pro.One.f.mm = 5))
#'
#' # Finally, previously constrained parameters can be freed
#' `#c5` <- ~ `#c3` + pro:all:{p0 <- ..free..}
#'
#' @rdname constraints
#' @name constraints
NULL

#' @rdname constraints
#' @export
cstrlist <- function (constraints = NULL, start = NULL,
                      lower = NULL, upper = NULL) {
  if ((!is.null(constraints) && !is.language(constraints)) ||
      (!is.null(start) && !is.language(start)) ||
      (!is.null(lower) && !is.language(lower)) ||
      (!is.null(upper) && !is.language(upper)))
    stop("parameters must be valid language elements ",
       "(name, call, expression)")
  ans <- structure(list(constraints = constraints,
           start = start, lower = lower, upper = upper),
        class = "cstrlist")
  if (!is.cstr(ans))
    stop("some constraints are not valid")
  ans
}

#' @rdname constraints
#' @export
fixcstr <- function (cstr, fit, before=TRUE) {
  if (!is.cstr(cstr))
    stop("'cstr' should be a valid constraint")
  if (!is.vector(fit)) {
    if (is.null(coef(fit)))
      stop("'fit' is not a vector nor an object with coefficients")
    fit <- coef(fit)
  }
  if (is.null(names(fit)))
    stop("'fit' does not have named elements")
  added <- paste0(
    "~ {",
    paste(
      paste0(gsub("\\.", "$", names(fit)), " <- ", fit),
      collapse=";\n"),
    "}")
  if (inherits(cstr, "cstrlist")) {
    constraints <- cstr$constraints
  } else {
    constraints <- cstr
  }
  if (before) constraints <- catcstr(as.formula(added), constraints)
  else constraints <- catcstr(constraints, as.formula(added))
  if (inherits(cstr, "cstrlist")) {
    cstr$constraints <- constraints
    cstr
  } else {
    constraints
  }
}

#' @rdname constraints
#' @export
catcstr <- function (..., drop = TRUE) {
  arglst <- list(...)
  arglst <- arglst[!sapply(arglst, is.null)]
  if (length(arglst) == 0) return (NULL)
  if (is.language(arglst[[length(arglst)]]))
    arglst[[length(arglst)]] <- list(arglst[[length(arglst)]])
  if (any(!sapply(arglst[[length(arglst)]], is.language)))
    stop("some arguments are not constraints")
  ans <- Reduce(
    function (x, y) {
      if (is.language(x)) x <- list(x)
      if (any(!sapply(x, is.language)))
        stop("some arguments are not constraints")
      unlist(lapply(x, function (xx)
        lapply(y, function (yy)
          call("+", xx, yy))))
    },
    arglst, right=TRUE
  )
  if (drop && length(ans) == 1) ans[[1L]]
  else ans
}

#' @keywords internal
#' @export
is.cstr <- function (cstr) {
  is.null(cstr) || (inherits(cstr, "cstrlist") && is.list(cstr) &&
    all(names(cstr) %in% c("constraints", "start", "lower", "upper")) &&
    all(sapply(cstr, is.language) | sapply(cstr, is.null))) ||
    is.language(cstr)
}

#' @keywords internal
#' @export
expandcstr <- function (expr, sub = NULL) {
  if (!is.null(sub) && !is.language(sub)) stop("'sub' must be language")
  if (!is.cstr(expr)) stop("'expr' must be a constraint")
  expandcstr_ <- function (expr, sub = NULL) {
    dollar <- FALSE
    if (is.symbol(expr) || !is.language(expr) ||
        (dollar <- length(expr) > 0L && expr[[1L]] == "$")) {
      expr_s <- deparse(expr)
      matched <- expr_s == names(sub)
      if (sum(matched) > 1) stop("More than one substitution possible")
      else if (any(matched)) {
        expr <- sub[[matched]]
        expr_s <- deparse(expr)
      }
      if (stringr::str_sub(expr_s, 1L, 1L) == "#" || dollar) {
        Recall(eval(expr), sub)
      } else expr
    } else if (length(expr) == 0 || expr[[1L]] == "{") {
      expr
    } else {
      for (i in 1:length(expr)) {
        expr[[i]] <- Recall(expr[[i]], sub)
      }
      expr
    }
  }
  expandcstr_(expr, sub)
}

# Internal ----------------------------------------------------------------

unbox.param <- function (param, lower, upper, finite=10) {
  if (length(lower) != length(upper) || length(param) != length(upper))
    stop("'param', 'lower' and 'upper' must have the same length")
  if (any(lower > upper)) stop("lower > upper")
  if (any(param < lower | param > upper))
    stop("'param' is outside boundary")

  param <- (param - lower) / (upper - lower)
  param <- log(param / (1 - param))
  if (!is.na(finite)) {
    # For boundary reaching: give an arbitrary (high) bound instead
    # of Inf.  Graphical investigation of abs(tranh(x)) shows that
    # 4 is probably good (tanh(4) > 0.999) but let's go with 10
    param[!is.finite(param)] <- sign(param[!is.finite(param)]) * finite
  }
  return (param)
}

box.param <- function (param, lower, upper) {
  if (length(lower) != length(upper) || length(param) != length(upper))
    stop("'param', 'lower' and 'upper' must have the same length")
  if (any(lower > upper)) stop("lower > upper")
  param <- exp(param) / (exp(param) + 1)
  param  * (upper - lower) + lower
}

parsecstr <- local({
  parse_prefix <- function (cstr) {
    if (is.symbol(cstr)) {
      deparse(cstr)
    } else if (cstr[[1L]] == "~") {
      parse_prefix(cstr[[2L]])
    } else if (cstr[[1L]] == "/") {
      c(parse_prefix(cstr[[2L]]), parse_prefix(cstr[[3L]]))
    } else if (cstr[[1L]] == "(") {
      parse_prefix(cstr[[2L]])
    } else if (cstr[[1L]] == ":") {
      outer(parse_prefix(cstr[[2L]]), parse_prefix(cstr[[3L]]),
          paste, sep="$")
    } else {
      print(cstr)
      stop("wrong prefix format")
    }
  }
  parsecstr_ <- function (cstr) {
    if (is.null(cstr)) {
      return (list(prefix = NULL, value = NULL))
    } else if (is.expression(cstr)) {
      if (length(cstr) > 1)
        stop("'parsecstr' does not allow expressions with ",
           "multiple terms")
      if (length(cstr) == 0)
        stop("empty expression")
      Recall(cstr[[1L]])
    } else if (is.symbol(cstr) || cstr[[1L]] == "$") {
      Recall(eval(cstr))
    } else if (cstr[[1L]] == "~") {
      Recall(cstr[[2L]])
    } else if (cstr[[1L]] == ":") {
      prefix <- parse_prefix(cstr[[2L]])
      ret <- Recall(cstr[[3L]])
      list(
        prefix = unlist(
          lapply(prefix,
               function (cur) paste0(cur, "$", ret$prefix))),
        value = rep(ret$value, length(prefix))
      )
    } else if (cstr[[1L]] == "{") {
      if (length(cstr) == 1)
        stop("empty expression")
      ret <- sapply(as.list(cstr[-1]), function (cur) {
        if (cur[[1L]] != "<-") {
          print(cstr)
          stop("terms must be of the form 'a <- b'")
        }
        c(prefix = deparse(cur[[2L]]), value = deparse(cur[[3L]]))
      })
      list(prefix = as.vector(ret[1, ]), value = as.vector(ret[2, ]))
    } else if (cstr[[1L]] == "+") {
      ret1 <- Recall(cstr[[2L]])
      if (length(cstr) > 2) {
        ret2 <- Recall(cstr[[3L]])
        list(prefix = c(ret1$prefix, ret2$prefix),
           value = c(ret1$value, ret2$value))
      } else {
        ret1
      }
    } else if (cstr[[1L]] == "(") {
      Recall(cstr[[2L]])
    } else {
      print(cstr)
      stop("unknown syntax")
    }
  }

  function (cstr) {
    ans <- parsecstr_(cstr)
    for (free in which(ans$value == "..free..")) {
      rem <- which(ans$prefix[1:free] == ans$prefix[free])
      ans$prefix[rem] <- NA
      ans$value[rem] <- NA
    }
    if (!is.null(ans$prefix)) {
      ans$prefix <- na.omit(ans$prefix)
      ans$value <- na.omit(ans$value)
    }
    ans
  }
})


applycstr <- function (object, cstr) {
  eval(cstr %>%
    parsecstr() %>%
    assigncstr("object", .) %>%
    parse(text=.))
  return (object)
}

assigncstr <- function (object, parsed) {
  if (length(parsed$prefix) > 0L) {
    L1p <- "\\$[^$]+$"
    L1 <- paste0("object$", sub(L1p, "", parsed$prefix))
    L1[L1 == parsed$prefix] <- "NULL"
    L2p <- "\\$[^$]+\\$[^$]+$"
    L2 <- paste0("object$", sub(L2p, "", parsed$prefix))
    L2[L2 == parsed$prefix] <- "NULL"
    L3p <- "\\$[^$]+\\$[^$]+\\$[^$]+$"
    L3 <- paste0("object$", sub(L3p, "", parsed$prefix))
    L3[L3 == parsed$prefix] <- "NULL"
    with(parsed, paste(
      paste0(
        object, "$", prefix,
        " <- with(list(", object,
        ", .L1 = ", L1, ", .L2 = ", L2, ", .L3 = ", L3, "),",
        value, ")"
      ),
      collapse = "\n"))
  } else character(0)
}