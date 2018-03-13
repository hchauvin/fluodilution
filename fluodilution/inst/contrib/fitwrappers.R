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

assign_model <- function (fit, data, addcstr, mgen) {
  if (!is.null(addcstr)) {
    mdl_ <- model(fit)
    if (is.null(mdl_)) mdl_ <- model(fit[[1L]])
    if (is.null(mdl_)) stop("'fit' does not have any model that goes along!")
    mdl <<- update(mdl_, data, addcstr=addcstr, start = fit)
    start_coefs <- start(mdl)
  } else {
    mdl <<- update(model(fit), data)
    start_coefs <- coef(fit)
  }
  if (is.null(mgen)) mgen <- fit$.mgen
  form <- fd_formula("mdl", mgen = mgen)
  list(form = form, start = start_coefs)
}

#' Title
#'
#' @param fit
#' @param data
#' @param random
#' @param addcstr
#' @param mgen
#' @param verbose
#' @param control
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
fd_nlme <- function (fit, data, random, addcstr = NULL, mgen = NULL,
                     verbose = TRUE,
                     control = NULL,
                     ...) {
  parm <- assign_model(fit, data, addcstr, mgen)
  fixed <- as.formula(paste0(paste(names(parm$start), collapse = " + "), " ~ 1"))
  ans <- nlme(parm$form,
              data,
              fixed = fixed,
              random = random,
              start = parm$start,
              verbose = verbose,
              control = control,
              ...)
  ans$.model <- mdl
  ans$.modelname <- "mdl"
  class(ans) <- c("fd_nls", class(ans))
  ans
}

#' Title
#'
#' @param fit
#' @param data
#' @param random
#' @param addcstr
#' @param mgen
#' @param control
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
fd_nlsList <- function (fit, data, random, addcstr = NULL, mgen = NULL,
            control = NULL, ...,
            nlsListfun = nlsList) {
  parm <- assign_model(fit, data, addcstr, mgen)
  ans <- nlsListfun(parm$form, data,
                    lower = lower(mdl),
                    upper = upper(mdl),
                    start = parm$start,
                    control = control,
                    ...)
  for (i in seq_along(ans)) {
    ans[[i]]$.model <- mdl
    ans[[i]]$.modelname <- "mdl"
  }
  class(ans) <- c("fd_nls", class(ans))
  ans
}
