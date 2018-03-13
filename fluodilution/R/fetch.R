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

# Internal ----------------------------------------------------------------

get_fmm <- function (object, data = NULL) {
  if (is.character(object)) {
    s <- paste0("fd_fmm_", object)
    object <- match.fun(s)
  }
  if (is.function(object)) {
    arglst <- attr(data, "fmm")
    if (is.null(arglst)) {
      message("A character string for an FMM object has been given, ",
          "but 'data' is 'NULL': ",
          "using default values for sd0/m0")
      object <- object(m0 = 10, sd0 = 0.3)
    } else object <- do.call(object, arglst)
  }
  if (!inherits(object, "fd_fmm"))
    stop("FMM object does not inherit \"fd_fmm\"")
  object
}

get_proliferation <- function (object, data = NULL) {
  if (is.character(object)) {
    s <- paste0("fd_proliferation_", object)
    object <- match.fun(s)
  }
  if (is.function(object)) {
    if (is.null(data)) {
      message("A character string for a proliferation object has ",
          "been given, but 'data' ",
          "is 'NULL': using default values")
      object <- object()
    } else {
      if (NROW(data) == 0)
        stop("a character string for a proliferation object has ",
           "been given, but 'data' ",
           "is empty")
      object <- object(categories = levels(data$Category))
    }
  }
  if (!inherits(object, "fd_proliferation"))
    stop("proliferation model object does not inherit ",
       "\"fd_proliferation\"")
  object
}

get_name <- function (object) {
  UseMethod("get_name")
}

get_name.default <- function (object) {
  if (is.character(object)) make.names(object)
  else make.names(deparse(substitute(object)))
}

get_name.fd_proliferation <- function (object) {
  object$name
}

get_name.fd_fmm <- function (object) {
  object$name
}