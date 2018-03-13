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

#' @export
#' @rdname finite-mixture
#' @include fmm-.R
fd_fmm_gaussian <- fd_fmm(
  name = "gaussian",
  new = function (sd0, m0, ...) {
    if (is.null(sd0) || is.null(m0))
      stop("either 'sd0' or 'm0' is NULL")
    list(
      sd0 = sd0,
      m0 = m0,
      ...
    )
  },
  diff_cdf = function (self, psi, a, b, gen) {
    if (!is.integer(gen) || any(gen < 0))
      stop("'gen' must specify valid generations")
    if (length(self$m0) == 1) {
      cur_m0 <- self$m0
      cur_sd0 <- self$sd0
    } else {
      if (is.null(psi$i))
        stop("'psi$i' must be specified")
      cur_m0 <- self$m0[psi$i]
      cur_sd0 <- self$sd0[psi$i]
    }
    if (is.na(cur_m0) || is.na(cur_sd0))
      stop("'psi$i' is not a valid inoculum")
    t(t(plnorm(outer(b, 2 ^ gen), cur_m0, cur_sd0) -
        plnorm(outer(a, 2 ^ gen), cur_m0, cur_sd0)))
  },
  constrain = function (self, object, type) {
    # nolint start
    catcstr(
       switch(type,
         start = ~ fmm:{c0 <- 5},
         lower = ~ fmm:{c0 <- 0},
         upper = ~ fmm:{c0 <- 10}), object)
    # nolint end
  }
)
