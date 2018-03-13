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

# AF FMM ------------------------------------------------------------------

#' @export
#' @rdname finite-mixture
#' @include fmm-.R
fd_fmm_af <- fd_fmm(
  name = "af",
  new = function (sd0, m0, ...) {
    if (is.null(sd0) || is.null(m0))
      stop("either 'sd0' or 'm0' is NULL")
    if (length(m0) != 1 || length(sd0) != 1)
      stop("multiple m0/sd0 not supported with fd_fmm_af")

    list(
      sd0 = sd0,
      m0 = m0,
      ...
    )
  },
  diff_cdf = function (self, psi, a, b, gen) {
    if (!is.integer(gen) || any(gen < 0))
      stop("'gen' must specify valid generations")
    if (is.null(psi$sdaf))
        stop("psi$sdaf' must be specified")

    if (psi$sdaf == 0)
      return ((t(t(plnorm(outer(b, 2 ^ gen), self$m0, self$sd0) -
        plnorm(outer(a, 2 ^ gen), self$m0, self$sd0)))))

    result <- matrix(0, ncol = length(gen), nrow = length(a))
    for (j in 1:NCOL(result)) {
      for (i in 1:NROW(result)) {
        int <- try(
          integrate(function (x) {
            plnorm(
              (b[i] - x) *
                2 ^ gen[j],
              self$m0, self$sd0) * dnorm(x, 0, psi$sdaf) -
            plnorm(
              (a[i] - x) *
                2 ^ gen[j],
              self$m0, self$sd0) * dnorm(x, 0, psi$sdaf)
          },
          -3 * psi$sdaf, min(a[i], 3 * psi$sdaf), subdivisions=1000L),
          silent=TRUE)
        if (inherits(int, "try-error")) {
          warning(as.vector(int))
          result[i, j] <- 0
        } else {
          result[i, j] <- int$value
        }
      }
    }
    result
  },
  constrain = function (self, object, type) {
    # nolint start
    catcstr(
       switch(type,
         start = ~ fmm:{c0 <- 5; sdaf <- 100},
         lower = ~ fmm:{c0 <- 0; sdaf <- 10},
         upper = ~ fmm:{c0 <- 10; sdaf <- 400}), object)
    # nolint end
  }
)
