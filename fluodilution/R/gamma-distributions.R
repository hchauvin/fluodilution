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

#' Reparametrization and evaluation of shifted gamma distributions.
#'
#' Our proliferation models use gamma distributions for times to division/death.
#' They are parametrized by their mean, coefficient of variation and a
#' parametrization of the skewness. These functions help in the manipulation of
#' such distributions.
#'
#' Gamma distributions are usually parametrized using scale and shape
#' parameters. Additionally, a shifted gamma distribution makes use of a
#' location parameter. It is possible to use a more "natural" parametrization
#' from the mean, coefficient of variation and an alternative reparametrization
#' of the location parameter (we could think of it as a \emph{second-order}
#' shape).
#'
#' @return See \emph{details}.
#'
#' @details \code{fd_gamma_dist} takes the parameters of a shifted gamma
#' distribution in their non-reparametrized (natural) form and returns a one-row
#' matrix with those parameters. \code{fd_reparam_gamma_dist} just returns a
#' list made with its parameters, in the reparametrized form.
#' \code{fd_pack_dist} takes a repametrized form and returns a natural form in
#' the format of \code{fd_gamma_dist}.
#'
#' \code{fd_pdist} evaluates the cumulative distribution function of the gamma
#' distribution \code{g} (natural form) at times \code{t}.  \code{fd_ddist}
#' evaluates the density function \eqn{g(t)} when \code{dt == NULL}.  When
#' \code{dt} is a numeric vector, it evaluates the small increments \eqn{g(t) \
#' dt}{g(t) dt} (so \code{dt} must have the size of \code{t}).
#'
#' @section (Shifted) gamma distributions: Times to division/death can be
#'   parametrized by any continuous distribution.  However, as observed by
#'   Hawkins et al. (2007) and Miao et al. (2012), a gamma
#'   distribution is often flexible enough.  A gamma distribution, contrary to a
#'   normal distribution, has a nonnegative support and is skewed, two
#'   properties that are expected of distributions of times to division/death
#'   (Hyrien et al. (2008)).  To ascertain whether finer parametrization
#'   can lead to better results, we allow the user to set independently not only
#'   the mean and variance of the distribution, but also its skewness (third
#'   moment) through a distribution \emph{shift}.
#'
#'   More specifically, is possible to extend the gamma distribution by shifting
#'   it by a location parameter \eqn{\lambda \ge 0} (type III distribution in
#'   the Pearson Distribution System): \eqn{\Gamma_{shifted}(t; \delta, \bar
#'   \tau, \lambda) = \Gamma(t - \lambda; \delta, \bar \tau) }{\Gamma/shifted(t;
#'   \delta, mean \tau, \lambda) = \Gamma(t - \lambda; \delta, mean \tau)} (with
#'   the convention that \eqn{\Gamma(t) = 0} for \eqn{t < 0}). (The gamma
#'   distribution is scale-invariant, so it would not make sense to add a
#'   fourth, scale parameter as is usually done to get location-scale
#'   families.)   To get relevant bounds and help convergence, we can
#'   reparametrize the distribution as follows.  Noting \eqn{\gamma > 0} the
#'   skewness, \eqn{\Delta > 0} the coefficient of variation and \eqn{\mu > 0}
#'   the mean, we have \deqn{\alpha = \frac{4}{\gamma^2} \quad\quad \beta =
#'   \frac{2}{\mu \gamma \Delta} \quad\quad \lambda = \mu (1 - 2 \Delta /
#'   \gamma)}{ \alpha = \gamma^2 / 4; \beta = \mu \gamma \Delta / 2.0; \lambda =
#'   \mu (1 - 2 \Delta / \gamma)}
#'
#'   We clearly see with this last equality that we must have \eqn{\Delta /
#'   \gamma \le 1/2}: instead of parametrizing with \eqn{\gamma}, we parametrize
#'   with the ratio \eqn{0 < u = \Delta / \gamma \le 1/2} of the coefficient of
#'   variation over the skewness.
#'
#' @section Special values: All the proliferation models in this package allow
#'   the special, "degenerate" value \code{delta = 0}.  In this case, the
#'   distribution of times to division/death is a Dirac function with the mean
#'   still given by \code{mm} and \code{ss} is irrelevant. Moreover, because the
#'   distribution does not have a density anymore, \code{fd_ddist} must be
#'   called with \code{dt} non NULL (it gives a numeric vector filled with
#'   \code{0} except for one element). Also, \code{fd_pdist} reduces to a step
#'   function.
#'
#' @param a The shape \eqn{\alpha} of the gamma distribution.
#' @param b The rate \eqn{\beta} of the gamma distribution.
#' @param loc The location parameter \eqn{\lambda}.
#' @param mul A multiplier factor.  Used internally for representing mixtures.
#' @param mm The mean of the gamma distribution, in hours.
#' @param delta The coefficient of variation (i.e., standard deviation over
#'   mean) of the gamma distribution.  It is dimensionless and should lie
#'   between 0 and 1.
#' @param ss The reparametrization of the location parameter.  It is
#'   dimensionless and should lie between 0, excluded (high skewness as compared
#'   to a non-shifted gamma distribution), and 0.5, included (same skewness as
#'   compared to a shifted gamma distribution).
#' @param t A numeric vector of times in hours at which to evaluate the gamma
#'   distribution.
#' @param dt A list of small time increments for evaluating the gamma
#'   distribution (see details).
#' @param dst A list of the parameters of a gamma distribution in their
#'   reparametrized form: \code{dst$mm}, \code{dst$delta}, \code{dst$ss}.
#' @param g A gamma distribution in the natural form (output of
#'   \code{fd_gamma_dist}).
#'
#' @rdname gamma-distributions
#' @name gamma-distributions
#' @examples
#' mdl <- fd_model()
#' dst <- fd_pack_dist(relist(start(mdl), mdl)$pro$One$f)
#' curve(fd_pdist(x, dst), from=0, to=20)
#' curve(fd_ddist(x, dst), from=0, to=20, add=TRUE, col="red")
#'
#' @references Hawkins ED, Turner ML, Dowling MR, van Gend C, Hodgkin PD (2007).
#' A model of immune regulation as a consequence of randomized lymphocyte
#' division and death times. \emph{Proc Natl Acad Sci USA} \strong{104} (12):
#' 5032-5037.
#'
#' Hyrien O, Zand MS (2008). A Mixture Model with Dependent Observations for the
#' Analysis of CSFE-Labeling Experiments. \emph{Journal of the American
#' Statistical Association} \strong{103} (481): 222-239.
#'
#' Miao H, Jin X, Perelson AS, Wu H (2012). Evaluation of multitype mathematical
#' models for CFSE-labeling experiment data. \emph{Bull Math Biol} \strong{74}
#' (2): 300-326.
NULL

#' @export
#' @rdname gamma-distributions
fd_gamma_dist <- function (a, b, loc=0, mul=1) {
  matrix(c(a, b, loc, mul), nrow=1,
      dimnames = list(NULL, c("a", "b", "loc", "mul")))
}

#' @export
#' @rdname gamma-distributions
fd_reparam_gamma_dist <- function (mm, delta = 1, ss = 0.5) {
  list(mm = mm, delta = delta, ss = ss)
}

#' @export
#' @rdname gamma-distributions
fd_pack_dist <- function (dst, mul = 1) {
  if (dst$mm <= 0)
    stop("'mm' (=", dst$mm, ") is negative or 0")
  if (dst$delta == 0) {
    if (dst$ss != 0.5)
      stop("'ss' (=", dst$ss, ") must be fixed to 0.5 when 'delta == 0'")
    fd_gamma_dist(
      a = Inf,
      b = Inf,
      loc = dst$mm,
      mul = mul
    )
  } else {
    if (dst$ss <= 0 || dst$ss > 0.5)
      stop("'ss' (=", dst$ss, ") must be between 0 (excluded) and 0.5")
    if (dst$delta < 0)
      stop("'delta' (=", dst$delta,
         ") must be greater than or equal to 0")
    fd_gamma_dist(
      a = 4 * dst$ss ^ 2 / dst$delta ^ 2,
      b = 2 * dst$ss / (dst$mm * dst$delta ^ 2),
      loc = dst$mm * (1 - 2 * dst$ss),
      mul = mul
    )
  }
}

#' @export
#' @rdname gamma-distributions
fd_pdist <- function (t, g) {
  # Mind that t must be sorted!!!
  if (!is.finite(g[1, "b"])) {
    ans <- rep(1, length(t))
    ans[1:findInterval(g[1, "loc"], t)] <- 0
  } else {
    pgamma(t - g[1, "loc"], shape = g[1, "a"], rate = g[1, "b"])
  }
}

#' @export
#' @rdname gamma-distributions
fd_ddist <- function (t, g, dt = NULL) {
  if (!is.finite(g[1, "b"])) {
    if (!is.vector(t))
      stop("fd_ddist with delta=0 can only be called on a vector!!!")
    ans <- rep(0, length(t))
    if (is.finite(g[1, "loc"]))
      ans[findInterval(g[1, "loc"], t)] <- 1
    ans
  } else {
    ans <- dgamma(t - g[1, "loc"], shape = g[1, "a"], rate = g[1, "b"])
    if (!is.null(dt)) ans * dt
    else ans
  }
}
