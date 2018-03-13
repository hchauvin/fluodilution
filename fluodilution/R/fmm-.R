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

#' Finite mixture models.
#'
#' Finite mixture models relate the fluorescence to the numbers of generations.
#' They are based on the decomposition, using Bayes' law, of the fluorescence
#' histograms into various cohorts whose individual shape is derived from
#' properties of the inoculum and hypotheses concerning the impact of
#' autofluorescence (AF) and binomial partitioning (BP).  By extension, our
#' implementations also contain variables concerning the behaviour of cell
#' counts and act more generally as an interface between the dataset and the
#' proliferation model.
#'
#' The simplest hypothesis to make is that the successive cohorts are no more
#' than translations, in logarithmic units, of the inoculum.  If the inoculum is
#' furthermore considered log-normally distributed, the resulting model is the
#' Gaussian FMM (\code{fd_fmm_gaussian}). \code{fd_fmm_af} also takes
#' autofluorescence into account and \code{fd_fmm_af_bp} takes both
#' autofluorescence and binomial partitioning into account (Chauvin \emph{et
#' al.} 2016). Autofluorescence impacts the shape of the cohorts: as the "real"
#' fluorescence dilutes, autofluorescence becomes more important relative to
#' "real" fluorescence. Binomial partitioning, or how fluorescent molecules are
#' partitioned at division between daughter cells, creates an asymmetry whose
#' magnitude depends on the number of molecules in the mother cell.
#' \code{fd_fmm_af} and \code{fd_fmm_af_bp} are significantly slower than
#' \code{fd_fmm_gaussian}. An alternative to mitigate the effect of AF and BP is
#' to gate away the low fluorescence, where it is the most distorted, using
#' \code{\link{cutoff}}.
#'
#' @section Parallel computing: \code{fd_fmm_af_bp} makes use of
#'   \code{\link[parallel]{mclapply}}. The number of clusters to use can be set
#'   using, e.g., \code{options(mc.cores = 4L)}.
#'
#' @section Required packages: For \code{fd_fmm_af_bp}, please install
#'   \pkg{parallel}.
#'
#' @param ... Not used by \code{summary}.  For the other functions: \describe{
#'   \item{\code{m0, sd0}}{Log-normal mean and standard deviation of
#'   fluorescence in the inoculum/of progenitors.  Either a scalar or a named
#'   numeric vector, with the names being the levels of the \code{"Inoculum"}
#'   column of the \code{\link{fd_data}} dataset, in the exact same order.}
#'   \item{\code{cctrans}}{The transformation that has been applied to the cell
#'   counts of the dataset, if any. \code{"log10"} by default.  See package
#'   \pkg{scales} for other possible values.} \item{\code{htrans}}{The
#'   transformation that has been applied to the proportions found in the
#'   fluorescence histograms, if any.  \code{"identity"} by default. See package
#'   \pkg{scales} for other possible values.} }
#'
#' @return An \code{fd_fmm} object, suitable for use with
#'   \code{\link{fd_model}}, \code{\link{fd_unif}}, ...
#'
#' @section Parametrization: \describe{ \item{\code{c0}}{Make the link between
#'   the relative cell counts \code{Nrel} given by the model and the
#'   experimental cell counts \code{Nexp}: \code{Nexp = Nrel * 10^c0}.}
#'   \item{\code{sdaf}}{Standard deviation of fluorescence, in linear units, for
#'   \code{fd_fmm_af} and \code{fd_fmm_af_bp}.} \item{\code{ftor}}{Log-specific
#'   fluorescence: that is, the number of molecules for a given fluorescence
#'   level \code{x} (in linear units) is \code{x / exp(ftor)}.} }
#'
#' @section Reparametrization: The natural parametrization (see
#'   \code{\link{fd_model}}) does not differ from the reparametrized form.
#'
#' @section Word of caution: In all those models, the mean autofluorescence is
#'   assumed to be 0.  This supposes that the flow cytometer was correctly
#'   calibrated.  If this was not the case, autofluorescence needs to be
#'   measured independently and the original fluorescence histograms shifted.
#'
#' @seealso \code{\link{fd_gaussian_fmm_solve}}.
#' @name finite-mixture
#' @examples
#' fd_model(fmm="gaussian")
#' fd_model(fmm=fd_fmm_gaussian(m0 = 10, sd0 = 0.3,
#'                              cctrans="log1p", htrans="identity"))
#' fd_fmm_af(m0 = 10, sd0 = 0.3)
#' fd_fmm_af_bp(m0 = 10, sd0 = 0.3)
#'
#' @references McLachlan JG, Peel D (2000). \emph{Finite Mixture Models.} John
#' Wiley & Sons.  ISBN 978-0-471-00626-8.
NULL

#' @export
#' @keywords internal
fd_fmm <- function (name, new, diff_cdf, constrain,
          trans = function (self, params) params,
          trans_inverse = function (self, params) params,
          ...) {
  structure(function (...) {
    self <- new(...)
    if (is.null(self$cctrans)) self$cctrans <- "log10"
    if (is.null(self$htrans)) self$htrans <- "identity"

    if (is.character(self$cctrans))
      self$cctrans <- match.fun(paste0(self$cctrans, "_trans"))()
    if (is.character(self$htrans))
      self$htrans <- match.fun(paste0(self$htrans, "_trans"))()
    if (any(self$sd0 > self$m0))
      stop("'sd0 > m0' for some inoculums")

    self$name <- name

    structure(as.environment(self),
          class=c("fd_fmm", "fd_class"),
          expanded_name = paste0("fd_fmm_", name))
  },
  methods = list(
    diff_cdf = diff_cdf,
    constrain = constrain,
    trans = trans,
    trans_inverse = trans_inverse,
    ...
  ))
}

#' @export
#' @keywords internal
`$.fd_class` <- function (self, member) {
  obj <- attr(get(attr(self, "expanded_name"),
          globalenv()), "methods")[[member]]
  if (is.null(obj)) self[[member]]
  else function (...) obj(self, ...)
}

#' @export
#' @keywords internal
print.fd_fmm <- function (x, ...) {
  print(summary(x))
}

#' @export
#' @keywords internal
summary.fd_fmm <- function (object, ...) {
  structure(list(name = object$name,
           ninoc = length(object$sd0),
           mm0 = mean(object$m0),
           msd0 = mean(object$sd0)),
        class = "summary.fd_fmm")
}

#' @export
#' @keywords internal
print.summary.fd_fmm <- function (x, ...) {
  cat(sep="", "'fd_fmm' object \"",
    x$name, "\" with ", x$ninoc, " inoculums\n")
  cat(sep="", "Mean (log): m0: ", signif(x$mm0, 2),
    "; sd0: ", signif(x$msd0, 2), "\n")
}