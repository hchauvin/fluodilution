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

#' Display a histogram with unit area in \pkg{ggplot2}.
#'
#' This function requires the \pkg{numDeriv} package (this package is not
#' automatically installed when the \code{fluodilution} package is installed).
#'
#' In \pkg{ggplot2}, densities displayed with
#' \code{\link[ggplot2]{stat_density}} automatically have unit area when the
#' \emph{x} axis is transformed.  If one wants to use \code{geom_line} to
#' display histograms from an \code{\link{fd_data}} object this transformation
#' is not made.  For instance, as an \code{\link{fd_data}} object stores its
#' breakpoints in linear units of fluorescence, when one wants to apply a
#' log-transformation on the \emph{x} axis then the area under the curve does
#' not sum to one anymore.  \code{stat_unitarea} takes those transformations
#' automatically into account.  For the arguments, see
#' \code{\link[ggplot2]{stat_density}}.  Additionally, \code{transform_hist} can
#' be used to apply a unit area transformation and return a numeric vector of
#' transformed proportions.
#'
#' @param y A numeric vector of histogram proportions to be transformed
#'   conserving unit area.
#' @param a The left breakpoints of the cells of the histogram.
#' @param b The right breakpoints.
#' @param trans Either a character string or a transformation object directly,
#'   see package \pkg{scales} for possible transformers.
#' @inheritParams ggplot2::stat_density
#'
#' @section Aesthetics:
#'
#'   \code{stat_unitarea} understands the following required aesthetics:
#'   \code{x} and \code{yu} (the \emph{y} component of the histograms for which
#'   we want to ensure a unit area).
#'
#' @section Computed variables:
#'
#'   \describe{
#'
#'   \item{\code{unitarea}}{transformation of \code{yu} to ensure unit area (the
#'   default \code{y})}. \item{\code{scaled}}{transformation of \code{yu} to
#'   ensure unit area, scaled to a maximum of 1}.
#'
#'   }
#'
#' @return An object that can be chained to a \code{\link{ggplot}} call.
#'
#' @export
#'
#' @rdname stat_unitarea
#' @examples
#' data(FdClearPeaks)
#' ggplot(subset(FdClearPeaks, Type == "hists"),
#'        aes(x = (a + b) / 2.0, yu = y))+
#'   stat_unitarea()+
#'   facet_wrap(~ Timepoint, scales="free")+
#'   labs(y = "Unit area")+
#'   scale_x_log10()
#'
#' ggplot(subset(FdClearPeaks, Type == "hists"),
#'        aes(x = (a + b) / 2.0, yu = y, y=..scaled..))+
#'   stat_unitarea()+
#'   labs(y = "Unit area, scaled to mode")+
#'   facet_wrap(~ Timepoint, scales="free")+
#'   scale_x_log10()
#'
#' ggplot(transform(subset(FdClearPeaks, Type == "hists"),
#'                  y = transform_hist(y, a, b, "log10")),
#'        aes(x = (a + b) / 2.0, y = y))+
#'   geom_line()+
#'   labs(y = "Unit area")+
#'   facet_wrap(~ Timepoint, scales="free")+
#'   scale_x_log10()
stat_unitarea <- function (mapping = NULL, data = NULL,
              geom = "line", position = "identity",
              ...,
             na.rm = FALSE,
             show.legend = NA,
             inherit.aes = TRUE) {
  if (!requireNamespace("numDeriv", quietly=TRUE)) {
    stop("package 'numDeriv' has to be installed for stat_unitarea to work")
  }

   ggplot2::layer(
  data = data,
  mapping = mapping,
  stat = StatUnitarea,
  geom = geom,
  position = position,
  show.legend = show.legend,
  inherit.aes = inherit.aes,
  params = list(
    na.rm = na.rm,
    ...
  )
  )
}

StatUnitarea <- ggproto("StatUnitarea", Stat,
  required_aes = c("x", "yu"),
  # nolint start
  default_aes = aes(y = ..unitarea.., fill = NA),
  # nolint end

  compute_group = function(data, scales, bw = "nrd0", adjust = 1,
               kernel = "gaussian",
               trim = FALSE, na.rm = FALSE) {
    grad_trans <- function (...) numDeriv::grad(scales$x$transform, ...)
    ret <- data.frame(
      x = data$x,
      xo = scales$x$trans$inverse(data$x))
    ret$dx <- c(diff(ret$xo)[1], diff(ret$xo))
    num <- data$yu
    den <- grad_trans(ret$xo) * ret$dx
    ret$unitarea <- num / den
    ret$unitarea <- ifelse(ret$unitarea < 0, NA, ret$unitarea)
    ret$scaled <- ret$unitarea / max(ret$unitarea)
    return (ret)
  }
)

#' @rdname stat_unitarea
#' @export
transform_hist <- function (y, a, b, trans) {
  if (!requireNamespace("numDeriv", quietly=TRUE)) {
    stop("package 'numDeriv' has to be installed for transform_hist to work")
  }
  if (is.character(trans)) {
    trans <- paste0(trans, "_trans")
    trans <- match.fun(trans)()
  }
  grad_trans <- function (...) numDeriv::grad(trans$transform, ...)
  y / (grad_trans((a + b) / 2.0) * (b - a))
}