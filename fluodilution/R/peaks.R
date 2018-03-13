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

#' @title Locate the distinct peaks in fluorescence dilution histograms.
#'
#' @description Find the peaks, using a smoothing method, of fluorescence
#' histograms.  This could be useful in order to eschew the \emph{a priori}
#' assumptions concerning the evolution over time and generations of the space
#' between those peaks that finite mixture models make (Hyrien et al.
#' (2008) have showed that the "naive", Gaussian FMM is sometimes
#' inappropriate).  However, clear peaks are not discernible in microbiology
#' applications so far, restricting this technique to immunological
#' applications.
#'
#' @param data An \code{\link{fd_data}} dataset.
#' @param plot Whether a diagnostic plot of the peaks should be displayed.
#' @param w Window parameter: how large (in number of histogram bins) should the
#'   window used to calculate the maximum be?
#' @param span The parameter that controls the degree of smoothing in
#'   \code{\link[stats]{loess}}.
#' @param peaks The return value of an \code{fd_findpeaks} call, with an
#'   additional \code{Generation} column mapping the peaks to numbers of
#'   generations.
#' @param rm Whether the rows with \code{Weight=="hist"} should be removed in
#'   the returned \code{\link{fd_data}} object.
#'
#' @return \code{fd_findpeaks} and \code{fd_peaks} both return a
#'   \code{data.frame} (see \emph{details}).
#'
#' @details \code{fd_findpeaks} returns a \code{data.frame} giving the
#' \code{"Timepoint"}, position of the peak \code{"x"}, proportion of cells in
#' the peak \code{"prop"} and position of the peaks from right to left
#' \code{"i"}, starting from 1. Notice that the number of generation is not
#' directly given as sometimes very low proportions could mean that a generation
#' is skipped and not recorded as a peak.  Since this should not pose any
#' problem, apart from labelling, the labelling should be done by the user (see
#' the examples section for a straightforward simple mapping).
#'
#' \code{fd_peaks} takes such a mapping and appends to the dataset the
#' proportions, with weight \code{"prop"} and types \code{"props"}
#' (corresponding to \code{"hists"} histograms) and \code{"props_lost"}
#' (corresponding to \code{"hists_lost"} histograms).
#'
#' @section Required packages: Please install \pkg{zoo}.
#'
#' @name fd_peaks
#' @seealso \code{\link{finite-mixture}}, \code{\link{fd_gaussian_fmm_solve}}.
#' @examples
#' data(FdClearPeaks)
#' dat <- FdClearPeaks
#'
#' # The above dataset has a low noise level, therefore
#' # no transformation is really needed to map peaks to
#' # generation numbers
#' peaks <- transform(fd_findpeaks(dat, span=0.05, w=6),
#'                    Generation = i - 1)
#' fd_peaks(dat, peaks, rm=TRUE)
NULL

#' @rdname fd_peaks
#' @export
fd_findpeaks <- function (data, plot=TRUE, w=5, span=0.1) {
  if (!requireNamespace("zoo", quietly=TRUE)) {
    stop("'zoo' package is required for 'fd_findpeaks'")
  }

  # http://stats.stackexchange.com/questions/36309/
  # how-do-i-find-peaks-in-a-dataset (whuber)
  argmax <- function(x, y, w=1, ...) {
    n <- length(y)
    y.smooth <- loess(y ~ x, ...)$fitted
    y.max <- zoo::rollapply(zoo::zoo(y.smooth), 2 * w + 1,
                max, align="center")
    delta <- y.max - y.smooth[-c(1:w, n + 1 - 1:w)]
    i.max <- which(delta <= 0) + w
    list(x=x[i.max], i=i.max, y.hat=y.smooth)
  }

  # First, get the position of the peaks
  df_peaks <- plyr::ddply(subset(data, Weight == "hist"),
              "Timepoint", function (df) {
    xt <- with(df, (asinh(a) + asinh(b)) / 2.0)
    tr <- with(df, sqrt(1 + ((a + b) / 2.0) ^ 2) / (b - a))
    yt <- with(df, y * tr)
    yt[!is.finite(yt)] <- 0
    am <- argmax(xt, yt, w=w, span=span)

    prop <- yt[am$i] / (asinh(data$b[am$i]) - asinh(data$a[am$i]))
    prop <- prop / sum(prop)

    data.frame(x = sinh(am$x), Proportion = prop, i = length(prop):1L)
  })

  if (plot) {
    print(ggplot(subset(data, Weight == "hist"),
           aes(x = (a + b) / 2.0))+
      stat_unitarea(aes(yu = y), geom ="line")+
      geom_vline(data = df_peaks, aes(xintercept = x))+
      geom_label(data = df_peaks, aes(x = x, y = 0.8, label = i))+
      scale_x_log10()+
      facet_wrap(~ Timepoint)+
      labs(x = "Fluorescence", y = "Unit area")+
      theme_classic())
  }

  df_peaks
}

#' @rdname fd_peaks
#' @export
fd_peaks <- function (data, peaks, rm=TRUE) {
  meta <- subset.data.frame(data, Weight == "hist") %>%
    dplyr::select(-a, -y, -b)
  meta <- meta[!duplicated(meta$Timepoint), ]
  df_props <- peaks %>%
    dplyr::select(Timepoint, a = Generation, y = Proportion) %>%
    merge(meta)
  df_props$b <- 0
  df_props$Weight <- "prop"
  df_props$Type <- ifelse(df_props$Type == "hists", "props", "props_lost")
  if (rm) data <- subset(data, Weight != "hist")
  fd_data(suppressWarnings(rbind.data.frame(
    data, df_props,
    make.row.names = FALSE)),
    categories = levels(data$Category))
}
