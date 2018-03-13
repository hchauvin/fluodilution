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

#' Get the proportions of cells in any given number of generations.
#'
#' The result of this function can be fed to \code{\link{fd_nls}} (it returns a
#' \code{data.frame} that can be \code{rbind}-ed to form an
#' \code{\link{fd_data}}).  Notice that this intermediary step is usually NOT
#' performed by the algorithm that gets the biological constants in this manner:
#' the finite mixture model layer (see \code{\link{finite-mixture}}) and the
#' proliferation model are better combined in a hierarchical fashion and fitted
#' together, as advocated by Hyrien (2008).
#'
#'
#' @details This function gives the proportion of cells in any given number of
#' generations from fluorescence histograms.  This function is similar to
#' implementations in the bioconductor package \pkg{flowFit} but is integrated
#' with \code{fd_data} dataset thinking and uses a constrained least square
#' approach (quadratic programming) instead of a Levenberg-Marquardt least
#' square (this makes more sense from a mathematical point of view).
#'
#' @param data An \code{\link{fd_data}} object that represents the population to
#'   fit.
#' @param fmm The finite mixture model to use for the fitting (see
#'   \code{\link{finite-mixture}}).
#' @param mgen The maximum number of generations to fit.  Notice that by
#'   convention the first generation has number 0, so there will be
#'   \code{mgen+1} generations fitted.
#' @param rm Whether to remove the \code{"hists"} and \code{"hists_lost"} types
#'   in the return value.
#'
#' @return The \code{\link{fd_data}} object passed as argument augmented by two
#'   new types: \code{"props"} and \code{"props_lost"} (if any
#'   \code{"hists_lost"} data).
#' @seealso \code{\link{fd_fmm_gaussian}}
#'
#' @references Hyrien, Ollivier and Zand, Martin S (2008).  A Mixture Model With
#' Dependent Observations for the Analysis of CSFE-Labeling Experiments.
#' \emph{Journal of the American Statistical Association} \strong{103} (481):
#' 222-239.
#' @family optimization-related entries
#' @examples
#' # Load an artificial data set with clear peaks
#' data(FdClearPeaks)
#' plot(FdClearPeaks)
#'
#' # Solve for proportions
#' fd_gaussian_fmm_solve(FdClearPeaks,
#'                       model(attr(FdClearPeaks, "params"))$fmm)
#'
#' @export
#' @section Required packages: Please install \pkg{limSolve}.
fd_gaussian_fmm_solve <- function (data, fmm = "gaussian", mgen = 10L,
                                   rm = TRUE) {
  fmm <- get_fmm(fmm, data)
  if (!requireNamespace("limSolve", quietly=TRUE))
    stop("'limSolve' is required for 'fd_gaussian_fmm_solve'")
  ans <- plyr::ddply(
    subset(data, Weight == "hist"),
    "Timepoint",
    function (df) {
      hst <- df$y

      gen <- 0:mgen
      phi <- fmm$diff_cdf(list(i=df$Inoculum[1]),
                 df$a, df$b, gen)
      sel <- colSums(phi) > .Machine$double.eps
      gen <- gen[sel]
      phi <- phi[, sel, drop=FALSE]

      slv <- limSolve::lsei(
        A = phi, B = hst,
        E = rbind(rep(1, NCOL(phi))), F = matrix(1, ncol = 1),
        G = diag(rep(1, NCOL(phi))), H = cbind(rep(0, NCOL(phi))),
        type = 2
      )

      if (!(abs(sum(slv$X) - 1) < 1e-5 & all(slv$X >= 0))) {
        stop("FD_gaussian_fmm_solve: internal error")
      }

      meta <- df[1, ] %>% dplyr::select(-y, -a, -b)
      suppressWarnings(transform(meta,
        a = gen,
        b = 0,
        y = slv$X
      ))
    }
  )
  if (NROW(ans) > 0) {
    ans$Weight <- "prop"
    ans$Type <- ifelse(ans$Type == "hists", "props", "props_lost")
    ans$Timepoint <- paste0(ans$Timepoint, "_prop")
    ans <- fd_data(ans)
  }
  if (rm) {
    data <- subset(data, Weight != "hist")
  }
  rbind(data, ans)
}
