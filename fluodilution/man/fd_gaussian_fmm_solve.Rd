% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/proportions.R
\name{fd_gaussian_fmm_solve}
\alias{fd_gaussian_fmm_solve}
\title{Get the proportions of cells in any given number of generations.}
\usage{
fd_gaussian_fmm_solve(data, fmm = "gaussian", mgen = 10L, rm = TRUE)
}
\arguments{
\item{data}{An \code{\link{fd_data}} object that represents the population to
fit.}

\item{fmm}{The finite mixture model to use for the fitting (see
\code{\link{finite-mixture}}).}

\item{mgen}{The maximum number of generations to fit.  Notice that by
convention the first generation has number 0, so there will be
\code{mgen+1} generations fitted.}

\item{rm}{Whether to remove the \code{"hists"} and \code{"hists_lost"} types
in the return value.}
}
\value{
The \code{\link{fd_data}} object passed as argument augmented by two
  new types: \code{"props"} and \code{"props_lost"} (if any
  \code{"hists_lost"} data).
}
\description{
The result of this function can be fed to \code{\link{fd_nls}} (it returns a
\code{data.frame} that can be \code{rbind}-ed to form an
\code{\link{fd_data}}).  Notice that this intermediary step is usually NOT
performed by the algorithm that gets the biological constants in this manner:
the finite mixture model layer (see \code{\link{finite-mixture}}) and the
proliferation model are better combined in a hierarchical fashion and fitted
together, as advocated by Hyrien (2008).
}
\details{
This function gives the proportion of cells in any given number of
generations from fluorescence histograms.  This function is similar to
implementations in the bioconductor package \pkg{flowFit} but is integrated
with \code{fd_data} dataset thinking and uses a constrained least square
approach (quadratic programming) instead of a Levenberg-Marquardt least
square (this makes more sense from a mathematical point of view).
}
\section{Required packages}{
 Please install \pkg{limSolve}.
}

\examples{
# Load an artificial data set with clear peaks
data(FdClearPeaks)
plot(FdClearPeaks)

# Solve for proportions
fd_gaussian_fmm_solve(FdClearPeaks,
                      model(attr(FdClearPeaks, "params"))$fmm)

}
\references{
Hyrien, Ollivier and Zand, Martin S (2008).  A Mixture Model With
Dependent Observations for the Analysis of CSFE-Labeling Experiments.
\emph{Journal of the American Statistical Association} \strong{103} (481):
222-239.
}
\seealso{
\code{\link{fd_fmm_gaussian}}

Other optimization-related entries: \code{\link{fd_comb}},
  \code{\link{fd_minuslogl}}, \code{\link{fd_nls}}
}
