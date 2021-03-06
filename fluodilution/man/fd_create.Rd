% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create.R
\name{fd_create}
\alias{fd_create}
\alias{fd_moments}
\title{Create an \code{\link{fd_data}} object from separate cell counts, flow
cytometry raw data, ...}
\usage{
fd_create(value, fmm = NULL, categories = NULL, breaks = "Sturges",
  momentControl = NULL)

fd_moments(x, channel, clean = TRUE, varnames = NULL)
}
\arguments{
\item{value}{A named list.  Each element represents the data for a different
weight class (the names of the list give the weight classes).}

\item{fmm}{An optional list of finite mixture model arguments to annotate the
result. Alternatively, could be \code{"use0"} for timepoints at 0h to be
used as inoculums, or \code{"usemin"} for earliest timepoints to be used
as inoculums (when \code{value} contains a \code{"flow"} element). In
this case, the inoculum timepoints are not present as \code{fd_data}
rows.  Instead, the \code{fmm} attribute is updated with inoculum
information.}

\item{categories}{A character vector of category names to give their precise
order (see \code{\link{proliferation}} for why it is important).}

\item{breaks}{For flow cytometry raw data, the breaks to use when forming
histograms. If it is a list, each element of the list is used separately
for the corresponding flow cytometry timepoint and has the same format as
what \code{\link[graphics]{hist}} allows. Otherwise, the same breaking
rules are used across the board.  The rules are the same as with
\code{\link[graphics]{hist}}, but they apply on an an \code{asinh}
transformation of the linear fluorescence.}

\item{momentControl}{A list of arguments that override the default arguments
of \code{fd_moments} when it is called (that is, when \code{fmm = "use0"}).}

\item{x}{A \code{\link[flowCore]{flowSet}} or
\code{\link[flowCore]{flowFrame}}}

\item{channel}{The flow cytometry channel featuring the fluorescence
dilution.}

\item{clean}{Whether to clean data by using model-based clustering
(`flowClust` package).}

\item{varnames}{If cleaning is enabled, channels on which to do the
model-based clustering on.}
}
\value{
An \code{\link{fd_data}} object.
}
\description{
\code{fd_create} combines together the various datasets that
could be used for fluorescence dilution and converts flow cytometry raw data
into histograms. \code{fd_moments} can be used to get a "robust" estimate of
the geometric mean and geometric standard deviation of flow cytometry
samples.
}
\details{
If provided, the element \code{value$flow} is either a
\code{\link[flowWorkspace]{GatingSet}} (bioconductor package
\pkg{flowWorkspace}) or a \code{\link[flowCore]{flowSet}}: that is, a list of
flow cytometry raw timepoints.  Additional attributes instruct
\code{\link{fd_create}} how to interpret this raw data: \describe{
\item{\code{channel}}{Which flow cytometry channel features the fluorescence
dilution.} \item{\code{meta}}{A \code{data.frame}, each row giving for the
corresponding raw timepoint some meta information that will be
\code{cbind}-ed to the histograms.  Meta information include \code{Time},
\code{Timepoint}, \code{Category}, \code{Inoculum}, \code{Type}.  If any but
\code{Time} is missing, the corresponding columns are filled with default
values.} }

The other elements of \code{value} are simply \code{rbind}-ed together and
the column \code{Weight} filled in with the name of those elements.

\code{fd_moments} returns either a vector of two elements (\code{x} is a
\code{flowFrame}) or a matrix with two rows (\code{x} is a \code{flowSet})
giving respectively the log-mean and the log-standard deviation of the
samples provided.  Additional arguments allow the removal of outliers.
}
\section{Required packages}{
 Please install the bioconductor packages
  \pkg{flowWorkspace} for \code{fd_create} and \pkg{flowClust} for
  \code{fd_moments}.  Moreover when calling \code{fd_moments}, because of a
  bug in \pkg{flowClust}, \pkg{parallel} must be explicitly attached to the
  search path using \code{library(parallel)}. Therefore, \pkg{parallel} must
  be attached as well when calling \code{fd_create} with \code{fmm = "use0"}.
}

\examples{
# Unzip the archive with the FCS and workspace files
unzip(system.file("extdata", "Archive.zip", package="fluodilution"),
      exdir = tempdir())

library(flowWorkspace)

## Load a flowJo workspace

ws <- openWorkspace(paste0(tempdir(), "/20150325.wsp"))
print(ws)
flow <- suppressWarnings(parseWorkspace(
  ws, "All Samples", path = tempdir(), isNcdf = FALSE,
  cleanup = FALSE, keep.indices = TRUE,
  requiregates = FALSE
))

# Meta information
channel <- "Comp-488-530_30-A"
times <- sapply(
  seq_along(flow),
  function (i) {
    # Just an example of what is possible
    as.numeric(getKeywords(ws, sub("_[0-9]*$", "",
               flow[[i]]@name))$Time)
  }
)
print(getNodes(flow))

# Go for it
meta <- data.frame(Time = times)
ans <- fd_create(
  value = list(flow = structure(flow, meta = meta, channel = channel))
)

# Display result
plot(ans, main="Workspace (1)")

# One can use "use0"
library(parallel)
ans2 <- fd_create(
  value = list(flow = structure(flow, meta = meta, channel = channel)),
  fmm = "use0"
)
print(attr(ans, "fmm"))
print(attr(ans2, "fmm"))

## Alternatively, a flowSet can be provided

channel <- "FL2-H"
gate <- rectangleGate(`FL1-H` = c(18, Inf), filterId="Bacteria")
src <- read.flowSet(paste0(tempdir(), "/",
                           c("WT-t12.023", "WT-t2.016", "WT-t6.014")))
flow <- Subset(src, filter(src, gate))

meta <- data.frame(Category = "One", Time = c(12, 2, 6), Type = "hists",
                   Inoculum = "inoc_1",
                   Timepoint = c("tp_1", "tp_2", "tp_3"))
value <- list(flow = structure(flow, meta = meta, channel = channel))
ans3 <- fd_create(value = value)
plot(ans3, main="Workspace (2)")

## fd_moments can be used to get, e.g., m0 and sd0 directly
library(parallel)
print(fd_moments(flow, channel))

}
