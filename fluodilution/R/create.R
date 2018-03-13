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

#' @title Create an \code{\link{fd_data}} object from separate cell counts, flow
#' cytometry raw data, ...
#'
#' @description \code{fd_create} combines together the various datasets that
#' could be used for fluorescence dilution and converts flow cytometry raw data
#' into histograms. \code{fd_moments} can be used to get a "robust" estimate of
#' the geometric mean and geometric standard deviation of flow cytometry
#' samples.
#'
#' @param value A named list.  Each element represents the data for a different
#'   weight class (the names of the list give the weight classes).
#' @param fmm An optional list of finite mixture model arguments to annotate the
#'   result. Alternatively, could be \code{"use0"} for timepoints at 0h to be
#'   used as inoculums, or \code{"usemin"} for earliest timepoints to be used
#'   as inoculums (when \code{value} contains a \code{"flow"} element). In
#'   this case, the inoculum timepoints are not present as \code{fd_data}
#'   rows.  Instead, the \code{fmm} attribute is updated with inoculum
#'   information.
#' @param categories A character vector of category names to give their precise
#'   order (see \code{\link{proliferation}} for why it is important).
#' @param breaks For flow cytometry raw data, the breaks to use when forming
#'   histograms. If it is a list, each element of the list is used separately
#'   for the corresponding flow cytometry timepoint and has the same format as
#'   what \code{\link[graphics]{hist}} allows. Otherwise, the same breaking
#'   rules are used across the board.  The rules are the same as with
#'   \code{\link[graphics]{hist}}, but they apply on an an \code{asinh}
#'   transformation of the linear fluorescence.
#' @param momentControl A list of arguments that override the default arguments
#'   of \code{fd_moments} when it is called (that is, when \code{fmm = "use0"}).
#' @param x A \code{\link[flowCore]{flowSet}} or
#'   \code{\link[flowCore]{flowFrame}}
#' @param channel The flow cytometry channel featuring the fluorescence
#'   dilution.
#' @param clean Whether to clean data by using model-based clustering
#' (`flowClust` package).
#' @param varnames If cleaning is enabled, channels on which to do the
#' model-based clustering on.
#'
#' @return An \code{\link{fd_data}} object.
#' @export
#'
#' @details If provided, the element \code{value$flow} is either a
#' \code{\link[flowWorkspace]{GatingSet}} (bioconductor package
#' \pkg{flowWorkspace}) or a \code{\link[flowCore]{flowSet}}: that is, a list of
#' flow cytometry raw timepoints.  Additional attributes instruct
#' \code{\link{fd_create}} how to interpret this raw data: \describe{
#' \item{\code{channel}}{Which flow cytometry channel features the fluorescence
#' dilution.} \item{\code{meta}}{A \code{data.frame}, each row giving for the
#' corresponding raw timepoint some meta information that will be
#' \code{cbind}-ed to the histograms.  Meta information include \code{Time},
#' \code{Timepoint}, \code{Category}, \code{Inoculum}, \code{Type}.  If any but
#' \code{Time} is missing, the corresponding columns are filled with default
#' values.} }
#'
#' The other elements of \code{value} are simply \code{rbind}-ed together and
#' the column \code{Weight} filled in with the name of those elements.
#'
#' \code{fd_moments} returns either a vector of two elements (\code{x} is a
#' \code{flowFrame}) or a matrix with two rows (\code{x} is a \code{flowSet})
#' giving respectively the log-mean and the log-standard deviation of the
#' samples provided.  Additional arguments allow the removal of outliers.
#'
#' @examples
#' # Unzip the archive with the FCS and workspace files
#' unzip(system.file("extdata", "Archive.zip", package="fluodilution"),
#'       exdir = tempdir())
#'
#' library(flowWorkspace)
#'
#' ## Load a flowJo workspace
#'
#' ws <- openWorkspace(paste0(tempdir(), "/20150325.wsp"))
#' print(ws)
#' flow <- suppressWarnings(parseWorkspace(
#'   ws, "All Samples", path = tempdir(), isNcdf = FALSE,
#'   cleanup = FALSE, keep.indices = TRUE,
#'   requiregates = FALSE
#' ))
#'
#' # Meta information
#' channel <- "Comp-488-530_30-A"
#' times <- sapply(
#'   seq_along(flow),
#'   function (i) {
#'     # Just an example of what is possible
#'     as.numeric(getKeywords(ws, sub("_[0-9]*$", "",
#'                flow[[i]]@name))$Time)
#'   }
#' )
#' print(getNodes(flow))
#'
#' # Go for it
#' meta <- data.frame(Time = times)
#' ans <- fd_create(
#'   value = list(flow = structure(flow, meta = meta, channel = channel))
#' )
#'
#' # Display result
#' plot(ans, main="Workspace (1)")
#'
#' # One can use "use0"
#' library(parallel)
#' ans2 <- fd_create(
#'   value = list(flow = structure(flow, meta = meta, channel = channel)),
#'   fmm = "use0"
#' )
#' print(attr(ans, "fmm"))
#' print(attr(ans2, "fmm"))
#'
#' ## Alternatively, a flowSet can be provided
#'
#' channel <- "FL2-H"
#' gate <- rectangleGate(`FL1-H` = c(18, Inf), filterId="Bacteria")
#' src <- read.flowSet(paste0(tempdir(), "/",
#'                            c("WT-t12.023", "WT-t2.016", "WT-t6.014")))
#' flow <- Subset(src, filter(src, gate))
#'
#' meta <- data.frame(Category = "One", Time = c(12, 2, 6), Type = "hists",
#'                    Inoculum = "inoc_1",
#'                    Timepoint = c("tp_1", "tp_2", "tp_3"))
#' value <- list(flow = structure(flow, meta = meta, channel = channel))
#' ans3 <- fd_create(value = value)
#' plot(ans3, main="Workspace (2)")
#'
#' ## fd_moments can be used to get, e.g., m0 and sd0 directly
#' library(parallel)
#' print(fd_moments(flow, channel))
#'
#' @section Required packages: Please install the bioconductor packages
#'   \pkg{flowWorkspace} for \code{fd_create} and \pkg{flowClust} for
#'   \code{fd_moments}.  Moreover when calling \code{fd_moments}, because of a
#'   bug in \pkg{flowClust}, \pkg{parallel} must be explicitly attached to the
#'   search path using \code{library(parallel)}. Therefore, \pkg{parallel} must
#'   be attached as well when calling \code{fd_create} with \code{fmm = "use0"}.
#'
fd_create <- function (value, fmm = NULL, categories = NULL,
                       breaks = "Sturges",
                       momentControl = NULL) {
  if (!is.null(value$flow)) {
    meta <- attr(value$flow, "meta")
    channel <- attr(value$flow, "channel")

    if (is.null(meta)) {
      stop("'meta' attribute is required")
    }
    if (is.null(channel)) {
      stop("'channel' attribute is required")
    }

    if (is.null(meta$Type)) {
      meta$Type <- factor("hists")
    }
    if (is.null(meta$Inoculum)) {
      meta$Inoculum <- factor("inoc_1")
    }
    if (is.null(meta$Category)) {
      meta$Category <- factor("One")
    }
    if (is.null(meta$Timepoint)) {
      meta$Timepoint <-
        with(meta,
          factor(paste0("tp_", Inoculum, ".", Category, ".", Time)))
    }

    isWorkspace <- inherits(value$flow, "GatingSet")
    if (!isWorkspace && !inherits(value$flow, "flowSet"))
      stop("'value$flow' must be either a 'GatingSet' (flowWorkspace) ",
        "or a 'flowSet' (flowCore)")

    if (!requireNamespace(c("flowCore"), quietly=TRUE)) {
      stop("package 'flowCore' is required")
    }
    if (isWorkspace &&
        !requireNamespace(c("flowWorkspace"), quietly=TRUE)) {
      stop("package 'flowWorkspace' is required")
    }

    old <- value$hist
    fmm0 <- list(m0 = c(), sd0 = c())
    value$hist <- plyr::rbind.fill(lapply(1:NROW(meta), function (i) {
      if (isWorkspace) {
        # Read raw data in the proper gating set
        set <- meta$Set[i]
        G <- value$flow
        if (!is.null(set) && set > 0)
          G <- G[[as.numeric(as.character(set))]]

        sname <- meta$Sample[i]
        if (is.null(sname)) {
          if (!is.null(set) && set > 0)
            stop("'Sample' not in meta, but 'Set' is")
          sname <- i
        }
        trInv <- flowWorkspace::getTransformations(G[[sname]],
          inverse=TRUE)
        if (is.null(meta$Gate[i]) || meta$Gate[i] == "") {
          d <- flowWorkspace::getData(G[[sname]])
        } else {
          d <- flowWorkspace::getData(
            G, as.character(meta$Gate[i]))[[sname]]
        }
        names(trInv) <- gsub("^\\s+|\\s+$", "", names(trInv))
        rawData <- flowCore::`%on%`(
          flowCore::transformList(channel, trInv[channel]), d)
      } else {
        # Warning: we do not check for linearization
        rawData <- value$flow[[i]]
      }

      # Use this timepoint as inoculum, if required
      as_inoculum <- FALSE
      if (length(fmm) == 1) {
        if (fmm == "use0") as_inoculum <- meta$Time[i] == 0
        else if (fmm == "usemin") {
          tmin <- min(subset(meta, Inoculum == meta$Inoculum[i])$Time)
          as_inoculum <- meta$Time[i] == tmin
        }
      }
      if (as_inoculum) {
        if (sum(meta$Time == 0 & meta$Inoculum == meta$Inoculum[i]) > 1) {
          print(subset(meta, Time == 0 & Inoculum == meta$Inoculum[i]))
          warning(
            "more than one 't=min/0h' timepoint for the same inoculum ",
            "(see above); the last one is used"
          )
        }
        m <- do.call(fd_moments, append(list(rawData, channel),
          momentControl))
        fmm0$m0[as.character(meta$Inoculum[i])] <<- m["mean"]
        fmm0$sd0[as.character(meta$Inoculum[i])] <<- m["sd"]
        NULL
      }
      else {
        # Read specific channel
        if (NROW(flowCore::exprs(rawData)) == 0 ||
            NROW(d <- flowCore::exprs(rawData)) == 0)
          return (NULL)
        s <- d[, channel]

        # Create histogram
        if (is.list(breaks)) brk <- breaks[[i]]
        else brk <- breaks
        if (is.numeric(brk) & length(brk) > 1) {
          # Breaks explicitly given
          h <- hist(s, breaks = brk, plot = FALSE)
        } else {
          h <- hist(asinh(s), breaks = brk, plot = FALSE)
          brk <- sinh(h$breaks)
        }
        suppressWarnings(data.frame(
          meta[i, , drop=FALSE],
          data.frame(a = brk[-length(brk)],
            b = brk[-1],
            y = h$counts / sum(h$counts))
        ))
      }
    }))
    value$hist$Set <- NULL
    value$hist$Gate <- NULL
    value$hist$Sample <- NULL
    value$hist <- rbind(old, value$hist)
    value$flow <- NULL
  }

  if (length(fmm) == 1 && fmm %in% c("use0", "usemin")) fmm <- fmm0

  fd_data(structure(
    plyr::rbind.fill(lapply(
      names(value),
      function (nm) {
        if (NROW(value[[nm]]) == 0) data.frame()
        else transform(value[[nm]], Weight = nm)
      })),
    fmm = fmm
  ), categories = categories)
}

#' @rdname fd_create
#' @export
fd_moments <- function (x, channel, clean = TRUE, varnames = NULL) {
  process <- function (xx) {
    if (clean) {
      if (!requireNamespace(c("flowClust"), quietly=TRUE)) {
        stop("package 'flowClust' is required")
      }
      if (!("package:parallel" %in% search()))
        stop("package 'parallel' must be added to the search path ",
          "using: 'library(parallel)' ",
          "for 'flowClust' to work")
      res <- suppressMessages(flowClust::flowClust(
        xx,
        varNames=varnames,
        K=1, B=100))
      xx <- xx[flowClust::`%in%`(xx, res), ]
    }

    e <- flowCore::exprs(xx)[, channel]
    e <- log(e[e > 0])
    c(mean = mean(e), sd = sd(e))
  }

  if (inherits(x, "flowSet")) {
    ans <- sapply(seq_along(x), function (i) process(x[[i]]))
    data.frame(mean = ans["mean", ], sd = ans["sd", ])
  } else if (inherits(x, "flowFrame"))
    process(x)
  else stop("class not handled")
}
