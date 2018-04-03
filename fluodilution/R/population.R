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

#' Manage a fluorescence dilution dataset.
#'
#' A fluorescence dilution dataset is a \code{data.frame} of fluorescence
#' histograms and Colony Forming Units (CFUs)/cell counts for various
#' timepoints.  It can contain additional attributes that can be used to
#' simplify the construction of the finite mixture model (FMM) or proliferation
#' model when calling \code{\link{fd_model}} or \code{\link{fd_unif}}.
#' \code{fd_data} ensures the conformity of such a dataset by performing various
#' checks.  Subsetting and \code{rbind}-ing are also available for an
#' \code{fd_data}, as well as various transformation and display mechanisms.
#'
#' \code{fd_data} is used to construct and check the validity of a fluorescence
#' dilution dataset.  \code{fd_transform} applies transformations to histograms
#' and cell counts and change the attributes \code{"fmm"} accordingly.
#' Transformations are not chained: the dataset is back transformed using the
#' current (or default!) transformation, then the new transformation is applied.
#' \code{cutoff} returns a subset of the dataset where the points below the
#' cutoff have been discarded.  Notice that \code{subset} cannot be used instead
#' of \code{cutoff} as the implementation for \code{fd_data} checks that the
#' histograms sum to one.  \code{fd_expand} uses \code{predict} to fill in more
#' predicted timepoints between the experimental timepoints (interpolation),
#' usually for plotting purposes.  \code{plot} can be used for a quick look at
#' the dataset in order to get an overall "feel" and spot potential mistakes.
#'
#' @param data The original \code{data.frame} (or perhaps a
#'   \code{\link[nlme]{groupedData}}).
#' @param x An \code{fd_data} object.
#' @param categories The specific order for the category levels.  If \code{NULL}
#'   and the \code{Category} field is not yet a factor, the levels will be in
#'   lexicographical order.
#' @param inoculums Same with the \code{Inoculum} field.
#' @param timepoints Same with the \code{Timepoint} field.
#' @param na.action An \code{\link[stats]{na.action}} object to apply to
#'   \code{data}.
#' @param hist A transformation object from the (graphical) package \pkg{scales}
#'   (or a character string converted to a transformation object) and to be
#'   applied to the \code{y} field when \code{Weight == "hist"}.
#' @param N A transformation object from the (graphical) package \pkg{scales}
#'   (or a character string converted to a transformation object) and to be
#'   applied to the \code{y} field when \code{Weight == "N"}.
#' @param cutoff If specified, overrides the behaviour of the attribute
#'   \code{cutoff} of \code{data}.  If a named vector of multiple entries, it
#'   contains the cutoff points in linear units of fluorescence (the names are
#'   the timepoint identifiers).  If it is a scalar, gives the maximum number of
#'   generations to consider: the cut is then made at \verb{
#'   } \code{exp(m0) /
#'   2^cutoff}, with \code{m0} the log-mean of the sample.
#' @param threshold Below this value, proportions in fluorescence histograms are
#'   discarded. Can help with egregious heteroskedasticity problems.
#' @param length.N The new number of cell counts per group.
#' @param seq Alternatively, the times at which to calculate the cell counts can
#'   be explicitly given.
#' @param separate If \code{TRUE} (default), the new timepoints to consider are
#'   found independently for each group.  This is useful to avoid extrapolation
#'   when data is unbalanced.
#' @param by The additional columns with which to do the grouping (Grouping is
#'   already made by \code{Category}). If \code{NULL}, \code{separate = TRUE}
#'   and \code{data} inherits \code{groupedData}, this is taken from the
#'   grouping information of \code{data}. Otherwise, no additional grouping is
#'   used.
#' @param type A character vector containing the types of plots to show. If
#'   \code{"overview"}, a selection of plots are shown on the same page.
#'   Otherwise, the plots are shown sequentially.  Plots include \code{"hist"}
#'   (fluorescence histograms), \code{"N"} (cell count), \code{"range"} (an
#'   overview of cutoff points and initial log-mean fluorescence and
#'   log-standard deviation), \code{"balancing"} (data balancing across times
#'   and categories), \code{"coverage"} (in terms of number of cells observed
#'   under flow cytometry, if available) and \code{"cutoff"} (the percentage of
#'   the population below the cutoff point).  Additionally, \code{"all"} shows
#'   all these plots at once.
#' @param ... Additional arguments to the generics \code{summary} and
#'   \code{plot}, not used here.
#' @param main Plot title when \code{type="overview"}.
#'
#' @return An \code{fd_data} object.
#'
#' @section Format: \code{data} should be a \code{data.frame} (or a
#'   \code{groupedData}) with at least the following columns:
#'
#'   \describe{ \item{\code{Type}}{The type of the measurement.  Current allowed
#'   types are \code{"hists"} (histogram of living cells), \code{"hists_lost"}
#'   (histogram of dying/dead cells), \code{"Ns"} (cell counts/CFUs),
#'   \code{"Ns_lost"} (cell counts of dying/dead cells), \code{"props"}
#'   (proportions of living cells in any given number of generations),
#'   \code{"props_lost"} (same for dying/dead cells). In the case of histograms,
#'   one row represents a bin.  However, other types can added by the user,
#'   provided that a new processing function is given to \code{\link{fd_model}}
#'   through the argument \code{process}.}
#'
#'   \item{\code{Weight}}{Either \code{"hist"} (for types \code{"hists"} and
#'   \code{"hists_lost"}), \code{"N"}, \code{"prop"} or any other weight defined
#'   by the user (with the right processing function, see argument
#'   \code{process} of \code{\link{fd_model}}.}
#'
#'   \item{\code{Inoculum}}{An identifier for the inoculum/progenitors used.
#'   This allows to associate the histograms with the geometric mean and
#'   standard deviation of the initial population. Non-histogram based types can
#'   use \code{"none"} for this field.}
#'
#'   \item{\code{Timepoint}}{A unique identifier to group together all the
#'   measurements made together on the same sample.  Notice that if you do
#'   technical replicates in flow cytometry, the resulting histograms must be
#'   listed as different timepoints as within the same timepoint there would be
#'   no way to know which bin pertains to which histogram.}
#'
#'   \item{\code{Time}}{The time, in hours, at which the timepoint were taken.}
#'
#'   \item{\code{Category}}{The category/compartment (a factor) for the current
#'   row.}
#'
#'   \item{\code{a}, \code{b}}{The left (excluded) and right limits (included)
#'   of the bin if \code{Type \%in\% c("hists", "hists_lost")}, in linear units
#'   of fluorescence.  If \code{Type \%in\% c("props", "props_lost") }, \code{a}
#'   is the number of generations (starting from 0) and \code{b} is ignored. In
#'   the case of a different type, should NOT be \code{NA} as the whole row can
#'   be interpreted as missing data by the fitting algorithm or the fitting
#'   fails if the \code{na.action} is \code{\link[stats]{na.fail}}.  We
#'   recommend in this case to set \code{a} and \code{b} to 0, as they won't be
#'   read out by the algorithms anyway.}
#'
#'   \item{\code{y}}{If \code{Type \%in\% c("hists", "hists_lost") }, this is
#'   the proportion of cells in bin \code{(a, b]} (right-closed interval, the
#'   default of \code{\link[graphics]{hist}}). This is a probability and not a
#'   density and all the rows for a given histogram should sum to one (this is
#'   checked by this function).  Alternatively, if \code{Type \%in\% c("Ns",
#'   "Ns_lost") }, this is the cell count/CFU count and if \code{Type \%in\%
#'   c("prop", "prop_lost") }, this is the proportion of cells in a given number
#'   of generations. The \code{y} can be transformed, see below the
#'   \code{cctrans} and \code{htrans} parameters passed to the \code{FMM}. } }
#'
#'   The following \code{\link[base]{attributes}} can also be set (but setting
#'   them is not mandatory, default values are provided): \describe{
#'
#'   \item{\code{fmm}}{List the parameters that will be passed to the finite
#'   mixture model (FMM), see \code{\link{finite-mixture}}.}
#'
#'   \item{\code{proliferation}}{List the parameters that will be passed to the
#'   proliferation model, see \code{\link{proliferation}}.}
#'
#'   \item{\code{cutoff}}{A named numeric vector containing the cutoff points in
#'   linear units of fluorescence (the names are the timepoint identifiers).
#'   See \code{\link{cutoff}}.}
#'
#'   \item{\code{counts}}{If the user wants to store the number of events per
#'   histogram, they are advised to use this field.  It should then be a named
#'   numeric vector (the names would be the timepoint identifiers).}
#'
#'   }
#'
#' @export
#'
#' @examples
#' # FdSTyphimuriumWTC57 is already an fd_data, so the following
#' # does not do much
#' data(FdSTyphimuriumWTC57)
#' dat <- fd_data(FdSTyphimuriumWTC57)
#'
#' # Visually assess the data
#' plot(dat)
#'
#' # It is possible to perform standard operations on an fd_data,
#' # as you would do with a data.frame, but with additional checks
#' datsub <- subset(dat, Individual == "160408.WT.C57")
#' datbind <- rbind(datsub, subset(dat, Individual == "010708.WT.C57"))
#'
#' # To disable checks, use e.g. subset.data.frame
#' invisible(subset(dat, y > 0.1))
#'   # Warning: some histograms do not have proportions
#'   # that sum to 1
#' invisible(subset.data.frame(dat, y > 0.1))
#'   # Silently returns a data.frame (not an fd_data anymore)
#'
#' # Cutoff returns an fd_data for which histograms do not have
#' # proportions that sum to 1, without issuing any warning
#' invisible(cutoff(dat))
#'
#' # Sometimes, it is useful to extend the time domain for, e.g.,
#' # plotting
#' length(unique(dat$Time))
#' length(unique(fd_expand(dat, length.N = 50))$Time)
#' # This expanded dataset can then be used as the 'newdata' argument
#' # of 'predict'.
#'
#' # fd_transform can be used for arbitrary transformation
#' # Notice that here Ns are already log10-transformed
#' # (this is the default transform as defined in FMM)
#' dat_tr <- fd_transform(fd_simulate(fd_unif(), c(2, 12, 24)),
#'                        hist = "log1p")
fd_data <- function (data, categories = NULL, inoculums = NULL,
                     timepoints = NULL,
                     na.action = na.pass) {
  if (!is.data.frame(data))
    stop("'data' should be a 'data.frame' (or a 'groupedData')")
  if (NROW(data) == 0)
    return (data)

  reqcolnames <- c("Category", "Time", "Timepoint",
           "Inoculum", "a", "b", "y")
  if (!all(reqcolnames %in% colnames(data)))
    stop("'data' does not have the required column names ",
       paste(reqcolnames, collapse=", "))

  data$Time <- as.numeric(as.character(data$Time))
  data$a <- as.numeric(as.character(data$a))
  data$b <- as.numeric(as.character(data$b))

  data$Type <- ordered(data$Type)
  data$Weight <- factor(data$Weight)

  if (!is.null(categories))
    data$Category <- factor(data$Category, levels = categories)
  else if (!is.factor(data$Category))
    data$Category <- factor(data$Category)
  if (!is.null(inoculums))
    data$Inoculum <- factor(data$Inoculum, levels = inoculums)
  else if (!is.factor(data$Inoculum))
    data$Inoculum <- factor(data$Inoculum)
  if (!is.null(timepoints))
    data$Timepoint <- factor(data$Timepoint, levels = timepoints)
  else if (!is.factor(data$Timepoint))
    data$Timepoint <- factor(data$Timepoint)

  # Perform some checks
  if (any(!is.finite(data$y[data$Weight == "N"]))) {
    warning("Some 'Ns' are -Inf: you probably want to use a ",
        "'log1p' transform.")
  }
  if (any(!is.finite(data$y[data$Weight == "hist"]))) {
    warning("Some 'y' are not finite for histograms.")
  }
  if (any(!is.finite(data$y[data$Weight == "prop"]))) {
    warning("Some 'y' are not finite for proportions.")
  }
  if (any(is.na(data))) {
    warning("Some fields are not specified (NA), this could ",
        "lead to errors.")
  }
  if (suppressWarnings(
      any(!is.na(as.numeric(as.character(data$Category)))))) {
    warning("Some categories are convertible to 'numeric'; we advise ",
        "against this convertibility.  Perhaps you should prefix ",
        "current names with a character string such as ",
        "\"cate_\".")
  }
  if (suppressWarnings(
      any(!is.na(as.numeric(as.character(data$Individual)))))) {
    warning("Some individuals are convertible to 'numeric'; we advise ",
        "against this convertibility.  Perhaps you should prefix ",
        "current names with a character string such as ",
        "\"indiv_\".")
  }
  if (suppressWarnings(
      any(!is.na(as.numeric(as.character(data$Inoculum)))))) {
    warning("Some inoculums are convertible to 'numeric'; we advise ",
        "against this convertibility.  Perhaps you should prefix ",
        "current names with a character string such as ",
        "\"inoc_\".")
  }
  if (suppressWarnings(
      any(!is.na(as.numeric(as.character(data$Timepoint)))))) {
    warning("Some timepoints are convertible to 'numeric'; we advise ",
        "against this convertibility.  Perhaps you should prefix ",
        "current names with a character string such as ",
        "\"tp_\".")
  }
  inocs <- as.character(levels(
    subset.data.frame(data, Weight == "hist")$Inoculum))
  if (length(inocs) > 0L) {
    m0 <- attr(data, "fmm")$m0
    if (!is.null(m0) && (length(m0) != length(inocs) ||
               is.null(names(m0)) ||
               any(names(m0) != inocs))) {
      warning("'fmm$m0' does not have matching names with ",
          "'data$Inoculum' (levels must be ",
          "in the exact same order): ",
          "fmm$m0 names = ", dput(names(m0)),
          "; data$Inoculum names = ", dput(inocs), ".")
    }
    sd0 <- attr(data, "fmm")$sd0
    if (!is.null(sd0) && (length(sd0) != length(inocs) ||
                is.null(names(sd0)) ||
               any(names(sd0) != inocs))) {
      warning("'fmm$sd0' does not have matching names with ",
          "'data$Inoculum' (levels must be ",
          "in the exact same order): ",
          "fmm$sd0 names = ", dput(names(sd0)),
          "; data$Inoculum names = ", dput(inocs), ".")
    }
  }
  if (is.null(attr(data, ".was_cutoff"))) {
    df_hists_sum <- plyr::ddply(subset.data.frame(data, Weight == "hist"),
                                "Timepoint",
                                plyr::summarise,
                                y = sum(y))
    # nolint start
    if (any(abs(df_hists_sum$y - 1) > .Machine$double.neg.eps)) {
      warning("Some histograms do not have proportions that sum to 1.")
    }
    # nolint end
  }
  df_props_sum <- plyr::ddply(subset.data.frame(data, Weight == "prop"),
                              "Timepoint",
                              plyr::summarise,
                              y = sum(y))
  if (any(abs(df_props_sum$y - 1) > 1e-5)) {
    warning("Some timepoints do not have proportions that sum to 1.")
  }

  structure(na.action(data),
        class = c("fd_data", setdiff(class(data), "fd_data")))
}

#' @export
#' @keywords internal
subset.fd_data <- function (x, ...) {
  ans <- NextMethod()
  fd_data_refactor(ans, x)
}

#' @export
#' @keywords internal
`[.fd_data` <- function (x, i, j, ...) {
  ans <- NextMethod()
  if (!missing(j)) {
    ans
  } else {
    fd_data_refactor(ans, x)
  }
}

#' @export
#' @keywords internal
rbind.fd_data <- function (..., deparse.level = 1) {
  ans <- rbind.data.frame(..., deparse.level = 1)
  arglst <- list(...)
  attr(ans, "fmm") <- list(
    sd0 = unlist(lapply(arglst, function (cur) attr(cur, "fmm")$sd0)),
    m0 = unlist(lapply(arglst, function (cur) attr(cur, "fmm")$m0))
  )
  attr(ans, "proliferation") <- attr(arglst[[1L]], "proliferation")
  attr(ans, "cutoff") <-
    unlist(lapply(arglst, function (cur) attr(cur, "cutoff")))
  return (ans)
}

#' @export
#' @rdname fd_data
fd_transform <- function (x, hist = "identity", N = "log10") {
  data <- x

  if (is.character(hist))
    hist <- match.fun(paste0(hist, "_trans"))()
  if (is.character(N))
    N <- match.fun(paste0(N, "_trans"))()

  # Inverse current transformation
  if (is.null(old <- attr(data, "fmm")$htrans))
    old <- "identity"
  if (is.character(old)) old <- match.fun(paste0(old, "_trans"))()
  data[data$Weight == "hist", ]$y <-
    old$inverse(data[data$Weight == "hist", ]$y)

  if (is.null(old <- attr(data, "fmm")$cctrans))
    old <- "log10"
  if (is.character(old)) old <- match.fun(paste0(old, "_trans"))()
  data[data$Weight == "N", ]$y <- old$inverse(data[data$Weight == "N", ]$y)

  # Apply new transformation
  data[data$Weight == "hist", ]$y <-
    hist$transform(data[data$Weight == "hist", ]$y)
  data[data$Weight == "N", ]$y <-
    N$transform(data[data$Weight == "N", ]$y)

  attr(data, "fmm")$htrans <- hist
  attr(data, "fmm")$cctrans <- N

  data
}

#' @rdname fd_data
#' @export
cutoff <- function (x, cutoff = NULL, threshold = 0.01) {
  data <- x

  if (!is.null(cutoff)) {
    if (is.null(names(cutoff))) {
      if (length(cutoff) > 1) {
        stop("'cutoff' needs to be named")
      }
      m0s <- attr(data, "fmm")$m0
      if (is.null(m0s))
        stop("no 'fmm$m0' to calculate cutoff")
      cutoff <- m0s / 2.0 ^ cutoff
    }
  } else {
    cutoff <- attr(data, "cutoff")
    if (is.null(attr(data, "cutoff")))
      stop("'cutoff' called but the dataset presents no ",
         "'cutoff' attribute")
    if (any(names(cutoff) !=
        levels(subset.data.frame(data, Weight == "hist")$Timepoint)))
      stop("'cutoff' names and 'Timepoint' levels are not the same ",
         "or exactly in the same order")
  }
  get_cutoff <- if (length(cutoff) == 1) {
    function (tp) cutoff
  } else {
    function (tp) {
      ans <- cutoff[as.character(tp)]
      if (is.null(ans)) {
        stop(
          "Timepoint ",
          dput(tp),
          " was not found in 'cutoff'; available timepoints: ",
          paste0(collapse = ", ", names(cutoff)))
      }
      ans
    }
  }
  attr(data, ".was_cutoff") <- TRUE
  subset(data, Weight != "hist" |
         (a > get_cutoff(Timepoint) &
          is.finite(a) & y >= threshold))
}

#' @rdname fd_data
#' @export
fd_expand <- function (x, length.N = 50, seq = NULL,
                       separate = TRUE, by=NULL) {
  data <- x

  if (is.null(by) && inherits(data, "groupedData"))
    by <- deparse(getGroupsFormula(data)[[2L]])
  if (!separate) {
    if (is.null(seq))
      seq <- match.fun("seq")(min(data$Time),
                  max(data$Time), length.out = length.N)
    lst <- list(Category = unique(data$Category),
          Time = seq)
    if (!is.null(by))
      lst[[by]] <- unique(data[[by]])
    dat_Ns <- transform(
        do.call(expand.grid, lst),
        a = 0, b = 0, Inoculum = "none",
        Weight = "N", Type = "Ns", y = 0)
  } else {
    if (is.null(by))
      stop("nothing to separate")
    dat_Ns <-
      plyr::ddply(data, by, function (cur) {
        if (is.null(seq))
          seq <- match.fun("seq")(min(cur$Time),
                      max(cur$Time),
                      length.out = length.N)

        lst <- list(Category = unique(cur$Category),
              Time = seq)
        lst[[by]] <- unique(cur[[by]])
        transform(do.call(expand.grid, lst),
              a = 0, b = 0, Inoculum = "none",
              Weight = "N", Type = "Ns", y = 0)
      })
  }

  dat_Ns$Timepoint <- paste0("Ns_", 1:NROW(dat_Ns))
  ans <- rbind(
    subset.data.frame(data, Weight == "hist")[,
      c(
        "Category", "Time", "Individual", "a", "b",
        "Inoculum", "Weight", "Type", "y", "Timepoint"
      )],
    dat_Ns
  )
  ans$Category <- factor(ans$Category, levels=levels(data$Category))
  fd_data_refactor(ans, data)
}


# Internal ----------------------------------------------------------------


fd_data_refactor <- function (ans, x, refactor=TRUE, categories=FALSE) {
  for (nm in setdiff(names(attributes(x)),
             c(names(attributes(ans)), "class"))) {
    attr(ans, nm) <- attr(x, nm)
  }
  if (refactor) {
    ans$Inoculum <- factor(ans$Inoculum)
    inoculums <- setdiff(levels(ans$Inoculum), "none")
    if (!is.null(attr(ans, "fmm")$sd0)) {
      attr(ans, "fmm")$sd0 <-
        attr(ans, "fmm")$sd0[as.character(inoculums)]
    }
    if (!is.null(attr(ans, "fmm")$m0))
      attr(ans, "fmm")$m0 <- attr(ans, "fmm")$m0[as.character(inoculums)]
    if (categories && !is.null(attr(ans, "pro")$categories))
      attr(ans, "pro")$categories <-
        intersect(attr(ans, "pro")$categories,
              unique(ans$Category))
  }
  fd_data(ans)
}