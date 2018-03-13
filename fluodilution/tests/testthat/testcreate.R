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

context("create")

unzip(system.file("extdata", "Archive.zip", package="fluodilution"),
      exdir = tempdir())

test_that("list of flows can be loaded", {
  library(flowCore)

  channel <- "FL2-H"
  gate <- rectangleGate(`FL1-H` = c(18, Inf), filterId="Bacteria")
  src <- read.flowSet(paste0(tempdir(), "/",
                 c("WT-t12.023", "WT-t2.016", "WT-t6.014")))
  flow <- Subset(src, flowCore::filter(src, gate))

  meta <- data.frame(Category = "One", Time = c(12, 2, 6), Type = "hists",
             Inoculum = "inoc_1", Timepoint = c("tp_1", "tp_2", "tp_3"))
  value <- list(flow = structure(flow, meta = meta, channel = channel))
  ans <- fd_create(value = value)
  expect_silent(plot(ans))

  # Stress tests
  expect_warning(fd_create(value, fmm = list(m0 = 1, sd0 = 1)))
  expect_silent(
    fd_create(value, fmm = list(m0 = c(inoc_1 = 1), sd0 = c(inoc_1 = 1))))
  expect_warning(fd_create(value, categories="klqsd"))
  expect_equal(
    levels(fd_create(append(value, list(N = data.frame())))$Weight),
    "hist")
  df_N <- data.frame(a = 0, b = 0, y = c(1, 2, 3),
                     Time = c(4, 5, 6),
                     Inoculum = "none",
                     Type = "Ns",
                     Category = "One",
                     Timepoint = "tp_Ns")
  expect_silent(
    ans3 <- fd_create(append(value,
                      list(N = df_N))))
  expect_equal(
    subset.data.frame(ans3, Weight == "N") %>%
      dplyr::select(-Weight) %>% as.character(),
    df_N %>% as.character())
})

test_that("workspaces can be loaded", {
  library(flowWorkspace)

  # Create the Gating Set
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
      as.numeric(getKeywords(ws, sub("_[0-9]*$", "", flow[[i]]@name))$Time)
    }
  )
  print(getNodes(flow))
  meta <- data.frame(Time = times)

  # Go for it
  ans <- fd_create(
    value = list(flow = structure(flow, meta = meta, channel = channel))
  )
  expect_silent(plot(ans))

  # Stress tests
  ans2 <- fd_create(list(
    flow = structure(flow,
             meta = transform(meta, Gate = "cells"),
             channel = channel)
  ))
  expect_true(all(all.equal(ans$y, ans2$y) != TRUE))
  ans3 <- fd_create(list(
    flow = structure(flow,
             meta = transform(meta, Gate = "root"),
             channel = channel)
  ))
  expect_equal(ans3, ans)

  # "use0" works with one inoculum
  library(parallel)
  expect_silent(ans2 <- fd_create(
    value = list(flow = structure(flow, meta = meta, channel = channel)),
    fmm = "use0"
  ))
  expect_true(names(attr(ans2, "fmm")$m0) == "inoc_1" &&
    names(attr(ans2, "fmm")$sd0) == "inoc_1")

  expect_silent(ans3 <- fd_create(
    value = list(flow = structure(flow, meta = meta, channel = channel)),
    fmm = "use0",
    momentControl = list(clean=FALSE)
  ))
  expect_true(names(attr(ans3, "fmm")$m0) == "inoc_1" &&
    names(attr(ans3, "fmm")$sd0) == "inoc_1")
  expect_true(
    all(all.equal(
      unlist(attr(ans2, "fmm")),
      unlist(attr(ans3, "fmm"))) != TRUE))
  expect_equal(structure(ans2, fmm=NULL), structure(ans3, fmm=NULL))

  # It works with two inoculums as well
  expect_warning(
    fd_create(
      value =
        list(flow = structure(flow,
                    meta = data.frame(Time = c(0, 1, 0)),
                    channel = channel)),
      fmm = "use0"
    ),
    "more than one 't=min/0h' timepoint for the same inoculum"
  )
  expect_silent(ans2 <- fd_create(
    value =
      list(flow = structure(flow,
           meta = data.frame(Time = c(0, 1, 0),
                             Inoculum = c("A", "A", "B")),
           channel = channel)),
    fmm = "use0"
  ))
  expect_equal(levels(ans2$Inoculum), c("A", "B"))
  expect_equal(as.character(unique(ans2$Inoculum)), "A")
  fmm <- attr(ans2, "fmm")
  expect_true(length(fmm$m0) == 2 && names(fmm$m0) == c("A", "B") &&
              fmm$m0[1] != fmm$m0[2])
  expect_true(length(fmm$sd0) == 2 && names(fmm$sd0) == c("A", "B") &&
              fmm$sd0[1] != fmm$sd0[2])
})

test_that("twisted workspaces can be loaded as well", {
  library(flowWorkspace)

  # Create the Gating Set
  ws <- openWorkspace(paste0(tempdir(), "/20150429.wsp"))
  print(ws)
  flow <- suppressWarnings(parseWorkspace(
    ws, "All Samples", path = tempdir(),
    isNcdf = FALSE,
    cleanup = FALSE, keep.indices = TRUE,
    requiregates = FALSE
  ))

  # Meta
  channel <- "525_50 B-A"
  library(parallel)
  expect_silent(x <- fd_moments(flowWorkspace::getData(flow[1]), channel))
  expect_silent(x2 <- fd_moments(flowWorkspace::getData(flow), channel))
  expect_equal(colnames(x2), c("mean", "sd"))
  expect_equal(dim(x2), c(3, 2))
  expect_silent(x3 <- fd_moments(flowWorkspace::getData(flow[[1L]]), channel))
  expect_true(is.vector(x3))
  expect_equal(names(x3), c("mean", "sd"))
  meta <- data.frame(Time = c(24, 2),
             Category = c("Spleen", "MLN"))
  expect_silent(ans <- fd_create(
    list(
      flow = structure(flow[-1],
                       meta = meta,
                       channel = channel)
    ),
    fmm = list(m0 = setNames(x$mean, "inoc_1"),
               sd0 = setNames(x$sd, "inoc_1"))
  ))
  expect_silent(plot(ans))
})

