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

context("fd_simulate")

test_that("structure is well-formed", {
  un <- fd_unif("gaussian", "branching", n=5)
  expect_silent(ans <- fd_simulate(
    un,
    c(2, 4, 6, 12),
    select=c("hists", "hists_lost", "Ns", "Ns_lost"),
    mean=FALSE))
  expect_is(ans, "groupedData")
  expect_is(ans, "fd_data")
  expect_equal(levels(ans$Category), "One")
  expect_equal(unique(ans$Time), c(2, 4, 6, 12))
  expect_equal(levels(ans$Individual), paste0("indiv_", 1:5))
  expect_equal(levels(ans$Type), c("hists", "hists_lost", "Ns", "Ns_lost"))
  expect_equal(levels(ans$Weight), c("hist", "N"))
  expect_equal(levels(ans$Inoculum), c("inoc_1", "none"))
  expect_equal(colnames(ans), c("Time", "Category", "Individual", "a", "Type",
                                "Weight", "Inoculum", "y", "b", "Timepoint"))

  expect_silent(ans2 <- fd_simulate(
    un,
    c(2, 4, 6, 12),
    select=c("hists", "hists_lost", "Ns", "Ns_lost"),
    mean=TRUE))
  expect_equal(levels(ans2$Individual), "indiv")
})

test_that("different param formats can work", {
  un <- fd_unif("gaussian", "branching", n=1)

  # No model specified
  expect_error(fd_simulate(un[1, ], 2),
                "no model specified, either through 'params' or 'model'")

  # As vector
  invisible(fd_simulate(un[1, ], 2, model = model(un)))

  # As fit with coef()
  un2 <- list(coefficients = un[1, ], .model = model(un))
  invisible(fd_simulate(un2, 2))

  # As matrix
  un3 <- fd_unif("gaussian", "branching", n=5)
  invisible(fd_simulate(as.matrix(un3), 2, model = model(un3)))

  # Error
  expect_error(
    fd_simulate("ldfsdfsd", 2),
    "no model specified, either through 'params' or 'model'")
  expect_error(
    fd_simulate(list(.model = "ldfsdfsd"), 2),
    "'params' must either be the result of a fit, a vector or a matrix")
})

test_that("various ranges can be specified", {
  un <- fd_unif("gaussian", "branching", n=1)

  # Automatic breaks
  expect_silent(fd_simulate(un, 24, breaks=2))
  expect_silent(fd_simulate(un, 24, breaks=1000))

  # Manual breaks
  ans <- fd_simulate(un, c(12, 24), breaks=c(13, 33, 211, 8500))
  expect_equal(sort(unique(ans$a)), c(0, 13, 33, 211))
  expect_equal(sort(unique(ans$b)), c(0, 33, 211, 8500))

  # Range specified
  ans <- fd_simulate(un, c(12, 24), breaks=10, range=c(13, 8500))
  expect_equal(
    sort(unique(ans$a)),
    c(0, 13, 26.7418857331649, 54.9625952030159, 112.941629721289,
      232.070501075511, 476.848996135572, 979.807210941429, 2013.26114488653,
      4136.75239121832))
  expect_equal(
    sort(unique(ans$b)),
    c(0, 26.7418857331649, 54.9625952030159, 112.941629721289,
      232.070501075511,
      476.848996135572, 979.807210941429, 2013.26114488653, 4136.75239121832,
      8500))
})

test_that("'noise' works", {
  un <- fd_unif("gaussian", "branching", n=1)
  expect_silent(fd_simulate(un, 24, noise=c(0.05, 0.01)))
})

test_that("it can cope with multiple inoculums", {
  data(FdSTyphimuriumWTC57)
  mdl <- fd_model(cutoff(FdSTyphimuriumWTC57),
                  fmm="gaussian",
                  proliferation="branching")
  ans <- fd_simulate(start(mdl), 24, model=mdl)
  expect_equal(levels(ans$Inoculum), c("inoc_1", "none"))
  expect_setequal(names(attr(ans, "fmm")), c("m0", "sd0", "cctrans", "htrans"))
  expect_equal(
    attr(ans, "fmm")$m0,
    structure(6.32986737172484, .Names = "inoc_1"))
  expect_equal(
    attr(ans, "fmm")$sd0,
    structure(0.524521260581569, .Names = "inoc_1"))
  expect_equal(class(attr(ans, "fmm")$cctrans), "trans")
  expect_equal(attr(ans, "fmm")$cctrans$name, "identity")
  expect_equal(attr(ans, "fmm")$cctrans$transform(c(10, 20)), c(10, 20))
  expect_equal(attr(ans, "fmm")$cctrans$inverse(c(10, 20)), c(10, 20))
  expect_equal(class(attr(ans, "fmm")$htrans), "trans")
  expect_equal(attr(ans, "fmm")$htrans$name, "identity")
  expect_equal(attr(ans, "fmm")$htrans$transform(c(10, 20)), c(10, 20))
  expect_equal(attr(ans, "fmm")$htrans$inverse(c(10, 20)), c(10, 20))
})

test_that("'fd_proportions' works", {
  un <- fd_unif("gaussian", "branching", n=1)
  expect_silent(fd_proportions(un, 24))
})

test_that("'fd_norm' works", {
  un <- fd_unif("gaussian", "branching", n=1)
  expect_silent(fd_norm(10, un))
})
