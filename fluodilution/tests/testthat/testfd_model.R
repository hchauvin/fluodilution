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

context("fd_model")

test_that("weird parameters get worked out", {
  expect_silent(fd_model(FdSTyphimuriumWTC57,
                         fmm="gaussian",
                         proliferation="cyton"))
  expect_error(fd_model(data=cos,
                        fmm="gaussian",
                        proliferation="cyton"))
  expect_error(fd_model(FdSTyphimuriumWTC57,
                        fmm="kldsfsdfs",
                        proliferation="cyton"))
  expect_error(fd_model(FdSTyphimuriumWTC57,
                        fmm="gaussian",
                        proliferation="ldsfkdsfsd"))
  expect_error(fd_model(FdSTyphimuriumWTC57,
                        fmm="gaussian",
                        proliferation="cyton",
                        constraints="jdsfksdf"))
  expect_error(fd_model(FdSTyphimuriumWTC57,
                        fmm="gaussian",
                        proliferation="cyton",
                        boxed="ezrzef"))
  expect_error(fd_model(FdSTyphimuriumWTC57,
                        fmm="gaussian",
                        proliferation="cyton",
                        partial="ezrzef"))
  expect_error(fd_model(FdSTyphimuriumWTC57,
                        fmm="gaussian",
                        proliferation="cyton",
                        partial=c(1,6)))
  expect_error(fd_model(FdSTyphimuriumWTC57,
                        fmm="gaussian",
                        proliferation="cyton",
                        partial=c(1, 1)))
})

test_that("additional checks on constraints", {
    expect_error(
      fd_model(FdSTyphimuriumWTC57,
               fmm="gaussian",
               proliferation="branching",
      constraints = cstrlist(start = ~ pro:all:{p0 <- -1})),
      "'start' is not within model boundary")
    expect_error(
      fd_model(FdSTyphimuriumWTC57,
               fmm="gaussian",
               proliferation="branching",
               constraints = cstrlist(upper = ~ pro:all:{p0 <- -1},
                                      lower = ~ pro:all:{p0 <- 0})),
      "'start' is not within model boundary")
    expect_error(
      fd_model(FdSTyphimuriumWTC57,
               fmm="gaussian",
               proliferation="branching",
               constraints = cstrlist(upper = ~ pro:all:{dflsdf <- -1})),
      "lower, upper and start have different structure")
})

test_that("'fd_clean' works", {
    mdl <- fd_model(FdSTyphimuriumWTC57,
                    fmm="gaussian",
                    proliferation="branching")
    expect_true(is.null(mdl$assigncstr) && is.null(mdl$mem_diff_cdf))
    s <- fd_assigncstr(mdl)
    expect_equal(deparse(s), deparse(mdl$assigncstr))
    fd_clean(mdl)
    expect_true(is.null(mdl$assigncstr) && is.null(mdl$mem_diff_cdf) &&
                is.null(mdl$mem_diff_cdf))
})

test_that("transformation works", {
    mdl <- fd_model(FdSTyphimuriumWTC57, fmm="gaussian", proliferation="cyton")
    expect_equal(mdl$start, relist(start(mdl), mdl), tolerance=1e-6)
    expect_equal(mdl$lower, relist(lower(mdl), mdl), tolerance=1e-6)
    expect_equal(mdl$upper, relist(upper(mdl), mdl), tolerance=1e-6)
})

test_that("'fd_predict_mat' works", {
    dat <- FdSTyphimuriumWTC57
    mdl <- fd_model(dat,
                    fmm="gaussian",
                    proliferation="branching")
    expect_error(fd_predict_mat(mdl)(dat, start(mdl)),
                 "'params' must be a matrix")
    ans <- fd_predict_mat(mdl)(dat, rbind(start(mdl)))
    expect_is(ans, "numeric")
    expect_true(length(ans) == NROW(dat))
})

test_that("'fd_clone' works", {
    mdl <- fd_model(FdSTyphimuriumWTC57,
                    fmm="gaussian",
                    proliferation="branching")
    expect_silent(mdl2 <- fd_clone(mdl))
})

test_that("'update.fd_model' works", {
    mdl <- fd_model(FdSTyphimuriumWTC57,
                    fmm="gaussian",
                    proliferation="branching")

    mdl1 <- update(mdl, FdSTyphimuriumWTC57)
    expect_equal(mdl1$pro, mdl$pro)
    expect_equal(mdl1$fmm, mdl$fmm)
    expect_equal(mdl1$partial, mdl$partial)
    expect_equal(mdl1$boxed, mdl$boxed)
    expect_equal(mdl1$process, mdl$process)
    expect_equal(mdl1$start, mdl$start)
    expect_equal(mdl1$lower, mdl$lower)
    expect_equal(mdl1$upper, mdl$upper)

    mdl2 <- update(mdl,
                   subset(FdSTyphimuriumWTC57, Individual == "010708.WT.C57"),
                   fmm="af_bp",
                   proliferation="cyton",
                   process = "cos",
                   boxed = FALSE)
    expect_true(attr(mdl2$pro, "expanded_name") == "fd_proliferation_cyton")
    expect_true(attr(mdl2$fmm, "expanded_name") == "fd_fmm_af_bp")
    expect_equal(mdl2$process, "cos")
    expect_equal(mdl2$processname, "cos")
    expect_equal(mdl2$boxed, FALSE)

    mdl3 <- update(mdl, FdSTyphimuriumWTC57)
    expect_equal(start(mdl3), start(mdl))
    expect_equal(upper(mdl3), upper(mdl))
    expect_equal(lower(mdl3), lower(mdl))
    expect_equal(relist(start(mdl3), mdl3), relist(start(mdl), mdl))
    expect_equal(relist(lower(mdl3), mdl3), relist(lower(mdl), mdl))
    expect_equal(relist(upper(mdl3), mdl3), relist(upper(mdl), mdl))
})

test_that("'fd_predict' works", {
    mdl <- fd_model(cutoff(FdSTyphimuriumWTC57),
                    fmm="gaussian",
                    proliferation="branching")

    dat <- cutoff(FdSTyphimuriumWTC57)
    dat$Id <- 1:NROW(dat)
    ans <- fd_predict(mdl, start(mdl), dat)
    expect_equal(
        dplyr::select(as.data.frame(ans[order(ans$Id), ]), -y),
        dplyr::select(as.data.frame(dat), -y))
})

test_that("'start', etc. work", {
    mdl <- fd_model()
    expect_equal(relist(start(mdl), mdl), start(mdl, free=FALSE))
    expect_equal(relist(lower(mdl), mdl), lower(mdl, free=FALSE))
    expect_equal(relist(upper(mdl), mdl), upper(mdl, free=FALSE))
})
