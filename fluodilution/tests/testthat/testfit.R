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

context("fit")

source(system.file("contrib", "nlsSA.R", package="fluodilution"))

# Base constraints for Salmonella Typhimurium experiments
`#macr` <<- ~ FdCommonConstraints$`#noss` + pro:all:{res <- 0} + fmm:{c0 <- 0}
`#invivo` <<- ~ FdCommonConstraints$`#noss` + ~ pro:all:{res <- 0}
CC <<- FdCommonConstraints

test_that("'fd_minuslogl' works", {
  objective <- -21.33836

  # CC <- FdCommonConstraints
  cstr <- ~ `#macr` + CC$`#delta_1100` + CC$`#mm_xyxz`
  dat <- subset(cutoff(FdSTyphimuriumWTC57), Individual == "140508.WT.C57")
  mdl <- fd_model(dat, fmm="gaussian", proliferation="branching",
                    constraints = cstr)
  fun <- fd_minuslogl(mdl, dat, verbose=FALSE)
  expect_equal(do.call(fun, as.list(start(mdl))), objective, tolerance = 1e-5)

  fun <- fd_minuslogl(mdl, dat, verbose=TRUE)
  stats4::mle(fun, method="Nelder-Mead", control=list(maxit=1))

  ans <- bbmle::mle2(fun, start=formals(fun),
          optimfun=function (method, hessian, gr, ...) GenSA::GenSA(...),
          optimizer="user",
          lower=lower(mdl), upper=upper(mdl),
          control = list(maxit=1),
          skip.hessian = TRUE)

  # FIXME: bug in bbmle as of now
  ans@details$convergence <- 0

  print(ans)

  # Stress test: about memoise
  expect_silent(mdl <- fd_model(dat, fmm="gaussian", proliferation="branching",
                    constraints = cstr, memoise = FALSE))
  fun <- fd_minuslogl(mdl, dat, verbose=FALSE)
  expect_equal(do.call(fun, as.list(start(mdl))), objective, tolerance = 1e-5)

  expect_silent(mdl <- fd_model(dat, fmm="gaussian", proliferation="branching",
                    constraints = cstr, memoise = TRUE))
  fun <- fd_minuslogl(mdl, dat, verbose=FALSE)
  expect_equal(do.call(fun, as.list(start(mdl))), objective, tolerance = 1e-5)
})

test_that("'fd_nls' works", {
  # CC <- FdCommonConstraints
  cstr <- ~ `#macr` + CC$`#delta_1100` + CC$`#mm_xyxz`
  dat <- subset(cutoff(FdSTyphimuriumWTC57), Individual == "140508.WT.C57")
  mdl <<- fd_model(dat, fmm="gaussian", proliferation="branching",
                    constraints = cstr)
  ans <- fd_nls("mdl", dat, mgen=16, loop=0, trace=TRUE,
          nlsfun=nlsSA,
          control=list(maxit=1))

  # Test not enough mgen
  expect_warning(fd_nls("mdl", dat, mgen=3, loop=1, trace=TRUE,
          nlsfun=nlsSA,
          control=list(maxit=1)))

  # Nota: Other functions (i.e. guess, etc.) are well-tested by the article
  # expect_equal(guess_mGen(cutoff(FdSTyphimuriumWTC57)), 17)
  expect_equal(guess_mGen(cutoff(FdSTyphimuriumWTC57)), 28)
  expect_silent(invisible(fd_rates(cutoff(FdSTyphimuriumWTC57))))
  expect_silent(invisible(fd_rates(cutoff(FdSTyphimuriumWTC57),
                                   cut_first = TRUE, group_hists = FALSE)))
})

test_that("'fd_comb' works", {
  cstr <- list(~ `#macr` + CC$`#delta_1100` + CC$`#mm_xyxz`,
                ~ `#macr` + CC$`#delta_1100` + CC$`#mm_xxxx`)
  if (.Platform$OS.type == "unix") {
      cores <- 2
  } else {
      cores <- 1 # no forking on Windows!
  }
  invisible(fd_comb(cstr, data=cutoff(FdSTyphimuriumWTC57),
          continue = TRUE,
          mc.cores=cores,
          file=NULL,
          nlsfun=nlsSA,
          control=list(maxit=1)))
  expect_null(
    suppressWarnings(
      dplyr::failwith(default=NULL, fd_comb)(
        cstr,
        data=cutoff(FdSTyphimuriumWTC57),
        continue = FALSE,
        mc.cores=1,
        file=NULL,
        nlsfun=nls,  # nls, not nlsSA, for the thing to fail
        control=list(maxit=1))))
  monitor <- tempfile("monitor")
  fd_comb(cstr,
          data=cutoff(FdSTyphimuriumWTC57),
          continue = TRUE,
          mc.cores=cores,
          file=monitor,
          nlsfun=nlsSA,
          control=list(maxit=1))
  expect_true(file.exists(monitor))
  file.remove(monitor)
})

test_that("'relisted_coef' works", {
  # CC <- FdCommonConstraints
  cstr <- ~ `#macr` + CC$`#delta_1100` + CC$`#mm_xyxz`
  dat <- subset(cutoff(FdSTyphimuriumWTC57), Individual == "140508.WT.C57")
  mdl <<- fd_model(dat, fmm="gaussian", proliferation="branching",
                   constraints = cstr)
  ans <- fd_nls("mdl", dat, mgen=16, loop=0, trace=TRUE,
                nlsfun=nlsSA,
                control=list(maxit=1))
  expect_silent(relisted_coef(ans))
  expect_silent(relisted_coef(ans, drop=FALSE))
})

test_that("it works for _in vivo_ as well", {
  load(system.file("extdata", "FdSTyphimuriumMice.rda",
                   package="fluodilution"))
  dat <- fd_transform(cutoff(FdSTyphimuriumMice), N = "log1p")

  plot(dat)

  # Normal stuff
  mdl <- fd_model(data = dat, constraints = `#invivo`)
  do.call(fd_minuslogl(mdl, dat, stop_boundary=TRUE), as.list(start(mdl)))

  # With partial on
  get_categories <- function (co) {
    stringr::str_match(names(co), "^pro\\.([^.]+)")[, 2]
  }
  partial <- 1:2
  dat2 <- subset(dat, as.numeric(Category) %in% partial)
  mdlp <- fd_model(data = dat2, constraints = `#invivo`,
                   partial = partial)
  expect_true(all(get_categories(start(mdlp)) %in% c(NA, "SI", "PP")))
  do.call(fd_minuslogl(mdlp, dat2, stop_boundary=TRUE), as.list(start(mdlp)))

  # Weird partial
  partial <- 2L:4L
  expect_error(fd_model(
    data = subset(dat, as.numeric(Category) %in% partial),
    constraints = `#invivo`,
    partial = partial))

  # fixcstr
  expect_silent(fixed <- fixcstr(`#invivo`, start(mdlp)))
  expect_equal(
    fixcstr(cstrlist(`#invivo`, start=~{test <- 5}), start(mdlp)),
    cstrlist(fixed, start=~{test <- 5}))
  expect_false(
    paste(deparse(fixed), collapse="\n") ==
      paste(deparse(
        fixcstr(`#invivo`, start(mdlp), before=FALSE)),
        collapse="\n"))
  mdl2 <- fd_model(data = dat,
                   constraints = fixcstr(`#invivo`, start(mdlp)))
  expect_true(all(get_categories(start(mdl2)) %in% c("MLN", "Spleen")))
  expect_equal(relist(lower(mdl2), mdl2)$pro[c("SI", "PP")],
               relist(start(mdl), mdl)$pro[c("SI", "PP")])
  do.call(fd_minuslogl(mdl2, dat2, stop_boundary=TRUE), as.list(start(mdl2)))

  # Another way to call it
  expect_silent(fixcstr(`#invivo`, list(coefficients = start(mdl))))

  # Adding a transport compartment
  dat_tr <- fd_data(subset(dat, Category %in% c("SI", "PP")),
                    categories = c("SI", "Trans0", "PP"))
  mdl_tr <- fd_model(data = dat_tr,
                     constraints = `#invivo`,
                     partial = 1:3)
  expect_equal(
    unique(get_categories(start(mdl_tr))),
    c(NA, "SI", "Trans0", "PP"))
  do.call(
    fd_minuslogl(mdl_tr, dat_tr, stop_boundary=TRUE),
    as.list(start(mdl_tr)))

  # Precalculate: initialize
  mdl3 <- fd_model(data = dat,
                   constraints = fixcstr(`#invivo`, start(mdlp)))
  freecategories <- na.omit(unique(stringr::str_match(names(start(mdl3)),
                                              "^pro\\.([^.]+)")[, 2]))
  expect_equal(freecategories, c("MLN", "Spleen"))
  mdl3$pro$precalculate(mdl2$pro$trans_inverse(mdl3$start$pro),
                        min(match(freecategories, levels(dat$Category))),
                        sort(unique(dat$Time)),
                        mgen = 16L)
  expect_true(!is.null(mdl3$pro$precalc))
  expect_true(is.null(mdl2$pro$precalc))

  # Precalculate: ensure we get the same results
  for (i in 1:3) {
    rand <- fd_draw_unif(mdl2, 1)[1, ]
    pred1 <- fd_predict(
      mdl2, param = rand, data = dat, mgen = 16L)
    pred2 <- fd_predict(
      mdl3, param = rand[grepl("^pro\\.(MLN|Spleen)\\.", names(rand))],
      data = dat, mgen = 16L)
    expect_equal(pred1, pred2)
  }

  # Precalculate: stress test
  expect_error(fd_predict(mdl3, param = start(mdl3), data = dat, mgen = 13L))
  expect_silent(
      mdl3$pro$model(relist(start(mdl), mdl)$pro,
                     times = c(2, 6, 12, 24, 48), mgen = 16L))
  expect_error(
      mdl3$pro$model(relist(lower(mdl), mdl)$pro,
                     times = c(2, 6, 12, 24, 48), mgen = 16L)
  )

  # Other stress tests
  expect_silent(fd_predict(mdl, param = start(mdl), data = dat, mgen = NULL))
  expect_error(fd_predict(mdl, param = start(mdl), data = dat, mgen = "dsfsdf"))
  expect_error(fd_predict(mdl, param = start(mdl), data = dat, mgen = -1))
  expect_error(fd_predict(mdl, param = start(mdl), data = NULL, mgen = NULL))
  expect_error(fd_predict("mdl", param = start(mdl), data = dat))
})

test_that("various FMM and proliferation models work", {
  cstr <- ~ `#macr` + CC$`#delta_1100` + CC$`#mm_xyxz`
  dat <- subset(cutoff(FdSTyphimuriumWTC57), Individual == "140508.WT.C57")
  for (fmm in c("gaussian", "af", "af_bp")) {
    for (pro in c("cyton", "branching")) {
      mdl <- fd_model(dat, fmm=fmm, proliferation=pro,
                      constraints = cstr)
      expect_silent(suppressWarnings(
        do.call(fd_minuslogl(mdl, dat, verbose=FALSE), as.list(start(mdl)))))
    }
  }

  # log10_scale
  mdl <- fd_model(data=dat,
                  proliferation=fd_proliferation_branching(
                    mgen = 10L, log10_scale = TRUE),
                  constraints = cstr)
  expect_equal(
    start(mdl),
    structure(
      c(-0.698970004336019, -0.301029995663981, -0.309803919971486, 5, 5, 5),
      .Names = c(
        "pro.One.res0", "pro.One.p", "pro.One.p0",
        "pro.One.g.mm", "pro.One.g0.mm", "pro.One.f.mm")))
  expect_equal(relist(start(mdl), mdl), start(mdl, free=FALSE))
})

test_that("crazy parameters get turned around", {
  expect_error(get_fmm("dsfsdfsdf"))
  for (fmm in c("gaussian", "af", "af_bp")) {
    expect_message(get_fmm("gaussian"))
    expect_error(
      get_fmm("gaussian", structure(data.frame(), fmm="ldsfsdfds")))
    expect_error(
      get_fmm("gaussian", structure(data.frame(), fmm=list(abc = "lqsd"))))
  }

  expect_error(get_proliferation("dsfsdfsdf"))
  for (pro in c("cyton", "branching")) {
    expect_message(get_proliferation(pro))
    expect_error(get_proliferation(
      pro,
      structure(data.frame(), pro="ldsfsdfds")))
    expect_silent(get_proliferation(
      pro,
      structure(
        data.frame(Category="One"),
        pro=list(abc = "lqsd"))))
  }
})

test_that("decay for Ns_lost works", {
  load(system.file("extdata", "FdSTyphimuriumMice.rda",
                   package="fluodilution"))
  dat <- fd_transform(cutoff(FdSTyphimuriumMice), N = "log1p")

  mdl <- fd_model(
    data = dat,
    constraints = cstrlist(
      constraints = ~ `#invivo` + pro:all:(h):{delta <- 1; ss <- 0.5},
      start = ~ pro:all:(h):{mm <- 5},
      lower = ~ pro:all:(h):{mm <- 0},
      upper = ~ pro:all:(h):{mm <- 20}))
  ans <- do.call(
    fd_minuslogl(mdl, dat, stop_boundary=TRUE, mgen=4),
    as.list(start(mdl)))

  mdl2 <- fd_model(data = dat, constraints = ~ `#invivo`)
  ans2 <- do.call(
    fd_minuslogl(mdl2, dat, stop_boundary=TRUE, mgen=4),
    as.list(start(mdl2)))

  expect_true(ans != ans2)
})
