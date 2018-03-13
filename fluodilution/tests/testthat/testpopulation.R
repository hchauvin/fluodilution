context("population")

source(system.file("contrib", "nlsSA.R", package="fluodilution"))

# nolint start

# Base constraints for Salmonella Typhimurium experiments
`#macr` <<- ~ FdCommonConstraints$`#noss` + pro:all:{res <- 0} + fmm:{c0 <- 0}
`#invivo` <<- ~ FdCommonConstraints$`#noss` + ~ pro:all:{res <- 0}

# nolint end

CC <<- FdCommonConstraints

test_that("'fd_data' & co work and put up the right warnings", {
  expect_error(fd_data("ldsflksdfds"))
  data(FdSTyphimuriumWTC57)
  dat <- FdSTyphimuriumWTC57
  expect_warning(subset(dat, a < 100),
                 "Some histograms do not have proportions that sum to 1.")
  expect_silent(fd_data(data.frame()))

  # Numeric error
  expect_warning(fd_data(transform(dat, Category = factor(as.numeric(Category)))))
  expect_warning(fd_data(transform(dat, Individual = factor(as.numeric(Individual)))))
  expect_warning(fd_data(transform(dat, Timepoint = factor(as.numeric(Timepoint)))))

  # Inoculum errors
  dat2 <- dat
  attr(dat2, "fmm")$m0 <- attr(dat2, "fmm")$m0[1:5]
  expect_warning(fd_data(dat2))

  # Also test '['
  expect_equal(class(dat[, c("Individual", "Timepoint")]), "data.frame")
  expect_is(dat[as.numeric(dat$Individual) == 1, ], "fd_data")

  # rbind
  b1 <- subset(FdSTyphimuriumWTC57, Individual == "160408.WT.C57")
  b2 <- subset(FdSTyphimuriumWTC57, Individual == "290708.WT.C57")
  bound <- rbind(b1, b2)
  expect_equal(class(bound),
               c("fd_data", "nffGroupedData", "nfGroupedData", "groupedData",
               "data.frame"))
  expect_equal(
    attr(bound, "fmm"),
    structure(
      list(
        sd0 = structure(
          c(0.848338664172271, 0.464521581665721),
          .Names = c("inoc_160408.WT.C57", "inoc_290708.WT.C57")),
        m0 = structure(
          c(6.97558484327899, 6.37922906123026),
          .Names = c("inoc_160408.WT.C57", "inoc_290708.WT.C57"))),
      .Names = c("sd0", "m0")))
  expect_equal(levels(bound$Individual), c("160408.WT.C57", "290708.WT.C57"))
  expect_equal(NROW(bound), NROW(b1) + NROW(b2))

  # cutoff
  expect_equal(cutoff(FdSTyphimuriumWTC57),
               cutoff(FdSTyphimuriumWTC57, attr(FdSTyphimuriumWTC57, "cutoff")))
  expect_silent(cutoff(FdSTyphimuriumWTC57, 2^7))
})

test_that("'fd_residuals' and 'fd_predict' works (and also test 'pretty' at the same time)", {
  dat <- subset(cutoff(FdSTyphimuriumWTC57), Individual == "140508.WT.C57")
  print(summary(dat))

  cstr <- ~ `#macr` + CC$`#delta_1100` + CC$`#mm_xyxz`
  mdl <<- fd_model(dat, fmm="gaussian", proliferation="branching",
                   constraints = cstr)
  print(mdl)
  print(summary(mdl))
  print(mdl$pro)
  print(summary(mdl$pro))

  ans <- fd_nls("mdl", dat, mgen=16, loop=0, trace=TRUE,
                nlsfun=nlsSA,
                control=list(maxit=1))
  print(ans)
  print(summary(ans))

  # Check that 'fd_residuals' and 'predict' give coherent results
  res <- fd_residuals(dat, ans)
  pred <- predict(ans, dat)
  expect_equal(dat$y - pred, as.vector(res$resid))
  expect_equal(pred, as.vector(res$fitted))

  # Try different ways of calling 'fd_residuals'
  expect_equal(levels(fd_residuals(dat, ans)$algo), "mod_1")
  expect_equal(levels(fd_residuals(dat, `1` = ans)$algo), "1")
  expect_equal(levels(fd_residuals(dat, `1` = ans,
                                   .list = list(`2` = ans, `3` = ans))$algo),
               c("1", "2", "3"))

  invisible(fd_residuals(
    list(`1` = dat, `2` = dat),
    `1` = ans, `2` = ans))
  expect_error(fd_residuals(
    list(dat, dat),
    `1` = ans, `2` = ans))

  # With nlsList
  dat2 <- subset(cutoff(FdSTyphimuriumWTC57),
                 Individual %in% c("140508.WT.C57", "160408.WT.C57"))
  source(system.file("contrib", "fitwrappers.R", package="fluodilution"))
  ans2 <- fd_nlsList(ans, dat2, mgen=16, trace=TRUE,
                     nlsListfun=nlsListSA,
                     control=list(maxit=1))
  print(summary(ans2))

  expect_equal(names(ans2), c("140508.WT.C57", "160408.WT.C57"))
  res2 <- fd_residuals(dat2, ans2)
  expect_equal(levels(res2$Individual), c("140508.WT.C57", "160408.WT.C57"))

  pred2 <- predict(ans2, dat2)
  expect_equal(dat2$y - pred2, as.vector(res2$resid))
  expect_equal(pred2, as.vector(res2$fitted))

  # With nlme
  expect_warning(ans3 <- fd_nlme(ans, dat2,
                  pdDiag(pro.One.p + pro.One.p0 + pro.One.g.mm ~ 1),
                  control = nlmeControl(maxIter = 1, returnObject=TRUE)))
  print(ans3)
  print(summary(ans3))
  res3 <- fd_residuals(dat2, ans3)

  pred3 <- predict(ans3, dat2)
})

test_that("'fd_expand' works", {
    dat <- subset(cutoff(FdSTyphimuriumWTC57), Individual == "140508.WT.C57")
    expect_silent(fd_expand(dat, seq = NULL, by = NULL))
    expect_silent(fd_expand(dat, seq = NULL, by = "Individual"))
    expect_error(fd_expand(dat, seq = NULL, by = "dsfsdf"))
    expect_silent(fd_expand(dat, seq = NULL, separate = FALSE))
})

test_that("'plot.fd_data' works", {
  data(FdSTyphimuriumWTC57)

  # Plot of only one individual
  dat <- subset(cutoff(FdSTyphimuriumWTC57), Individual == "140508.WT.C57")
  expect_silent(plot(dat))
  expect_silent(plot(dat, c("hist", "N", "range", "balancing", "cutoff")))

  # Plot of multiple individuals
  datall <- cutoff(FdSTyphimuriumWTC57)
  expect_silent(plot(datall))
  expect_silent(plot(datall, c("hist", "N", "range", "balancing", "cutoff")))
  # Note: cannot test "coverage" as we don't have the "count" attribute

  # In vivo
  load(system.file("extdata", "FdSTyphimuriumMice.rda",
                    package="fluodilution"))
  datall_invivo <- fd_transform(cutoff(FdSTyphimuriumMice), N = "log1p")
  plot(datall_invivo, c("hist", "N", "range", "balancing", "cutoff"))
})

test_that("'fd_transform' works", {
  load(system.file("extdata", "FdSTyphimuriumMice.rda",
                    package="fluodilution"))
  dat <- cutoff(FdSTyphimuriumMice)
  dat_tr <- fd_transform(dat, hist = "log", N = "log1p")

  # Check that the previous transformation is unwind
  dat2 <- fd_transform(dat, hist = "identity", N = "identity")
  dat_tr2 <- fd_transform(dat_tr, hist = "identity", N = "identity")
  expect_equal(dat2, dat_tr2)

  dat_tr3 <- fd_transform(dat_tr, hist = "log1p", N = "log1p")
  dat_tr4 <- fd_transform(dat, hist = "log1p", N = "log1p")
  expect_equal(dat_tr3, dat_tr4)
})
