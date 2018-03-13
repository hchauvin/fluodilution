context("modsel")

source(system.file("contrib", "nlsSA.R", package="fluodilution"))

# nolint start

# Base constraints for Salmonella Typhimurium experiments
`#macr` <<- ~ FdCommonConstraints$`#noss` + pro:all:{res <- 0} + fmm:{c0 <- 0}
`#invivo` <<- ~ FdCommonConstraints$`#noss` + ~ pro:all:{res <- 0}

# nolint end

CC <<- FdCommonConstraints

test_that("modsel works", {
  cstr <- list(~ `#macr` + CC$`#delta_1100` + CC$`#mm_xyxz`,
                ~ `#macr` + CC$`#delta_1100` + CC$`#mm_xxxx`)
  data(FdSTyphimuriumWTC57)
  dat <- cutoff(FdSTyphimuriumWTC57)
  comb <- fd_comb(
    cstr,
    data=dat,
    continue=TRUE,
    mc.cores=1,
    file=NULL,
    nlsfun=nlsSA,
    control=list(maxit=1, simple.function=FALSE))
  expect_silent(invisible(fd_aictab(data=dat, .list=comb)))
  expect_silent(invisible(fd_aictab(data=dat, .list=comb, second.ord=TRUE)))
  ans <- fd_aictab(data=dat, comb1=comb[[1L]], comb2=comb[[2L]])
  expect_silent(fd_weights(ans))
  avg <- fd_modavg(ans, vcov. = sandwich::sandwich)
  expect_error(relisted_fit(avg))
  rel <- relisted_fit(comb[[1L]])
  expect_silent(vcov(rel))
  expect_silent(vcov(avg))
  expect_silent(weights(avg))
  expect_equal(coef(best(ans, Q == "")), coef(comb[[2L]]))
  expect_equal(coef(best(ans, pos=2)), coef(comb[[1L]]))
  expect_silent(best(ans, verbose=FALSE))
})
