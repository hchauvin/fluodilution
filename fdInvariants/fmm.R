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

# Test some invariants pertaining to FMM

library(fluodilution)
library(assertthat)
library(magrittr)

data(FdClearPeaks, package = "fluodilution")

dat <- FdClearPeaks

PLOT <- TRUE

params <- list(
  m0 = 10,
  sd0 = 0.1,
  psi = list(
    i = 1,
    ftor = 5.0,
    sdaf = 100.0
  )
)

fmm_gaussian <- fd_fmm_gaussian(m0 = params$m0, sd0 = params$sd0)
fmm_af <- fd_fmm_af(m0 = params$m0, sd0 = params$sd0)
fmm_af_bp <- fd_fmm_af_bp(m0 = params$m0, sd0 = params$sd0)

res_af <- fmm_af$diff_cdf(
  psi = params$psi, a = dat$a, b = dat$b, gen = 1:10)
res_gaussian <- fmm_gaussian$diff_cdf(
  psi = params$psi, a = dat$a, b = dat$b, gen = 1:10)
res_af_bp <- fmm_af_bp$diff_cdf(
  psi = params$psi, a = dat$a, b = dat$b, gen = 1:10)

if (PLOT) {
  library(reshape2)

  molten_res <- transform(
    rbind(
      melt(res_gaussian, varnames= c("index", "gen")) %>%
        transform(fmm = "gaussian"),
      melt(res_af, varnames= c("index", "gen")) %>%
        transform(fmm = "af"),
      melt(res_af_bp, varnames= c("index", "gen")) %>%
        transform(fmm = "af_bp")
    ),
    x = ((dat$a + dat$b) / 2.0)[index],
    density = value / (dat$b - dat$a)[index])

  print(ggplot(molten_res, aes(x = x, color = fmm, y = density))+
    facet_wrap(~ gen)+
    geom_line()+scale_x_log10())
}

# Kullback–Leibler divergence on the "discretized" probability measure,
# by generation
kl_div <- function (p, q) sapply(1:NCOL(p), function (i) {
  pp <- p[, i]
  qq <- q[, i]
  sum((pp * log(pp / qq))[qq > 0 & pp > 0])
})

# Jensen-Shannon divergence on the "discretized" probability measure,
# by generation
js_div <- function (p, q) 0.5 * (kl_div(p, q) + kl_div(q, p))

invariant <- function (description, expr) {
  tryCatch(expr, error = function (e) {
    stop("failed invariant: ", description, ": ", e)
  })
}

js_div_af_bp_for_ftor <- function (ftor) {
  psi <- params$psi
  psi$ftor <- ftor
  res_af_bp <- fmm_af_bp$diff_cdf(psi = psi, a = dat$a, b = dat$b, gen = 1:10)
  js_div(res_af, res_af_bp)
}

invariant("AF+BP: the greater the ftor, the greater the average divergence from AF only", {
  m1 <- mean(js_div_af_bp_for_ftor(4))
  m2 <- mean(js_div_af_bp_for_ftor(5))
  m3 <- mean(js_div_af_bp_for_ftor(6))

  assert_that(m1 < m2)
  assert_that(m2 < m3)
})

invariant("AF+BP: the divergence is always greater for lower generations", {
  # NOTE: we reduce here because the divergence is only computed on the intersection
  # of the supports of the probability measures, so with this "pseudo" metric
  # the invariant does not always hold.

  assert_that(all(diff(js_div_af_bp_for_ftor(5))[5L:9L] < 0))
})

num_err <- function (a, b) {
  mean(abs(a - b)) / mean(a + b)
}

psi <- params$psi
psi$ftor <- 1
res_af_bp_low_ftor <- fmm_af_bp$diff_cdf(
  psi = psi, a = dat$a, b = dat$b, gen = 1:10)

invariant("AF+BP: for very low ftor, the result is indistinguishable from autofluorescence only", {
  err <- num_err(res_af_bp_low_ftor, res_af)

  cat("err", err)

  assert_that(err < 0.05)
})

js_div_af_for_sdaf <- function (sdaf) {
  psi <- params$psi
  psi$sdaf <- sdaf
  res_af <- fmm_af$diff_cdf(psi = psi, a = dat$a, b = dat$b, gen = 1:10)
  js_div(res_af, res_gaussian)
}

invariant("AF: the greater the autofluorescence, the greater the average divergence from the Gaussian FMM", {
  m1 <- mean(js_div_af_for_sdaf(100.0))
  m2 <- mean(js_div_af_for_sdaf(200.0))
  m3 <- mean(js_div_af_for_sdaf(300.0))

  cat("m1", m1, "m2", m2, "m3", m3)

  assert_that(m1 < m2)
  assert_that(m2 < m3)
})

js_div_af_bp_for_sdaf <- function (sdaf) {
  psi <- params$psi
  psi$sdaf <- sdaf
  res_af_bp <- fmm_af_bp$diff_cdf(psi = psi, a = dat$a, b = dat$b, gen = 1:10)
  js_div(res_af, res_af_bp)
}

invariant("AF: the greater the autofluorescence, the greater the average divergence from the Gaussian FMM", {
  m1 <- mean(js_div_af_bp_for_sdaf(100.0))
  m2 <- mean(js_div_af_bp_for_sdaf(200.0))
  m3 <- mean(js_div_af_bp_for_sdaf(300.0))

  cat("m1", m1, "m2", m2, "m3", m3)

  assert_that(m1 < m2)
  assert_that(m2 < m3)
})

invariant("AF and AF+BP: the result is indistinguishable from the Gaussian FMM at the start", {
  err_af <- num_err(res_af[, 1], res_gaussian[, 1])
  err_af_bp <- num_err(res_af_bp_low_ftor[, 1], res_gaussian[, 1])

  cat("err_af", err_af, "err_af_bp", err_af_bp)

  assert_that(err_af < 0.005)
  assert_that(err_af_bp < 0.05)
})