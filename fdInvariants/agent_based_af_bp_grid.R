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

# Validate the AF+BP FMM with an agent-based simulation.

library(fluodilution)
library(ggplot2)
library(assertthat)
library(magrittr)
library(reshape2)
library(parallel)

agent_af_bp <- function (mgen, m0, sd0, ftor, sdaf, Ni = 100000) {
  ans <- list()

  get_ans <- function (nmol) {
    list(
      nmol = nmol,
      fluo = nmol * exp(ftor) + rnorm(length(nmol), 0, sdaf)
    )
  }

  ans[[1L]] <- get_ans(nmol = round(rlnorm(Ni, m0 - ftor, sd0)))

  for (m in 2L:(mgen + 1)) {
    ans[[m]] <- get_ans(
      nmol = sapply(ans[[m-1]]$nmol, function (n) rbinom(1, n, 0.5)))
  }

  return(ans)
}

agent_af_bp_hist <- function (a, b, params) {
  min_break <- min(a)
  max_break <- max(b)

  breaks <- unique(c(a, b))

  rand <- do.call(agent_af_bp, params)

  sapply(rand, function (cur) {
    hist(cur$fluo %>% .[. >= min_break & . <= max_break],
         breaks = breaks, plot = FALSE) %>%
        .$density * (b - a)
  })
}

kl_div <- function (p, q) sapply(1:NCOL(p), function (i) {
  pp <- p[, i]
  qq <- q[, i]
  sum((pp * log(pp / qq))[qq > 0 & pp > 0])
})

js_div <- function (p, q) 0.5 * (kl_div(p, q) + kl_div(q, p))

process <- function (params, plot=TRUE) {
  a <- unique(FdClearPeaks$a)
  b <- unique(FdClearPeaks$b)

  fmm_af_bp <- fd_fmm_af_bp(m0 = params$m0, sd0 = params$sd0)
  res_af_bp <- fmm_af_bp$diff_cdf(
    psi = list(ftor = params$ftor, sdaf = params$sdaf),
    a = a,
    b = b,
    gen = 0:params$mgen
  ) %>% apply(2, function (col) col / sum(col))

  res_agent_af_bp <- agent_af_bp_hist(a, b, params)

  df <- rbind(
    melt(res_af_bp, varnames = c("Idx", "Generation")) %>%
      transform(Algo = "Model: AF+BP"),
    melt(res_agent_af_bp, varnames = c("Idx", "Generation")) %>%
      transform(Algo = "Simulation", value = value)
  ) %>% transform(
    x = unique((a + b) / 2.0)[Idx]
  )

  if (plot) {
    print(ggplot(df, aes(x = x, y = value, color = Algo))+
      facet_wrap(~ Generation, scales = "free")+
      geom_line()+
      scale_x_log10())
  }

  rel_diff <- function (fn) {
    sapply(0:params$mgen, function (i) {
      model <- fn(res_af_bp[, i+1])
      agent <- fn(res_agent_af_bp[, i+1])
      abs((model - agent) / agent)
    })
  }

  x <- (asinh(a) + asinh(b)) / 2.0
  moment_1_diff <- rel_diff(function (d) sum(d * x))
  moment_2_diff <- rel_diff(function (d) sum(d * (x - sum(d * x)) ** 2))
  moment_3_diff <- rel_diff(function (d) sum(d * (x - sum(d * x)) ** 3))
  moment_4_diff <- rel_diff(function (d) sum(d * (x - sum(d * x)) ** 4))

  return (c(
    mean_diff = mean(abs(res_af_bp - res_agent_af_bp)
      [res_af_bp > 0 & res_agent_af_bp > 0]),
    max_diff = max(abs(res_af_bp - res_agent_af_bp)),
    js_div = mean(js_div(res_af_bp, res_agent_af_bp)),
    m1 = mean(moment_1_diff),
    m2 = mean(moment_2_diff),
    m3 = mean(moment_3_diff),
    m4 = mean(moment_4_diff)))
}

params_grid <- expand.grid(
  mgen = 10L,
  m0 = seq(7, 11, length.out = 2),
  sd0 = seq(0.20, 0.50, length.out = 2),
  ftor = seq(0.0, 3.0, length.out = 2),
  sdaf = seq(50.0, 400.0, length.out = 2)
)

result <- mclapply(
  1:NROW(params_grid),
  function (i) process(as.list(params_grid[i,]), plot = FALSE),
  mc.cores = 2)

df_result <- cbind(
  params_grid,
  transform(
    as.data.frame(do.call(rbind, result)),
    ok = mean_diff < 0.0005 & max_diff < 0.25 & js_div < 0.9))

print(df_result)

assert_that(all(df_result$ok))

