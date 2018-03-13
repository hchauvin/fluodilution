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

params <- list(
  mgen = 10L,
  m0 = 10,
  sd0 = 0.1,
  ftor = 5.0,
  sdaf = 100.0
)

data(FdClearPeaks, package = "fluodilution")

dat <- FdClearPeaks

fmm_af_bp <- fd_fmm_af_bp(m0 = params$m0, sd0 = params$sd0)
res_af_bp <- fmm_af_bp$diff_cdf(
  psi = list(ftor = params$ftor, sdaf = params$sdaf),
  a = unique(dat$a),
  b = unique(dat$b),
  gen = 0:params$mgen)

res_agent_af_bp <- agent_af_bp_hist(unique(dat$a), unique(dat$b), params)

df <- rbind(
  melt(res_af_bp, varnames = c("Idx", "Generation")) %>%
    transform(Algo = "Model: AF+BP"),
  melt(res_agent_af_bp, varnames = c("Idx", "Generation")) %>%
    transform(Algo = "Simulation", value = value)
) %>% transform(
  x = unique((dat$a + dat$b) / 2.0)[Idx]
)

print(ggplot(df, aes(x = x, y = value, color = Algo))+
  geom_line()+
  facet_wrap(~ Generation, scales = "free")+
  scale_x_log10())

assert_that(
  mean(abs(res_af_bp - res_agent_af_bp)
    [res_af_bp > 0 & res_agent_af_bp > 0]) < 0.001)
assert_that(max(abs(res_af_bp - res_agent_af_bp)) < 0.005)