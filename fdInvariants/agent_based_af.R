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

# Validate the AF FMM with an agent-based simulation.

library(fluodilution)
library(ggplot2)
library(assertthat)
library(magrittr)
library(reshape2)

agent_af <- function (mgen, m0, sd0, sdaf, Ni = 10000) {
  ans <- list()

  add_af <- function (fluo) {
    list(
      fluo = fluo,
      fluo_af = fluo + rnorm(length(fluo), 0, sdaf)
    )
  }

  ans[[1L]] <- add_af(fluo = rlnorm(Ni, m0, sd0))

  for (m in 2L:(mgen + 1)) {
    ans[[m]] <- add_af(fluo = ans[[m-1]]$fluo / 2.0)
  }

  return(ans)
}

agent_af_hist <- function (a, b, params) {
  min_break <- min(a)
  max_break <- max(b)

  breaks <- unique(c(a, b))

  rand <- do.call(agent_af, params)

  sapply(rand, function (cur) {
    hist(cur$fluo_af %>% .[. >= min_break & . <= max_break],
         breaks = breaks, plot = FALSE) %>%
      .$density * (b - a)
  })
}

params <- list(
  mgen = 10L,
  m0 = 10,
  sd0 = 0.1,
  sdaf = 100.0
)

data(FdClearPeaks, package = "fluodilution")

dat <- FdClearPeaks

fmm_af <- fd_fmm_af(m0 = params$m0, sd0 = params$sd0)
res_af <- fmm_af$diff_cdf(
  psi = list(sdaf = params$sdaf),
  a = unique(dat$a),
  b = unique(dat$b),
  gen = 0:params$mgen)

res_agent_af <- agent_af_hist(unique(dat$a), unique(dat$b), params)

df <- rbind(
  melt(res_af, varnames = c("Idx", "Generation")) %>%
    transform(Algo = "Model: AF+BP"),
  melt(res_agent_af, varnames = c("Idx", "Generation")) %>%
    transform(Algo = "Simulation", value = value)
) %>% transform(
  x = unique((dat$a + dat$b) / 2.0)[Idx]
)

print(ggplot(df, aes(x = x, y = value, color = Algo))+
  geom_line()+
  facet_wrap(~ Generation, scales = "free")+
  scale_x_log10())

assert_that(
  mean(abs(res_af - res_agent_af)[res_af > 0 & res_agent_af > 0]) < 0.001)
assert_that(max(abs(res_af - res_agent_af)) < 0.01)