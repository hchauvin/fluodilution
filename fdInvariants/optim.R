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

# Test optimization on simulated data, to prove that the objective function
# can be used to get back to the biological constants

library(fluodilution)
library(assertthat)

fmm_params <- list(
  m0 = 10,
  sd0 = 0.1
)

coefs <- c(pro.One.f.mm = 1,
           pro.One.f0.mm = 3,
           pro.One.f0.delta = 0.3)

times <- seq(0, 12, by = 2)

attach(FdCommonConstraints)
constraints <- fluodilution::cstrlist(
  # nolint start
  ~ `#nodeath` + `#noss` +
      {fmm$c0 <- 0} + pro:all:{res0 <- 0; res <- 0; f$delta <- 0},
  start = ~ pro:all:(f/f0):{mm <- 1},
  lower = ~ pro:all:(f/f0):{mm <- 0.5},
  upper = ~ pro:all:(f/f0):{mm <- 5}
  # nolint end
)

mdl <- fd_model(fmm = "gaussian",
                proliferation = "branching",
                constraints = constraints,
                boxed = T)

dat <- fd_simulate(coefs,
                   times = times,
                   mgen = 10L,
                   breaks = 50L,
                   range = c(-100.0, exp(fmm_params$m0) * 2.0),
                   select = c("hists"),
                   model = mdl)

fit <- nls(
  fd_formula("mdl", mgen = 10L),
  data = dat,
  start = start(mdl),
  lower = lower(mdl),
  upper = upper(mdl),
  algorithm = "port"
)

rbind(
  start = start(mdl),
  fit = coef(fit),
  actual = coefs
)

assert_that(max(abs(coef(fit) - coefs)) < 1e-12)