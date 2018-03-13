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

# Test the branching model against an agent-based model for the in vitro case
# (that is, without death, and within only one compartment).

library(plyr)
library(fluodilution)
library(assertthat)
library(magrittr)
library(reshape2)
library(ggplot2)

source(system.file("contrib", "agent.R", package="fluodilution"))

attach(FdCommonConstraints)
invitro_constraints <- fluodilution::cstrlist(
  # nolint start
  ~ `#nodeath` + `#noss` +
  {fmm$c0 <- 0} + pro:all:{res0 <- 0; res <- 0; f$delta <- 0},
  start = ~ pro:all:(f/f0):{mm <- 1},
  lower = ~ pro:all:(f/f0):{mm <- 0.5},
  upper = ~ pro:all:(f/f0):{mm <- 2}
  # nolint end
)

data(FdClearPeaks, package = "fluodilution")

dat <- FdClearPeaks

mdl_branching <- fd_model(dat,
  fmm="gaussian",
  proliferation="branching",
  constraints = invitro_constraints)

params <- structure(
  relist(start(mdl_branching), mdl_branching),
  model = mdl_branching)

mgen <- 10
times <- seq(0, 10, by = 2)
tree_bank <- agent_tree_bank(
  n =1000,
  params = params,
  initial = rbind(c(1, rep(0, mgen))),
  times = times)

N0 <- 3000
tt <- sample_agent_trees(tree_bank, Inf, N0)

tt2 <- mdl_branching$pro$model(params$pro, times=times, mgen=mgen)

df_tt2 <- transform(
  melt(tt2$live_pop, varnames = c("Idx", "Generation")),
  Category = L1,
  Time = times[Idx]
)

print(ggplot(subset(tt, !Lost),
             aes(x = factor(Generation), y = ..count.. / .GlobalEnv$N0))+
  geom_bar()+
  geom_point(data = df_tt2, aes(x = Generation, y = value))+
  facet_wrap(~ Time))

merged <- merge(
  ddply(
    subset(tt, !Lost),
    c("Time", "Generation", "Category"),
    summarize,
    value = length(Generation) / N0),
  df_tt2[, c("Time", "Generation", "Category", "value")] %>%
    transform(Generation = Generation - 1L),
  by = c("Time", "Generation", "Category"),
  all = TRUE)
merged$value.x[is.na(merged$value.x)] <- 0
merged <- transform(merged, diff = abs(value.x - value.y))

print(merged)

print(max(merged$diff))
assert_that(max(merged$diff) < 0.05)