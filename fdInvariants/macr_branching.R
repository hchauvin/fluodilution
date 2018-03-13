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

# Test the branching model against an agent-based model for the macrophage
# case (that is, with death, and within only one compartment).

library(plyr)
library(fluodilution)
library(assertthat)
library(magrittr)
library(reshape2)
library(ggplot2)

source(system.file("contrib", "agent.R", package="fluodilution"))

attach(FdCommonConstraints)
macr_constraints <- fluodilution::cstrlist(
  # nolint start
  ~ `#noss` +
  {fmm$c0 <- 0} + pro:all:{res0 <- 0; res <- 0; f$delta <- 0},
  start = ~ pro:all:(f/f0):{mm <- 1},
  lower = ~ pro:all:(f/f0):{mm <- 0.5},
  upper = ~ pro:all:(f/f0):{mm <- 2}
  # nolint end
)

data(FdClearPeaks, package = "fluodilution")

dat <- FdClearPeaks

mgen <- 10
times <- seq(0, 10, by = 2)

mdl_branching <- fd_model(dat,
  fmm="gaussian",
  proliferation="branching",
  constraints = macr_constraints)

params <- structure(
  relist(start(mdl_branching), mdl_branching),
  model = mdl_branching)

tree_bank <- agent_tree_bank(
  n = 1000,
  params = params,
  initial = rbind(c(1, rep(0, mgen))),
  times = times)

Ni <- 3000
agent_based <- sample_agent_trees(tree_bank, Inf, Ni)

branching_based <- mdl_branching$pro$model(params$pro, times=times, mgen=mgen)

df_props_branching_based <- transform(
  melt(branching_based$live_pop, varnames = c("Idx", "Generation")),
  Category = L1,
  Time = times[Idx]
)

print(ggplot(subset(agent_based, !Lost),
             aes(x = factor(Generation), y = ..count.. / .GlobalEnv$Ni))+
  geom_bar()+
  geom_point(data = df_props_branching_based, aes(x = Generation, y = value))+
  facet_grid(Category ~ Time)+
  ggtitle("Live pop"))

df_props_branching_based_lost <- transform(
  melt(branching_based$lost_pop, varnames = c("Idx", "Generation")),
  Category = L1,
  Time = times[Idx]
)

print(ggplot(subset(agent_based, Lost),
             aes(x = factor(Generation), y = ..count.. / .GlobalEnv$Ni))+
    geom_bar()+
    geom_point(data = df_props_branching_based_lost,
               aes(x = Generation, y = value))+
    facet_grid(Category ~ Time)+
    ggtitle("Lost pop"))

df_Ns_branching_based <- rbind(
  melt(branching_based$Ns, varnames = c("Idx", "Category")) %>%
    transform(Lost = FALSE),
  melt(branching_based$Ns_lost, varnames = c("Idx", "Category")) %>%
    transform(Lost = TRUE)
) %>% transform(Time = times[Idx])

df_Ns_agent_based <- ddply(
  tree_bank$trees,
  c("Time", "Category", "Lost"),
  summarize,
  y = length(Time) / tree_bank$n)

print(ggplot(df_Ns_branching_based, aes(x = Time, y = value))+
    geom_line()+
    geom_point(data = df_Ns_agent_based, aes(x = Time, y = y))+
    facet_wrap(Lost ~ Category)+
    ggtitle("Ns"))

merged <- merge(
  rbind(
    df_props_branching_based %>%
      .[, c("Time", "Generation", "Category", "value")] %>%
      transform(Lost = FALSE, Generation = Generation - 1L),
    df_props_branching_based_lost %>%
      .[, c("Time", "Generation", "Category", "value")] %>%
      transform(Lost = TRUE, Generation = Generation - 1L)
  ),
  ddply(
    agent_based,
    c("Time", "Generation", "Category", "Lost"),
    summarize,
    value = length(Generation) / Ni),
  by = c("Time", "Generation", "Category", "Lost"),
  all = TRUE)
merged$value.x[is.na(merged$value.x)] <- 0
merged$value.y[is.na(merged$value.y)] <- 0
merged <- transform(merged, diff = abs(value.x - value.y))

print(merged)

# About this threshold: it is a model, after all
assert_that(max(merged$diff) < 0.15)

assert_that(mean(merged$diff[merged$diff > 0]) < 0.02)
