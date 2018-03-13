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

# Test the Cyton model against an agent-based model for the macrophage case (that is,
# with death, and within only one compartment).  The Cyton model is a special case of a
# branching process with the probability of dividing and rates adapted
# (see Chauvin, Watson et al., forthcoming).

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
  ~ `#noss` + `#delta_1111` +
  {fmm$c0 <- 0} + pro:all:{res0 <- 0; res <- 0},
  start = ~ pro:all:(f/f0):{mm <- 1; delta <- 1} +
            pro:all:(g/g0):{mm <- 5; delta <- 1},
  lower = ~ pro:all:(f/f0):{mm <- 0.5},
  upper = ~ pro:all:(f/f0):{mm <- 2}
  # nolint end
)

data(FdClearPeaks, package = "fluodilution")

dat <- FdClearPeaks

times <- seq(0, 10, by = 2)

mdl_cyton <- fd_model(dat,
  fmm="gaussian",
  proliferation="cyton",
  constraints = macr_constraints)
mdl_branching <- update(mdl_cyton, data=dat, proliferation="branching")

params <- structure(relist(start(mdl_cyton), mdl_cyton), model = mdl_cyton)

params_branching <- params
params_branching$pro$One$p0 <-
  (1 / params_branching$pro$One$f0$mm) /
    (1 / params_branching$pro$One$f0$mm + 1 / params_branching$pro$One$g0$mm)
params_branching$pro$One$p <-
  (1 / params_branching$pro$One$f$mm) /
    (1 / params_branching$pro$One$f$mm + 1 / params_branching$pro$One$g$mm)
common_rate0 <-
  1 / params_branching$pro$One$f0$mm +
  1 / params_branching$pro$One$g0$mm
common_rate <-
  1 / params_branching$pro$One$f$mm +
  1 / params_branching$pro$One$g$mm
params_branching$pro$One$f0$mm <- 1 / common_rate0
params_branching$pro$One$g0$mm <- 1 / common_rate0
params_branching$pro$One$f$mm <- 1 / common_rate
params_branching$pro$One$g$mm <- 1 / common_rate

mgen <- 30

tree_bank <- agent_tree_bank(
  n = 1000,
  params = params_branching,
  initial = rbind(c(1, rep(0, mgen))),
  times = times,
  cyton = TRUE)

Ni <- 1000
agent_based <- sample_agent_trees(tree_bank, Inf, Ni)

branching_based <- mdl_branching$pro$model(
  params_branching$pro, times=times, mgen=mgen)
cyton_based <- mdl_cyton$pro$model(
  params$pro, times=times, mgen=mgen)

df_props_cyton_based <- rbind(
  melt(branching_based$live_pop, varnames = c("Idx", "Generation")) %>%
    transform(Model = "branching"),
  melt(cyton_based$live_pop, varnames = c("Idx", "Generation")) %>%
    transform(Model = "cyton")
) %>% transform(
  Category = L1,
  Time = times[Idx]
)

print(ggplot(subset(agent_based, !Lost),
             aes(x = factor(Generation), y = ..count.. / .GlobalEnv$Ni))+
  geom_bar()+
  geom_point(data = df_props_cyton_based,
             aes(x = Generation, y = value, color = Model))+
  facet_grid(Category ~ Time)+
  ggtitle("Live pop"))

df_props_cyton_based_lost <- rbind(
  melt(branching_based$lost_pop, varnames = c("Idx", "Generation")) %>%
    transform(Model = "branching"),
  melt(cyton_based$lost_pop, varnames = c("Idx", "Generation")) %>%
    transform(Model = "cyton")
) %>% transform(
  Category = L1,
  Time = times[Idx]
)

print(ggplot(subset(agent_based, Lost),
             aes(x = factor(Generation), y = ..count.. / .GlobalEnv$Ni))+
    geom_bar()+
    geom_point(data = df_props_cyton_based_lost,
               aes(x = Generation, y = value, color = Model))+
    facet_grid(Category ~ Time)+
    ggtitle("Lost pop"))

df_Ns_cyton_based <- rbind(
  melt(branching_based$Ns, varnames = c("Idx", "Category")) %>%
    transform(Lost = FALSE, Model = "branching"),
  melt(branching_based$Ns_lost, varnames = c("Idx", "Category")) %>%
    transform(Lost = TRUE, Model = "branching"),
  melt(cyton_based$Ns, varnames = c("Idx", "Category")) %>%
    transform(Lost = FALSE, Model = "cyton"),
  melt(cyton_based$Ns_lost, varnames = c("Idx", "Category")) %>%
    transform(Lost = TRUE, Model = "cyton")
) %>% transform(Time = times[Idx])

df_Ns_agent_based <- ddply(
  tree_bank$trees,
  c("Time", "Category", "Lost"),
  summarize,
  y = length(Time) / tree_bank$n)

print(ggplot(df_Ns_cyton_based, aes(x = Time, y = value, color = Model))+
  geom_line()+
  geom_point(data = df_Ns_agent_based, aes(x = Time, y = y, color = "agent"))+
  facet_wrap(Lost ~ Category)+
  ggtitle("Ns"))

merged <- lapply(c("branching", "cyton"), function (model) {
  ans <- merge(
    rbind(
      df_props_cyton_based %>%
        subset(Model == model) %>%
        .[, c("Time", "Generation", "Category", "value")] %>%
        transform(Lost = FALSE, Generation = Generation - 1L),
      df_props_cyton_based_lost %>%
        subset(Model == model) %>%
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

  ans$value.x[is.na(ans$value.x)] <- 0
  ans$value.y[is.na(ans$value.y)] <- 0

  return (transform(ans, diff = abs(value.x - value.y)))
}) %>% setNames(c("branching", "cyton"))

merged_Ns <- lapply(c("branching", "cyton"), function (model) {
  ans <- merge(
    df_Ns_cyton_based %>%
      subset(Model == model) %>%
      .[, c("Time", "Category", "Lost", "value")],
    df_Ns_agent_based %>%
      transform(value = y),
    by = c("Time", "Category", "Lost"),
    all = TRUE
  )

  ans$value.x[is.na(ans$value.x)] <- 0
  ans$value.y[is.na(ans$value.y)] <- 0

  ans <- transform(ans, diff = abs(log10(value.x) - log10(value.y)))
  ans$diff[!is.finite(ans$diff)] <- NA

  return (ans)
}) %>% setNames(c("branching", "cyton"))

cat("------ Summary ------\n")
print(rbind(
  branching = c(
    max = max(merged$branching$diff),
    mean = mean(merged$branching$diff %>% .[. > 0])),
  cyton = c(
    max = max(merged$cyton$diff),
    mean = mean(merged$cyton$diff %>% .[. > 0])),
  branching_Ns = c(
    max = max(merged_Ns$branching$diff, na.rm = TRUE),
    mean = mean(merged_Ns$branching$diff, na.rm = TRUE)),
  cyton_Ns = c(
    max = max(merged_Ns$cyton$diff, na.rm = TRUE),
    mean = mean(merged_Ns$cyton$diff, na.rm = TRUE))
))

assert_that(max(merged$branching$diff) < 0.06)
assert_that(max(merged$cyton$diff) < 0.05)

assert_that(mean(merged$branching$diff %>% .[. > 0]) < 0.005)
assert_that(mean(merged$cyton$diff %>% .[. > 0]) < 0.005)

assert_that(max(merged_Ns$branching$diff, na.rm = TRUE) < 0.1)
assert_that(max(merged_Ns$cyton$diff, na.rm = TRUE) < 0.15)

assert_that(mean(merged_Ns$branching$diff, na.rm = TRUE) < 0.1)
assert_that(mean(merged_Ns$cyton$diff, na.rm = TRUE) < 0.15)
