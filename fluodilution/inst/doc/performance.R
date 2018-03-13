## ----message=F, warning=F-------------------------------------------
options(width=70)
library(fluodilution)
library(microbenchmark)
library(plyr)

source(system.file("contrib", "nlsSA.R", package="fluodilution"))

# nolint start
`#macr` <- ~ FdCommonConstraints$`#noss` + pro:all:{
        res <- 0
    } + fmm:{
        c0 <- 0
    }
`#invivo` <- ~ FdCommonConstraints$`#noss` + ~ pro:all:{res <- 0}
# nolint end

attach(FdCommonConstraints)

sessionInfo()

## -------------------------------------------------------------------
TIMES <- 50L
TIMES2 <- 10L

LL <- function (mdl, mgen) {
    do.call(fd_minuslogl(mdl, dat, verbose=FALSE, mgen = mgen),
                as.list(start(mdl)))
}

PRO <- function (mdl, mgen, times = seq(0, 22, 2)) {
    fd_proportions(start(mdl), times, mgen=mgen, model = mdl)
}

PRED <- function (mdl, mgen) {
    fd_predict(model=mdl, param=start(mdl), data=dat, mgen = mgen)
}

## -------------------------------------------------------------------
cstr <- ~ `#macr`
dat <- cutoff(FdSTyphimuriumWTC57)
mdl_branching <- fd_model(dat, 
                          fmm="gaussian", 
                          proliferation="branching",
                          constraints = cstr)
mdl_cyton <- fd_model(dat, 
                          fmm="gaussian", 
                          proliferation="cyton",
                          constraints = cstr)

## -------------------------------------------------------------------
(macr_first <- microbenchmark(
    LL(mdl_branching, mgen=10),
    LL(mdl_branching, mgen=20),
    LL(mdl_cyton, mgen=10),
    LL(mdl_cyton, mgen=20),
    times = TIMES
))

## -------------------------------------------------------------------
(macr_second <- microbenchmark(
    LL(mdl_branching, mgen=10),
    LL(mdl_branching, mgen=20),
    LL(mdl_cyton, mgen=10),
    LL(mdl_cyton, mgen=20),
    times = TIMES
))

## -------------------------------------------------------------------
(macr_pro <- microbenchmark(
    PRO(mdl_branching, mgen=10),
    PRO(mdl_branching, mgen=20),
    PRO(mdl_cyton, mgen=10),
    PRO(mdl_cyton, mgen=20),
    times = TIMES
))

## -------------------------------------------------------------------
(macr_pred <- microbenchmark(
    PRED(mdl_branching, mgen=10),
    PRED(mdl_branching, mgen=20),
    PRED(mdl_cyton, mgen=10),
    PRED(mdl_cyton, mgen=20),
    times = TIMES
))

## -------------------------------------------------------------------
pro <- ddply(macr_pro, "expr", summarise, time = mean(time))
total <- ddply(macr_first, "expr", summarise, time = mean(time))
wo_fmm <- ddply(macr_second, "expr", summarise, time = mean(time))
pred <- ddply(macr_pred, "expr", summarise, time = mean(time))

data.frame(expr = total$expr, 
           pro = pro$time / total$time, 
           fmm = (total$time - wo_fmm$time) / total$time,
           loglik = (wo_fmm$time - pred$time) / total$time,
           other = (pred$time - pro$time) / total$time,
           total = 1)

## -------------------------------------------------------------------
load(system.file("extdata", "FdSTyphimuriumMice.rda", package="fluodilution"))
dat <- cutoff(FdSTyphimuriumMice)

mdl_all <- fd_model(data = dat, constraints = `#invivo`)
mdl_partial <- fd_model(
  data = subset(dat, Category %in% c("SI", "PP")),
  constraints = `#invivo`, partial = 1:2)
mdl_fixed_10 <- fd_model(
  data = dat,
  constraints = fixcstr(`#invivo`, start(mdl_partial)))
mdl_fixed_10$pro$precalculate(
  mdl_fixed_10$pro$trans_inverse(mdl_fixed_10$start$pro),
  3,
  sort(unique(dat$Time)),
  mgen = 10)
mdl_fixed_20 <- fd_model(
  data = dat,
  constraints = fixcstr(`#invivo`, start(mdl_partial)))
mdl_fixed_20$pro$precalculate(
  mdl_fixed_20$pro$trans_inverse(mdl_fixed_20$start$pro),
  3,
  sort(unique(dat$Time)),
  mgen = 20)
times <- c(2, 6, 12, 24, 48)

microbenchmark(
    PRO(mdl_all, mgen = 10),
    PRO(mdl_all, mgen = 20),
    PRO(mdl_partial, mgen = 10),
    PRO(mdl_partial, mgen = 20),
    PRO(mdl_fixed_10, mgen = 10, times = times),
    PRO(mdl_fixed_20, mgen = 20, times = times),
    times = TIMES2
)

