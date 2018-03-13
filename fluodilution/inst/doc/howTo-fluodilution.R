## ----message=F,warning=F,echo=F---------------------------
options(width=60)

## ----message=F, warning=F---------------------------------
library(fluodilution)
data(FdSTyphimuriumWTC57)
invisible(summary(FdSTyphimuriumWTC57))

## ----eval=F-----------------------------------------------
#  plot(FdSTyphimuriumWTC57, type=c("hist", "N"))

## ---------------------------------------------------------
# plot(FdSTyphimuriumWTC57, type="cutoff")
dat <- cutoff(FdSTyphimuriumWTC57)

## ---------------------------------------------------------
mdl <- fd_model(fmm="gaussian", proliferation="branching", data=dat)

## ---------------------------------------------------------
cbind(Start = start(mdl), Lower = lower(mdl), Upper = upper(mdl))

## ---------------------------------------------------------
# nolint start
cstr <- ~
    pro:all:(f0/f/g0/g):{
        ss <- 0.5
    } + pro:all:{
        res <- 0
    } + fmm:{
        c0 <- 0
    } + ~pro:all:((f0/g0):{
        delta <- 1
    } + (f/g):{
        delta <- 0.01
    }) + pro:all:f0:{
        mm <- .L2$f$mm
    }
# nolint end
mdl <- update(mdl, data=dat, addcstr=cstr)

## ---------------------------------------------------------
fd_formula("mdl", mgen=10)

## ---------------------------------------------------------
source(system.file("contrib", "nlsSA.R", package="fluodilution"))
# Put maxit=100 to actually run the optimization
fit <- nlsSA(
    fd_formula("mdl", mgen = 10),
    data = dat,
    start = start(mdl),
    lower = lower(mdl),
    upper = upper(mdl),
    trace = T,
    control = list(simple.function = T, maxit = 1)
)

## ---------------------------------------------------------
coef(fit)
unlist(relist(coef(fit), mdl))

