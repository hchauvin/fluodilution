context("peaks")

# nolint start

# Base constraints for Salmonella Typhimurium experiments
`#macr` <<- ~ FdCommonConstraints$`#noss` + pro:all:{res <- 0} + fmm:{c0 <- 0}
`#invivo` <<- ~ FdCommonConstraints$`#noss` + ~ pro:all:{res <- 0}

# nolint end

CC <<- FdCommonConstraints

test_that("peaks", {
  params <- fd_unif(fd_fmm_gaussian(m0 = 10.0, sd0 = 0.1), "branching",
                  constraints = ~ `#macr` + CC$`#delta_1111` + CC$`#mm_xyxz`)
  times <- c(2.0, 12.0, 24.0, 48.0)
  sim <- fd_simulate(params,
          times = times, breaks = 1000L,
          select=c("hists", "hists_lost", "Ns", "Ns_lost"))

  peaks <- transform(fd_findpeaks(sim, span=0.05, w=6L),
                     Generation = i - 1L)
  ans1 <- fd_peaks(sim, peaks, rm=TRUE)
  props <- fd_proportions(params, times = times)
  ans2 <- rbind(
    props$live_pop$One %>%
      reshape2::melt(varnames=c("Time", "a"), value.name="y2") %>%
      transform(Type = "props"),
    props$lost_pop$One %>%
      reshape2::melt(varnames=c("Time", "a"), value.name="y2") %>%
      transform(Type = "props_lost")
  ) %>%
    dplyr::mutate(Time = times[Time], a = a - 1)
  ansm <- merge(subset(ans1, Weight=="prop"), ans2)
  expect_equal(ansm$y, ansm$y2, tol=1e-3)

  # Check rm=TRUE/FALSE
  expect_equal(as.character(unique(fd_peaks(sim, peaks, rm=TRUE)$Type)),
               c("Ns", "Ns_lost", "props", "props_lost"))
  expect_equal(as.character(unique(fd_peaks(sim, peaks, rm=FALSE)$Type)),
               c("hists", "hists_lost", "Ns", "Ns_lost",
                 "props", "props_lost"))

  # Check silent when plot=FALSE
  expect_silent(fd_findpeaks(sim, span=0.05, w=6L, plot=FALSE))
})

test_that("proportions", {
    C <- FdCommonConstraints
    params <- fd_unif(fd_fmm_gaussian(m0 = 10.0, sd0 = 0.1),
                      "branching",
                      constraints = ~ `#macr` + CC$`#delta_1111` + CC$`#mm_xyxz`)
    times <- c(2.0, 12.0, 24.0, 48.0)
    sim <- fd_simulate(params,
                       times = times, breaks = 1000L,
                       select=c("hists", "hists_lost", "Ns", "Ns_lost"))

    expect_silent(fd_gaussian_fmm_solve(sim))
    ans1 <- fd_gaussian_fmm_solve(sim, model(sim)$fmm, mgen=8)

    props <- fd_proportions(params, times = times, mgen=8)
    ans2 <- rbind(
      props$live_pop$One %>%
        reshape2::melt(varnames=c("Time", "a"), value.name="y2") %>%
        transform(Type = "props"),
      props$lost_pop$One %>%
        reshape2::melt(varnames=c("Time", "a"), value.name="y2") %>%
        transform(Type = "props_lost")
    ) %>%
      dplyr::mutate(Time = times[Time], a = a - 1)

    ansm <- merge(subset(ans1, Weight=="prop"), ans2)
    expect_equal(ansm$y, ansm$y2, tol=1e-6)

    # Optimization for prop works
    expect_warning(subset(ans1, Weight == "prop" & a < 4),
                   "Some timepoints do not have proportions that sum to 1.")
    expect_error(fd_minuslogl(model(params), ans1, mgen = 7L),
                 "'mgen' not big enough for 'prop': should be at least 8")
    expect_silent(fd_minuslogl(model(params), ans1, mgen = 8L, verbose=FALSE)())
})