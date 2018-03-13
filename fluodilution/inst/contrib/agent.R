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

#' Agent-based simulation: generate a cell tree from a single progenitor, using a branching
#' process.
#'
#' @param params Parameters, in relisted form
#' @param times
#' @param initial Number of evaluations
#'
#' @return A cell tree (data.frame)
agent_tree <- function (params, times, initial, cyton = FALSE, length.out = 500L) {
  timestep <- max(times) / length.out

  age <- rep(0, sum(initial > 0))
  # Age, from 0 to ...
  generation <- (which(initial > 0) - 1) %/% (NCOL(initial) + 1)
  # Number of generation, from 0 to ...
  fate <- rep(NA, sum(initial > 0))
  # If 0: non-replicator, = 1L divides with time-to-divide fate,
  # = -1L lost with time-to-loss -fate
  time_to_evolve <- rep(NA, sum(initial > 0))
  # Time till next event: Inf if non-replicator, or time-to-loss or
  # time-to-divide.
  status <- (which(initial > 0) - 1) %% (NCOL(initial) + 1) + 1
  # 0 for dead, 1 for alive in compartment 1, 2 for compartment 2, etc.
  # (as given by the order of params).

  snapshots <- data.frame()

  cur.snapshot <- 1

  for (cur.t in seq(0, max(times) + timestep, timestep)) {
    #print(data.frame(age = age, generation = generation, fate = fate, time_to_evolve = time_to_evolve, status = status)[status == 2, ])
    while (any(is.na(time_to_evolve) |
           (status > 0 & time_to_evolve <= age))) {
      # deals with cells that evolve "immediately" (that's why loop)

      # Select fate of cells of age 0
      for (ct in 1:length(params)) {
        just_born <- is.na(time_to_evolve) & status == ct
        if (cur.t == 0) {
          p <- params[[ct]]$p0
          res <- params[[ct]]$res0
          div_dist <- fd_pack_dist(params[[ct]]$f0)
          lost_dist <- fd_pack_dist(params[[ct]]$g0)
        } else {
          p <- params[[ct]]$p
          res <- params[[ct]]$res
          div_dist <- fd_pack_dist(params[[ct]]$f)
          lost_dist <- fd_pack_dist(params[[ct]]$g)
        }

        fate[just_born] <- sample(
          c(1, 0, -1),
          sum(just_born),
          prob = c(p, res, max(0, 1 - p - res)),
          replace = TRUE
        )

        should_divide <- just_born & fate == 1
        time_to_evolve[should_divide] <-
          rgamma(sum(should_divide), shape = div_dist[, "a"],
               rate = div_dist[, "b"])+
          div_dist[, "loc"]

        should_loss = just_born & fate == -1
        time_to_evolve[should_loss] <-
          rgamma(sum(should_loss), shape = lost_dist[, "a"],
               rate = lost_dist[, "b"])+
          lost_dist[, "loc"]

        should_rest = just_born & fate == 0
        time_to_evolve[should_rest] <- Inf

        # Evolve cells
        lost_now = (status == ct & time_to_evolve <= age & fate == -1)
        assert_that(ct > 0 && ct <= length(params))
        if (ct < length(params)) {
          mrates <- sapply(
            params,
            function (cur) {
              nm <- names(params)[ct]
              ans <- if (cur.t == 0) cur$mrates0[[nm]] else cur$mrates[[nm]]
              if (is.null(ans)) 0 else ans
            }
          )
          status[lost_now] <- -ct
          evol <- sample(
            c(-ct, 1:length(params)),
            sum(lost_now),
            prob = c(1 - sum(mrates), mrates),
            replace = TRUE)
          status <- c(status, evol[evol > 0])
          time_to_evolve <- c(time_to_evolve, rep(NA, sum(evol > 0)))
          fate <- c(fate, rep(NA, sum(evol > 0)))
          age <- c(age, rep(0, sum(evol > 0)))
          generation <- c(generation, generation[lost_now][evol > 0])
        } else
          status[lost_now] <- -ct

        divide_now <- (status == ct & time_to_evolve <= age & fate == 1)
        new_cells <- sum(divide_now)
        age[divide_now] <- 0
        age <- c(age, rep(0, new_cells))
        generation[divide_now] <- generation[divide_now] + 1
        generation <- c(generation, generation[divide_now])
        fate[divide_now] <- NA
        fate <- c(fate, rep(NA, new_cells))
        time_to_evolve[divide_now] <- NA
        time_to_evolve <- c(time_to_evolve, rep(NA, new_cells))
        status <- c(status, rep(ct, new_cells))
      }
    }

    # Age cells
    age <- age + timestep

    # Take snapshot
    if (times[cur.snapshot] <= cur.t) {
      snapshots <- rbind(
        snapshots,
        data.frame(
          Time = times[cur.snapshot],
          Idx = cur.snapshot,
          Generation = generation,
          Category = names(params)[abs(status)],
          Lost = status < 0
        )
      )
      if (cur.snapshot == length(times)) break
      cur.snapshot <- cur.snapshot + 1
    }
  }

  return (snapshots)
}

#' Generate a tree bank (multiple calls to \code{\link{sim_generate_tree}}).
#' @export
agent_tree_bank <- function (n, params, times, initial, cyton=FALSE,
                             length.out = 500L) {
  trees <- new.env()

  trees$n <- n
  trees$size_trees <- vector(mode="numeric", length=n)
  trees$trees <- do.call("rbind", lapply(1:n, function (i) {
    cat(sprintf("Tree %s...\n", i))
    ret <- agent_tree(
      params$pro, times, initial, cyton=cyton, length.out = 500L)
    if (NROW(ret) > 0) ret$Tree <- i
    trees$size_trees[i] <- NROW(ret)
    return (ret)
  }))
  trees$times <- times
  trees$categories <- names(params$pro)
  trees$params <- params
  trees$mgen <- max(trees$trees$Generation)

  return (trees)
}

#' Sample fluorescence for agent-based simulation.
#'
#' @param trees
#' @param B1 Number of trees in the tree bank
#' @param N0 Number of initial bacteria; can be Inf, in this case sampling with replacement.
#' @param Ni Number of bacteria sampled at each time point
#'
#' @return The sample fluorescence (a \code{data.frame})
sample_agent_trees <- function (bank, N0, Ni) {
  # Current implementation is motivated by NOT exploding your memory,
  # so you cannot just copy whole trees and think, hey, that might
  # be good enough!

  # N0 = Inf: sample with replacement from trees
  if (is.finite(N0)) {
    if (N0/bank$n > 1e6)
      bank_idx <- ceiling(runif(N0, 0, bank$n))
    else bank_idx <- table(sample.int(bank$n, N0, replace = TRUE))
  }

  do.call("rbind", lapply(bank$times, function (curTime) {
    do.call("rbind", lapply(bank$categories, function (curCategory) {
      do.call("rbind", lapply(c(TRUE, FALSE), function (lost) {
        if (is.finite(N0)) {
          tally = unlist(lapply(1:bank$n, function (i) {
            rep(bank_idx[names(bank_idx) == i],
              sum(bank$trees$Time == curTime &
                  bank$trees$Category == curCategory &
                  bank$trees$Tree == i &
                  bank$trees$Lost == lost))
          }))
          if (sum(tally) < Ni) {
            warning("not enough bacteria to sample from at time ", curTime,
                " in category ", curCategory, ": ", sum(tally), " given, ",
                Ni, " required; sampling with replacement instead")
          }
          if (length(tally) == 0) return (data.frame())

          # Using sampling *with* replacement, the probability of choosing
          # the same event twice, when you choose Ni events out of a list
          # of N0 events, is (Ni/N0)^2.  If we want this to be less than 1%,
          # we must have Ni < N0 * 1e-4.  In this case, we perform a
          # sampling *with* replacement, which is much more efficient.
          sampled <- vector(mode = "numeric", length = Ni)
          for (i in 1:Ni) {
            sampled[i] <- sample(which(tally > 0), 1)
            tally[sampled[i]] <- tally[sampled[i]] - 1
          }
          ret = bank$trees[
            bank$trees$Time == curTime &
              bank$trees$Category == curCategory &
              bank$trees$Lost == lost, ][sampled, ]
        } else {  # N0 = Inf: sampling with replacement
          ret = bank$trees[
            bank$trees$Time == curTime &
              bank$trees$Category == curCategory &
              bank$trees$Lost == lost, ]
          if (NROW(ret) > 0) {
            sampled <- sample.int(NROW(ret), Ni, replace=TRUE)
            ret <- ret[sampled, ]
          }
        }
        return (ret)
      }))
    }))
  }))
}
