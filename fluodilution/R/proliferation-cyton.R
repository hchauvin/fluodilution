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

#' @export
#' @rdname proliferation
#'
#' @include proliferation-.R
fd_proliferation_cyton <- fd_proliferation(
  name = "cyton",
  new = function (mgen = 8L, length.out = 500L, ...) {
    list(
      categories="One",
      mgen=mgen,
      length.out = length.out,
      ...
    )
  },
  model = function (self, theta, times, mgen=NULL) {
    if (is.null(mgen)) mgen <- self$mgen
    value <- theta[[1L]]
    ret <- model_cyton(value, times, mgen,
              length.out = self$length.out)
    ret$Ns <- cbind(ret$Ns)
    colnames(ret$Ns) <- self$categories
    ret$Ns_lost <- cbind(ret$Ns_lost)
    colnames(ret$Ns_lost) <- self$categories
    ret$live_pop <- setNames(list(ret$live_pop), self$categories)
    ret$lost_pop <- setNames(list(ret$lost_pop), self$categories)
    return (ret)
  },
  constrain = function (self, object, type) {
    # nolint start
    catcstr(
      switch(type,
        start = ~ pro:all:({res0 <- 0.2; res <- 0.2} +
                  (g/g0):{mm <- 5; delta <- 1; ss <- 0.5} +
                  (f/f0):{mm <- 5; delta <- 0.2; ss <- 0.5}),
        lower = ~ pro:all:({res0 <- 0; res <- 0} +
                  (g/g0/f/f0):{mm <- 0.5; delta <- 0.1; ss <- 0.1} +
                    (f/g):{delta <- 0.01}),
        upper = ~ pro:all:({res0 <- 1; res <- 1} +
                  (g/g0/f/f0):{mm <- 5; delta <- 1; ss <- 0.5})
      ),
      object) %>% expand_categories(self$categories)
    # nolint end
  }
)

# Models are implemented on the side for easy customisation
model_cyton <- function (params, times = c(2, 4, 6, 8, 12, 24), mgen = 15L,
                         length.out = 500L) {
  etimes <- sort(unique(c(times, seq(0, max(times), length.out = length.out))))
  outtimes <- which(etimes %in% times)

  # Precalculate F_g and G_f
  dt <- c(etimes[-1], etimes[length(etimes)]) - etimes

  F0dist <- fd_pack_dist(params$f0)
  G0dist <- fd_pack_dist(params$g0)
  F0_g0 <- (1 - fd_pdist(etimes, F0dist)) * fd_ddist(etimes, G0dist)
  G0_f0 <- (1 - fd_pdist(etimes, G0dist)) * fd_ddist(etimes, F0dist)
  F0_g0[1] <- 0
  G0_f0[1] <- 0

  Fdist <- fd_pack_dist(params$f)
  Gdist <- fd_pack_dist(params$g)
  F_g_dt <- (1 - fd_pdist(etimes, Fdist)) * fd_ddist(etimes, Gdist, dt)
  G_f_dt <- (1 - fd_pdist(etimes, Gdist)) * fd_ddist(etimes, Fdist, dt)
  F_g_dt[1] <- 0
  G_f_dt[1] <- 0

  # Rates (r: division rate; d: death rate)
  r <- matrix(NA, ncol=length(etimes), nrow=mgen + 1)
  d <- matrix(NA, ncol=length(etimes), nrow=mgen + 1)

  r[1, ] <- (1 - params$res0) * G0_f0
  d[1, ] <- (1 - params$res0) * F0_g0

  # Ns
  N <- matrix(NA, ncol=length(etimes), nrow=mgen + 1)
  N_lost <- matrix(NA, ncol=length(etimes), nrow=mgen + 1)
  N_lost[1, ] <- cumsum(d[1, ] * dt)
  N[1, ] <- pmax(1 - cumsum(r[1, ] * dt) - N_lost[1, ], 0)

  for (i in 1:mgen) {
    r[i + 1, ] <- 2 * (1 - params$res) * cytonConv(r[i, ], G_f_dt)
    d[i + 1, ] <- 2 * (1 - params$res) * cytonConv(r[i, ], F_g_dt)

    N_lost[i + 1, ] <- cumsum(d[i + 1, ] * dt)
    N[i + 1, ] <-
      pmax(cumsum((2 * r[i, ] - r[i + 1, ]) * dt) - N_lost[i + 1, ], 0)
  }

  N_out <- N[, outtimes]
  N_lost_out <- N_lost[, outtimes]

  Ns <- colSums(N_out)
  Ns_lost <- colSums(N_lost_out)

  live_pop <- t(N_out) / Ns
  live_pop[!is.finite(live_pop)] <- 0

  lost_pop <- t(N_lost_out) / Ns_lost
  lost_pop[!is.finite(lost_pop)] <- 0

  list(
    Ns = Ns,
    Ns_lost = Ns_lost,
    live_pop = live_pop,
    lost_pop = lost_pop
  )
}
