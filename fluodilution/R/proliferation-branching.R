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

#' @rdname proliferation
#' @export
#'
#' @include proliferation-.R
fd_proliferation_branching <- fd_proliferation(
  name = "branching",
  new = function (categories="One", mgen=8, log10_scale=FALSE,
          initial=NULL, ...) {
    if (is.null(initial)) {
      initial <- matrix(0, ncol=mgen + 1, nrow=length(categories))
      initial[1, 1] <- 1
    } else if (!is.matrix(initial) ||
           dim(initial) != c(length(categories), mgen + 1)) {
      stop("'initial' must be a matrix of dimensions ",
         "(number of categories) x (mgen+1)")
    }
    list(
      categories = categories,
      mgen = mgen,
      log10_scale = log10_scale,
      initial = initial,
      ...
    )
  },
  adapt = function (self, theta, mgen) {
    if (is.null(mgen)) {
      mgen <- self$mgen
      initial_ <- self$initial
    } else {
      initial_ <- matrix(0, nrow=NROW(self$initial), ncol=mgen + 1)
      initial_[, 1:min(NCOL(self$initial), NCOL(initial_))] <-
        self$initial[, 1:min(NCOL(self$initial), NCOL(initial_))]
    }
    nms <- names(theta)
    if (any(all.equal(self$categories, nms) != TRUE))
      stop("wrong categories or in the wrong order")
    for (nm in nms) {
      theta[[nm]]$f <- unlist(theta[[nm]]$f)
      theta[[nm]]$f0 <- unlist(theta[[nm]]$f0)
      theta[[nm]]$g <- unlist(theta[[nm]]$g)
      theta[[nm]]$g0 <- unlist(theta[[nm]]$g0)
      if (!is.null(theta[[nm]]$h))
        theta[[nm]]$h <- unlist(theta[[nm]]$h)
      theta[[nm]]$mrates0 <- unlist(theta[[nm]]$mrates0)
      if (is.null(theta[[nm]]$mrates0))
        theta[[nm]]$mrates0 <- numeric(0)
      theta[[nm]]$mrates <- unlist(theta[[nm]]$mrates)
      if (is.null(theta[[nm]]$mrates))
        theta[[nm]]$mrates <- numeric(0)
    }
    list(mgen = mgen, theta = theta, initial = initial_)
  },
  precalculate = function (self, theta, startpop, times, mgen = NULL) {
    if (!is.null(self$obj)) stop("'precalculate' already called")
    if (startpop > 1) {
      a <- self$adapt(theta, mgen)
      self[["precalc"]] <- list(
        obj = modelBranchingInit(a$theta[1L:(startpop - 1)],
                     startpop - 1,
                     a$mgen,
                     times,
                     a$initial),
        times = times,
        mgen = a$mgen,
        theta = a$theta[1L:(startpop - 1)],
        startpop = startpop
      )
    }
  },
  model = function (self, theta, times, mgen = NULL) {
    a <- self$adapt(theta, mgen)
    if (is.null(self$precalc)) {
      obj <- modelBranchingInit(a$theta, 0, a$mgen, times, a$initial)
    } else {
      if (any(all.equal(self$precalc$times, times) != TRUE) ||
          a$mgen != self$precalc$mgen ||
          any(all.equal(a$theta[1L:(self$precalc$startpop - 1)],
                 self$precalc$theta) != TRUE)) {
        stop("'precalculate' was called with different arguments")
      }
      obj <- self$precalc$obj
    }

    ret <- tryCatch(
      modelBranching(a$theta, obj),
      finally = {
        # Release immediately
        if (is.null(self$precalc)) modelBranchingRelease(obj)
      })

    ret$Ns <- do.call(cbind, ret$Ns)
    ret$Ns_lost <- do.call(cbind, ret$Ns_lost)
    names(ret$lost_pop) <- self$categories
    colnames(ret$Ns) <- self$categories
    colnames(ret$Ns_lost) <- self$categories
    names(ret$live_pop) <- self$categories
    names(ret$cum_influx) <- self$categories
    return (ret)
  },
  constrain = function (self, object, type) {
    mrates <- switch(type, start = 0.5, lower = 0, upper = 1)
    if (!is.null(mrates) & length(self$categories) > 1) {
      if (self$log10_scale) mrates <- log10(mrates)
      mrates_cstr <- sapply(
        2L:length(self$categories),
        function (i) {
          sq <- paste(paste0(self$categories[1L:(i - 1)],
                     " <- ", mrates),
                collapse="; ")
          paste0("pro:",
               self$categories[i],
               ":(mrates0:{", sq, "} + mrates:{", sq, "})")
        }
      ) %>%
        paste(collapse=" + ") %>%
        parse(text=.) %>%
        .[[1L]]
    } else {
      mrates_cstr <- NULL
    }
    if (self$log10_scale) {
      # nolint start
      catcstr(
        switch(type,
          start = ~ pro:all:({res0 <- log10(0.2);
                      res <- log10(0.2); p <- log10(0.5);
                      p0 <- log10(0.49)} +
                    (g/g0):{mm <- 5; delta <- 1; ss <- 0.5} +
                    (f/f0):{mm <- 5; delta <- 0.2;
                      ss <- 0.5}),
          lower = ~ pro:all:({res0 <- log10(0); res <- log10(0);
                      p <- log10(0.1);
                      p0 <- log10(0.1)} +
                    (g/g0):{mm <- 1; delta <- 0.01;
                      ss <- 0.1} +
                    (f/f0):{mm <- 1; delta <- 0.01;
                      ss <- 0.1}),
          upper = ~ pro:all:({res0 <- log10(1); res <- log10(1);
                      p <- log10(0.9);
                      p0 <- log10(0.9)} +
                    (g/g0):{mm <- 20; delta <- 10; ss <- 0.5} +
                    (f/f0):{mm <- 20; delta <- 10; ss <- 0.5})
        ),
        mrates_cstr,
        object
      ) %>%
        expand_categories(self$categories)
      # nolint end
    } else {
      # nolint start
      catcstr(
        switch(type,
          start = ~ pro:all:({res0 <- 0.2; res <- 0.2; p <- 0.5; p0 <- 0.49} +
                    (g/g0):{mm <- 5; delta <- 1; ss <- 0.5} +
                    (f/f0):{mm <- 5; delta <- 0.2; ss <- 0.5}),
          lower = ~ pro:all:({res0 <- 0; res <- 0; p <- 0.1; p0 <- 0.1} +
                    (g/g0):{mm <- 1; delta <- 0.01; ss <- 0.1} +
                    (f/f0):{mm <- 1; delta <- 0.01; ss <- 0.1}),
          upper = ~ pro:all:({res0 <- 1; res <- 1; p <- 0.9; p0 <- 0.9} +
                    (g/g0):{mm <- 20; delta <- 10; ss <- 0.5} +
                    (f/f0):{mm <- 20; delta <- 10; ss <- 0.5})
        ),
        mrates_cstr, object
      ) %>%
        expand_categories(self$categories)
      # nolint end
    }
  },
  trans = function (self, params) {
    if (length(params) > 2) {
      for (i in 1L:(length(params) - 2)) {
        from <- names(params)[i]
        cur_sum <- 0
        cur_sum0 <- 0
        for (j in (i + 1):length(params)) {
          cur_rate <- params[[j]]$mrates[[from]]
          params[[j]]$mrates[[from]] <- ifelse(
            cur_sum >= 1, 0, cur_rate / (1 - cur_sum))
          cur_sum <- cur_sum + cur_rate

          cur_rate0 <- params[[j]]$mrates0[[from]]
          params[[j]]$mrates0[[from]] <- ifelse(
            cur_sum0 >= 1, 0, cur_rate0 / (1 - cur_sum0))
          cur_sum0 <- cur_sum0 + cur_rate0
        }
      }
    }

    for (i in 1:length(params)) {
      params[[i]]$p0 <- params[[i]]$p0 + params[[i]]$res0
      params[[i]]$res0 <- params[[i]]$res0 / params[[i]]$p0
      if (is.na(params[[i]]$res0)) params[[i]]$res0 <- 0

      params[[i]]$p <- params[[i]]$p + params[[i]]$res
      params[[i]]$res <- params[[i]]$res / params[[i]]$p
      if (is.na(params[[i]]$res)) params[[i]]$res <- 0
    }

    # From linear scale to log10 scale
    if (self$log10_scale) {
      for (i in 1:length(params)) {
        params[[i]]$res0 <- log10(params[[i]]$res0)
        params[[i]]$p0 <- log10(params[[i]]$p0)

        params[[i]]$res <- log10(params[[i]]$res)
        params[[i]]$p <- log10(params[[i]]$p)

        params[[i]]$mrates <- log10(params[[i]]$mrates)
        params[[i]]$mrates0 <- log10(params[[i]]$mrates0)
      }
    }

    return (params)
  },
  trans_inverse = function (self, params) {
    # From log10 scale to linear scale
    if (self$log10_scale) {
      for (i in 1:length(params)) {
        params[[i]]$res0 <- 10^params[[i]]$res0
        params[[i]]$p0 <- 10^params[[i]]$p0

        params[[i]]$res <- 10^params[[i]]$res
        params[[i]]$p <- 10^params[[i]]$p

        params[[i]]$mrates <- 10^params[[i]]$mrates
        params[[i]]$mrates0 <- 10^params[[i]]$mrates0
      }
    }

    # Reparametrization of rates
    if (length(params) > 2) {
      for (i in 1L:(length(params) - 2)) {
        from <- names(params)[i]
        cur_sum <- 0
        cur_sum0 <- 0
        for (j in (i + 1):length(params)) {
          params[[j]]$mrates[[from]] <-
            (1 - cur_sum) * params[[j]]$mrates[[from]]
          cur_sum <- cur_sum + params[[j]]$mrates[[from]]

          params[[j]]$mrates0[[from]] <-
            (1 - cur_sum0) * params[[j]]$mrates0[[from]]
          cur_sum0 <- cur_sum0 + params[[j]]$mrates0[[from]]
        }
      }
    }

    # Reparametrization of probabilities
    for (i in 1:length(params)) {
      params[[i]]$res0 <- params[[i]]$p0 * params[[i]]$res0
      params[[i]]$p0 <- params[[i]]$p0 - params[[i]]$res0

      params[[i]]$res <- params[[i]]$p * params[[i]]$res
      params[[i]]$p <- params[[i]]$p - params[[i]]$res
    }

    return (params)
  }
)
