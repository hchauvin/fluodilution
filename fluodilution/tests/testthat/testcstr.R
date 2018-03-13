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

context("cstr")

test_that("cstr are constructed correctly", {
  expect_true(is.cstr(cstrlist()))
  expect_true(is.cstr(~ "jdfsdfs"))
  expect_false(is.cstr("kdsfsdfds"))
  expect_error(cstrlist("kdsfsf"))
  expect_error(cstrlist(constraints = "ldsfdsfds"))
  expect_false(
    is.cstr(structure(list(constraints = "lsfsdf"), class="cstrlist")))
  expect_true(
    is.cstr(structure(list(constraints = ~ "lsfsdf"), class="cstrlist")))
})

test_that("cstr are parsed correctly", {
  expect_error(parsecstr(~ "ldsfsdfs"), "unknown syntax")
  expect_error(parsecstr("ldsfdsfsd"), "unknown syntax")
  expect_equal(parsecstr(~ {a <- 10}), list(prefix = "a", value = "10"))
  expect_equal(parsecstr(~ NULL), list(prefix = NULL, value = NULL))
  expect_equal(
    parsecstr(expression({a <- 10})),
    list(prefix = "a", value = "10"))
  expect_error(parsecstr(expression()), "empty expression")
  expect_error(parsecstr(expression({})), "empty expression")
  expect_equal(parsecstr(~ ~ {a <- 10}), list(prefix = "a", value = "10"))
  expect_error(parsecstr(~ {a = 10}), "terms must be of the form 'a <- b'")
  expect_equal(parsecstr(~ {a <- 10; b <- 20}),
               list(prefix = c("a", "b"), value = c("10", "20")))
  expect_error(parsecstr(~ {a <- 10; {kdsfds}}),
               "terms must be of the form 'a <- b'")
  expect_equal(parsecstr(~ aa:bb:{c <- 30}),
               list(prefix = "aa$bb$c", value = "30"))
  expect_equal(parsecstr(~ NULL + {a <- 10}), list(prefix = "a", value = "10"))
  expect_equal(parsecstr(~ {a <- 10} + {b <- 20}),
               list(prefix = c("a", "b"), value = c("10", "20")))
  expect_equal(parsecstr(~ (AA/aa):(BB/bb):{c <- 10}),
               list(prefix = c("AA$BB$c", "aa$BB$c", "AA$bb$c", "aa$bb$c"),
                    value = rep("10", 4)))
  expect_equal(
    parsecstr(~ pro:(all:f0:{delta <- 1; ss <- 0.5} +
                     One:(f0:{delta <- ..free..} +
                          f0:{delta <- 0.5}))),
    structure(
      list(
        prefix = structure(
          c("pro$all$f0$delta", "pro$all$f0$ss", "pro$One$f0$delta"),
          na.action = structure(3L, class = "omit")),
        value = structure(c("1", "0.5", "0.5"),
        na.action = structure(3L, class = "omit"))),
      .Names = c("prefix", "value")))
})

test_that("boxing/unboxing works as intended", {
  expect_error(box.param(c(1, 2, 3), 1, c(4, 5)),
               "'param', 'lower' and 'upper' must have the same length")
  expect_error(unbox.param(c(1, 2, 3), 1, c(4, 5)),
               "'param', 'lower' and 'upper' must have the same length")
  expect_true(is.na(box.param(NA, 1, 2)))
  expect_error(unbox.param(-1, 0, 1),
               "'param' is outside boundary")
  for (i in 1:100) {
    lower <- runif(10, -100, 100)
    upper <- runif(10, 500, 1000)
    param <- runif(10, 100, 500)
    expect_equal(box.param(unbox.param(param, lower, upper),
                           lower, upper),
                 param)
  }
})

test_that("'..free..' works as intended", {
  cstr1 <- ~ pro:all:{p0 <- .L1$p}
  mdl1 <- fd_model(data = NULL, constraints = cstr1)
  expect_true(!("pro.One.p0" %in% names(start(mdl1))))
  drawn <- fd_draw_unif(mdl1, 10)
  for (i in 1:NROW(drawn)) {
    rel <- relist(drawn[i, ], mdl1)$pro$One
    expect_equal(rel$p + rel$res, rel$p0 + rel$res0)
  }

  cstr2 <- catcstr(cstr1, ~ pro:all:{p0 <- ..free..})
  mdl2 <- fd_model(data = NULL, constraints = cstr2)
  expect_true("pro.One.p0" %in% names(start(mdl2)))
  drawn <- fd_draw_unif(mdl2, 10)
  for (i in 1:NROW(drawn)) {
    rel <- relist(drawn[i, ], mdl2)$pro$One
    expect_true(rel$p + rel$res != rel$p0 + rel$res0)
  }
})

test_that("'catcstr' works as intended", {
  CC <<- FdCommonConstraints
  expect_equal(
    catcstr(CC$`#noss`, CC$`#mm_xyxz`, CC$`#mm_xyxy`),
    catcstr(CC$`#noss`, CC$`#mm_xyxz`, CC$`#mm_xyxy`, drop=F)[[1L]]
  )
  ans <- catcstr(
    CC$`#noss`,
    list(CC$`#mm_xyxz`, CC$`#mm_xyxy`, CC$`#mm_xxxx`),
    list(CC$`#delta_1100`, CC$`#delta_1111`))
  expect_equal(
    ans,
    list(substitute((~pro:all:(f0/f/g0/g):{
      ss <- 0.5
    }) + ((~pro:all:(f0:{
      mm <- .L2$f$mm
    })) + ~pro:all:((f0/g0):{
      delta <- 1
    } + (f/g):{
      delta <- 0.01
    }))), substitute((~pro:all:(f0/f/g0/g):{
      ss <- 0.5
    }) + ((~pro:all:(f0:{
      mm <- .L2$f$mm
    })) + ~pro:all:(f0/f/g0/g):{
      delta <- 1
    })), substitute((~pro:all:(f0/f/g0/g):{
      ss <- 0.5
    }) + ((~pro:all:(f0:{
      mm <- .L2$f$mm
    } + g0:{
      mm <- .L2$g$mm
    })) + ~pro:all:((f0/g0):{
      delta <- 1
    } + (f/g):{
      delta <- 0.01
    }))), substitute((~pro:all:(f0/f/g0/g):{
      ss <- 0.5
    }) + ((~pro:all:(f0:{
      mm <- .L2$f$mm
    } + g0:{
      mm <- .L2$g$mm
    })) + ~pro:all:(f0/f/g0/g):{
      delta <- 1
    })), substitute((~pro:all:(f0/f/g0/g):{
      ss <- 0.5
    }) + ((~pro:all:(f0/g0/g):{
      mm <- .L2$f$mm
    }) + ~pro:all:((f0/g0):{
      delta <- 1
    } + (f/g):{
      delta <- 0.01
    }))), substitute((~pro:all:(f0/f/g0/g):{
      ss <- 0.5
    }) + ((~pro:all:(f0/g0/g):{
      mm <- .L2$f$mm
    }) + ~pro:all:(f0/f/g0/g):{
      delta <- 1
    }))))
})
