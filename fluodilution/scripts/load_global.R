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

# Load package in .GlobalEnv
############################

for (fname in list.files("src", pattern="\\.cpp$")) {
  if (fname != "RcppExports.cpp") {
    fname_exp <- paste0("src/", fname)
    cat(sep="", "CppSourcing ", fname_exp, "...\n")
    Rcpp::sourceCpp(fname_exp)
  }
}

# Collation
con <- file("DESCRIPTION")
L <- readLines(con)
close(con)
collate <- stringr::str_match(L, "^Collate:")
if (all(is.na(collate))) {
  fnames <- list.files("R/", pattern="\\.R$")
} else {
  fnames <- stringr::str_match(
    L[which(!is.na(collate)):length(L)], "'([^']+)'$")[, 2]
  end <- which(is.na(fnames[-1]))
  if (length(end) == 0) end <- length(fnames) else end <- end[1]
  fnames <- fnames[1:end]
  fnames <- fnames[!is.na(fnames)]
}

for (fname in fnames) {
  fname_exp <- paste0("R/", fname)
  cat(sep="", "Sourcing ", fname_exp, "...\n")
  source(fname_exp)
}

