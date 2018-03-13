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

release_questions <- function() {
  c(
    "Have you cleaned the package, especially of private matters?  You can do so by executing 'R CMD build --no-build-vignettes --no-manual --no-resave-data fluodilution' and looking at the '.tgz' result",
    "Have you checked for spelling errors (systematically do it, using, e.g., RStudio spell checking)?",
    "Have you checked coherence of words between vignette and manual?",
    "Have you checked coherence of parameter/function names in documentation?",
    "Have you compiled the vignettes?",
    "Have you generated the documentation?",
    "Have you checked that README.md and NEWS.md display properly?",
    "Have you run the linters (C++ and R)?",
    "Have you run the integration tests (Bazel tests)?"
  )
}
