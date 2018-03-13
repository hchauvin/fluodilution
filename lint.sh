#!/bin/bash

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

set -eou pipefail

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

function header() {
  echo "========================================================================"
  echo "$1"
  echo "========================================================================"
}

function test_copyright() {
  header "COPYRIGHT"
  files=$(grep -rL "Copyright" "$DIR" \
    --include '*.go' \
    --include '*.R' \
    --include '*.py' \
    --include '*.bzl' \
    --include '*.h' \
    --include '*.hpp' \
    --include '*.cc' \
    --include '*.c' \
    --include '*.cpp' \
    --include '*.sh' \
    --include 'BUILD.bazel' \
    --include 'BUILD' \
    --exclude-dir 'bazel-*' \
    --exclude-dir 'inst' \
    --exclude-dir 'man-*' \
    --exclude 'RcppExports.R' \
    --exclude 'RcppExports.cpp' \
    --exclude-dir 'vendor')
  ok=1
  for i in $files; do
    if [ -s "$i" ]; then
      # File is not empty
      echo "$i"
      ok=0
    fi
  done
  if [ $ok == 0 ]; then
    echo "==> Those files lack a copyright notice"
    return 1
  fi
}

function suggest_go() {
  header "GO (SUGGEST)"
  # This is not enforced, just given here for information.
  if hash golint 2>/dev/null; then
    golint "$DIR/..."
  else
    echo "(golint not installed: skipped)"
  fi
}

function test_bazel_config() {
  header "BAZEL CONFIG"
  for i in BUILD WORKSPACE '*.bzl'; do
    files=$(find "$DIR" -name "$i" -type f -not -path '**/bazel-*')
    if [ -n "$files" ]; then
      buildifier -mode=check $files
    fi
  done
}

test_copyright
suggest_go
test_bazel_config