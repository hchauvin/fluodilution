# fdInvariantScripts

Various `r_test` Bazel rules that provide some additional tests for the `fluodilution` package.
Think of the tests executed with `R CMD check` as unit tests, and tests in this directory as
integration tests.

One way to test the validity of our models is to compare the results to a completely different
model and see if there is a match.  Here, we compare the proliferation and FMM models to
agent-based, or particle models where we randomly generate single-cell fluorescence and division
behavior and compare the populations to their "limit" representation.