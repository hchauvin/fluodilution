# fluodilution

[![Build status](https://badge.buildkite.com/1045d20ee21e2f2e21bd186af81829676e08593f8735001379.svg)](https://buildkite.com/hchauvin/fluodilution)

In biology, fluorescence Dilution (FD) experiments / proliferation
assays are used to assess growth and migration at the single-cell level.
Mathematical models can be used to extract parameters such as the division,
death and migration rates and their variability.
This package provides generic functions to manipulate FD data and
solve for these parameters.  It can be used both in immunology
and microbiology.

## Documentation

Full documentation is available [here](fluodilution/inst/doc/fluodilution.pdf).  An example 
workflow can be found [here](fluodilution/inst/doc/howTo-fluodilution.pdf) and a performance
benchmark [here](fluodilution/inst/doc/performance.pdf).  Details are in 
our upcoming article (Chauvin et al. 2018).

## Installation

The development version can be installed from github:

```R
if (!require(devtools))
    install.packages("devtools")

devtools::install_github("hchauvin/fluodilution", subdir="fluodilution")

# To install all the suggested packages as well:
devtools::install_github("hchauvin/fluodilution", 
                         subdir="fluodilution",
                         dependencies=TRUE)
```

A release version will be available soon on 
[bioconductor](http://www.bioconductor.org).

## Development

The `fluodilution` package can be consumed either using `devtools`, `bioconductor` or
[bazel/rules_r](https://github.com/grailbio/rules_r), which is a set of
[Bazel](https://bazel.build) rules for the R statistical language.

With Bazel, the external dependencies of the project are described in [`WORKSPACE`](./WORKSPACE)
and [`cran`](./cran).  Additional test scripts are available in [`fdInvariants`](./fdInvariants).
These scripts are run with `bazel test` and, among other things, test the validity of our
distribution models against agent-based simulations.  Because they are more expensive to run
than your usual R package tests, they are not run when `R CMD check` is invoked.

