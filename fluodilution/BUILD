load("@com_grail_rules_r//R:defs.bzl", "r_pkg", "r_unit_test", "r_pkg_test")

package(
    default_visibility = ["//visibility:public"],
    licenses = ["notice"],  # MIT
)

r_pkg(
    name = "fluodilution",
    srcs = glob(
        ["**"],
        exclude = ["BUILD"],
    ),
    lazy_data = True,
    deps = [
        "@R_Rcpp//:Rcpp",
        "@R_RcppArmadillo//:RcppArmadillo",
        "@R_RcppNumerical//:RcppNumerical",
        "@R_RcppParallel//:RcppParallel",
        "@R_dplyr//:dplyr",
        "@R_ggplot2//:ggplot2",
        "@R_magrittr//:magrittr",
        "@R_memoise//:memoise",
        "@R_plyr//:plyr",
        "@R_reshape2//:reshape2",
        "@R_scales//:scales",
        "@R_stringr//:stringr",
    ],
)

SUGGESTED_DEPS = [
    "@R_zoo//:zoo",
    "@R_lmtest//:lmtest",
    "@R_sandwich//:sandwich",
    "@R_testthat//:testthat",
    "@R_GenSA//:GenSA",
    "@R_truncnorm//:truncnorm",
    "@R_mvtnorm//:mvtnorm",
    "@R_limSolve//:limSolve",
    "@R_flowCore//:flowCore",
    "@R_flowWorkspace//:flowWorkspace",
    "@R_flowClust//:flowClust",
    #"@R_devtools//:devtools",
    "@R_AICcmodavg//:AICcmodavg",
    "@R_knitr//:knitr",
    "@R_gridExtra//:gridExtra",
    "@R_RColorBrewer//:RColorBrewer",
    "@R_numDeriv//:numDeriv",
    "@R_bbmle//:bbmle",
    "@R_microbenchmark//:microbenchmark",
    #"@R_roxygen2//:roxygen2",
]

r_unit_test(
    name = "test",
    timeout = "long",
    pkg = "fluodilution",
    suggested_deps = SUGGESTED_DEPS,
)

r_pkg_test(
    name = "check",
    timeout = "long",
    check_args = "--as-cran",
    pkg = "fluodilution",
    suggested_deps = SUGGESTED_DEPS,
)
