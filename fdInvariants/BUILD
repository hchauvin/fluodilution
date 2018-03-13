load("@com_grail_rules_r//R:defs.bzl", "r_test")

package(licenses = ["notice"])  # MIT

r_test(
    name = "invitro_branching",
    src = "invitro_branching.R",
    deps = [
        "@R_assertthat//:assertthat",
        "@R_ggplot2//:ggplot2",
        "@R_magrittr//:magrittr",
        "@R_plyr//:plyr",
        "@R_reshape2//:reshape2",
        "@fluodilution//fluodilution",
    ],
)

r_test(
    name = "macr_branching",
    src = "macr_branching.R",
    flaky = True,
    deps = [
        "@R_assertthat//:assertthat",
        "@R_ggplot2//:ggplot2",
        "@R_magrittr//:magrittr",
        "@R_plyr//:plyr",
        "@R_reshape2//:reshape2",
        "@fluodilution//fluodilution",
    ],
)

r_test(
    name = "macr_cyton",
    size = "large",
    src = "macr_cyton.R",
    flaky = True,
    deps = [
        "@R_assertthat//:assertthat",
        "@R_ggplot2//:ggplot2",
        "@R_magrittr//:magrittr",
        "@R_plyr//:plyr",
        "@R_reshape2//:reshape2",
        "@fluodilution//fluodilution",
    ],
)

r_test(
    name = "invivo_branching",
    src = "invivo_branching.R",
    deps = [
        "@R_assertthat//:assertthat",
        "@R_ggplot2//:ggplot2",
        "@R_magrittr//:magrittr",
        "@R_plyr//:plyr",
        "@R_reshape2//:reshape2",
        "@fluodilution//fluodilution",
    ],
)

r_test(
    name = "fmm",
    src = "fmm.R",
    deps = [
        "@R_assertthat//:assertthat",
        "@R_ggplot2//:ggplot2",
        "@R_magrittr//:magrittr",
        "@R_reshape2//:reshape2",
        "@fluodilution//fluodilution",
    ],
)

r_test(
    name = "agent_based_af_bp",
    src = "agent_based_af_bp.R",
    flaky = True,
    deps = [
        "@R_assertthat//:assertthat",
        "@R_ggplot2//:ggplot2",
        "@R_magrittr//:magrittr",
        "@R_reshape2//:reshape2",
        "@fluodilution//fluodilution",
    ],
)

r_test(
    name = "agent_based_af_bp_grid",
    src = "agent_based_af_bp_grid.R",
    flaky = True,
    deps = [
        "@R_assertthat//:assertthat",
        "@R_ggplot2//:ggplot2",
        "@R_magrittr//:magrittr",
        "@R_reshape2//:reshape2",
        "@fluodilution//fluodilution",
    ],
)

r_test(
    name = "agent_based_af",
    src = "agent_based_af.R",
    deps = [
        "@R_assertthat//:assertthat",
        "@R_ggplot2//:ggplot2",
        "@R_magrittr//:magrittr",
        "@R_reshape2//:reshape2",
        "@fluodilution//fluodilution",
    ],
)

r_test(
    name = "optim",
    src = "optim.R",
    deps = [
        "@R_assertthat//:assertthat",
        "@fluodilution//fluodilution",
    ],
)

r_test(
    name = "optim_methods",
    src = "optim_methods.R",
    deps = [
        "@R_GenSA//:GenSA",
        "@R_assertthat//:assertthat",
        "@R_minpack_lm//:minpack.lm",
        "@R_nls2//:nls2",
        "@R_optimx//:optimx",
        "@fluodilution//fluodilution",
    ],
)
