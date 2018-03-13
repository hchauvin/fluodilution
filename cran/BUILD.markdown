load("@com_grail_rules_r//R:defs.bzl", "r_pkg")

package(default_visibility = ["//visibility:public"])

r_pkg(
    name = "markdown",
    srcs = glob(
        ["**"],
        exclude = [],
    ),
    deps = [
      "@R_mime//:mime",
    ],
)
