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

# Names used in non-standard evaluation:
# just so that R CMD CHECK does not complain
globalVariables(
  c(
    "Category", "Individual", "Label", "Size", "Count", "Replicates",
    "Time", "Cutoff", "Weight", "Inoculum", "y", "a", "b", "Proportion",
    "Generation", "Timepoint", "Weight", "Type", "i", "x", "mdl",
    ".", "Effect", "Modnames", "gen",
    "xmean"
  )
)

#' @title Generic functions for analysing fluorescence dilution experiments /
#'   proliferation assays.
#'
#' @description In Chauvin et al. (2016), we devise a framework to
#'   quantitatively analyse fluorescence dilution (Helaine et al. 2010,
#'   2014), a
#'   technique aimed at understanding the history of bacteria in terms of
#'   division, migration and death in an \emph{in vitro} setting, in macrophage
#'   cultures or in systemic infection of mice. This package builds upon and
#'   translates into open-source code previous efforts in immunology (Hyrien
#'   et al. 2008, Hyrien et al. 2010, Hawkins et al. 2006)
#'   and adapts those analytical techniques to the novel use case of
#'   microbiology.  Therefore, this package can be used by microbiologists and
#'   immunologists alike.
#'
#' @section Background: Proliferation assays make use of a fluorescent reporter,
#'   such as a dye like CFSE (Lyons and Parish 1994) or a fluorescence protein
#'   (e.g., GFP) expressed by a plasmid (Helaine 2010, Roostalu 2008). In both
#'   cases, the fluorophore is present in fixed quantities in progenitor
#'   cells/bacteria. Whenever a cell/bacterium divides, the fluorescence is
#'   approximately halved.  The distribution of fluorescence in the population
#'   can then be recorded over time using flow cytometry.  Augmented by cell
#'   counts / Colony Forming Units (Chauvin et al.
#'   2016) or
#'   the staining of dying/dead
#'   cells (Hyrien et al. 2010), the resulting dataset can be used to
#'   infer a
#'   set of \emph{biological constants} such as the growth/death rates, their
#'   variability, ...
#'
#' @section Aim: We took great care to cover a number of particular cases that
#'   might arise in fluorescence dilution and extensively documented our
#'   methods, both concerning their use and more broadly their philosophy.  Our
#'   aim was not only to provide a robust package that could be reused many
#'   times in the future and extended for new usages but also to transmit some
#'   important practical knowledge concerning fluorescence dilution and cell
#'   growth.  We believe that a package is the best way to embody this
#'   knowledge.
#'
#' @section A two-step hierarchical approach: An intermediary step before the
#'   estimate of the biological constants involves
#'   extracting the proportion of cells in any given \emph{number of
#'   generations} (a cell in generation 0 has not divided since the beginning of
#'   the experiment, a cell in generation 1 or 2 has divided once or twice,
#'   ...).  This is done through a Finite Mixture Model (FMM).  In immunology,
#'   with the exception of Hyrien et al. (2008) and Hyrien et al. (2010), this
#'   step has usually been glossed over as fluorescence dilution in monocytes
#'   produces distinct peaks in the distribution of fluorescence, one peak for
#'   each generational cohort.  In microbiology, the generational cohorts are
#'   not as clearly delineated, probably owing to more variability of the
#'   fluorescence in the bacterial inoculum than in the monocytes (Chauvin
#'   et al. 2016).  Thus, in microbiology, this first step requires a
#'   more thorough treatment.  In this respect, we implemented FMMs taking into
#'   account autofluorescence and the stochastic partitioning of fluorescent
#'   molecules upon division (binomial partitioning).  Our system is flexible
#'   enough for other FMMs to be easily "plugged-in".
#'
#'   In a second step, the proportions inferred by the FMM are then fed to a
#'   \emph{proliferation model} that links those proportions to the underlying
#'   proliferation/death mechanisms through \emph{biological constants}.  With
#'   this package, a large variety of proliferation models can be implemented.
#'   They all rest on two broad sets of assumptions. On one hand, a generalized
#'   Cyton model and its nested, simpler models, consider that division and
#'   death occur concurrently, in competition (Hawkins et al. 2007). On
#'   the other hand, a branching process (Hyrien et al. 2010) considers
#'   that birth and death have both a given probability of occurring in a cell
#'   (probability of transition) and a given probability of occurring after a
#'   certain time (time to transition).  Despite these apparently distinct
#'   behaviour, Hyrien et al. (2010) have shown that a Cyton model is
#'   actually a special case of a branching process with very specific
#'   probabilities of transition (see also Chauvin et al. (2016)).
#'   However, as those probabilities are difficult to calculate, we decided to
#'   implement the Cyton model independently.
#'
#'   As we have seen, a variety of models in this two-step hierarchical approach
#'   can be used.  Finite mixture models and proliferation models can be used in
#'   any combination imaginable.  Indeed, we believe that not all models are
#'   good for all datasets.  After 15 years of literature and countless models
#'   produced, the question of model selection is therefore central.  Our
#'   integrated approach allows to test novel and proven models alike in a quick
#'   way so as to find an acceptable fit.  Criteria researchers should look for
#'   include some information criterion such as the Akaike Information Criterion
#'   (AIC) to balance out better model fit with parcimony, whether the model is
#'   practically identifiable and whether heroskedasticity and systematic bias,
#'   across time or experiments, persist in an overly problematic way.
#'
#' @section Overview of the implementation: Loading of experimental data from a
#'   flowJo workspace (\url{www.flowjo.org}) can be carried out with
#'   \code{\link{fd_create}}.
#'   Alternatively, an \code{\link{fd_data}} object, which is the form the
#'   datasets take in this package, can be generated directly.  A hierarchical
#'   model can be specified through \code{\link{fd_model}}.  The resulting model
#'   can be fed to generic nonlinear optimization algorithms such as
#'   \code{\link[stats]{nls}}, \code{\link[nlme]{nlsList}},
#'   \code{\link[nlme]{nlme}} or \code{\link[nlme]{gnls}}.  We also implemented
#'   the wrapper \code{\link{fd_nls}} around a generic nonlinear optimizer.
#'   This wrapper automatically sets the boundary and starting values from an
#'   \code{\link{fd_model}} object and performs additional checks, especially on
#'   the maximum number of generations considered by the algorithm (even if the
#'   number of generations one can follow using fluorescence dilution is
#'   limited, setting a maximum number of generations too low can have
#'   consequences on the fitting of the cell counts and it is absolutely
#'   necessary to check for that).
#'
#' @section Strategy of nonlinear optimization: We found out, echoing Miao
#'   et al. (2012), that traditional fitting algorithms perform poorly on
#'   fluorescence dilution problems (Chauvin et al. 2016).  However,
#'   generalized simulated annealing (a stochastic global search approach),
#'   implemented in package
#'   \pkg{GenSA} (Xiang 2013), is promising. However, \code{\link[GenSA]{GenSA}}
#'   is a general optimizer and cannot be directly fed a nonlinear formula of
#'   the form \code{y ~ f(x)} as provided by \code{\link{fd_formula}}.
#'   Therefore, we provide the function \code{\link{fd_minuslogl}} that, instead
#'   of returning a formula, returns a minus log-likelihood function that can be
#'   directly fed to \code{GenSA}. To get confidence intervals and other such
#'   niceties, the user is invited to use \code{GenSA} through the package
#'   \pkg{bbmle}, an extension of package \pkg{stats4} that allows an arbitrary
#'   function as a nonlinear optimizer.  Alternatively, the user can use the
#'   nonlinear square wrapper \code{nlsSA}, not part of this package (the focus
#'   here is on fluorescence dilution, not nonlinear optimization) but included
#'   on the side in \code{contrib/nlsSA.R} (see examples in
#'   \code{\link{fd_nls}}). \code{nlsSA} borrows code from \code{nls} but
#'   extends it so that any nonlinear optimizer, whether specifically tailored
#'   to least-square problems (e.g., Levenberg-Marquardt) or general (e.g.,
#'   BFGS), can be used.  (The \pkg{stats} package, from which \code{nls}
#'   originates, has been authored by the R Core Team and contributors
#'   worldwide.)  For generalized nonlinear square (\link[nlme]{gnls})
#'   and nonlinear mixed-effects fitting (\link[nlme]{nlme}), a global
#'   optimization strategy is (as of now) still not available.  Therefore, the
#'   user is invited to use the result of \code{nlsSA} as starting value for a
#'   round of \code{gnls} or \code{nlme}.  Wrappers are given in
#'   \code{contrib/fitwrappers.R}.
#'
#'   We believe that generalized simulated annealing provides a very flexible
#'   way to deal with nonlinear regressions with overparametrization or poor
#'   choices of starting values.  Other possibilities in R include self-starting
#'   models and initial grid search using brute force (package \pkg{nls2}).
#'   Although there is no straightforward way to self-start an FD model, brute
#'   force could be a viable option, but generalized simulated annealing offers
#'   in our opinion a more elegant way to go about.  Moreover, we found out that
#'   \pkg{nls2} performed poorly in our benchmarking (Chauvin et al.
#'   2016). In any case, neither self-starting nor \pkg{nls2} deal with
#'   overparametrization (in the worst case, they fail after encountering a
#'   singular gradient). Now, overparametrization is not good either for
#'   \code{GenSA} that will struggle to find an optimum (slow convergence, if
#'   any), but at least the results of \code{GenSA} can be used for further
#'   improvements to the model.
#'
#' @section Generalized Estimating Equations: Notice that as a fluorescence
#'   dilution is not an ordinary nonlinear regression but a generalized
#'   estimating equation (GEE) with potential clustering (repeat experiments)
#'   and "space" and time autocorrelations, the degree of freedom is difficult
#'   to assess.  Therefore, correcting for small-sample bias using the
#'   \code{\link[AICcmodavg]{AICc}} instead of the \code{\link[stats]{AIC}}
#'   information criterion (or QIC, see Pan [2001]) should be done with caution.
#'   If the correction can be done with the number of clusters, as one usually
#'   does with logistic growth curves in the context of nonlinear mixed-effects
#'   models, the AIC will be overestimated.
#'   If the correction is made taking the number of residuals, the AIC
#'   will be underestimated.
#'
#'   Moreover, if a precise determination of the confidence intervals is
#'   required, the user is advised, when possible, to perform a parametric
#'   boostrap for one-time experiments (Hyrien 2008) or a clustered bootstrap or
#'   jack-knife when the data is clustered in the context of repeat experiments
#'   (Chauvin et al. 2016). Sometimes, a sandwich correction is enough
#'   (see the \pkg{sandwich} package) and we accommodate the use of such
#'   "nonstandard" variance-covariance matrices.
#'
#' @section Representing, loading and simulating datasets: \describe{
#'   \item{\code{\link{fd_data}}}{Represent a fluorescence dilution dataset.}
#'   \item{\code{\link{fd_create}}}{Create a histogram-based dataset from row
#'   fluorescence data.} \item{\code{\link{fd_simulate}}}{Simulate a dataset.} }
#'
#' @section Available finite mixture models (FMM): \describe{
#'   \item{\code{\link{fd_fmm_gaussian}}}{Log-normal cohorts.}
#'   \item{\code{\link{fd_fmm_af}}}{Correction for autofluorescence only.}
#'   \item{\code{\link{fd_fmm_af_bp}}}{Correction for autofluorescence and
#'   binomial partitioning.} }
#'
#' @section Available proliferation models: \describe{
#'   \item{\code{\link{fd_proliferation_cyton}}}{A fairly general Cyton model.}
#'   \item{\code{\link{fd_proliferation_branching}}}{A fairly general branching
#'   process.} }
#'
#' @section Optimization: \describe{ \item{\code{\link{fd_model}}}{Two-tier
#'   hierarchical fluorescence dilution model.}
#'   \item{\code{\link{fd_nls}}}{Nonlinear least square optimizer.}
#'   \item{\code{\link{fd_minuslogl}}}{Minus log-likelihood for more general
#'   optimization.} }
#'
#' @section Agent-based model: In \code{contrib/agent.R}, we implemented an
#' agent-based model.  It can be used to explore the validity of our various
#' assumptions used for the proliferation models and to factor
#' stochasticity in.  However, this model should not be used for "routine"
#' fitting and as a consequence is not formally included in the package.
#'
#' @references Hawkins ED, Turner ML, Dowling MR, van Gend C, Hodgkin PD (2007).
#' A model of immune regulation as a consequence of randomized lymphocyte
#' division and death times. \emph{Proc Natl Acad Sci USA} \strong{104} (12):
#' 5032-5037.
#'
#' Hawkins ED, Markham JF, McGuinness LP, Hodgkin PD (2006). A single-cell
#' pedigree analysis of alternative stochastic lymphocyte fates. \emph{Proc Natl
#' Acad Sci USA} \strong{106} (32): 13457-13462.
#'
#' Helaine S, Thompson JA, Watson KG, Liu M, Boyle C, and Holden DW (2010).
#' Dynamics of intracellular bacterial replication at the single cell level.
#' \emph{Proc Natl Acad Sci USA} \strong{107} (8): 3746-3751.
#'
#' Helaine S, Chverton AM, Watson KG, Faure LM, Matthews SA, Holden DW (2014).
#' Internalization of Salmonella by macrophages induces formation of
#' nonreplicating persisters. \emph{Science} \strong{343}: 204-8.
#'
#' Hyrien O, Zand MS (2008). A Mixture Model with Dependent Observations for the
#' Analysis of CSFE-Labeling Experiments. \emph{Journal of the American
#' Statistical Association} \strong{103}: 222-239.
#'
#' Hyrien O, Chen R, Zand MS (2010). An age-dependent branching process model
#' for the analysis of CFSE-labeling experiments. \emph{Biology Direct}
#' \strong{5} (41).
#'
#' Lyons AB, Parish CR (1994). Determination of lymphocyte division by flow
#' cytometry. \emph{Journal of Immunological Methods} \strong{171} (1): 131-137.
#'
#' Miao H, Jin X, Perelson AS, Wu H (2012). Evaluation of multitype mathematical
#' models for CFSE-labeling experiment data. \emph{Bull Math Biol} \strong{74}
#' (2): 300-326.
#'
#' Pan W (2001). Akaike's Information Criterion in Generalized Estimating
#' Equations. \emph{Biometrics} \strong{57} (1): 120-125.
#'
#' Roostalu J, Joers A, Luidalepp H, Kaldalu N, Tenson T (2008). Cell division
#' in Escherichia coli cultures monitored at single cell resolution. \emph{BMC
#' Microbiology} \strong{8} (68).
#'
#' Xiang Y, Gubian S, Suomela B, Hoeng J (2013). Generalized Simulated Annealing
#' for Global Optimization: The GenSA Package. \emph{The R Journal} \strong{5}
#' (1): 13-28.
#'
#' @docType package
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom utils relist
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom scales percent_format identity_trans log10_trans log1p_trans
#' @importFrom nlme getGroupsFormula groupedData
#' @import ggplot2
#' @importFrom graphics hist
#' @import stats
#' @importFrom plyr rbind.fill summarise ddply .
#' @importFrom dplyr select failwith
#' @importFrom memoise memoise
#' @importFrom stringr str_sub
#' @importFrom reshape2 melt
#' @useDynLib fluodilution
"_PACKAGE"

#' @template constraints
#' @export
FdCommonConstraints <-
  eval(parse("inst/contrib/FdCommonConstraints.R"))

#' Example dataset for one-compartment fluorescence dilution.
#'
#' In this experiment, primary BMM extracted from C57 B1/6 mice were infected by
#' wild-type 12023s \emph{Salmonella} Typhimurium harbouring the pFCcGi plasmid,
#' a plasmid with a constitutive reporter helping in bacterial recognition and
#' an inducible reporter used for fluorescence dilution.  The data was acquired
#' on a LSRFortessa cytometer (DB).  Raw data is available upon request.
#'
#' This dataset can be best fitted with a Gaussian finite mixture model
#' (\code{\link{fd_fmm_gaussian}}) and a branching model
#' (\code{\link{fd_proliferation_branching}}) with the common constraints\verb{
#' }\code{~ `#noss` + `#nodeath` + `#delta_1100` + `#mm_xyxz`} (see
#' \code{\link{FdCommonConstraints}}).  As this was a repeat experiment (the
#' field \code{Individual} gives an identifier for the experiment), it is also
#' possible to explore the use of mixed-effects models \emph{via}
#' \code{\link[nlme]{nlme}}.
#'
#' @usage data(FdSTyphimuriumWTC57)
#' @format This is an \code{\link{fd_data}} with the FMM and proliferation
#'   parameters preset using attributes.
#' @source Sophie Helaine, Jessica A. Thompson, Kathryn G. Watson, Mei Liu,
#'   Cliona Boyle, and David W. Holden (2010).  Dynamics of intracellular
#'   bacterial replication at the single cell level.  \emph{Proc Natl Acad Sci
#'   USA.} \strong{107} (8):3746-3751.
#' @seealso \code{\link{fd_data}} for data format, \code{\link{fd_model}} for an
#'   example of use.
#' @examples
#' data(FdSTyphimuriumWTC57)
#' plot(FdSTyphimuriumWTC57)
"FdSTyphimuriumWTC57"

#' Artificial example dataset for one-compartment fluorescence dilution.
#'
#' This dataset was generated by \code{\link{fd_simulate}} and features clear,
#' distinct peaks such as those found in immunological applications.
#'
#' This dataset can be best fitted with a Gaussian finite mixture model
#' (\code{\link{fd_fmm_gaussian}}) and a branching model
#' (\code{\link{fd_proliferation_branching}}).
#'
#' @usage data(FdClearPeaks)
#' @format This is an \code{\link{fd_data}} with the FMM and proliferation
#'   parameters preset using attributes.
#' @seealso \code{\link{fd_data}} for data format, \code{\link{fd_model}} for an
#'   example of use.
#' @examples
#' data(FdClearPeaks)
#' plot(FdClearPeaks)
"FdClearPeaks"