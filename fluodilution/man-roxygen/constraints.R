#' Common constraints.
#' 
#' Those constraints are a useful starting point to constrain a model to make it more
#' identifiable or to capture salient features.
#' 
#' Those various constraints can be 
#' "combinatorially"
#' combined and their relative merits assessed with \code{\link[stats]{AIC}} (or quasi-AIC, see 
#' Pan [2001]).  The package \pkg{AICcmodavg} can be used to this end.  Moreover, practical 
#' identifiability can be 
#' assessed by looking at the 
#' \emph{effect}.  See \code{\link{fd_aictab}} for a wrapper around those two concepts.
#' 
#' The constraints can be combined together as such: \code{~ `#noss` + `#nodeath` + `#delta_1111`}
#' or using \code{\link{catcstr}}: \code{catcstr(`#noss`, `#nodeath`, `#delta_1111`)}.
#' 
#' They are fed to the hierarchical model through function \code{\link{fd_model}}.
#' 
#' @seealso \code{\link{constraints}}, \code{\link{fd_model}}, \code{\link{fd_aictab}}
#' 
#' @examples
#' attach(FdCommonConstraints)
#' @format The common constraints are usually loaded in the global environment using 
#' \code{\link[base]{attach}}.  Their names start by
#' convention with a '#', as such: `#constraint_name`.  See \code{\link{constraints}} for
#' the format of constraints in general.
#' 
#' @return
#' The value of \code{FdCommonConstraints} is:
#' 
#' \preformatted{
#' <% cat(paste(readLines("inst/contrib/FdCommonConstraints.R"), collapse="\r")) %>
#' }


