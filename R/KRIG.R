# KRIG ---------------------------------------------------------------------------------------------
#' @title Spatial statistics with kriging.
#' @description Implements different methods for spatial statistics, in particular focused with 
#' kriging based models. We count with different implemented models, simple, ordinary and 
#' universal forms of kriging, cokriging and regression kriging models. Includes, multivariate 
#' sensitivity analysis under an approximation designed over Reproducing Kernel Hilbert Spaces and 
#' computation of Sobol indexes under this framework.
#' 
#' The linear algebra operations are supported by RcppArmadillo.
#' @examples 
#' library( KRIG )
#' vignette( topic = 'simple_kriging', package = 'KRIG' )
#' vignette( topic = 'ordinary_kriging', package = 'KRIG' )
#' vignette( topic = 'universal_kriging', package = 'KRIG' )
#' vignette( topic = 'copper_mining_2d', package = 'KRIG' )
#' 
#' @importFrom Rdpack reprompt
#' 
#' @references
#' \insertRef{MatheronTH1970}{KRIG}
#' 
#' \insertRef{Chiles2011Geos}{KRIG}
#' 
#' \insertRef{Kanova:2013}{KRIG}
#' 
#' \insertRef{Rasmussen2005GPML}{KRIG}
#' 
#' \insertRef{SteinInterSpatial}{KRIG}
#' 
#' \insertRef{Clark1979Geo}{KRIG}
#' 
"_PACKAGE"
