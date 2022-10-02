# KRIG ---------------------------------------------------------------------------------------------
#' @title Spatial Statistics with Kriging
#' @description Implements different methods for spatial statistics, in particular focused with 
#' Kriging based models. We count with different implemented models, simple, ordinary and 
#' universal forms of Kriging, co-Kriging and regression Kriging models. Includes, multivariate 
#' sensitivity analysis under an approximation designed over reproducing kernel Hilbert spaces and 
#' computation of Sobol indexes under this framework.
#' 
#' @examples 
#' library( KRIG )
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
#' \insertRef{Aronszajn:ThRK}{KRIG}
#' 
"_PACKAGE"
