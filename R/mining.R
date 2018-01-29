#' Copper mining data
#' @docType data
#' @description This is reproduced from the original description for the dataset. 
#' A simulation based on a stockpile of mined material in the former Soviet Union. Boreholes have 
#' been drilled into the dump. The drill core is cut every 5 metres and assayed for copper and 
#' cobalt content in percentage by weight. This is the only three dimensional set of tutorial data. 
#' Coordinates are in metres.
#' @keywords datasets, mining, copper
#' @references Dr. Isobel Clark
#' @source \url{www.edumine.com}.
#' @examples 
#' data( 'Copper', packages = 'KRIG' )
#' @format A data table with 442 and 7 columns.
#' \describe{
#'   \item{a}{cod id}
#'   \item{s}{sample id}
#'   \item{x1}{first coordinate}
#'   \item{x2}{second coordinate}
#'   \item{x3}{third coordinate}
#'   \item{Z}{copper grade}
#'   \item{C}{value}
#' }
"Copper"  
