#' @title Convert hot.deck output to Amelia format.
#'
#' @description
#' Converts the output from hot.deck to the format used by Amelia for use with the Zelig package.
#'
#' @usage hd2amelia(object)
#'
#' @param object Output from a run of the \code{hot.deck} function.
#'
#' @return An object of class \dQuote{amelia} that can be used with Zelig.
#'
#' @export
hd2amelia <-
function(object){
    class(object) <- "amelia"
    names(object)[which(names(object) == "data")] <- "imputations"
    return(object)
}
