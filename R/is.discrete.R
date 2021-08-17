#' @title Identify whether variables are discrete or continuous.
#'
#' @description
#' Variables are considered discrete if they have fewer unique, non-missing values than \code{cutoff} or they are factors.  Otherwise, variables are considered continuous.
#'
#' @usage
#' is.discrete(data, cutoff = 10)
#'
#' @param data A data frame, matrix or vector of values to be evaluated.
#' @param cutoff A numeric scalar identifying the cutoff relative to the number of unique, non-missing values for \sQuote{discreteness}.
#'
#' @return A logical vector indicating whether variables are discrete (\code{TRUE}) or continuous \code{FALSE}.
#'
#' @export
is.discrete <-
function(data, cutoff = 10){
    discrete <- vector()
    if(is.vector(data)){
    data <- na.omit(data)
    discrete <- length(unique(data)) <= cutoff | is.factor(data)
    }
    if(is.data.frame(data) | is.matrix(data)){
        for(j in 1:ncol(data)){
            tmp <- na.omit(data[,j])
            discrete[[j]] <- length(unique(tmp)) <= cutoff | is.factor(tmp)
        }
    }
    return(discrete)
}
