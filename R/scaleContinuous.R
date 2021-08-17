#' @title Standardize continuous variables.
#'
#' @description
#' Standardizes (centers and scales) continuous variable in a dataset, leaving discrete variables untouched.
#'
#' @usage
#' scaleContinuous(data, discrete, sdx = 1)
#'
#' @param data A data frame or matrix of variables to be scaled.
#' @param discrete Either a logical vector which is \code{TRUE} for discrete variables and \code{FALSE} for continuous ones or a vector of column numbers of discrete variables.
#' @param sdx The standard deviation of the columns for the continuous variables.
#'
#' @return A data frame with the same dimensions as \code{data} where the continuous variables are centered and scaled.
#'
#' @export
# DA 9/17/14: Added a pass-through path such that if the data are all discrete, the function just returns the original data.
scaleContinuous <-
function(data, discrete, sdx=1){
    allNum <- FALSE
    if(is.logical(discrete) & !any(discrete)){
        allNum <- TRUE
    }
    if(length(discrete) == 0){
        allNum <- TRUE
    }

    if(!allNum){
    	if(is.logical(discrete)){
    		discrete <- which(discrete)
    	}
        if(length(discrete) < ncol(data)){
        	num.ind <- (1:ncol(data))[-discrete]
        	nn.ind <- discrete
            if(length(num.ind) == 0){stop("No Numeric Variables in Data Frame")}
            num.dat <- data[,num.ind, drop=FALSE]
            nonnum.dat <- data[,nn.ind, drop=FALSE]
            if(is.null(dimnames(num.dat))){
            	nn <- names(data)[num.ind]
            	num.dat <- data.frame(num.dat)
            	names(num.dat) <- nn
            }
            if(is.null(dimnames(nonnum.dat))){
            	nnn <- names(data)[nn.ind]
            	nonnum.dat <- data.frame(nonnum.dat)
            	names(nonnum.dat) <- nnn
            }
            newdat <- cbind(as.data.frame(scale(num.dat)*sdx), as.data.frame(nonnum.dat))
            newdat <- newdat[, match(names(data), names(newdat))]
        }
        else{
            newdat <- data
        }
    }
    else{
        newdat <- scale(data)*sdx
    }
    class(newdat) <- c(class(newdat), "scaledDataFrame")
    return(newdat)
}
