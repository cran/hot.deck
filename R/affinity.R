#' @title Affinity Calculation.
#'
#' @description
#' Calculates affinity based on Cranmer and Gill (2013). The function performs the original method (as described in the article) and also
#' a method that takes into account the correlation structure of the observed data that increases efficiency in making matches.
#' Affinity is calculated by first identifying whether two observations are sufficiently \sQuote{close} on each variable.  Consider the target observation number 1.
#' If observation \emph{i} is close to the target observation on variable \emph{j}, then \code{A[i,j] = 1} otherwise, it equals zero.
#' Close for two discrete variables is defined by them taking on the same value.  Close for continuous variables is taking on a distance no greater than 1 from each other.
#' While this may seem restrictive and arbitrary, arguments exist in the main package function \code{hot.deck} that allows the user to set how many standard deviations equal a distance of 1 (with the \code{cutoffSD} argument.
#'
#' @usage affinity(data, index, column = NULL, R = NULL, weighted = FALSE)
#'
#' @param data A data frame or matrix of values for which affinity should be calculated.
#' @param index A row number identifying the target observation.  Affinity will be calculated between this observation and all others in the dataset.
#' @param column A column number identifying the variable with missing information.  This is only needed for the optional correlation-weighted affinity score. The correlation that is used is the correlation of all variables with the focus variable (i.e., the column).
#' @param R A correlation matrix for \code{data}.
#' @param weighted Logical indicating whether or not the correlation-weighted affinity measure should be used.
#'
#' @return A number of missing observation-variable combinations-by-number of observations in data matrix of affinity scores.
#'
#' @examples
#' data(D)
#' out <- hot.deck(D)
#'
#' @export
affinity <-
function(data, index, column=NULL, R = NULL, weighted=FALSE){
	tmp <- abs(data - matrix(data[index, ], nrow=nrow(data), ncol=ncol(data), byrow=TRUE)) < 1
	if(is.null(column) | is.null(R)){
		if(weighted)warning("Correlation Matrix and/or Column Number not provided, switching to Unweighted Affinity\n")
	}
	if(!weighted){
		affinity <- rowMeans(tmp, na.rm=TRUE)
	}
	if(weighted){
		tmp[which(is.na(tmp), arr.ind=T)] <- FALSE
		affinity <- tmp %*% R[ ,column]
	}
	affinity[index] <- 0
	affinity
}
