#' @title Multiple Hot Deck Imputation.
#'
#' @description
#' This function performs multiple hot deck imputation on an input data frame with missing observations using either the \dQuote{best cell} method (default) or
#' the \dQuote{probabilistic draw} method as described in Cranmer and Gill (2013). This technique is best suited for missingness in discrete variables, though it also performs well on continuous missing data.
#'
#' @usage hot.deck(data, m = 5, method = c("best.cell", "p.draw"), cutoff = 10, sdCutoff = 1,
#' optimizeSD = FALSE, optimStep = 0.1, optimStop = 5, weightedAffinity = FALSE,
#' impContinuous = c("HD", "mice"), IDvars = NULL, ...)
#'
#' @param data A data frame with missing values to be imputed using multiple hot deck imputation.
#' @param m Number of imputed datasets required.
#' @param method Method used to draw donors based on affinity either \dQuote{best.cell} (the default) or \dQuote{p.draw} for probabilistic draw.
#' @param cutoff A numeric scalar such that any variable with fewer than \code{cutoff} unique non-missing values will be considered discrete and necessarily imputed with hot deck imputation.
#' @param sdCutoff Number of standard deviations between observations such that observations fewer than \code{sdCutoff} standard deviations away from each other are considered sufficiently close to be a match, otherwise they are considered too far away to be a match.
#' @param optimizeSD Logical indicating whether the \code{sdCutoff} parameter should be optimized such that the smallest possible value is chosen that produces no thin cells from which to draw donors.  Thin cells are those where the number of donors is less than \code{m}.
#' @param optimStep The size of the steps in the optimization if \code{optimizeSD} is \code{TRUE}.
#' @param optimStop The value at which optimization should stop if it has not already found a value that produces no thin cells.  If this value is reached and thin cells still exist, a warning will be returned, though the routine will continue using \code{optimStop} as \code{sdCutoff}.
#' @param weightedAffinity Logical indicating whether a correlation-weighted affinity score should be used.
#' @param impContinuous Character string indicating how continuous missing data should be imputed.  Valid options are \dQuote{HD} (the default) in which case hot-deck imputation will be used, or \dQuote{mice} in which case multiple imputation by chained equations will be used.
#' @param IDvars A character vector of variable names not to be used in the imputation, but to be included in the final imputed datasets.
#' @param ... Optional additional arguments to be passed down to the \code{mice} routine.
#'
#' @return The output is a list with the following elements:
#' \itemize{
#'   \item{data}{An object of class \code{mi} which contains \code{m} imputed datasets.}
#'   \item{affinity}{A matrix of affinity scores see \code{\link{affinity}}.}
#'   \item{donors}{A list of donors for each missing observation based on the affinity score.}
#'   \item{draws}{The \code{m} observations drawn from donors that were used for the multiple imputations.}
#'   \item{max.emp.aff}{Normalization constant for each row of affinity scores; the maximum possible value of the affinity scores if correlation-weighting is used.}
#'   \item{max.the.aff}{Normalization constant for each row of affinity scores; the number of columns in the original data. }
#'   }
#'
#' @importFrom stats aggregate as.formula coef cor na.omit
#' @examples
#' data(D)
#' hot.deck(D)
#'
#' @export
# DA 9/10/14: Added argument for IDvars (ID variables that are not to be used in imputation, but should remain in data)
hot.deck <-
function(data, m = 5, method=c("best.cell", "p.draw"), cutoff=10, sdCutoff=1, optimizeSD = FALSE,
    optimStep = 0.1, optimStop = 5, weightedAffinity = FALSE, impContinuous = c("HD", "mice"),
    IDvars = NULL, ...){
	method <- match.arg(method)
    impContinuous <- match.arg(impContinuous)
# DA 9/15/14 Added warning about weighted affinity calculations and correlations among or with categorical variables.
    if(weightedAffinity){
        warning("Affinity calculations made as a function of pearson correlations among variables coerced to class 'numeric'\ntake care when using this on categorical, especially nominal variables")
    }
# DA 9/10/14: If IDvars is specified, remove them from the data and save in a different file.
    if(!is.null(IDvars)){
        IDdata <- data[, which(names(data) %in% IDvars), drop=FALSE]
        data <- data[,-which(names(data) %in% IDvars), drop=FALSE]
# DA 9/10/14: Added code to remove any observations that are always missing
        allNA <- apply(data, 1, function(x)all(is.na(x)))
        if(any(allNA)){
            IDdata <- IDdata[-which(allNA), , drop=FALSE]
            data <- data[-which(allNA), , drop=FALSE]
        }
    }
    else{
        allNA <- apply(data, 1, function(x)all(is.na(x)))
        if(any(allNA)){
            data <- data[-which(allNA), , drop=FALSE]
        }
    }
    if(any(allNA)){
        warning(paste(sum(allNA), " observations with no observed data.  These observations were removed\n", sep="") )
    }
	facs <- sapply(1:ncol(data), function(x)is.factor(data[,x]))
	disc.miss <- which(is.discrete(data, cutoff) & apply(data, 2, function(x)any(is.na(x))))
	alldisc <- is.discrete(data, cutoff)
	allmiss <- which(is.na(data), arr.ind=TRUE)
	cont.miss <- allmiss[-which(allmiss[,2] %in% disc.miss), ]
# DA 10/2/14: moved the warning here in response to Toby's problem where the error was getting tripped
# even though there was no continuous data.  This goes a step further and doesn't trip the warning unless
# there is any continuous data with missing observations.
    if(impContinuous == "HD" & method == "p.draw" & length(cont.miss) > 0){
        stop("Hot Deck imputation of continuous values can only be used with the best cell method\n")
    }
	whichna <- which(is.na(data), arr.ind=TRUE)
    if(impContinuous == "mice"){
	    whichna <- whichna[which(whichna[,2] %in% disc.miss), ]
    }
# DA 9/5/14: added condition any(!alldisc) so that optimization only happens when there are continuous variables on which to optimize
    if(optimizeSD & any(!alldisc)){
        mm <- 0
        while(sdCutoff <= optimStop & mm < m){
            tmp <- scaleContinuous(data, alldisc, sdx=1/sdCutoff)
            numdata <- sapply(1:ncol(tmp), function(i)as.numeric(tmp[,i]))
            R <- abs(cor(numdata, use="pairwise"))
            diag(R) <- 0
# DA 9/5/14: Commented out 2 lines below because I changed the looping mechanism for the affinity calculation
            # aff <- t(apply(whichna, 1, function(x)affinity(numdata, x[1], x[2], R, weightedAffinity)))
            # aff[which(!is.finite(aff), arr.ind=TRUE)] <- 0

# DA 9/5/14: changed the way affinity is calculated so no duplicates are calculated.  This should speed
# up computation particularly when there are many observations missing on many variables.  This also makes
# the weighted measure potentially much slower
            unnaobs <- unique(whichna[,1])
            if(!weightedAffinity){
        	    aff <- t(sapply(unnaobs, function(x)affinity(numdata, x, weighted=FALSE)))
                aff <- aff[match(whichna[,1], unnaobs), ]
            }
            if(weightedAffinity){
            	aff <- t(apply(whichna, 1, function(x)affinity(numdata, x[1], x[2], R, weightedAffinity)))
            }
        	if(any(!is.finite(aff))){aff[which(!is.finite(aff), arr.ind=TRUE)] <- 0}

# DA 9/5/14: added the following 4 lines to ensure that only valid donors (i.e., those with observed values) have
# non-zero affinity scores.
            wnadat <- matrix(1, nrow=nrow(data), ncol=ncol(data))
            wnadat[which(is.na(data), arr.ind=TRUE)] <- 0
            wnadat <- t(wnadat[, whichna[,2]])
            aff <- aff*wnadat
            w <- apply(aff, 1, function(x)which(x == max(x)))
            donors <- lapply(1:nrow(whichna), function(x)na.omit(data[w[[x]], whichna[x,2]]))
        	matches <- sapply(donors, length)
            mm <- min(matches)
            cat("SD Cutoff = ", sprintf("%.2f", sdCutoff), ", # Thin Cells = ", sum(matches < m), "\n", sep="")
            if(mm < m & sdCutoff == optimStop){warning(paste("Optimization unsuccessful, ", sum(matches < m), " thin cells remain with SD cutoff of ", sdCutoff, "\n", sep=""))}
            if(sdCutoff < optimStop){sdCutoff <- sdCutoff + optimStep}
        }

    }
# DA 9/10/14: changed result of scaleContinuous here to tmp from data so that the draws for the donors will not come from the scaled, but from the unscaled data.
    tmp <- scaleContinuous(data, alldisc, sdx=1/sdCutoff)
	numdata <- sapply(1:ncol(tmp), function(i)as.numeric(tmp[,i]))
	R <- abs(cor(numdata, use="pairwise"))
	diag(R) <- 0
	max.emp.aff <- 	apply(R, 2, sum)[whichna[,2]] # new
	max.the.aff <- rep(dim(R)[2] - 1, nrow(whichna)) # new

# DA 9/5/14: Commented out 2 lines below because I changed the looping mechanism for the affinity calculation
    # aff <- t(apply(whichna, 1, function(x)affinity(numdata, x[1], x[2], R, weightedAffinity)))
    # aff[which(!is.finite(aff), arr.ind=TRUE)] <- 0

# DA 9/5/14: changed the way affinity is calculated so no duplicates are calculated.  This should speed
# up computation particularly when there are many observations missing on many variables.  This also makes
# the weighted measure potentially much slower
    unnaobs <- unique(whichna[,1])
    if(!weightedAffinity){
	    aff <- t(sapply(unnaobs, function(x)affinity(numdata, x, weighted=FALSE)))
        aff <- aff[match(whichna[,1], unnaobs), ]
    }
    if(weightedAffinity){
    	aff <- t(apply(whichna, 1, function(x)affinity(numdata, x[1], x[2], R, weightedAffinity)))
    }
	if(any(!is.finite(aff))){aff[which(!is.finite(aff), arr.ind=TRUE)] <- 0}

# DA 9/5/14: added the following 4 lines to ensure that only valid donors (i.e., those with observed values) have
# non-zero affinity scores.
    wnadat <- matrix(1, nrow=nrow(data), ncol=ncol(data))
    wnadat[which(is.na(data), arr.ind=TRUE)] <- 0
    wnadat <- t(wnadat[, whichna[,2]])
    aff <- aff*wnadat
	if(method == "best.cell"){
		w <- apply(aff, 1, function(x)which(x == max(x)))
		donors <- lapply(1:nrow(whichna), function(x)na.omit(data[w[[x]], whichna[x,2]]))
		matches <- sapply(donors, length)
		if(any(matches < m)){warning(paste(sum(matches < m ), " of ", length(matches), " imputations with # donors < ", m, ", consider increasing sdCutoff or using method='p.draw'\n", sep=""))}
		repl <- ifelse(matches < m, TRUE, FALSE)
	    draws <- lapply(1:length(donors), function(x)sample(donors[[x]], m, replace=repl[x]))
	}
	if(method == "p.draw"){
		donors <- lapply(1:nrow(whichna), function(x)aggregate(aff[x, ], list(data[, whichna[x,2]]), mean, na.rm=TRUE))
	    draws <- lapply(1:length(donors), function(x)sample(donors[[x]][,1], m, replace=TRUE, prob=donors[[x]][,2]))
	}
	res <- vector(mode="list", length=m)
	inp.D <- lapply(1:m, function(x)data)
	for(md in 1:m){
		for(i in 1:nrow(whichna)){
			inp.D[[md]][whichna[i,1], whichna[i,2]] <- draws[[i]][md]
		}
		if(length(cont.miss) > 0 & impContinuous == "mice"){
			mice.D <- mice::mice(inp.D[[md]], m = 1, ...)
			res[[md]] <- tidyr::complete(mice.D)
		}
		else{
			res[[md]] <- inp.D[[md]]
		}
# DA 9/10/14: added three lines to put ID variables back in dataset
        if(!is.null(IDvars)){
            res[[md]] <- cbind(IDdata, res[[md]])
        }
	}
	class(res)  <- c("mi","list")
	return(list(data = res, affinity = aff, donors = donors, draws = draws, max.emp.aff = max.emp.aff, max.the.aff = max.the.aff))
}
