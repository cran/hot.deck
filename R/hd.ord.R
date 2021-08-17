#' @title Implement hot deck multiple imputation with ordinal variables.
#'
#' @description
#' This function adapts the \dQuote{hot.deck} function to impute data with missing observations by specifically accounting for ordinal variables.
#' The ordinal variable is regressed on specified meaningful explanatory variables with the \code{polr} ordered probit approach.
#' The approach assumes an underlying latent continuous variable and estimates the distances between ordinal variable categories.
#' Ordinal levels are replaced with mid-cutpoints of the newly estimated intercepts. Categories that are not supported by the data are dropped.
#' The resulting categories are used to impute the data with multiple hot deck imputation with either the \dQuote{best cell} method (default) or the \dQuote{probabilistic draw} method.
#' Any number of ordinal variables can be specified. The specified ordinal variables must not contain missing values.
#'
#' @usage hd.ord(data, ord, evs, m = 5, method=c("best.cell", "p.draw"),
#' cutoff=10, sdCutoff=1, optimizeSD = FALSE, optimStep = 0.1, optimStop = 5,
#' weightedAffinity = FALSE, impContinuous = c("HD", "mice"), IDvars = NULL, ...)
#'
#' @param data A data frame with missing values to be imputed using multiple hot deck imputation.
#' @param ord A vector of ordinal variables to be used on the LHS of the ordered probit regression. Variables must not contain missing values
#' @param evs A vector of explanatory variables to be used on the RHS of the ordered probit regression. Variables may contain missing values.
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
#'   \item{data.orig}{Original data fed into the function}
#'   \item{data.orig.na.omit}{Original data without missing values}
#'   \item{data.cut}{Data after cutpoint replacements}
#'   \item{plr.out}{Results \code{polr}}
#'   \item{plr.df}{Results of \code{polr} as a data frame}
#'   \item{int.dfs}{A list of intercepts as data frames}
#'   \item{ord.new.lev}{New ordinal variable levels}
#'   \item{ord.new.lev.num}{Numeric version of new ordinal levels}
#'   }
#'
#' @importFrom stats aggregate as.formula coef cor na.omit
#' @examples
#' data(ampData)
#' hd.ord(data = ampData,
#'       ord = c("Educ", "Interest"),
#'       evs = c("Dem", "Black", "Empl", "Male", "Inc", "Age"))
#'
#' @export
hd.ord <- function(data, ord, evs, m = 5, method=c("best.cell", "p.draw"), cutoff=10, sdCutoff=1, optimizeSD = FALSE,
                  optimStep = 0.1, optimStop = 5, weightedAffinity = FALSE, impContinuous = c("HD", "mice"),
                  IDvars = NULL, ...){
  # save original version of data
  data.orig <- data
  # use na.omit version of the data in the function
  data.orig.na.omit <- na.omit(data)
  ord.new.lev <- list()
  ord.new.lev.num <- list()
  int.dfs <- list()
  # turn column classes depending on original class
  for(i in 1:length(ord)){
    if(is.factor(data.orig.na.omit[, ord[i]]) == TRUE){
      # change it to numeric and save under new numeric name
      data.orig.na.omit[, paste0(ord[i], ".num")] <- as.numeric(data.orig.na.omit[, ord[i]])
    }
    if(is.numeric(data.orig.na.omit[, ord[i]]) == TRUE){
      # save it under new numeric name
      data.orig.na.omit[, paste0(ord[i], ".num")] <- data.orig.na.omit[,ord[i]]
      # change it to factor and save under current name
      data.orig.na.omit[, ord[i]] <- as.factor(data.orig.na.omit[, ord[i]])
    }
    if(is.character(data.orig.na.omit[, ord[i]]) == TRUE){
      # change it to numeric and save under new numeric name
      data.orig.na.omit[, paste0(ord[i], ".num")] <- as.numeric(data.orig.na.omit[,ord[i]])
      # change it to factor and save under current name
      data.orig.na.omit[, ord[i]] <- as.factor(data.orig.na.omit[, ord[i]])
    }
    # if(ord[i] %in% NAcols == TRUE){
    #   stop(paste0("The variable '", ord[i], "' cannot be the ordinal variable and a column with NAs at the same time"))
    # }
    variables <- c(ord[i], evs)
    a <- as.formula(data.orig.na.omit[,variables])
    plr.out <- MASS::polr(a, data = data.orig.na.omit, Hess=TRUE)
    # turn plr.out output into a data frame
    plr.df <- data.frame(coef(summary(plr.out)))
    # empty vector for storage
    storage <- c()
    # loop to fill empty vector with intercept names
    for (x in 1:(length(levels(data.orig.na.omit[,ord[i]]))-1)){
      storage[x] <- paste0(levels(data.orig.na.omit[,ord[i]])[x], "|", levels(data.orig.na.omit[,ord[i]])[x+1])
    }
    # select only rows with intercept names
    int.df <- plr.df[storage,]
    # turn row names into a column
    int.df <- data.table::setDT(int.df, keep.rownames = TRUE)[]
    # removes class "data.table", which was added by setDT
    int.df <- data.frame(int.df)
    colnames(int.df) <- c("Intercepts", "Values", "SE", "t-values")
    # factorize intercepts with correctly ordered levels
    int.df$Intercepts <- factor(int.df$Intercepts, levels = int.df$Intercepts)
    # empty df to fill with assigned new cases
    df.cases <- data.frame(matrix(NA, nrow(data.orig.na.omit), length(levels(data.orig.na.omit[,ord[i]]))))
    # assign cases that fall underneath the lowest intercept with the respective education category, put results in first column
    df.cases[,1] <- ifelse(plr.out$lp <= int.df$Values[1], levels(data.orig.na.omit[,ord[i]])[1], NA)
    # assign cases that fall between all intercepts except the lowest and the highest, put results in all but first and last columns
    for(q in 1:(length(levels(data.orig.na.omit[,ord[i]]))-2)){
      df.cases[,q+1] <- ifelse(plr.out$lp > int.df$Values[q] &
                                 plr.out$lp <= int.df$Values[q+1], levels(data.orig.na.omit[,ord[i]])[q+1], NA)
    }
    # assign all cases that fall above the highest intercept, put in last column
    df.cases[,length(levels(data.orig.na.omit[,ord[i]]))] <- ifelse(plr.out$lp > int.df$Values[length(levels(data.orig.na.omit[,ord[i]]))-1],
                                                                    levels(data.orig.na.omit[,ord[i]])[length(levels(data.orig.na.omit[,ord[i]]))], NA)
    # combine all columns, omit NAs, add column to df
    data.orig.na.omit[,paste0(ord[i], ".new")] <- factor(c(na.omit(c(t(df.cases)))))
    # empty small df to extract numbers for the respective re-estimated categories
    df.factors.1 <- data.frame(matrix(NA, length(levels(data.orig.na.omit[,paste0(ord[i], ".new")])), 2))
    # for each remaining category, extract its name and match it with its respective number, add by rows
    for (w in 1:length(levels(data.orig.na.omit[,paste0(ord[i], ".new")]))){
      df.factors.1[w,] <- c(unique(data.orig.na.omit[,paste0(ord[i], ".num")][data.orig.na.omit[,ord[i]] == levels(data.orig.na.omit[,paste0(ord[i], ".new")])[w]]),
                            levels(data.orig.na.omit[,paste0(ord[i], ".new")])[w])
    }
    # empty long df to set up numbers for re-estimated categories for entire data set
    df.factors.2 <- data.frame(matrix(NA, nrow(data.orig.na.omit), nrow(df.factors.1)))
    # assign numbers, add by columns
    for (n in 1:nrow(df.factors.1)){
      df.factors.2[,n] <- ifelse(data.orig.na.omit[,paste0(ord[i], ".new")] == df.factors.1[n,2], df.factors.1[n,1], NA)
    }
    # combine all columns, omit NAs, add column to df
    data.orig.na.omit[,paste0(ord[i], ".new.num")] <- as.numeric(na.omit(c(t(df.factors.2))))
    # refactor levels in the order of the numbers
    data.orig.na.omit[,paste0(ord[i], ".new")] <- factor(data.orig.na.omit[,paste0(ord[i], ".new")],
                                                         levels = unique(data.orig.na.omit[,paste0(ord[i], ".new")][order(data.orig.na.omit[,paste0(ord[i], ".new.num")])]))
    # make numbers for .new.num column start at 1, based on newly refactored .new levels
    data.orig.na.omit[,paste0(ord[i], ".new.num")] <- as.numeric(data.orig.na.omit[,paste0(ord[i], ".new")])

    int.dfs[[i]] <- int.df
    ord.new.lev[[i]] <- levels(data.orig.na.omit[,paste0(ord[[i]], ".new")])
    ord.new.lev.num[[i]] = sort(unique(data.orig.na.omit[,paste0(ord[[i]], ".new.num")]))
  }
  names(ord.new.lev) <- names(ord.new.lev.num) <- names(int.dfs) <- ord
  # turn column classes depending on original class
  for(i in 1:length(ord)){
    # if(ord[i] %in% NAcols == TRUE){
    #   stop(paste0("The variable '", ord[i], "' cannot be the ordinal variable and a column with NAs at the same time"))
    # }
    if(is.factor(data[, ord[i]]) == FALSE){
      # change it to factor and save under current name
      data[, ord[i]] <- as.factor(data[, ord[i]])
    }
    lev <- levels(data[, ord[i]])
    penult <- length(lev)-1
    # category span between first and second level
    cat_second_span <- abs(int.dfs[[ord[i]]][1,2] - int.dfs[[ord[i]]][2,2])

    # category span between penultimate and last level
    cat_penult_span <- abs(int.dfs[[ord[i]]][(penult-1),2] - int.dfs[[ord[i]]][penult,2])
    # beginning cutpoint for first level
    cat_first_cut <- int.dfs[[ord[i]]][1,2] - cat_second_span
    # end cutpoint for last level
    cat_last_cut <- int.dfs[[ord[i]]][penult,2] + cat_penult_span
    # empty vector to store middle between cutpoints for each level
    midpoints <- c()
    # mid-cutpoint for first level
    midpoints[1] <- (cat_first_cut + int.dfs[[ord[i]]][1,2])/2
    # mid-cutpoints for all except first and last levels

    for (x in 1:(penult-1)){
      midpoints[x+1] <- (int.dfs[[ord[i]]][x,2] + int.dfs[[ord[i]]][(x+1),2])/2
    }
    # mid-cutpoint for last level
    midpoints[length(lev)] <- (int.dfs[[ord[i]]][penult,2] + cat_last_cut)/2
    # empty df to store ord level replacements
    replace.lev.df <- data.frame(matrix(NA, nrow(data), length(lev)))
    # replace ord levels with mid-cutpoints
    for (x in 1:length(lev)){
      replace.lev.df[, x] <- ifelse(data[, ord[i]] == lev[x], midpoints[x], NA)
    }
    # overwrite ord with mid-cutpoints replacement
    data[, ord[i]] <- c(na.omit(c(t(replace.lev.df))))
    # save cut version of data
  }
  data.cut <- data
  method <- match.arg(method)
  impContinuous <- match.arg(impContinuous)
  if(weightedAffinity){
    warning("Affinity calculations made as a function of pearson correlations among variables coerced to class 'numeric'\ntake care when using this on categorical, especially nominal variables")
  }
  if(!is.null(IDvars)){
    IDdata <- data[, which(names(data) %in% IDvars), drop=FALSE]
    data <- data[,-which(names(data) %in% IDvars), drop=FALSE]
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
  if(impContinuous == "HD" & method == "p.draw" & length(cont.miss) > 0){
    stop("Hot Deck imputation of continuous values can only be used with the best cell method\n")
  }
  whichna <- which(is.na(data), arr.ind=TRUE)
  if(impContinuous == "mice"){
    whichna <- whichna[which(whichna[,2] %in% disc.miss), ]
  }
  if(optimizeSD & any(!alldisc)){
    mm <- 0
    while(sdCutoff <= optimStop & mm < m){
      tmp <- scaleContinuous(data, alldisc, sdx=1/sdCutoff)
      numdata <- sapply(1:ncol(tmp), function(i)as.numeric(tmp[,i]))
      R <- abs(cor(numdata, use="pairwise"))
      diag(R) <- 0
      unnaobs <- unique(whichna[,1])
      if(!weightedAffinity){
        aff <- t(sapply(unnaobs, function(x)affinity(numdata, x, weighted=FALSE)))
        aff <- aff[match(whichna[,1], unnaobs), ]
      }
      if(weightedAffinity){
        aff <- t(apply(whichna, 1, function(x)affinity(numdata, x[1], x[2], R, weightedAffinity)))
      }
      if(any(!is.finite(aff))){
        aff[which(!is.finite(aff), arr.ind=TRUE)] <- 0
      }
      wnadat <- matrix(1, nrow=nrow(data), ncol=ncol(data))
      wnadat[which(is.na(data), arr.ind=TRUE)] <- 0
      wnadat <- t(wnadat[, whichna[,2]])
      aff <- aff*wnadat
      w <- apply(aff, 1, function(x)which(x == max(x)))
      donors <- lapply(1:nrow(whichna), function(x)na.omit(data[w[[x]], whichna[x,2]]))
      matches <- sapply(donors, length)
      mm <- min(matches)
      cat("SD Cutoff = ", sprintf("%.2f", sdCutoff), ", # Thin Cells = ", sum(matches < m), "\n", sep="")
      if(mm < m & sdCutoff == optimStop){
        warning(paste("Optimization unsuccessful, ", sum(matches < m), " thin cells remain with SD cutoff of ", sdCutoff, "\n", sep=""))
      }
      if(sdCutoff < optimStop){
        sdCutoff <- sdCutoff + optimStep
      }
    }
  }
  tmp <- scaleContinuous(data, alldisc, sdx=1/sdCutoff)
  ord.disc <- FALSE  # added for scaleContinuous for all ordinal variables
  for(i in 1:length(ord)){ # added for scaleContinuous for all ordinal variables
    tmp[, ord[i]] <- scaleContinuous(data[, ord[i]], ord.disc, sdx = 1/sdCutoff)[,1]
  }
  # runs scaleContinuous only on the ordinal variables and replaces the unscaled versions of the variables in tmp
  numdata <- sapply(1:ncol(tmp), function(i)as.numeric(tmp[,i]))
  R <- abs(cor(numdata, use="pairwise"))
  diag(R) <- 0
  max.emp.aff <- 	apply(R, 2, sum)[whichna[,2]]
  max.the.aff <- rep(dim(R)[2] - 1, nrow(whichna))
  unnaobs <- unique(whichna[,1])
  if(!weightedAffinity){
    aff <- t(sapply(unnaobs, function(x)affinity(numdata, x, weighted=FALSE)))
    aff <- aff[match(whichna[,1], unnaobs), ]
  }
  if(weightedAffinity){
    aff <- t(apply(whichna, 1, function(x)affinity(numdata, x[1], x[2], R, weightedAffinity)))
  }
  if(any(!is.finite(aff))){
    aff[which(!is.finite(aff), arr.ind=TRUE)] <- 0
  }

  wnadat <- matrix(1, nrow=nrow(data), ncol=ncol(data))
  wnadat[which(is.na(data), arr.ind=TRUE)] <- 0
  wnadat <- t(wnadat[, whichna[,2]])
  aff <- aff*wnadat
  if(method == "best.cell"){
    w <- apply(aff, 1, function(x)which(x == max(x)))
    donors <- lapply(1:nrow(whichna), function(x)na.omit(data[w[[x]], whichna[x,2]]))
    matches <- sapply(donors, length)
    if(any(matches < m)){
      warning(paste(sum(matches < m ), " of ", length(matches), " imputations with # donors < ", m, ", consider increasing sdCutoff or using method='p.draw'\n", sep=""))
    }
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
    if(!is.null(IDvars)){
      res[[md]] <- cbind(IDdata, res[[md]])
    }
  }
  class(res)  <- c("mi","list")
  output <- list(data = res, # data after applying hot decking
                 affinity = aff,
                 donors = donors,
                 draws = draws,
                 max.emp.aff = max.emp.aff,
                 max.the.aff = max.the.aff,
                 data.orig = data.orig, # original data fed into the function
                 data.orig.na.omit = data.orig.na.omit, # na.omit() of original data
                 data.cut = data.cut, # data after cutpoint replacements
                 plr.out = plr.out, # results of polr()
                 plr.df = plr.df, # results of polr() as a data frame
                 int.dfs = int.dfs, # intercepts as data frames (in a list)
                 ord.new.lev = ord.new.lev, # new ordinal levels
                 ord.new.lev.num = ord.new.lev.num) # numeric version of new ordinal levels
  return(output)
}
