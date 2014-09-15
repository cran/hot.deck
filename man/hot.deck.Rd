%% DA, 9/5/14: Changed order of arguments to impContinuous so that HD is the default
\name{hot.deck}
\alias{hot.deck}
\title{Multiple Hot Deck Imputation}
\description{
This function performs multiple hot deck imputation on an input data frame with missing observations using either the \dQuote{best cell} method (default) or the \dQuote{probabilistic draw} method as described in Cranmer and Gill (2013). This technique is best suited for missingness in discrete variables, though it also performs well on continuous missing data.
}
\usage{
hot.deck(data, m = 5, method = c("best.cell", "p.draw"), cutoff = 10, sdCutoff = 1, 
    optimizeSD = FALSE, optimStep = 0.1, optimStop = 5, weightedAffinity = FALSE, 
    impContinuous = c("HD", "mice"), IDvars = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A data frame or matrix with missing values to be imputed using multiple hot deck imputation.}
  \item{m}{Number of imputed datasets required.}
  \item{method}{Method used to draw donors based on affinity either \dQuote{best.cell} (the default) or \dQuote{p.draw} for probabilistic draw}
  \item{cutoff}{A numeric scalar such that any variable with fewer than \code{cutoff} unique non-missing values will be considered discrete and necessarily imputed with hot deck imputation.}
  \item{sdCutoff}{Number of standard deviations between observations such that observations fewer than \code{sdCutoff} standard deviations away from each other are considered sufficiently close to be a match, otherwise they are considered too far away to be a match.}
  \item{optimizeSD}{Logical indicating whether the \code{sdCutoff} parameter should be optimized such that the smallest possible value is chosen that produces no thin cells from which to draw donors.  Thin cells are those where the number of donors is less than \code{m}.}
  \item{optimStep}{The size of the steps in the optimization if \code{optimizeSD} is \code{TRUE}.}
  \item{optimStop}{The value at which optimization should stop if it has not already found a value that produces no thin cells.  If this value is reached and thin cells still exist, a warning will be returned, though the routine will continue using \code{optimStop} as \code{sdCutoff}. }
  \item{weightedAffinity}{Logical indicating whether a correlation-weighted affinity score should be used.}
  \item{impContinuous}{Character string indicating how continuous missing data should be imputed.  Valid options are \dQuote{HD} (the default) in which case hot-deck imputation will be used, or \dQuote{mice} in which case multiple imputation by chained equations will be used.}
  \item{IDvars}{A character vector of variable names not to be used in the imputation, but to be included in the final imputed datasets.}
  \item{...}{Optional additional arguments to be passed down to the \code{mice} routine.}
}
\value{
A list with the following elements: 
\item{data}{An object of class \code{mi} which contains \code{m} imputed datasets.}
\item{affinity}{A matrix of affinity scores see \code{\link{affinity}}.}
\item{donors}{A list of donors for each missing observation based on the affinity score.}
\item{draws}{The \code{m} observations drawn from donors that were used for the multiple imputations.}
\item{max.emp.aff}{Normalization constant for each row of affinity scores; the maximum possible value of the affinity scores if correlation-weighting is used.}
\item{max.the.aff}{Normalization constant for each row of affinity scores; the number of columns in the original data. }
}
\references{
Cranmer, S.J. and Gill, J.M.. (2013) \dQuote{We Have to Be Discrete About This: A Non-Parametric Imputation Technique for Missing Categorical Data.} \emph{British Journal of Political Science} 43:2 (425-449).  

van Buuren, S. and Karin Groothuis-Oudshoorn (2011). \dQuote{mice: Multivariate Imputation by Chained Equations in R.} \emph{Journal of Statistical Software}, 45:3 (1-67).
}
\author{
Skyler Cranmer, Jeff Gill, Natalie Jackson, Andreas Murr and Dave Armstrong
}
\seealso{
\code{\link{mice}}, \code{\link{affinity}}
}
\examples{
data(D)
hot.deck(D)
}
\keyword{multiple imputation}
