\name{affinity}
\alias{affinity}
\title{Affinity Calculation}
\description{Calculates affinity based on Cranmer and Gill (2013).  The function performs the original method (as described in the article) and also a method that takes into account the correlation structure of the observed data that increases efficiency in making matches.}
\usage{
affinity(data, index, column = NULL, R = NULL, weighted = FALSE)
}
\arguments{
  \item{data}{A data frame or matrix of values for which affinity should be calculated}
  \item{index}{A row number identifying the target observation.  Affinity will be calculated between this observation and all others in the dataset.}
  \item{column}{A column number identifying the variable with missing information.  This is only needed for the optional correlation-weighted affinity score.  The correlation that is used is the correlation of all variables with the focus variable (i.e., the column).}
  \item{R}{A correlation matrix for \code{data}.}
  \item{weighted}{Logical indicating whether or not the correlation-weighted affinity measure should be used.}
}
\details{
Affinity is calculated by first identifying whether two observations are sufficiently \sQuote{close} on each variable.  Consider the target observation number 1.  If observation \emph{i} is close to the target observation on variable \emph{j}, then \code{A[i,j] = 1} otherwise, it equals zero.  Close for two discrete variables is defined by them taking on the same value.  Close for continuous variables is taking on a distance no greater than 1 from each other.  While this may seem restrictive and arbitrary, arguments exist in the main package function \code{hot.deck} that allows the user to set how many standard deviations equal a distance of 1 (with the \code{cutoffSD} argument).
}
\value{
A number of missing observation-variable combinations-by-number of observations in data matrix of affinity scores.  
}
\references{
Cranmer, S.J. and Gill, J.M.. (2013) \dQuote{We Have to Be Discrete About This: A Non-Parametric Imputation Technique for Missing Categorical Data.} \emph{British Journal of Political Science} 43:2 (425-449).  
}
\author{
Skyler Cranmer, Jeff Gill, Natalie Jackson, Andreas Murr and Dave Armstrong
}
\seealso{
\code{\link{hot.deck}}
}
\examples{
data(D)
out <- hot.deck(D)
}
\keyword{multiple imputation}

