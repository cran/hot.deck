\name{hot.deck-package}
\alias{hot.deck-package}
\docType{package}
\title{Multiple Hot-Deck Imputation}
\description{This package contains all of the functions necessary to perform multiple hot deck imputation on an input data frame with missing observations using either the \dQuote{best cell} method (default) or the \dQuote{probabilistic draw} method as described in Cranmer and Gill (2013). This technique is best suited for missingness in discrete variables, though it also works well for continuous missing observations.}
\details{
\tabular{ll}{
Package: \tab hot.deck\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2015-11-19\cr
License: \tab GPL (>= 2)\cr
}
In multiple hot deck imputation, several observed values of the variable with missing observations are drawn conditional on the rest of the data and are used to impute each missing value. The advantage of this class of methods over multiple imputation is that the imputed values are actually draws from the observed data. As such, when discrete variables are imputed with a hot deck method, their discrete properties are maintained.

Two methods for weighting the imputations are provided in this package. The \dQuote{best cell} [called as \dQuote{best.cell}] technique uses the degree of affinity between the row with missing data and each potential donor row to generate weights such that rows more closely resembling the row with missingness are more likely to be drawn as donors. The probabilistic draw method is the default method. The \dQuote{probabilistic draw}  [called as \dQuote{p.draw}] technique is also available. The best cell method draws randomly from the cell of best matches to the row with a missing observation. }
\author{Skyler Cranmer, Jeff Gill, Natalie Jackson, Andreas Murr and Dave Armstrong
Maintainer: Dave Armstrong <dave@quantoid.net>
}
\references{
Cranmer, S.J. and Gill, J.M.. (2013) \dQuote{We Have to Be Discrete About This: A Non-Parametric Imputation Technique for Missing Categorical Data.} \emph{British Journal of Political Science} 43:2 (425-449).  
}
\keyword{ package }
\keyword{multiple imputation}
