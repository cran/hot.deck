\name{D}
\alias{D}
\docType{data}
\title{Example data for multiple hot deck imputation.}
\usage{data(D)}
\description{Simulated example data for multiple hot deck imputation.}

\format{
  A data frame with 20 observations on the following 5 variables.
  \describe{
    \item{\code{x1}}{a numeric vector}
    \item{\code{x2}}{a numeric vector}
    \item{\code{x3}}{a numeric vector}
    \item{\code{x4}}{a numeric vector}
    \item{\code{x5}}{a numeric vector}
  }
}
\examples{
data(D)
out <- hot.deck(D)
}
\keyword{datasets}
