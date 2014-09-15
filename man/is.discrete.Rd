\name{is.discrete}
\alias{is.discrete}
\title{Identify whether variables are discrete or continuous
}
\description{
Variables are considered discrete if they have fewer unique, non-missing values than \code{cutoff} or they are factors.  Otherwise, variables are considered continuous.
}
\usage{
is.discrete(data, cutoff = 10)
}
\arguments{
  \item{data}{A data frame, matrix or vector of values to be evaluated.}
  \item{cutoff}{A numeric scalar identifying the cutoff relative to the number of unique, non-missing values for \sQuote{discreteness}.}
}
\value{
A logical vector indicating whether variables are discrete (\code{TRUE}) or continuous \code{FALSE}.}
\author{
Skyler Cranmer, Jeff Gill, Natalie Jackson, Andreas Murr and Dave Armstrong
}
