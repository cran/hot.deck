\name{scaleContinuous}
\alias{scaleContinuous}
\title{Standardize continuous variables.}
\description{
Standardizes (centers and scales) continuous variable in a dataset, leaving discrete variables untouched.
}
\usage{
scaleContinuous(data, discrete, sdx = 1)
}
\arguments{
  \item{data}{
A data frame or matrix of variables to be scaled. 
}
  \item{discrete}{
Either a logical vector which is \code{TRUE} for discrete variables and \code{FALSE} for continuous ones or a vector of column numbers of discrete variables.
}
  \item{sdx}{
The standard deviation of the columns for the continuous variables.  
}
}
\value{
A data frame with the same dimensions as \code{data} where the continuous variables are centered and scaled.
}
\author{
Skyler Cranmer, Jeff Gill, Natalie Jackson, Andreas Murr and Dave Armstrong
}
