\name{ampData}
\alias{ampData}
\docType{data}
\title{Example data for multiple hot deck imputation with ordinal variables.}
\usage{data(ampData)}
\description{Simulated example data for multiple hot deck imputation with ordinal variables.}

\format{
  A data frame with 1000 observations on the following 20 variables.
  \describe{
    \item{\code{Ind}}{a numeric binary vector indicating Independent party ID}
    \item{\code{Black}}{a numeric binary vector indicating African-American ethnicity}
    \item{\code{Hisp}}{a numeric binary vector indicating Hispanic ethnicity}
    \item{\code{Asian}}{a numeric binary vector indicating Asian ethnicity}
    \item{\code{Empl}}{a numeric binary vector indicating employment}
    \item{\code{Stud}}{a numeric binary vector indicating students}
    \item{\code{Interest}}{a numeric vector indicating political interest}
    \item{\code{Educ}}{a numeric vector indicating education level}
    \item{\code{Religious}}{a numeric binary vector indicating religious affiliation}
    \item{\code{InternetHome}}{a numeric binary vector indicating the presence of internet at home}
    \item{\code{OwnHome}}{a numeric binary vector indicating home ownership}
    \item{\code{Rally}}{a numeric binary vector indicating attendance at political rallies}
    \item{\code{Donate}}{a numeric binary vector indicating donations}
    \item{\code{Moderate}}{a numeric binary vector indicating moderate political ideology}
    \item{\code{Married}}{a numeric binary vector indicating marriage}
    \item{\code{Separated}}{a numeric binary vector indicating separation}
    \item{\code{Dem}}{a numeric binary vector indicating Democratic party ID, contains missing values}
    \item{\code{Male}}{a numeric binary vector indicating men, contains missing values}
    \item{\code{Inc}}{a numeric vector indicating income, contains missing values}
    \item{\code{Age}}{a numeric vector indicating age, contains missing values}
  }
}
\examples{
data(ampData)
hd.ord(data = ampData,
      ord = c("Educ", "Interest"),
      evs = c("Dem", "Black", "Empl", "Male", "Inc", "Age"))
}
\keyword{datasets}
