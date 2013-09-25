\name{calculateCrowding}
\alias{calculateCrowding}
\title{calculateCrowding}
\usage{
  calculateCrowding(FdpData = list(bci.full6, bci.full7),
    r = seq(1, 50, 10), statistic = "count", ...)
}
\arguments{
  \item{FdpData}{A list of all census years on which the
  calculation must be conducted}

  \item{r}{a vector containing all distances at which the
  below statistic must be calculated}

  \item{statsistic}{character defining which method to use.
  Methods include: 1) "count", simple count of neighbours
  within r 2) "distance index", a distance weighted
  function 3) "gravitation", a size and distance weighted
  function}

  \item{...}{further arguments to be passed to statistic
  calculation}
}
\description{
  calculates Crowding statistics for each individual alive
  in a given census year for the BCI plot data for further
  analysis.
}
\author{
  Marco D. Visser
}

