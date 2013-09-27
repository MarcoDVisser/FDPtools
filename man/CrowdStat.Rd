\name{CrowdStat}
\alias{CrowdStat}
\title{CrowdStat}
\usage{
  CrowdStat(censusdata, r, ...)
}
\arguments{
  \item{censusdata}{a single census year, with only
  individuals alive within this census included.}

  \item{r}{a vector containing all distances at which the
  below statistic must be calculated}

  \item{...}{further arguments to be passed to statistic
  calculation}
}
\description{
  \code{CrowdStat} Generic function which calls specific
  crowding statistic to calculate crowding for each
  individual in a dataset for each radius.
}
\author{
  Marco D. Visser
}

