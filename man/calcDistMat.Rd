\name{calcDistMat}
\alias{calcDistMat}
\title{Calculate distances between points within and between
datasets}
\usage{
  calcDistMat(pointdata, referencedata = NULL,
    xy = c("x", "y"))
}
\arguments{
  \item{pointdata}{a dataset containing x and y
  coordinates.}

  \item{referencedata}{a reference data set, used when the
  distance between all points in the dataset pointdata and
  all points in the reference data set is needed. If NULL
  this is ignored and the distance between all points in
  pointdata is calculated.}

  \item{xy}{column names containing coordinates. Length
  should be 4 if using a reference dataset, with the first
  and second elements referering to pointdata and the
  second and third to the reference dataset.}
}
\description{
  \code{calcDistMat} - Generic optimized function for
  calculation of 1) distances between all points in a
  datatset containing x and y coordinates or 2) between all
  points in two datasets containing x and y coordinates.
}
\examples{
\dontrun{
x <- data.frame(x.cor=runif(10),y.cor=runif(10))
calcDistMat(x,xy=c("x.cor","y.cor"))
}
}
\author{
  Marco D. Visser
}

