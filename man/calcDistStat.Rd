\name{calcDistStat}
\alias{calcDistStat}
\title{Calculate distance-based statistics between
datasets or within datasets (with spatial data).}
\usage{
calcDistStat(pointdata, referencedata = NULL, xy = c("x", "y"),
  statistic = NULL, radii = NULL, annuli = NULL)
}
\arguments{
  \item{pointdata}{a dataset containing x and y
  coordinates. Example, fdp spatial data.}

  \item{referencedata}{a reference data set, used when the
  distance between all points in the dataset pointdata and
  all points in the reference data set is needed. If NULL
  this is ignored and the distance between all points in
  pointdata is calculated.  Example, reference data can be
  traplocations or seedling plots.}

  \item{xy}{column names containing coordinates. Length
  should be 4 if using a reference dataset, with the first
  and second elements referering to pointdata and the
  second and third to the reference dataset.}

  \item{statistic}{a vector of information about each point
  in pointdata to summerize. Example: if the statistic
  contains pi*(dbh/2)^2 values, the sum of basal area
  within each radii is returned. The values in "statistic"
  should correspond to pointdata.}

  \item{radii}{a vector of distances.}

  \item{annuli}{a matrix/array/data.frame of annuli. The
  object should contain two columns, containing the radii
  of each pair of two concentric circles that form the
  annuli. The first column must contain the smallest circle
  radius. When annuli are supplied, radii are ignored.}
}
\description{
\code{calcDistStat} - Generic optimized function for
calculation of distance based statistics as e.g. total
basal area with radii
}
\examples{
\dontrun{
x <- data.frame(x.cor=runif(10),y.cor=runif(10))
calcDistMat(x,xy=c("x.cor","y.cor"),statistic=runif(10))
}
}
\author{
Marco D. Visser
}

