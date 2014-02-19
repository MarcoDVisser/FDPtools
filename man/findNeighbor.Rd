\name{findNeighbor}
\alias{findNeighbor}
\title{findNeighbor a lower-level helper function}
\usage{
findNeighbor(x1, y1, x2, y2, radius, annulus = NULL)
}
\arguments{
  \item{x1}{x coordinates of points}

  \item{y1}{y coordinates of points}

  \item{x2}{reference point x coordinate}

  \item{y2}{reference point y coordinate}

  \item{radius}{radius with which to test distance}

  \item{annulus}{radii of the two concentric annulus
  circles. The smallest radius if the first element.}
}
\description{
Wrapper function which calls some underlying C code and
returns a matrix with 1 or 0 indicating whether the
distance between the points is smaller or equal to radius
or within an annulus.
}
\author{
Marco D. Visser
}

