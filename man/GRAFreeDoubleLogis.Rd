\name{GRAFreeDoubleLogis}
\alias{GRAFreeDoubleLogis}
\alias{GRAdFreeDoubleLogis}
\alias{GRAiFreeDoubleLogis}
\title{the GRAFreeDoubleLogis set

\code{GRAFreeDoubleLogis} - A family of double logistic functions used in
gradient modelling. These models are usually used to predict frequency,
 a probabily or fractions over x (distance, time, concentration).
The 4 parameter model is limited to the range [beta4,beta3].
The set includes an inverse and differential, which
give the inverse estimate (x at a given y) or
dy/dx rates over x.}
\usage{
GRAFreeDoubleLogis(x, beta0, beta1, beta3, beta4, betas = NULL)

GRAiFreeDoubleLogis(M, M0, r)

GRAdFreeDoubleLogis(M, r)
}
\arguments{
  \item{x}{the gradient (size, time, concentration)}

  \item{y}{prediction, used in the inverse form}

  \item{beta0}{intercept}

  \item{beta1}{rate of increase}

  \item{beta3}{second intercept}

  \item{beta4}{second rate of increase}

  \item{betas}{vector to input all parameters in one go}
}
\description{
the GRAFreeDoubleLogis set

\code{GRAFreeDoubleLogis} - A family of double logistic
functions used in gradient modelling. These models are
usually used to predict frequency, a probabily or fractions
over x (distance, time, concentration). The 4 parameter
model is limited to the range [beta4,beta3]. The set
includes an inverse and differential, which give the
inverse estimate (x at a given y) or dy/dx rates over x.
}
\author{
Marco D. Visser

Marco D. Visser

Marco D. Visser
}

