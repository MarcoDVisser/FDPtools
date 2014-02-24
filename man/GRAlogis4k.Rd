\name{GRAlogis4k}
\alias{GRAdlogis4k}
\alias{GRAilogis4k}
\alias{GRAlogis4k}
\title{the GRAlogis4k set

\code{GRAlogis4k} - A family of 4 parameter logistic functions for
gradient modelling. These models are usually used to predict frequency,
 a probabily or fractions over x (distance, time, concentration).
The 4 parameter model is limited to the range [beta4,beta3].
The set includes an inverse and differential, which
give the inverse estimate (x at a given y) or
dy/dx rates over x.}
\usage{
GRAlogis4k(x, beta0, beta1, beta3, beta4, betas = NULL)

GRAilogis4k(M, M0, r)

GRAdlogis4k(M, r)
}
\arguments{
  \item{x}{the gradient (size, time, concentration)}

  \item{y}{prediction, used in the inverse form}

  \item{beta0}{intercept}

  \item{beta1}{rate of increase}

  \item{beta3}{Asymtotic rate}

  \item{beta4}{lower bound}

  \item{betas}{vector to input all parameters in one go}
}
\description{
the GRAlogis4k set

\code{GRAlogis4k} - A family of 4 parameter logistic
functions for gradient modelling. These models are usually
used to predict frequency, a probabily or fractions over x
(distance, time, concentration). The 4 parameter model is
limited to the range [beta4,beta3]. The set includes an
inverse and differential, which give the inverse estimate
(x at a given y) or dy/dx rates over x.
}
\author{
Marco D. Visser

Marco D. Visser

Marco D. Visser
}

