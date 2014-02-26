\name{GROlogis4k}
\alias{GROdlogis4k}
\alias{GROilogis4k}
\alias{GROlogis4k}
\title{the GROlogis 4 parameter set

\code{GROlogis4k} - A family of 4 parameter logistic functions for size
growth modelling,
predicting size at time t, from initial size M0, and a set of
parameters. The set includes an inverse and differential, which
predict age at size M, and growth rates dM/dtime for a size M.}
\usage{
GROlogis4k(M0, K, r, L, P, time = 1, betas = NULL)

GROilogis4k(M, M0, K, r, L, P, betas = NULL)

GROdlogis4k(M, L, K, r)
}
\arguments{
  \item{time}{time since initial size}

  \item{M}{mass}

  \item{M0}{initial mass}

  \item{K}{asymptotic size}

  \item{r}{growth rate}

  \item{L}{lowerlimit}

  \item{P}{scaling parameter?}
}
\description{
the GROlogis 4 parameter set

\code{GROlogis4k} - A family of 4 parameter logistic
functions for size growth modelling, predicting size at
time t, from initial size M0, and a set of parameters. The
set includes an inverse and differential, which predict age
at size M, and growth rates dM/dtime for a size M.
}
\author{
Marco D. Visser

Marco D. Visser

Marco D. Visser
}

