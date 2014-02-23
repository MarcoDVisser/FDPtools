\name{EMlogis4k}
\alias{EMdlogis4k}
\alias{EMilogis4k}
\alias{EMlogis4k}
\title{the EMlogis 4 parameter set

\code{EMlogis4k} - A family of 4 parameter logistic functions
predicting size at time t, from initial size M0, and a set of
parameters. The set includes an inverse and differential, which
predict age at size M, and growth rates dM/dtime for a size M.}
\usage{
EMlogis4k(M0, K, L, P, r, time = 1)

EMilogis4k(M, M0, K, L, P, r)

EMdlogis4k(M, L, K, r)
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
the EMlogis 4 parameter set

\code{EMlogis4k} - A family of 4 parameter logistic
functions predicting size at time t, from initial size M0,
and a set of parameters. The set includes an inverse and
differential, which predict age at size M, and growth rates
dM/dtime for a size M.
}
\author{
Marco D. Visser

Marco D. Visser

Marco D. Visser
}

