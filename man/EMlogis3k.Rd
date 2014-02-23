\name{EMlogis3k}
\alias{EMdlogis3k}
\alias{EMilogis3k}
\alias{EMlogis3k}
\title{the EMlogis 3 parameter set

\code{EMlogis3k} - A family of 3 parameter logistic functions
predicting size at time t, from initial size M0, and a set of
parameters. The set includes an inverse and differential, which
predict age at size M, and growth rates dM/dtime for a size M.}
\usage{
EMlogis3k(M0, K, r, time = 1)

EMilogis3k(M, M0, K, r)

EMdlogis3k(M, r)
}
\arguments{
  \item{time}{time since initial size}

  \item{M}{mass}

  \item{M0}{initial mass}

  \item{K}{asymptotic size}

  \item{r}{growth rate}
}
\description{
the EMlogis 3 parameter set

\code{EMlogis3k} - A family of 3 parameter logistic
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

