\name{EMGompertz}
\alias{EMGompertz}
\alias{EMdGompertz}
\alias{EMiGompertz}
\title{the EMGompertz set

\code{EMGompertz} - A family of Gompertz functions
predicting size at time t, from initial size M0, and a set of
parameters. The set includes an inverse and differential, which
predict age at size M, and growth rates dM/dtime for a size M.}
\usage{
EMGompertz(M0, K, r, time = 1)

EMiGompertz(M, M0, K, r)

EMdGompertz(M, K, r)
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
the EMGompertz set

\code{EMGompertz} - A family of Gompertz functions
predicting size at time t, from initial size M0, and a set
of parameters. The set includes an inverse and
differential, which predict age at size M, and growth rates
dM/dtime for a size M.
}
\author{
Marco D. Visser

Marco D. Visser

Marco D. Visser
}

