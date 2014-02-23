\name{EMpow}
\alias{EMdpow}
\alias{EMipow}
\alias{EMpow}
\title{the EMpow set
\code{EMpow} - A family of power functions,
predicting size at time t, from initial size M0, and a set of
parameters. The set includes an inverse and differential, which
predict age at size M, and growth rates dM/dtime for a size M.}
\usage{
EMpow(M0, r, beta, time)

EMipow(M, r, beta, time)

EMdpow(M, r, beta)
}
\arguments{
  \item{time}{time since initial size}

  \item{M}{mass}

  \item{M0}{initial mass}

  \item{r}{growth rate}

  \item{beta}{scaling factor}
}
\description{
the EMpow set \code{EMpow} - A family of power functions,
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

