\name{EMmonomol}
\alias{EMdmonomol}
\alias{EMimonomol}
\alias{EMmonomol}
\title{the EMmonomol set

\code{EMmonomol} - A family of Monomolecular functions
predicting size at time t, from initial size M0, and a set of
parameters. The set includes an inverse and differential, which
predict age at size M, and growth rates dM/dtime for a size M.}
\usage{
EMmonomol(M0, K, r, time = 1)

EMimonomol(M, M0, K, r)

EMdmonomol(M, K, r)
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
the EMmonomol set

\code{EMmonomol} - A family of Monomolecular functions
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

