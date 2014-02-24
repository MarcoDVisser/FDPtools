\name{GROmonomol}
\alias{GROdmonomol}
\alias{GROimonomol}
\alias{GROmonomol}
\title{the GROmonomol set

\code{GROmonomol} - A family of Monomolecular functions for size
growth modelling,
predicting size at time t, from initial size M0, and a set of
parameters. The set includes an inverse and differential, which
predict age at size M, and growth rates dM/dtime for a size M.}
\usage{
GROmonomol(M0, K, r, time = 1)

GROimonomol(M, M0, K, r)

GROdmonomol(M, K, r)
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
the GROmonomol set

\code{GROmonomol} - A family of Monomolecular functions for
size growth modelling, predicting size at time t, from
initial size M0, and a set of parameters. The set includes
an inverse and differential, which predict age at size M,
and growth rates dM/dtime for a size M.
}
\author{
Marco D. Visser

Marco D. Visser

Marco D. Visser
}

