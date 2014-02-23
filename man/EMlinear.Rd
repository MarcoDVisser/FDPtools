\name{EMlinear}
\alias{EMdlinear}
\alias{EMilinear}
\alias{EMlinear}
\title{the EMlinear set

\code{EMlinear} - A family of linear functions
predicting size at time t, from initial size M0, and a set of
parameters. The set includes an inverse and differential, which
predict age at size M, and growth rates dM/dtime for a size M.}
\usage{
EMlinear(M0, K, r, time = 1)

EMilinear(M, M0, r)

EMdlinear(M, r)
}
\arguments{
  \item{time}{time since initial size}

  \item{M}{mass}

  \item{M0}{initial mass}

  \item{r}{growth rate}
}
\description{
the EMlinear set

\code{EMlinear} - A family of linear functions predicting
size at time t, from initial size M0, and a set of
parameters. The set includes an inverse and differential,
which predict age at size M, and growth rates dM/dtime for
a size M.
}
\author{
Marco D. Visser

Marco D. Visser

Marco D. Visser
}

