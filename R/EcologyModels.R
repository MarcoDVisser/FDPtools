############################################################
## A list of models used in biology and ecology
## with their inverse and differential forms
## Marco Visser, Gamboa, Panama, February 2014
############################################################
#' @name GROlinear
#' @aliases GROilinear
#' @aliases GROdlinear
#' @title the GROlinear set
#' 
#' \code{GROlinear} - A family of linear functions for size
#' growth modelling. These models predict size at time t,
#' from initial size M0, and a set of
#' parameters. The set includes an inverse and differential, which
#' predict age at size M, and growth rates dM/dtime for a size M.
#'
#' @param time time since initial size
#' @param M mass
#' @param M0 initial mass
#' @param r growth rate
#' 
#' @author Marco D. Visser
#' 
#' @export
GROlinear<-function(M0,r,time=1,betas=NULL) {
  if(!is.null(betas)) {
    r <- betas[1] }
  (M0+r*time) }

#' @rdname GROlinear
#' @author Marco D. Visser
#' 
#' @export
GROilinear <- function(M,M0,r,betas=NULL) {
  if(!is.null(betas)) {
    r <- betas[1]
    }
(M-M0)/r
}

#' @rdname GROlinear
#' @author Marco D. Visser
#' 
#' @export
GROdlinear <- function(M,r) {
r
}

#' @name GROexp
#' @aliases GROiexp
#' @aliases GROdexp
#' @title the GROexp set
#' \code{GROexp} - A family of exponential functions for size
#' growth modelling. These models predict size at time t,
#' from initial size M0, and a set of
#' parameters. The set includes an inverse and differential, which
#' predict age at size M, and growth rates dM/dtime for a size M.
#'
#' 
#' @param time time since initial size
#' @param M mass
#' @param M0 initial mass
#' @param r growth rate
#' 
#' @author Marco D. Visser
#' 
#' @export
GROexp <- function(M0,r,time,betas=NULL) {
  if(!is.null(betas)) {
    r <- betas[1]
    }
 M0*exp(r*time)
}

#' @rdname GROexp
#' @author Marco D. Visser
#' 
#' @export
GROiexp <- function(M0,r,time,betas=NULL) {
  if(!is.null(betas)) {
    r <- betas[1]
    }
 M0*exp(r*time)
}

#' @rdname GROexp
#' @author Marco D. Visser
#' 
#' @export
GROdexp <- function(M,r) {
 r*M
}

#' @name GROpow
#' @aliases GROipow
#' @aliases GROdpow
#' @title the GROpow set
#' \code{GROpow} - A family of power functions for size
#' growth modelling, predicting size at time t, from initial size M0,
#' and a set of
#' parameters. The set includes an inverse and differential, which
#' predict age at size M, and growth rates dM/dtime for a size M.
#'
#' 
#' @param time time since initial size
#' @param M mass
#' @param M0 initial mass
#' @param r growth rate
#' @param beta scaling factor
#' 
#' @author Marco D. Visser
#' 
#' @export
GROpow <- function(M0,r,beta,time,betas=NULL) {
  if(!is.null(betas)) {
    r <- betas[1]
    beta <- betas[2] }
((M0^(1-beta))+r*time*(1-beta))^(1/(1-beta))
}

#' @rdname GROpow
#' @author Marco D. Visser
#' 
#' @export
GROipow <- function(M,M0,r,beta,time,betas=NULL) {
  if(!is.null(betas)) {
    r <- betas[1]
    beta <- betas[2] }
M^(1-beta)/((beta-1)*r)
}

#' @rdname GROpow
#' @author Marco D. Visser
#' 
#' @export
GROdpow <- function(M,r,beta) {
 r*M^beta
}

#' @name GROlogis3k
#' @aliases GROilogis3k
#' @aliases GROdlogis3k
#' @title the GROlogis 3 parameter set
#' 
#' \code{GROlogis3k} - A family of 3 parameter logistic functions for size
#' growth modelling,
#' predicting size at time t, from initial size M0, and a set of
#' parameters. The set includes an inverse and differential, which
#' predict age at size M, and growth rates dM/dtime for a size M.
#'
#' @param time time since initial size
#' @param M mass
#' @param M0 initial mass
#' @param K asymptotic size
#' @param r growth rate
#' 
#' @author Marco D. Visser
#' 
#' @export
GROlogis3k<-function(M0,K,r,time=1,betas=NULL) {
  if(!is.null(betas)) {
    K <- betas[1]
    r <- betas[2] }
(M0*K)/(M0+(K-M0)*exp(-r*time))
}

#' @rdname GROlogis3k
#' @author Marco D. Visser
#' 
#' @export
GROilogis3k <- function(M,M0,K,r,betas=NULL) {
  if(!is.null(betas)) {
    K <- betas[1]
    r <- betas[2] }
log((K*M)/(K*M0-M*M0)-(M*M0)/(K*M0-M*M0))/r
}

#' @rdname GROlogis3k
#' @author Marco D. Visser
#' 
#' @export
GROdlogis3k <- function(M,r) {
 r*M
}

#' @name GROlogis4k
#' @aliases GROilogis4k
#' @aliases GROdlogis4k
#' @title the GROlogis 4 parameter set
#' 
#' \code{GROlogis4k} - A family of 4 parameter logistic functions for size
#' growth modelling,
#' predicting size at time t, from initial size M0, and a set of
#' parameters. The set includes an inverse and differential, which
#' predict age at size M, and growth rates dM/dtime for a size M.
#'
#' @param time time since initial size
#' @param M mass
#' @param M0 initial mass
#' @param K asymptotic size
#' @param r growth rate
#' @param L lowerlimit
#' @param P scaling parameter?
#' 
#' @author Marco D. Visser
#' 
#' @export
GROlogis4k<-function(M0,K,r,L,P,time=1,betas=NULL) {
  if(!is.null(betas)) {
    K <- betas[1]
    r <- betas[2]
    L <- betas[3]
    P <- betas[4]
  }
L+((M0*(K-L))/(M0+P*exp(-r*time)))
}
 

#' @rdname GROlogis4k
#' @author Marco D. Visser
#' 
#' @export
GROilogis4k <- function(M,M0,K,r,L,P,betas=NULL) {
  if(!is.null(betas)) {
    K <- betas[1]
    r <- betas[2]
    L <- betas[3]
    P <- betas[4]
      }
log((M*P)/(K*M0-M*M0)-(L*P)/(K*M0-M*M0))/r
}

#' @rdname GROlogis4k
#' @author Marco D. Visser
#' 
#' @export
GROdlogis4k <- function(M,L,K,r) {
r*(M-L)*((K-M)/(K-L))
}

#' @name GROmonomol
#' @aliases GROimonomol
#' @aliases GROdmonomol
#' @title the GROmonomol set
#' 
#' \code{GROmonomol} - A family of Monomolecular functions for size
#' growth modelling,
#' predicting size at time t, from initial size M0, and a set of
#' parameters. The set includes an inverse and differential, which
#' predict age at size M, and growth rates dM/dtime for a size M.
#'
#' @param time time since initial size
#' @param M mass
#' @param M0 initial mass
#' @param K asymptotic size
#' @param r growth rate
#' @param L lowerlimit
#' @param P scaling parameter?
#' 
#' @author Marco D. Visser
#' 
#' @export
GROmonomol<-function(M0,K,r,time=1,betas=NULL) {
  if(!is.null(betas)) {
    K <- betas[1]
    r <- betas[2] }
K - exp(r*time) * (K- M0)
}
 

#' @rdname GROmonomol
#' @author Marco D. Visser
#' 
#' @export
GROimonomol <- function(M,M0,K,r,betas=NULL) {
  if(!is.null(betas)) {
    K <- betas[1]
    r <- betas[2] }
log((K/(K-M0))-(M/(K-M0)))/r
}

#' @rdname GROmonomol
#' @author Marco D. Visser
#' 
#' @export
GROdmonomol <- function(M,K,r) {
r*(K-M)
}

#' @name GROGompertz
#' @aliases GROiGompertz
#' @aliases GROdGompertz
#' @title the GROGompertz set
#' 
#' \code{GROGompertz} - A family of Gompertz functions for size
#' growth modelling,
#' predicting size at time t, from initial size M0, and a set of
#' parameters. The set includes an inverse and differential, which
#' predict age at size M, and growth rates dM/dtime for a size M.
#'
#' @param time time since initial size
#' @param M mass
#' @param M0 initial mass
#' @param K asymptotic size
#' @param r growth rate
#' @param L lowerlimit
#' @param P scaling parameter?
#' 
#' @author Marco D. Visser
#' 
#' @export
GROGompertz<-function(M0,K,r,time=1,betas=NULL) {
  if(!is.null(betas)) {
    K <- betas[1]
    r <- betas[2] }
K*(M0/K)^exp(-r*time)
}
 
#' @rdname GROGompertz
#' @author Marco D. Visser
#' 
#' @export
GROiGompertz <- function(M,M0,K,r,betas=NULL) {
  if(!is.null(betas)) {
    K <- betas[1]
    r <- betas[2] }
log(log(M0/K)/log(M/K))/r
}

#' @rdname GROGompertz
#' @author Marco D. Visser
#' 
#' @export
GROdGompertz <- function(M,K,r) {
r*M*log(K/M)
}

#' @name GRAlogis2k
#' @aliases GRAilogis2k
#' @aliases GRAdlogis2k
#' @title the GRAlogis2k set
#' 
#' \code{GRAlogis2k} - A family of 2 parameter logistic functions for 
#' gradient modelling. These models are usually used to predict frequency,
#'  a probabily or fractions over x (distance, time, concentration). The 2k
#' version is bound to the range [0,1].
#' The set includes an inverse and differential, which
#' give the inverse estimate (x at a given y) or
#' dy/dx rates over x.
#'
#' @param x the gradient (size, time, concentration)
#' @param y prediction, used in the inverse form 
#' @param beta0 intercept
#' @param beta1 rate of increase
#' @param betas vector to input all parameters in one go
#' @author Marco D. Visser
#' 
#' @export
GRAlogis2k<-function(x,beta0,beta1,beta3,betas=NULL) {
if(is.null(betas)) {betas<-c(beta0,beta1)}
plogis(betas[1]+betas[2]*x)
}


#' @rdname GRAlogis2k
#' @author Marco D. Visser
#' 
#' @export
GRAilogis2k <- function(y,beta0,beta1) {
plogis((qlogis(y)-beta0)/beta[1])
}

#' @rdname GRAlogis2k
#' @author Marco D. Visser
#' 
#' @export
GRAdlogis2k <- function(M,r) {
print("NOT IMPLEMNETED")
}


#' @name GRAlogis3k
#' @aliases GRAilogis3k
#' @aliases GRAdlogis3k
#' @title the GRAlogis3k set
#' 
#' \code{GRAlogis3k} - A family of 3 parameter logistic functions for 
#' gradient modelling. These models are usually used to predict frequency,
#'  a probabily or fractions over x (distance, time, concentration) though
#' though they are not limited to the range [0,1].
#' The set includes an inverse and differential, which
#' give the inverse estimate (x at a given y) or
#' dy/dx rates over x.
#'
#' @param x the gradient (size, time, concentration)
#' @param y prediction, used in the inverse form 
#' @param beta0 intercept
#' @param beta1 rate of increase
#' @param beta3 Asymtotic rate
#' @param betas vector to input all parameters in one go
#' @author Marco D. Visser
#' 
#' @export
GRAlogis3k<-function(x,beta0,beta1,beta3,betas=NULL) {
if(is.null(betas)) {betas<-c(beta0,beta1,beta3)}
betas[3]/(1+exp((betas[1]-x)/betas[2]))
}


#' @rdname GRAlogis3k
#' @author Marco D. Visser
#' 
#' @export
GRAilogis3k <- function(M,M0,r) {
print("NOT IMPLEMNETED")
}

#' @rdname GRAlogis3k
#' @author Marco D. Visser
#' 
#' @export
GRAdlogis3k <- function(M,r) {
print("NOT IMPLEMNETED")
}

#' @name GRAlogis4k
#' @aliases GRAilogis4k
#' @aliases GRAdlogis4k
#' @title the GRAlogis4k set
#' 
#' \code{GRAlogis4k} - A family of 4 parameter logistic functions for 
#' gradient modelling. These models are usually used to predict frequency,
#'  a probabily or fractions over x (distance, time, concentration).
#' The 4 parameter model is limited to the range [beta4,beta3].
#' The set includes an inverse and differential, which
#' give the inverse estimate (x at a given y) or
#' dy/dx rates over x.
#'
#' @param x the gradient (size, time, concentration)
#' @param y prediction, used in the inverse form 
#' @param beta0 intercept
#' @param beta1 rate of increase
#' @param beta3 Asymtotic rate
#' @param beta4 lower bound
#' @param betas vector to input all parameters in one go
#' @author Marco D. Visser
#' 
#' @export
GRAlogis4k<-function(x,beta0,beta1,beta3,beta4,betas=NULL) {
if(is.null(betas)) {betas<-c(beta0,beta1,beta3,beta4)}
betas[3]+(betas[4]-betas[3])/(1+exp((betas[1]-x)/betas[2]))
}


#' @rdname GRAlogis4k
#' @author Marco D. Visser
#' 
#' @export
GRAilogis4k <- function(M,M0,r) {
print("NOT IMPLEMNETED")
}

#' @rdname GRAlogis4k
#' @author Marco D. Visser
#' 
#' @export
GRAdlogis4k <- function(M,r) {
print("NOT IMPLEMNETED")
}

#' @name GRASymDoubleLogis
#' @aliases GRAiSymDoubleLogis
#' @aliases GRAdSymDoubleLogis
#' @title the GRASymDoubleLogis set
#' 
#' \code{GRASymDoubleLogis} - A family of double symmetric logistic functions for 
#' gradient modelling. These models are usually used to predict frequency,
#'  a probabily or fractions over x (distance, time, concentration).
#' The 4 parameter model is limited to the range [beta4,beta3].
#' The set includes an inverse and differential, which
#' give the inverse estimate (x at a given y) or
#' dy/dx rates over x.
#'
#' @param x the gradient (size, time, concentration)
#' @param y prediction, used in the inverse form 
#' @param beta0 intercept
#' @param beta1 rate of increase
#' @param beta3 second intercept
#' @param betas vector to input all parameters in one go
#' @author Marco D. Visser
#' 
#' @export
GRASymDoubleLogis<-function(x,beta0,beta1,beta3,betas=NULL) {
if(is.null(betas)) {betas<-c(beta0,beta1,beta3)}
plogis(beta0+beta1*x)*plogis(beta3-beta1*x)
}


#' @rdname GRASymDoubleLogis
#' @author Marco D. Visser
#' 
#' @export
GRAiSymDoubleLogis <- function(M,M0,r) {
print("NOT IMPLEMNETED")
}

#' @rdname GRASymDoubleLogis
#' @author Marco D. Visser
#' 
#' @export
GRAdSymDoubleLogis <- function(M,r) {
print("NOT IMPLEMNETED")
}

#' @name GRAFreeDoubleLogis
#' @aliases GRAiFreeDoubleLogis
#' @aliases GRAdFreeDoubleLogis
#' @title the GRAFreeDoubleLogis set
#' 
#' \code{GRAFreeDoubleLogis} - A family of double logistic functions used in  
#' gradient modelling. These models are usually used to predict frequency,
#'  a probabily or fractions over x (distance, time, concentration).
#' The 4 parameter model is limited to the range [beta4,beta3].
#' The set includes an inverse and differential, which
#' give the inverse estimate (x at a given y) or
#' dy/dx rates over x.
#'
#' @param x the gradient (size, time, concentration)
#' @param y prediction, used in the inverse form 
#' @param beta0 intercept
#' @param beta1 rate of increase
#' @param beta3 second intercept
#' @param beta4 second rate of increase 
#' @param betas vector to input all parameters in one go
#' @author Marco D. Visser
#' 
#' @export
GRAFreeDoubleLogis<-function(x,beta0,beta1,beta3,beta4,betas=NULL) {
if(is.null(betas)) {betas<-c(beta0,beta1,beta3,beta4)}
plogis(beta0+beta1*x)*plogis(beta3+beta4*x)
}


#' @rdname GRAFreeDoubleLogis
#' @author Marco D. Visser
#' 
#' @export
GRAiFreeDoubleLogis <- function(M,M0,r) {
print("NOT IMPLEMNETED")
}

#' @rdname GRAFreeDoubleLogis
#' @author Marco D. Visser
#' 
#' @export
GRAdFreeDoubleLogis <- function(M,r) {
print("NOT IMPLEMNETED")
}
