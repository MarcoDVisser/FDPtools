############################################################
## A list of models general purpose jags 
## Marco Visser, Wageningen 09.06.2015
############################################################
#' load a standard Jags models
#' 
#' \code{loadModel} - load a predifined BUGS model
#' @name loadModel
#' @title Set of common Jags models
#' 
#' \code{loadModel} - Contains a set of jags models 
#' which include:
#'  - univariate probability density function (weibull,
#'  lognormal, normal, exponential).
#' - more to come
#' 
#'
#' @param model the selected model
#' @param class the model class (pdf,glm,lm)
#' 
#' @author Marco D. Visser
#' 
#' @export
loadModel <- function(model="weibull", class="pdf"){
  if(class=="pdf"){
    if(model=="weibull"){

modelstring <- 
"
model
{

## Hyperparameters and priors
v ~ dunif(1,10)
lambda ~ dunif(0,100)
## loop over N
for(i in 1:N){
 Y[i] ~ dweib(v,lambda)
	}

}

"
interestlist <- c("v","lambda")
jagsdata <- c("Y","N")

    }
    if(model=="lognormal"){

modelstring <- 
"
model
{

# Hyperparameters and priors
mu ~ dunif(0,10)
prec ~ dunif(0,100)
# loop over N
for(i in 1:N){
  Y[i] ~ dlnorm(mu,prec)
	}

}

"
interestlist <- c("mu","prec")
jagsdata <- c("Y","N")

    }
    if(model=="exponential"){

modelstring <- 
"
model
{

# Hyperparameters and priors
lambda ~ dunif(0,10)

# loop over N
for(i in 1:N){
Y[i] ~ dexp(lambda)
	}

}

"
interestlist <- c("lambda")
jagsdata <- c("Y","N")

    }


  }

return(list(modelstring,interestlist,jagsdata))

}
