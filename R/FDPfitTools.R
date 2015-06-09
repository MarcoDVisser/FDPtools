############################################################
# Helper functions for fitting models to FDP datasets
# Marco Visser, Gamboa, Panama, February 2014
############################################################

#' Fit and Converge a list of Jags models
#' 
#' \code{FitConv} - runs a Jags model, or list of models,
#' and checks for convergence using gelman stats. It then
#' returns a list of MCMC samples, a jags model and DIC values.
#' 
#' @param JAGSdata A pre-prepared jags data object
#' @param modellist A vector containing the location and names
#' of the jags model files. These files should be read in R,
#' as a modelstring. These files must also generate a vector of
#' model parameter names to be returned (in a object called interestlist)
#' @param gelmanthr Gelman threshold values to except fit
#' @param StepIt what is the standard number of MCMC samples to
#' take before checking for convergence?
#' @param MaxIt What is the maximum number of MCMC samples
#' @param FinalSampleLength How many posterior sample do you need
#' when the model is returned?
#' @param ... additional parameters to be passed to \code{jags.model}
#' @param hlist Use this when you have a list of higher parameters
#' that should converge first
#' before convergence of lower parameters (in interestlist) should be
#' checked.
#' @author Marco D. Visser
#' 
#' @export
FitConv <- function(JAGSdata=NULL,modellist=NULL,gelmanthr=1.05,
                    StepIt=1000, MaxIt=1e4, FinalSampleLength=500,
                    hlist=NULL,...) {
  if(is.null(JAGSdata)) {stop("JAGSdata is null, supply data?")}
  if(is.null(modellist)) {stop("modellist is null, supply model(s)?")}
  if(round(MaxIt-StepIt)<=0) {print("MaxIt is too small, changing to
                                     2 x StepIt")
                              MaxIt <-  2 * StepIt
                            }

  ## Number of candidate models?
  Ncandidates <- length(modellist)

  tempmodelsamples <-  vector("list",Ncandidates)
  tempmodelfitstats <-  vector("list",Ncandidates)
  tempmodels <- vector("list",Ncandidates)

  
    for(j in 1:Ncandidates){
                                          #Load Baysesian models
    source(modellist[j])

           model <- jags.model(data=JAGSdata,file
                               = textConnection(modelstring)
                              ,n.adapt=round(StepIt),...)
    gelman <- data.frame(mpsrf=gelmanthr+1)
    CurIt <- StepIt

    if(!is.null(hlist)){
    if(length(hlist)>1){
    while(gelman$mpsrf>=gelmanthr) {
    codasamp <- coda.samples(model,hlist,StepIt)
                                        #Smaller list for saving
    gelman <- gelman.diag(codasamp)
    CurIt <- CurIt + StepIt
    print(paste("Hyper parameters - gelman statistic:", gelman$mpsrf))
    if(CurIt>=MaxIt) {break}
  }
  } else {message("hlist only has one parameter, 
                   skipping convergence check")}

    }

    if(length(interestlist)>1){
    while(gelman$mpsrf>=gelmanthr) {
    codasamp <- coda.samples(model,interestlist,StepIt)
                                        #Smaller list for saving
    gelman <- gelman.diag(codasamp)
    CurIt <- CurIt + StepIt
    print(paste("Gelman statistic:", gelman$mpsrf))
    if(CurIt>=MaxIt) {break}
  }
  } else {message("interestlist only has one parameter, 
                   skipping convergence check")}
    
    codasamp <- coda.samples(model,interestlist,FinalSampleLength)
    DIC<-dic.samples(model,FinalSampleLength)
    tempmodelsamples[[j]] <- codasamp
    tempmodelfitstats[[j]] <- DIC
    tempmodels[[j]] <- model
  }
 
  
  FDPtoolsfitlist <- list(tempmodels,tempmodelsamples,
                          tempmodelfitstats,modellist)
  class(FDPtoolsfitlist) <- c(class(FDPtoolsfitlist),"FDPmodel")

  return(FDPtoolsfitlist)
}

#' Summarize models and fitstatistics 
#'
#' \code{summary.FDPmodel} - summarizes of a list of
#' previously fit models returned by \code{FitConv}
#' 
#' @param modelfit A list of FDPmodel objects returned by FitConv
#' (e.g. FitConv run on multiple species)
#' @param ... additional parameters
#' 
#' @author Marco D. Visser
#' 
#' @export
summary.FDPmodel <- function(modelfit,...) {
  
  if(!is.element("FDPmodel",class(modelfit))){
    stop("modelfit is not of the FDPmodel class")
  }
  
  if(class(modelfit[[1]][[3]][[1]])!='dic') {
    stop("modelfit of unexpected format")
  }
  
  DICs <- t(sapply(modelfit, function(X) {
  sapply(X[[3]], function(y) sum(y$deviance)+sum(y$penalty))
       }))

  colnames(DICs) <- modelfit[[4]]
  bestmod <- apply(DICs,1,function(X) which(X==min(X)))
  secondbest <- apply(DICs,1,function(X) which(X==sort(X)[2]))
  DeltaDic <- sapply(1:nrow(DICs), function(X)
                     DICs[X,secondbest[X]] - DICs[X,bestmod[X]])
  data.frame(DICs) -> DICs
  DICs <- cbind(DICs,bestmod=bestmod,Delta=DeltaDic)
  
  ##Summarize results (alpha + upper lower CI) per species
  
  MeanEst<-lapply(modelfit,function(X)
                  summary(modelfit[[2]][[X]][[bestmod[X]]])$statistics)
  
  CI<-lapply(1:length(sdlhghtsp),function(X)
             summary(modelfit[[X]][[bestmod[X]]])$quantiles)

  return(list('DIC'=DICs,'Parameters'=MeanEst,'CI'=CI))
}



#' Calculate fit statistics 
#' 
#' \code{FitStats} - Calculates various fit statistics, including
#' Delta AIC/DIC/BIC values from a set of stats, AIC weighst and ...
#' list to grow.
#' 
#' @param fitvalues a list of AIC, DIC or BIC values
#' @param stat do you wish to calculate 'delta' values, 'weights' or ...
#' list to grow.
#' when the model is returned?
#' @param ... additional parameters to be passed to \code{jags.model}
#' 
#' @author Marco D. Visser
#' 
#' @export
FitStats <- function(fitvalues,stat='delta'){

  deltas <- fitvalues - sort(fitvalues)[1]
  if(stat=='delta') return(deltas)
  weights <- exp(-.5*deltas)/sum(exp(-.5*deltas))
  if(stat=='weights') return(weights)
}
