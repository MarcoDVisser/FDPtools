############################################################
# Helper functions used to organise the FDP datasets
# Marco Visser, Gamboa, Panama, February 2014
############################################################

#' Prepare Growth and Survival datasets 
#' 
#' \code{prepareFDPdata} - prepares raw data for analysis in two
#' posible formats for each individual for the BCI
#' plot data. 
#' 
#' @param FDPobjects A vector of fdp census objects stored in
#'  the global environment
#' @param census Which census is included? The "tree" census (> 1cm dbh)
#' or the "seedling" census.
#' @param type Which format is desired? Can be either 'paired' or 'trajectory'
#' A paired format, is where each consecutive census is paired. This is the
#' standard growth analysis format. Or a trajectory format, where all data are
#' organized next to each other (this is for growth trajectory analysis).
#' Only trajectory format is allowed for seedling data.
#' @param SpFitList A vector of 6 char species codes, if null all species
#' are used.
#' 
#' @author Marco D. Visser
#' 
#' @export
prepareFDPdata <- function(FDPobjects = ls()[grep("bci.full",ls())][3:7],
                           census = "tree", type = 'paired',SpFitList=NULL) {

  if(census=='seedling'){
    if(type=='paired') {stop("paired format not available for seedlings")}

    ## combine Joe and Liza's census data in one format
    if(is.null(SpFitList)){SpFitList <- unique(get(FDPobjects[1])$SPP)}
  prepdata <- lapply(FDPobjects, function(X) subset(get(X),
                                                     SPP%in%toupper(SpFitList)))

      ## subset data and only include what is needed
  
   censusdata <- lapply(prepdata,function(X)
                     X[,c("SPP","TAG","plotnum","alt","dbh","status","date")])

    
  ## extract data
    dbhdata <- do.call(cbind,lapply(censusdata,function(X) X[,"dbh"]))
    altdata <- do.call(cbind,lapply(censusdata,function(X) X[,"alt"]))
    survivaldata <- do.call(cbind,lapply(censusdata,function(X)
                                         as.character(X[,"status"])))

    datedata <- do.call(cbind,lapply(censusdata,function(X) as.Date(X[,"date"])))

    ## rename columns
    colnames(dbhdata) = paste0("dbh",1:dim(dbhdata)[2])
    colnames(survivaldata) = paste0("survival",1:dim(survivaldata)[2])
    colnames(altdata) = paste0("alt",1:dim(altdata)[2])
    colnames(datedata) = paste0("date",1:dim(datedata)[2])
    ## get meta data
    metadata <- censusdata[[1]][,c("TAG","SPP")]
    growclean <- cbind(metadata,dbhdata,altdata,survivaldata,datedata)
    class(growclean) <- c(class(growclean),"fdpdata",type,census)
  
  }

  
 if(census=='tree'){
  if(is.null(SpFitList)){SpFitList <- unique(get(FDPobjects[1])$sp)}
   prepdata <- lapply(FDPobjects, function(X) subset(get(X),
                                                    sp%in%tolower(SpFitList)))

  ## subset data and only include what is needed
  
   censusdata <- lapply(prepdata,function(X)
                     X[,c("sp","tag","dbh","DFstatus","date","codes")])
  

 if(type=='paired')  {
    
  
  censuspairs <- vector("list",length(censusdata)-1)
  
  for(i in 1:length(censuspairs)){
    censuspairs[[i]] <- c((1:(length(FDPobjects)-1))[i],(2:length(FDPobjects))[i])
  }
  
  tempdatapairs <- lapply(censuspairs, function(X) {
    cbind(censusdata[[X[1]]],censusdata[[X[2]]])
  })
  
  censustracker <- lapply(censuspairs, function(X) {
    rep(paste(X,collapse=''), nrow(censusdata[[X[1]]]))
  })
  
  finalgrowthlong <- do.call(rbind,tempdatapairs)
  names(finalgrowthlong) <- c(paste(names(censusdata[[1]]),1,sep=''),
                              paste(names(censusdata[[1]]),2,sep=''))
  finalgrowthlong$census <- as.factor(unlist(censustracker))
  rm(tempdatapairs)
  
  growclean<-finalgrowthlong
  class(growclean) <- c(class(growclean),"fdpdata",type,census)
  
} else {

  ## extract data
  dbhdata <- do.call(cbind,lapply(censusdata,function(X) X[,"dbh"]))
  survivaldata <- do.call(cbind,lapply(censusdata,function(X) X[,"DFstatus"]))
  codedata <- do.call(cbind,lapply(censusdata,function(X) X[,"codes"]))
  datedata <- do.call(cbind,lapply(censusdata,function(X) X[,"date"]))
  ## rename columns
  colnames(dbhdata) = paste0("dbh",1:dim(dbhdata)[2])
  colnames(survivaldata) = paste0("survival",1:dim(survivaldata)[2])
  colnames(codedata) = paste0("code",1:dim(codedata)[2])
  colnames(datedata) = paste0("date",1:dim(datedata)[2])
  ## get meta data
  metadata <- censusdata[[1]][,c("tag","sp")]
  growclean <- cbind(metadata,dbhdata,survivaldata,datedata,codedata)
  class(growclean) <- c(class(growclean),"fdpdata",type,census)
}
}
  return(growclean)
}
# Compile function for speed
prepareFDPdata<-cmpfun(prepareFDPdata)

#' Extract Relavent Data
#' 
#' \code{extractData} - extracts growth, survival or a combined dataset from
#' from objects previously prepared by \code{prepareFDPdata}
#' 
#' @param FDPdata An object return by \code{prepareFDPdata}
#' @param type Which format is desired? Can be either 'growth', 'survival',
#' 'resprouter' or a paired format "paired" - in which each growth
#' and survival information is paired.
#' The paired format must be created with \code{prepareFDPdata} using the option
#' type='trajectory'.
#' @param datestan A date standardizer, when datestan = 21 the first cenus
#' dates are set to +- 0. This is only relavent for trajectory modelling.
#' 
#' @author Marco D. Visser
#' 
#' @export
extractData <- function(FDPdata=NULL,type='growth', datestan=21) {

  if(is.element("tree",class(FDPdata))){
    if(type=='growth'&is.element("paired",class(FDPdata))) {
      
      ##remove individuals that died in a second census
      
      finalgrowthlong <- subset(FDPdata,DFstatus2=="alive")
      
      finalgrowthlong$years <- (finalgrowthlong$date2-
                                finalgrowthlong$date1)/365.25
      finalgrowthlong$grw <- ((finalgrowthlong$dbh2-finalgrowthlong$dbh1)
                              /finalgrowthlong$years)
      
      ## remove  broken and dying individuals (resprouts)
      growclean<-subset(finalgrowthlong, !is.na(dbh1)&!is.na(dbh2))
      return(growclean)
    }

  if(type=='survival'&is.element("paired",class(FDPdata))) {
    survivallong <- subset(FDPdata,DFstatus1=='alive')
    survivallong$surv <- survivallong$DFstatus2=='alive'
    survivallong$T <- (survivallong$date2-survivallong$date1)/356.25
    survclean <- subset(survivallong,!is.na(dbh1))
    return(survclean)
  }

  if(type=='resprouter'&is.element("paired",class(FDPdata))) {
    ## Get reprouter transition

    ##remove individuals that died in a second census
    finallong <- subset(FDPdata,DFstatus2!="dead"&DFstatus2!='missing'&
                        codes1!='R')
    
    finallong$years <- (finalgrowthlong$date2-
                        finalgrowthlong$date1)/365.25
    ## id resprouters
    finallong$R <- ifelse(finallong$codes2=='R')
    
  }
  
  if(is.element("trajectory",class(FDPdata))) {
    ##How many census points?
    Ncen <- length(grep("dbh",names(FDPdata)))

    ##Kick out all trees dead on the first census if any
    deadonarrival <- FDPdata[,paste0("survival",1)]=="dead"
    FDPdata <- FDPdata[!deadonarrival,]
 
  }
  if(type=='resprouter'&is.element("trajectory",class(FDPdata))) {

    ## 1)  Resprouter cohorts
    ## find resprouters based on code (as this is most reliable
    codesub<-FDPdata[,paste0("code",1:5)]

    ## remove individuals that were already resprouting in the
    ## first census
    cen1resprout <- codesub[,1]=='R'
    cen1resprout[is.na(cen1resprout)] <- FALSE
    resprouters <- FDPdata[!cen1resprout,]

    ## limit data to the rest of the resprouters
    Rid <- apply(resprouters[,paste0("code",1:5)],1
                 ,function(x) is.element("R",x))

    resprouters <- resprouters[Rid,]

    ## Id start and end of resprouter status
    codesub<-resprouters[,paste0("code",1:5)]
    firstresprout <- suppressWarnings(sapply(1:dim(codesub)[1],function(X) 
                                             min(which(codesub[X,]=='R'))))
    resprouter$fr <- firstresprout

    ## finalize
    return(resprouter)
  }

  
  if(is.element(type,c('survival','growth','paired'))&
                is.element("trajectory",class(FDPdata))) {
    
    ## find and drop all reprouters, these are treated differently
    ## once an individual respouts,it is always considered a reprouter
    ## and it is treated in a seperate analysis
    resprouter <- FDPdata[,paste0("survival",1:Ncen)]=="lost_stem"
    ## when did it first respout?
    firstresprout <- suppressWarnings(sapply(1:dim(resprouter)[1],function(X)
                                            min(which(resprouter[X,],arr.ind=T))))

    firstresprout[!is.finite(firstresprout)]<-NA

    ## Find missing individuals that are actually resprouters
    missing <- FDPdata[,paste0("survival",1:Ncen)]=="missing"

    
    ## run through all respouters and fix
    correction <-resprouter[!is.na(firstresprout),]

    pb <- txtProgressBar(min = 0, max = dim(correction)[1], style = 3)    
    for(i in 1:dim(correction)[1]){
      correction[i,][na.omit(firstresprout)[i]:Ncen]<-TRUE
      setTxtProgressBar(pb, i)
    }
    close(pb)

    ## id missing

    resprouter[!is.na(firstresprout),]<-correction
    
    # remove respouter information from data
    FDPdata[,paste0("dbh",1:Ncen)][resprouter] <- NA
    ## which data points have information in them?
    infopoints <- !is.na(FDPdata[,paste0('dbh',1:Ncen)])

    if(type=='growth'){
    ## dicard those with only one or less points 
    growdatpoints <- apply(infopoints,1,sum)
    growclean <- subset(FDPdata,growdatpoints>1)
  } else {
    ## dicard those no points of meassurement (eg. missing)
    growdatpoints <- apply(infopoints,1,sum)
    growclean <- subset(FDPdata,growdatpoints>=1)
    
  }

    
    ## update infopoints
    infopoints <- !is.na(growclean[,paste0('dbh',1:Ncen)])
    ##id the first and last census with information
    growclean$start <- sapply(1:dim(infopoints)[1],function(X)
                              min(which(infopoints[X,],arr.ind=T)))
    growclean$end <- sapply(1:dim(infopoints)[1],function(X)
                            max(which(infopoints[X,],arr.ind=T)))
    ## stardardize dates
    dates <-  growclean[,paste0('date',1:Ncen)]/365.25
    growclean[,paste0('date',1:Ncen)]<-dates
    growclean$startdate <- sapply(1:nrow(growclean),
                                  function(X) dates[X,growclean$start[X]])

    ## remove those with missing dates, as dates cannot be sampled in the MCMC
    ## these are a handful of missing individuals
    nainc=apply(is.na(growclean[,paste0('date',1:5)]),1,sum)>1
    growclean <- growclean[!nainc,]
    
    if(is.element(type,c('paired','survival'))){
      ## transform survival columns into binary data
      survdata <- growclean[,paste0('survival',1:Ncen)]
      survdata <- ifelse(survdata!='dead',1,0)
      ## shift survival to correspond to DBH t-1
      survdata[,1:(Ncen-1)] <- survdata[,2:Ncen]
      ## no information on survival beyond the census
      survdata[,Ncen] <- NA
      ## finalize
      survdata -> growclean[,paste0('survival',1:Ncen)]

      if(type=='survival'){
        ## Make time since previous census
        survdates <- growclean[,paste0('date',1:Ncen)]
        ## Shift dates to correspond to t-1 and subtract from date t
        survdates[,1:Ncen] <- (survdates[,c(2:Ncen,5)] - survdates[,1:Ncen])
        ## no information on survival beyond the census
        survdates[,Ncen] <- 0
        ## repair names
        colnames(survdates) <- paste0('date',1:Ncen)
        ## finalize
        survdates -> growclean[,paste0('date',1:Ncen)]
        
        return(growclean)}

      if(type=='paired'){
        ## Make time since previous census
        survdates <- growclean[,paste0('date',1:Ncen)]
        ## Shift dates to correspond to t-1 and subtract from date t
        survdates[,1:Ncen] <- (survdates[,c(2:Ncen,5)] - survdates[,1:Ncen])
        ## no information on survival beyond the census
        print(1)
        survdates[,Ncen] <- 0
        ## rename to T
        colnames(survdates) <- paste("T",1:Ncen,sep='')
        print(1)
        ## finalize
        pairdata <- cbind(growclean,survdates)
        return(pairdata)}

    } else {return(growclean)}
  }


    if(type=='paired'&!is.element("trajectory",class(FDPdata))) {
    stop("paired option only valid on trajectory data")
  }

}
  else{
    if(!is.element("seedling",class(FDPdata))){
      stop('expected FDPdata of class seedling')}

    finalseedlings <- extractDataSeedling(FDPdata=FDPdata,type=type) 
    return(finalseedlings)

  }

}



#' Extract Relavent Data from seedling database
#' 
#' \code{extractDataSeeding} - extracts growth, survival or a combined dataset
#' from seedling data previously prepared by \code{prepareFDPdata}
#' 
#' @param FDPdata An object return by \code{prepareFDPdata}
#' @param type Which format is desired? Can be either 'growth', 'survival',
#' 'resprouter' or a paired format "paired" - in which each growth
#' and survival information is
#' The paired format must be created with \code{prepareFDPdata} using the option
#' type='trajectory'.
#' @param datestan A date standardizer, when datestan = 21 the first cenus
#' dates are set to +- 0. This is only relavent for trajectory modelling.
#' 
#' @author Marco D. Visser
#' 
#'@export
extractDataSeedling <- function(FDPdata=NULL,type='growth') {

  Ncen <- length(grep("alt",names(FDPdata)))

  ## make all -2 and -9 into NA
  heightdata <- FDPdata[,paste0('alt',1:Ncen)]
  makeNApoints <- (is.na(heightdata)|heightdata<=0)
  heightdata[makeNApoints] <- NA
  ## overwrite original
  heightdata -> FDPdata[,paste0('alt',1:Ncen)]
  
   ## which data points have information in them?
  heightdata <- FDPdata[,paste0('alt',1:Ncen)]
  infopoints <- !(is.na(heightdata))
  
    if(type=='growth'){
    ## dicard those with only one or less points 
    growdatpoints <- apply(infopoints,1,sum)
    growclean <- subset(FDPdata,growdatpoints>1)
  } else {
    ## dicard those no points of meassurement (eg. missing)
    growdatpoints <- apply(infopoints,1,sum)
    growclean <- subset(FDPdata,growdatpoints>=1)
     }

    
  ## update infopoints
  heightdata <- growclean[,paste0('alt',1:Ncen)]
  infopoints <- !(is.na(heightdata))
  ##id the first and last census with information
  growclean$start <- sapply(1:dim(infopoints)[1],function(X)
                            min(which(infopoints[X,],arr.ind=T)))
  growclean$end <- sapply(1:dim(infopoints)[1],function(X)
                          max(which(infopoints[X,],arr.ind=T)))
  ## stardardize dates
  dates <-  growclean[,paste0('date',1:Ncen)]/365.25
  ## replace those with missing dates, as dates cannot be sampled in the MCMC
  ## replace with average censusdate
  averagedates <- colMeans(dates,na.rm=T)
  dates <- sapply(1:dim(dates)[2],function(X)
                  ifelse(!is.na(dates[,X]),dates[,X],averagedates[X]))

  growclean[,paste0('date',1:Ncen)]<-dates
  growclean$startdate <- sapply(1:nrow(growclean),
                                function(X) dates[X,growclean$start[X]])

  
  
  if(is.element(type,c('paired','survival'))){
    ## transform survival columns into binary data
    survdata <- growclean[,paste0('survival',1:Ncen)]
    survdata <- ifelse(survdata=='A',1,0)
    ## shift survival to correspond to DBH t-1
    survdata[,1:(Ncen-1)] <- survdata[,2:Ncen]
    ## no information on survival beyond the census
    survdata[,Ncen] <- NA
    ## finalize
    survdata -> growclean[,paste0('survival',1:Ncen)]

      if(type=='survival'){
        ## Make time since previous census
        survdates <- growclean[,paste0('date',1:Ncen)]
        ## Shift dates to correspond to t-1 and subtract from date t
        survdates[,1:Ncen] <- (survdates[,c(2:Ncen,5)] - survdates[,1:Ncen])
        ## no information on survival beyond the census
        survdates[,Ncen] <- 0
        ## repair names
        colnames(survdates) <- paste0('date',1:Ncen)
        ## finalize
        survdates -> growclean[,paste0('date',1:Ncen)]
        
        return(growclean)}

      if(type=='paired'){
        ## Make time since previous census
        survdates <- growclean[,paste0('date',1:Ncen)]
        ## Shift dates to correspond to t-1 and subtract from date t
        survdates[,1:Ncen] <- (survdates[,c(2:Ncen,5)] - survdates[,1:Ncen])
        ## no information on survival beyond the census
        print(1)
        survdates[,Ncen] <- 0
        ## rename to T
        colnames(survdates) <- paste("T",1:Ncen,sep='')
        print(1)
        ## finalize
        pairdata <- cbind(growclean,survdates)
        return(pairdata)}

    } else {return(growclean)}

}

extractData<-cmpfun(extractData)
