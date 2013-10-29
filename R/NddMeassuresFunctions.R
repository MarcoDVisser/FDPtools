############################################################
# Helper functions used to evaluate different meassures of
# local conspecific density.
# Marco Visser, Nijmegen, September 2013
############################################################

#' Calculate crowding statistics
#' 
#' \code{calculateCrowding} - Calculates crowding statistics
#' for each individual alive, based on the surrounding
#' environment, in each included census for the BCI
#' plot data. Alternatively, if a second data set is included
#' (refdata) it calculates crowding statistics for all points
#' in the second dataset based on all points in each FDP census
#' included. 
#' 
#' @param FdpData A list of fdp census datasets on which the
#'  calculation must be conducted
#' @param refdata a second dataset including points for
#' which crowding statistics based on the FdpData must
#' be calculated. This option is ignored if NULL.
#' @param r a vector containing all distances at which the
#' below statistic must be calculated. Used if relevant.
#' @param statsistic character defining which method to use.
#' Methods include: 1) "count", simple count of neighbours
#' within r 2) "distance index", a distance weighted function
#' 3) "gravitation", a size and distance weighted function
#' @param ... further arguments to be passed to statistic
#' calculation    
#' 
#' @author Marco D. Visser
#' 
#' @export
calculateCrowding <- function(FdpData=list(bci.full6,bci.full7),
                              refdata=NULL,
                              r=seq(1,50,10),
                              statistic="count",
                              ...) {
  
# Prune data to include only individuals alive
FdpData <- lapply(FdpData, function(X)
                  subset(X,DFstatus=="alive"))

# Classify FdpData by statistic
 for(i in 1:length(FdpData))   {
class(FdpData[[i]]) <- c(paste("ndd",statistic,sep=""),"data.frame")
                 }

if(is.null(refdata)){                  
# create crowding metrices for each dataset
CrowdData<-lapply(FdpData,function(X) CrowdStat(X,r,...))
}

# return crowding data
return(CrowdData)
}

#' CrowdStat
#' 
#'  \code{CrowdStat} Generic function which calls specific
#' crowding statistic to calculate crowding for each individual
#' in a dataset for each radius.
#' 
#' @param censusdata a single census year, with only
#' individuals alive within this census included.
#' @param r a vector containing all distances at which the
#' below statistic must be calculated
#' @param ... further arguments to be passed to statistic
#' calculation    
#' 
#' @author Marco D. Visser
#' 
#' @export
CrowdStat <- function(censusdata,r,...){
UseMethod("CrowdStat")
}

##' CrowdStat.nddcount
##' 
##' primitive function which counts points (trees)
##' surrounding a given focal point (tree) in the census data
##' 
##' @param distancemat a matrix containing the distances
##' between all individuals
##' @param r a vector containing all distances at which the
##' below statistic must be calculated
##' 
##' @author Marco D. Visser
##' @export

CrowdStat.ndd<-function(FdpData,r){

    CrowdMat<-matrix(nrow=nrow(FdpData),ncol=length(r))


    for(i in 1:nrow(FdpData)){
   Distvec <- calcDistMat(FdpData,FdpData[i,],xy=c("gx","gy","gx","gy"))
   CrowdMat[i,] <- sapply(r,function(X) sum(X<=Distvec))
    }

    return(CrowdMat)
}
                             
