############################################################
# Helper functions used to evaluate different meassures of
# local conspecific density.
# Marco Visser, Nijmegen, September 2013
############################################################

#' calculateCrowding
#' 
#' calculates Crowding statistics for each individual alive
#' in a given census year for the BCI plot data for
#' further analysis.
#' 
#' @param FdpData A list of all census years on which the
#'  calculation must be conducted
#' @param r a vector containing all distances at which the
#' below statistic must be calculated
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
                              r=seq(1,50,10),
                              statistic="count",
                              ...) {
  
# Prune data to include only individuals alive
FdpData <- lapply(FdpData, function(X)
                  subset(X,DFstatus=="alive"))

 # calculate a distance matrix for each individual
 # to each other individual
DistMatList <- lapply(FdpData,function(X)
                    calcDistMat(X,xy=c("gx","gy")))

# Classify matrices by statistic
lapply(DistMatlist,function(X) class(X) <- statistic)

# create crowding matrices
CrowdMat<-lapply(FdpData,function(X)
                 matrix(-1.0,nrow=dim(X)[1],
                        ncol=length(r))
                 )

# Next loop over CrowdMat and calculate stats for each r

}


#' calcDistMat
#' 
#' Generic function to calculate distances between all
#' point is a datatset containing x and y coordinates
#' 
#' @param pointdata a dataset containing x and y
#' coordinates
#' @param xy column names containing coordinates
#' 
#' @author Marco D. Visser
#' 
calcDistMat <- function(pointdata,xy=c("x","y")){

  # Limit data to x and y
  pointdata <- pointdata[,is.element(xy,colnames(ponitdatat))]

  DistMat<-matrix(0,nrow=nrow(pointdata))

  #Vectorized over rows so should be pretty fast

  for(i in 1:nrow(pointdata)){
    DistMat[i,] <- sqrt( (pointdata[i,1]-pointdata[,1])^2 +
                        (pointdata[i,2] - pointdata[,2])^2)
  }

  return(DistMat)

}

#' CrowdStat
#' 
#' Generic function which calls specific crowding statistic
#' to calculate crowding for each individual in a dataset
#' for each radius.
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

##' CrowdStat.count
##' 
##' primitive function which counts individual
##' surrounding a given individual in the census data
##' 
##' @param distancemat a matrix containing the distances
##' between all individuals
##' @param r a vector containing all distances at which the
##' below statistic must be calculated
##' 
##' @author Marco D. Visser
##' not exported
CrowdStat.count<-function(distancemat,r){




}
                             
