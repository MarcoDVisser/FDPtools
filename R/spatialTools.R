############################################################
# A List of optimized spatial tools
# Marco Visser, Nijmegen, September 2013
############################################################

#' Calculate distance-based statistics between
#' datasets or within datasets (with spatial data).
#'
#' \code{calcDistStat} - Generic optimized function for
#' calculation of distance based statistics as
#' e.g. total basal area with radii
#' 
#' @param pointdata a dataset containing x and y
#' coordinates. Example, fdp spatial data.
#' @param referencedata a reference data set, used when
#' the distance between all points in the dataset pointdata
#' and all points in the reference data set is needed. If
#' NULL this is ignored and the distance between all
#' points in pointdata is calculated.  Example, reference
#' data can be traplocations or seedling plots.
#' 
#' @param xy column names containing coordinates. Length should be
#' 4 if using a reference dataset, with the first and second elements
#' referering to pointdata and the second and third to the
#' reference dataset.
#' 
#' @param statistic a vector of information about each point in
#' pointdata to summerize. Example: if the statistic contains
#' pi*(dbh/2)^2 values, the sum of basal area within each radii is
#' returned. The values in "statistic" should correspond to
#' pointdata.
#' 
#' @param radii a vector of distance 
#' 
#' @author Marco D. Visser
#' @examples
#' \dontrun{
#' x <- data.frame(x.cor=runif(10),y.cor=runif(10))
#' calcDistMat(x,xy=c("x.cor","y.cor"),statistic=runif(10))
#' } 
#' @export
calcDistStat <- function(pointdata,referencedata=NULL,xy=c("x","y"),
                        statistic=NULL,radii=NULL){

  if(is.null(radii)|is.null(statistic)){
  stop("statistic or radii cannot be NULL")}

    if(length(statistic)<dim(pointdata)[1]){
  stop("pointdata does not conform to statistic")
  }

  
  if(is.null(referencedata)){

  # Limit data to x and y
  pointdata <- as.matrix(pointdata[,is.element(xy,colnames(pointdata))])
 #make ready results
  results <- matrix(nrow=dim(pointdata)[1],ncol=length(radii))
 
   #run through all radii
  for(i in 1:length(radii)){   
   NeighborMat <- findNeighbor(pointdata[,1], pointdata[,2],
                         pointdata[,1], pointdata[,2],radius=radii[i])
  
  resulti <- apply(NeighborMat,1,function(X)
                   sum(statistic*X))
  results[,i] <- resulti
   
  } } else
    {
  # Limit data to x and y
  if(length(xy)!=4) {
    stop("supply the x and y coordinates col.names for both data sets")
                   }

   #make ready results
  results <- matrix(nrow=dim(referencedata)[1],ncol=length(radii))
 
  pointdata <- pointdata[,is.element(colnames(pointdata),xy[1:2])]
  referencedata <- referencedata[,is.element(colnames(referencedata),xy[3:4])]
  pointdata <- as.matrix(pointdata)
  referencedata <- as.matrix(referencedata)

                     #run through all radii
  for(i in 1:length(radii)){   
  
                    
    NeighborMat <-findNeighbor(referencedata[,1], referencedata[,2],
                               pointdata[,1], pointdata[,2],
                               radius=radii[i])
    
  resulti <- apply(NeighborMat,1,function(X)
                   sum(statistic*X))
  results[,i] <- resulti
  
 }  
}

  colnames(results) <- radii
  return(results)

}
                             
# Compile function for speed
calcDistStat<-cmpfun(calcDistStat)

##' calcDist a lower-level helper function
##'
##' Wrapper function which calls some
##' underlying C code.
##'
##' @param x1 x coordinates of points
##' @param y1 y coordinates of points
##' @param x2 reference point x coordinate
##' @param y2 reference point y coordinate 
##'    
##' @export
##' @useDynLib FDPtools    
calcDist <- function(x1,y1,x2,y2) {

  if(length(x2)<1|length(y2)<1){
  stop("Input of incorrect length")
  }

 out <- .Call("calcDist",
           x1 = as.double(x1),
           y1 = as.double(y1),
           x2 = as.double(x2),
           y2 = as.double(y2)
          )
return(out)
 }

##' findNeighbor a lower-level helper function
##'
##' Wrapper function which calls some
##' underlying C code and returns a
##' matrix with 1 or 0 indicating
##' whether the distance between the points
##' is smaller or equal to radius.
##'
##' @param x1 x coordinates of points
##' @param y1 y coordinates of points
##' @param x2 reference point x coordinate
##' @param y2 reference point y coordinate 
##' @param radius radius with which to test distance
##' @export 
findNeighbor <- function(x1,y1,x2,y2,radius) {

  if(length(x2)<1|length(y2)<1){
  stop("Input of incorrect length")
  }

 out <- .Call("Rneighbor",
           x1 = as.double(x1),
           y1 = as.double(y1),
           x2 = as.double(x2),
           y2 = as.double(y2),
           radius = as.double(radius)
          )
return(out)
 }
