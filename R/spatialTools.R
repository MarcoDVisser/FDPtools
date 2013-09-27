############################################################
# A List of optimized spatial tools
# Marco Visser, Nijmegen, September 2013
############################################################

#' Calculate distances between points within and between
#' datasets
#'
#' \code{calcDistMat} - Generic optimized function for
#' calculation of 1) distances between all points in a datatset
#' containing x and y coordinates or 2) between all points in
#' two datasets containing x and y coordinates.
#' 
#' @param pointdata a dataset containing x and y
#' coordinates.
#' @param referencedata a reference data set, used when
#' the distance between all points in the dataset pointdata
#' and all points in the reference data set is needed. If
#' NULL this is ignored and the distance between all
#' points in pointdata is calculated.
#' 
#' @param xy column names containing coordinates. Length should be
#' 4 if using a reference dataset, with the first and second elements
#' referering to pointdata and the second and third to the
#' reference dataset.
#' 
#' @author Marco D. Visser
#' @examples
#' \dontrun{
#' x <- data.frame(x.cor=runif(10),y.cor=runif(10))
#' calcDistMat(x,xy=c("x.cor","y.cor"))
#' } 
#' @export
calcDistMat <- function(pointdata,referencedata=NULL,xy=c("x","y")){


  if(is.null(referencedata)){
    # Limit data to x and y
  pointdata <- as.matrix(pointdata[,is.element(xy,colnames(pointdata))])

  DistMat<-matrix(0,nrow=nrow(pointdata),
                  ncol=nrow(pointdata))

  #Vectorized over rows so should be pretty fast
  for(i in 1:nrow(pointdata)){
    DistMat[i,] <- calcDist(pointdata[,1], pointdata[,2],
                            pointdata[i,1], pointdata[i,2])
                             }
} else {

  # Limit data to x and y
  if(length(xy)!=4) {
    stop("supply the x and y coordinates col.names for both data sets")
                   }
  
  pointdata <- pointdata[,is.element(colnames(pointdata),xy[1:2])]
  referencedata <- referencedata[,is.element(colnames(referencedata),xy[3:4])]
  pointdata <- as.matrix(pointdata)
  referencedata <- as.matrix(referencedata)

DistMat<-matrix(ncol=nrow(pointdata),
                nrow=nrow(referencedata))

    for(i in 1:nrow(pointdata)){

      DistMat[,i] <- sqrt((pointdata[i,1]-referencedata[,1])^2 +
                        (pointdata[i,2] - referencedata[,2])^2)   
    }
}
  
  return(DistMat)

}
                             
# Compile function for speed
calcDistMat<-cmpfun(calcDistMat)

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
calcDist <- function(x1,y1,x2,y2) {

  if(length(x2)>1|length(y2)>1){
  stop("Input of incorrect length")
  }

 out <- .C("calcDist",
           x1 = as.double(x1),
           y1 = as.double(y1),
           x2 = as.double(x2),
           y2 = as.double(y2),
           r  = as.double(x1),
           n  = as.integer(length(x1))	   		   
           )

  return(out$r)
}
