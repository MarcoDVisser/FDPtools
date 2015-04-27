FDPtools
========

Optimized tools for working with forest dynamics plot data. The package is an 
attempt to create many optimized analytical tools for calculating statistics
of interest from forest dynamics plot data. Tools range from spatail tools such
as distance calculation to Bayesian variable selection and inverse modelling.
The package is currently under development and contains source code in C 
that has only been tested under linux. Windows users should take care.

## Installation

You can download the [development version of FDPtools]
(https://github.com/MarcoDVisser/FDPtools) as [zip](https://github.com/MarcoDVisser/FDPtools/zipball/master) or [tar ball](https://github.com/MarcoDVisser/FDPtools/tarball/master). Then `R CMD INSTALL` on it, or use the **devtools** package to install the development version:

```r
## Make sure your current packages are up to date
update.packages()
## devtools is required
library(devtools)
install_github("FDPtools", "MarcoDVisser")
```
As the package depends on compiled code, Microsoft Windows users will need to install [Rtools](http://www.murdoch-sutherland.com/Rtools/) first.

## Examples 
```r
# Make spatial data
require(FDPtools)
x1<-runif(10000)
y1<-runif(10000)
x2<-runif(1000)
y2<-runif(1000)
#pure R way to calculate distance
a <- matrix(ncol=1000,nrow=1000)
system.time(for(i in 1:100) a[i,]<-sqrt((x1[i]-x2)^2 + (y1[i]-y2)^2))
#using FDPtools
system.time(b<-calcDist(x1,y1,x2,y2))

```

## Bugs
* bug-reports can be submitted through: <https://github.com/MarcoDVisser/FDPtools/issues>

##Contact
* emails can be sent to: <marco.d.visser@gmail.com>
