# Version: 30-11-2012, Daniel Fischer

plot.re <- function(x,...){
  rejectionPlot(x$X,crit=x$crit)
}