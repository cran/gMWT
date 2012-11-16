`print.re` <- function(x,...){
 X <- list()
 X$repMat <- x$repMat
 X$opt.alpha <- x$opt.alpha
 print(X,...)
}