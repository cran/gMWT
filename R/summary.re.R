# Version: 30-11-2012, Daniel Fischer

summary.re <- function(object, ...){
  cat("Rejection Summary\n")
  cat("---------------\n")
  cat("Alpha         :",object$alpha,"\n")
  cat("Optimal Alpha :",object$opt.alpha,"\n")
  cat("Obs. Sig.     :",object$obs,"\n")
  cat("Exp. Sig.     :",object$exp,"\n")
  invisible(object)
} 
