plot.estPI <- function(x,col="black",highlight=NULL,hlCol="red",pch=20,zoom=FALSE,...){
  if(x$type=="single")
  {
    estPlotSingle(x,col=col,highlight=highlight,hlCol=hlCol,pch=pch,zoom=zoom,...)
  } else if(x$type=="pair")
  {
    estPlotPair(x,col=col,highlight=highlight,hlCol=hlCol,pch=pch,zoom=zoom,...)
  } else if(x$type=="triple")
  {
    estPlotTriple(x,col=col,highlight=highlight,hlCol=hlCol,pch=pch,zoom=zoom,...)
  }
}