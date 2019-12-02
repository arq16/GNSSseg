#' Graph of the result
#'
#' plot the signal with the estimated average
#' @param Data a data frame, with size [n x 2], containing the signal (e.g. the daily GPS-ERAI series for GNSS) and the dates (in format yyyy-mm-dd of type "calendar time" (class POSIXct))
#' @param segmentation the estimated segmentation (result of the GNSSseg function)
#' @param functional the estimated functional (result of the GNSSseg function)
#'
#' @return a plot of the results with the signal
#'
#' @details
#' The function gives the plot of the results with the signal
#'
#' @examples
#' data(Data)
#' lyear=100
#' Kmax=8
#' lmin=1
#' result=GNSSseg(Data,lyear,Kmax=Kmax)
#' plot_GNSS(Data,result$segmentation,result$functional)
#' @export

plot_GNSS=function(Data,segmentation,functional){
  graphics::plot(Data$date,Data$signal,cex=0.7,type="l",col="#009999",xlab="date",ylab="signal")
  mean.est.t  = rep(segmentation$mean,diff(c(0,segmentation$end)))
  if (functional==FALSE) {functional=0}
  average=mean.est.t+functional
  graphics::lines(Data$date,average,col="red")
  graphics::abline(v=Data$date[segmentation$end],col="red")
}
