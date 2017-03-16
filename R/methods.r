#' PLotting Methods for DRtrace, doseResponse Objects
#'
#'
#' Plotting methods for \code{\link{doseResponse}} and \code{\link{DRtrace}} classes.
#'
#'
#' Generic methods for dose-response trajectory/trace (\code{\link{DRtrace}}), and dose-response summary  (\code{\link{doseResponse}}) class objects. 

#' The \code{\link{DRtrace}} plotting uses the typical convention of plotting dose-finding experimental trace, with dose levels (x) in the vertical axis and 1/0 responses (y) denoted via filled/empty circles, respectively. In other words, this generic plotting method is only relevant for binary 0/1 outcomes.

#' The \code{\link{doseResponse}} plotting has response rate on the y-axis and dose on the x-axis, and plots symbols whose area is proportional to the weights. 

#' @seealso \code{\link{doseResponse}}, \code{\link{DRtrace}}
#' @param x 	the object, whether DRtrace or doseResponse
#' @param xlab,ylab		x-axis and y-axis labels passed on to \code{\link{plot}}
#' @param pch	the plotting character (doseResponse only), the default being 'X' marks
#' @param varsize 	(doseResponse only) logical, should symbol size vary by sample size? Default \code{TRUE}
#' @param refsize 	(doseResponse only) a reference size by which the plotting sizes will be divided. Larger values make the symbols smaller. Default is \code{mean(dr$weight)}.
#' @param ...	Other arguments passed on to \code{\link{plot}}. 

##' @author Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}	  
#' @example inst/examples/classExamples.r
#' @export
#' @import graphics
#' 
plot.DRtrace<-function(x,xlab="Patient Order",ylab="Dose",...) {

n=dim(x)[1]
plot(x$x,pch=ifelse(x$y==1,16,1),type="b",xaxt="n",yaxt="n",xlab=xlab,ylab=ylab,...)
axis(1,at=1:n)
axis(2,at=sort(unique(x$x)))
}


#############
##' @rdname plot.DRtrace
#' @export
plot.doseResponse<-function(x,xlab="Dose",ylab="Response",pch='X',varsize=TRUE,refsize=mean(x$weight),...) {
cexy=refsize
if(varsize) cexy=sqrt(x$weight/refsize)
plot(y~x,data=x,pch=pch,xlab=xlab,ylab=ylab,cex=cexy,xaxt='n',...)
axis(1,at=sort(unique(x$x)))
}