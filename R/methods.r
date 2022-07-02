#' Plotting Methods for DRtrace, doseResponse Objects
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
#' @param shape the plotting shape (DRtrace only): 'circle' (default), 'square', or 'triangle'
#' @param mcol The color of the main plotting symbols and connecting lines. Default 1 (the current palette's first color). Note: if you change the color and inadvertently use \code{col} instead, there will be an error message.
#' @param connect logical: whether to connect the symbols (generic plotting type 'b'). Default \code{TRUE} for \code{\link{DRtrace}} and \code{FALSE} for \code{\link{doseResponse}}
#' @param varsize 	(\code{doseResponse} only) logical, should symbol size vary by sample size? Default \code{TRUE}
#' @param refsize 	(\code{doseResponse} only) a reference size by which the plotting sizes will be divided. Larger values make the symbols smaller. Default is \code{mean(dr$weight)}.
#' @param dosevals (\code{DRtrace} only) Dose values to be plotted along the y-axis. If \code{NULL} (default), those will be the doses in the dataset (i.e.,\code{sort(unique(x$x))}).
#' @param ...	Other arguments passed on to \code{\link{plot}}. 

##' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	  
#' @example inst/examples/classExamples.r
#' @export
#' @import graphics
#' 
plot.DRtrace <- function(x, xlab="Patient Order", ylab="Dose", shape='circle', connect=TRUE, mcol=1, dosevals=NULL, ...) {

n=dim(x)[1]
# Setting plotting symbol
ch1=16
if(shape=='square') ch1=15
if(shape=='triangle') ch1=17

plot(x$x, pch = ifelse(x$y==1, ch1, ch1-15), type=ifelse(connect,'b','p'), xaxt="n", yaxt="n", xlab=xlab, ylab=ylab, col=mcol, ...)
axis(1,at=1:n, ...)
# Y-axis options
if(is.null(dosevals[1]))  dosevals = sort(unique(x$x))
axis(2, at=dosevals, ...)
}


#############
##' @rdname plot.DRtrace
#' @export
plot.doseResponse<-function(x, xlab="Dose", ylab="Response", pch='X', varsize=TRUE,refsize=1/mean(x$weight), connect=FALSE, mcol=1, ...) {
cexy=refsize
if(varsize) cexy = sqrt(x$weight*refsize)
plot(y~x, data=x, pch=pch, xlab=xlab, ylab=ylab, cex=cexy,
	type=ifelse(connect,'b','p'), col=mcol, ...)
#axis(1,at=x$x,cex.axis=cex.axis)
}