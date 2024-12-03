#' Plotting Methods for DRtrace, doseResponse Objects
#'
#'
#' Plotting methods for \code{\link{doseResponse}} and \code{\link{DRtrace}} classes.
#'
#'
#' Generic methods for dose-response trajectory/trace (\code{\link{DRtrace}}), and dose-response summary  (\code{\link{doseResponse}}) class objects. 

#' The \code{\link{DRtrace}} plotting uses the typical convention of plotting dose-finding experimental trace, with dose levels (x) in the vertical axis and 1/0 responses (y) denoted via filled/empty circles, respectively. In other words, this generic plotting method is only relevant for binary 0/1 outcomes. If cohort information is provided via `x$cohort` (i.e., multiple observations considered as collected together rather than each data point sequentially), then the plotting will respect cohort structure.

#' The \code{\link{doseResponse}} plotting has response rate on the y-axis and dose on the x-axis, and plots symbols whose area is proportional to the weights. 

#' @seealso \code{\link{doseResponse}}, \code{\link{DRtrace}}
#' @param x 	the object, whether DRtrace or doseResponse
#' @param xlab,ylab		x-axis and y-axis labels passed on to \code{\link{plot}}
#' @param pch	the plotting character (doseResponse only), the default being 'X' marks
#' @param shape the plotting shape (DRtrace only): 'circle' (default), 'square', or 'triangle'
#' @param mcol The color of the main plotting symbols and connecting lines. Default 1 (the current palette's first color). Note: if you change the color and inadvertently use \code{col} instead, there will be an error message.
#' @param connect logical: whether to connect the symbols (generic plotting type 'b'). Default \code{TRUE} for \code{\link{DRtrace}} and \code{FALSE} for \code{\link{doseResponse}}
#' @param varsize 	(\code{doseResponse} only) logical, should symbol size vary by sample size? Default \code{TRUE}
#' @param refsize 	(\code{doseResponse} only) a reference size by which the plotting sizes will be multiplied. Default is \code{1/sqrt(mean(dr$weight))}, scaled so that if `varsize = TRUE` the weighted-average symbol size is 1. If `varsize = FALSE`, this argument is equivalent to `cex` in an ordinary x-y `plot()` call.
#' @param dosevals Dose values to be plotted along the x-axis (`plot.doseResponse`) or y-axis (`plot.DRtrace`) . If \code{NULL} (default), those will be the doses in the dataset (i.e.,\code{sort(unique(x$x))}).
#' @param offset (\code{DRtrace} only) In case of a cohort-based experiment, the relative vertical offset between symbols for outcomes within the same cohort (as fraction of dose spacing). Default 0.2.
#' @param ...	Other arguments passed on to \code{\link{plot}}. 
#' 
#' Conversely, putting values on a different scale into `dosevals`, or even text labels instead of numbers, won't work. For the former, change the scale at the source data (i.e., in the plotted object). For the latter, sorry but no solution at present. 

#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	  
#' @example inst/examples/classExamples.r
#' @export
#' @import graphics
#' 
plot.DRtrace <- function(x, xlab="Patient Order", ylab="Dose", shape='circle', connect=TRUE, 
                         mcol=1, dosevals=NULL, offset=0.2, ...) {

n=dim(x)[1]
# Setting plotting symbol
ch1=16
if(shape=='square') ch1=15
if(shape=='triangle') ch1=17

if(is.null(dosevals[1]))  dosevals = sort(unique(x$x))

if(length(unique(x$cohort)) < n) # Case with cohorts
{
  xlab = "Cohort"
  spacings = table( abs(diff(x$x)[diff(x$x) != 0]) )
  spacing = as.numeric(names(spacings))[which.max(spacings)]
  if(length(spacing)==0) spacing = 1

  x$within = unlist(tapply(x$cohort, x$cohort, seq_along))
  x$midpoint = rep(tapply(x$within, x$cohort, max), tapply(x$within, x$cohort, max)) / 2 + 0.5
  plot(x$cohort, x$x + spacing*offset*(x$within - x$midpoint), 
       pch = ifelse(x$y==1, ch1, ch1-15),
       xaxt='n', yaxt='n', xlab=xlab, ylab=ylab, col=mcol, ...)
  if(connect) points(unique(x$cohort), tapply(x$x, x$cohort, function(x) x[1]), 
                     type = 'b', cex=0, ...)
  axis(1, at = 1:max(x$cohort), ...)
    
}  else {
  
  plot(x$x, pch = ifelse(x$y==1, ch1, ch1-15), type=ifelse(connect,'b','p'), 
     xaxt='n', yaxt='n', xlab=xlab, ylab=ylab, col=mcol, ...)
#ylim = range(dosevals),
  axis(1, at = 1:n, ...)
}

axis(2, at=dosevals, ...)
}


#############
#' @rdname plot.DRtrace
#' @export
plot.doseResponse<-function(x, xlab="Dose", ylab="Response", pch='X', varsize=TRUE,
                            refsize=sqrt(1/mean(x$weight)), connect=FALSE, mcol=1, 
                            dosevals=NULL, ...) 
{
cexy = refsize
if(varsize) cexy = refsize * sqrt(x$weight)

if(is.null(dosevals[1]))  dosevals = sort(unique(x$x))

#if(is.null(xlim)) xlim = range(x$x)

plot(y~x, data=x, pch=pch, xlab=xlab, ylab=ylab, cex=cexy, xaxt="n",
	type=ifelse(connect, 'b', 'p'), col=mcol, ...)

axis(1, at=dosevals, ...)
}