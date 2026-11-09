#' Plotting Methods for DRtrace, doseResponse Objects
#'
#'
#' Plotting methods for \code{\link{doseResponse}} and \code{\link{DRtrace}} classes.
#'
#'
#' Generic methods for dose-response trajectory/trace (\code{\link{DRtrace}}), and dose-response summary  (\code{\link{doseResponse}}) class objects. 

#' The \code{\link{DRtrace}} plotting uses the typical convention of plotting dose-finding experimental trace, with dose levels (x) in the vertical axis and 1/0 responses (y) denoted via filled/empty symbols, respectively. In other words, this generic plotting method is only relevant for binary 0/1 outcomes. If cohort information is provided via `x$cohort` (i.e., multiple observations considered as collected together rather than each data point sequentially), then the plotting will respect cohort structure.
#'
#' The \code{\link{doseResponse}} plotting has response rate on the y-axis and dose on the x-axis, and plots symbols whose area is proportional to the weights. 

#' @seealso \code{\link{doseResponse}}, \code{\link{DRtrace}}
#' @param x 	the object, whether DRtrace or doseResponse
#' @param xlab,ylab		x-axis and y-axis labels passed on to \code{\link{plot}}
#' @param pch	the plotting character (doseResponse only), the default being 'X' marks
#' @param shape the plotting shape (`DRtrace` only): 'circle' (default), 'square', or 'triangle'
#' @param connect logical: whether to connect the symbols (generic plotting type 'b'). Default \code{TRUE} for \code{\link{DRtrace}} and \code{FALSE} for \code{\link{doseResponse}}
#' @param scol The color of the main plotting symbols. Default 1 (the current palette's first color). Note: if you inadvertently use \code{col} to set the color, there will be an error message.
#' @param lcol The color of the connecting lines, if relevant. Default same as `scol`.
#' @param size (\code{DRtrace} only) symbol size, analogous to `cex` (which you should not use, see note for `scol` above). A different name is used because there is some nuance related to enabling `scol` and `lcol` to differ, as well as accommodating larger symbols. Don't ask *(but you can look at the code of course)*.
#' @param varsize 	(\code{doseResponse} only) logical, should symbol size vary by sample size? Default \code{TRUE}
#' @param refsize 	(\code{doseResponse} only) a reference size by which the plotting sizes will be multiplied. Default is \code{1/sqrt(mean(dr$weight))}, scaled so that if `varsize = TRUE` the weighted-average symbol size is 1. If `varsize = FALSE`, this argument is equivalent to `cex` in an ordinary x-y `plot()` call.
#' @param obsvals (`DRtrace` only) values to be plotted along the x-axis denoting observation order. If \code{NULL} (default), the graphical environment will choose which values to plot automatically, per default R behavior.
#' @param dosevals Dose values to be plotted along the x-axis (`plot.doseResponse`) or y-axis (`plot.DRtrace`). Must be on the same scale as `x$x`. If \code{NULL} (default), those will be the doses in the dataset (i.e.,\code{sort(unique(x$x))}).
#' @param offset (\code{DRtrace} only) In case of a cohort-based experiment, the relative vertical offset between symbols for outcomes within the same cohort (as fraction of dose spacing). Default 0.2.
#' @param ylim (`DRtrace` only) Similar to the base `plot()` function's `ylim` argument. If kept `NULL` (default), will be determined by `range(dosevals)`.
#' @param xlim (`doseResponse` only) Analogous to the above. 
#' @param ...	Other arguments passed on to \code{\link{plot}}. 
#' 

#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	  
#' @example inst/examples/classExamples.r
#' @export
#' @import graphics
#' 
plot.DRtrace <- function(x, xlab="Patient Order", ylab="Dose", shape = 'circle', connect = TRUE,
                         scol=1, lcol = scol, obsvals = NULL, dosevals = NULL, offset = 0.1, 
                         ylim = NULL, size = 1, ...) 
{

n=dim(x)[1]
# Setting plotting symbols
symb1=16
if(shape=='square') symb1=15
if(shape=='triangle') symb1=17

if(is.null(dosevals[1]))  dosevals = sort(unique(x$x))
xaxis = ifelse(is.null(obsvals[1]), 's', 'n')

if(length(unique(x$cohort)) < n) # Case with cohorts
{
  xlab = "Cohort"
  spacings = table( abs(diff(x$x)[diff(x$x) != 0]) )
  spacing = as.numeric(names(spacings))[which.max(spacings)]
  if(length(spacing)==0) spacing = 1

  x$within = unlist(tapply(x$cohort, x$cohort, seq_along))
  x$midpoint = rep(tapply(x$within, x$cohort, max), tapply(x$within, x$cohort, max)) / 2 + 0.5
  
  if(is.null(ylim)) ylim = c(min(dosevals) - offset * max(x$within)/2, 
                             max(dosevals) + offset * max(x$within)/2)
  plot(x$cohort, x$x + spacing*offset*(x$within - x$midpoint), 
       pch = ifelse(x$y==1, symb1, symb1-15), cex = size,
       xaxt=xaxis, yaxt='n', xlab=xlab, ylab=ylab, col=scol, ylim = ylim, ...)
  if(connect) points(unique(x$cohort), tapply(x$x, x$cohort, function(x) x[1]), 
                     type = 'b', pch = NA, col = lcol, ...)
  if(xaxis == 'n') axis(1, at = obsvals, xaxt = 's', ...)
    
}  else {

  if(is.null(ylim)) ylim = range(dosevals)
  plot(x$x, pch = ifelse(x$y==1, symb1, symb1-15), xaxt=xaxis, yaxt='n', 
       xlab=xlab, ylab=ylab, col=scol, ylim = ylim, cex = size, ...)
  if(connect) {
    points (x$x, type = 'b', cex = 0, col = lcol, ...)
  # Now need to erase the connectors butting into the hollow symbols...
    points ((1:n)[x$y==0], x$x[x$y==0], pch = symb1, cex = size *0.8 , col = par()$bg) 
  }
  if(xaxis == 'n') axis(1, at = obsvals) # , xaxt = 's')
}

axis(2, at=dosevals)
}


#############
#' @rdname plot.DRtrace
#' @export
plot.doseResponse<-function(x, xlab = "Dose", ylab = "Response", pch = 'X', varsize = TRUE,
                            refsize = NULL, connect = FALSE, scol = 1, 
                            dosevals = NULL, xlim = NULL, ...) 
{
if(!is.doseResponse(x)) stop('x must be a doseResponse object.\n')
if(is.null(refsize)) refsize = sqrt(1/mean(x$weight) )
cexy = refsize
if(varsize) cexy = refsize * sqrt(x$weight)

if(is.null(dosevals[1]))  dosevals = sort(unique(x$x))
if(is.null(xlim)) xlim = range(dosevals)

plot(y~x, data=x, pch=pch, xlab=xlab, ylab=ylab, cex=cexy, xaxt="n",
	type=ifelse(connect, 'b', 'p'), col=scol, xlim = xlim, ...)

axis(1, at=dosevals, ...)
}
