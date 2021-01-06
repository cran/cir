##' Returns analytical interval estimates, given isotonic-regression (centered or not) point estimates
#'
#'
#' For confidence intervals at design points ($x$ values with obesrvations), this function calls \code{intfun} to do the work. In addition, CIs for any $x$ value are calculated using linear interpolation between design points (note that for CIR, this differs from the interpolation of point estimates which is carried out between shrinkage points, as explained in \code{\link{quickIsotone}})
#' 
#' @note All provided algorithm and formulae are for Binomial data only. For other data, write your own \code{intfun}, returning a two-column matrix. The interval estimation method is presented and discussed by Oron and Flournoy (2017).
#'
#' @seealso \code{\link{quickIsotone}},\code{\link{quickInverse}},\code{\link{morrisCI}},
#' @example inst/examples/fwdCiExamples.r
#'
##' @author Assaf P. Oron \code{<aoron.at.idmod.org>}
#' @export
#' @references Oron, A.P. and Flournoy, N., 2017. Centered Isotonic Regression: Point and Interval Estimation for Dose-Response Studies. Statistics in Biopharmaceutical Research, In Press (author's public version available on arxiv.org).

#' @return a data frame with two variables \code{ciLow, ciHigh} containing the estimated lower and upper confidence bounds, respectively.
#' 
#' @param isotPoint a \code{\link{doseResponse}} object with the x values and isotonic point estimates at design points.
#' @param outx vector of x values for which estimates will be made. If \code{NULL} (default), this will be set to the set of unique values in isotPoint$x argument (or equivalently in y$x).
#' @param conf numeric, the interval's confidence level as a fraction in (0,1). Default 0.9.
#' @param intfun the function to be used for interval estimation. Default \code{\link{morrisCI}} (see help on that function for additional options).
#' @param parabola logical: should the confidence-interval's interpolation between points with observations follow a parabola (\code{TRUE}) creating broader intervals between observations, or a straight line (\code{FALSE}, default)?
#' @param ... additional arguments passed on to \code{intfun}

isotInterval<-function(isotPoint,outx=isotPoint$x,conf=0.9,intfun=morrisCI,parabola=FALSE,...)
{
## Validation
if(conf<=0 || conf>=1) stop("Confidence must be between 0 and 1.\n")
if(!is.doseResponse(isotPoint)) stop("Point-estimate data must be in doseResponse format.\n")
if(min(outx)<min(isotPoint$x) || max(outx)>max(isotPoint$x)) stop("Cannot predict outside data boundaries.\n")

ycount=round(isotPoint$weight*isotPoint$y)
designInt=intfun(phat=isotPoint$y,n=isotPoint$weight,y=ycount,conf=conf,...)

#if(all(outx %in% isotPoint$x)) return(designInt[match(outx,isotPoint$x),])

if(parabola) lcl=parapolate(isotPoint$x,designInt[,1],xout=outx,upward=TRUE) else
	lcl=approx(isotPoint$x,designInt[,1],xout=outx)$y
if(parabola) ucl=parapolate(isotPoint$x,designInt[,2],xout=outx,upward=FALSE) else
	ucl=approx(isotPoint$x,designInt[,2],xout=outx)$y

return(data.frame(ciLow=lcl,ciHigh=ucl))
}

####################### Inverse

#' Calculate inverse (dose-finding) intervals, using local inversion and the Delta method
#'
#'
#' Calculate left-bound to right-bound intervals for the dose point estimates, using local slopes at design points (places where observations exist) to invert the forward lower-upper bounds.
#'
#'
#' The Delta method in this application boils down to dividing the distance to the forward (vertical) bounds, by the slope, to get the left/right interval width. Slope estimates are performed by \code{\link{slope}}. An alternative method (dubbed "global") is hard-coded into \code{\link{quickInverse}}. 

#' @return two-column matrix with the left and right bounds, respectively

#' @param y  can be either of the following: y values (response rates), a 2-column matrix with positive/negative response counts by dose, a \code{\link{DRtrace}} object or a \code{\link{doseResponse}} object. 
#' @param x dose levels (if not included in y). 
#' @param wt weights (if not included in y).
#' @param target A vector of target response rate(s), for which the interval is needed. If \code{NULL} (default), interval will be returned for the point estimates at design points (e.g., if the forward point estimate at $x_1$ is 0.2, then the first returned interval is for the 20th percentile).
#' @param estfun the function to be used for point estimation. Default \code{\link{cirPAVA}}.
#' @param intfun the function to be used for initial (forward) interval estimation. Default \code{\link{morrisCI}} (see help on that function for additional options).
#' @param conf numeric, the interval's confidence level as a fraction in (0,1). Default 0.9.
#' @param parabola logical: should the confidence-interval's interpolation between points with observations follow a parabola (\code{TRUE}) creating broader intervals between observations, or a straight line (\code{FALSE}, default)?
#' @param adaptiveShrink logical, should the y-values be pre-shrunk towards an experiment's target? Recommended if data were obtained via an adaptive dose-finding design. See \code{\link{DRshrink}}.
#' @param starget The shrinkage target. Defaults to \code{target[1]}.
#' @param ... additional arguments passed on to \code{\link{quickIsotone}}

#' @seealso \code{\link{quickIsotone}},\code{\link{quickInverse}},\code{\link{isotInterval}},
#' \code{\link{slope}}; \code{\link{DRshrink}} for the shrinkage fix.
#' @example inst/examples/invCiExamples.r


#' @export

deltaInverse<-function(y,x=NULL,wt=NULL,target=NULL,estfun=cirPAVA, intfun = morrisCI, conf = 0.9,adaptiveShrink=FALSE,starget=target[1],parabola=FALSE,...)
{
dr=doseResponse(y,x,wt)
# Optional pre-shrinking of y for adaptive designs
if(adaptiveShrink) dr=DRshrink(y=dr,target=starget,...)

k=length(target)
# We start by constructing inverse intervals based on design-point estimates
forward=quickIsotone(dr,outx=NULL,conf=conf,intfun=intfun,estfun=estfun,...)
forward$y=round(forward$y,10) ### avoid rounding errors from PAVA
yvals=sort(unique(forward$y))
#cat(yvals)
if(length(yvals)==1 || var(yvals)<.Machine$double.eps*1e3) return(cbind(rep(NA,k),rep(NA,k))) ## degenerate case, completely flat
fslopes=slope(dr$x,forward$y)

# inverse widths raw
rwidths=(forward$y-forward$lower)/fslopes
lwidths=(forward$y-forward$upper)/fslopes
# Adding the widths to the mean curve, self-consistently
rbounds=rev(cummin(rev(tapply(dr$x+rwidths,forward$y,max))))
lbounds=cummax(tapply(dr$x+lwidths,forward$y,min))
designCIs=cbind(lbounds,rbounds)[match(forward$y,yvals),]

if (is.null(target)) return(designCIs)

if(parabola) {  ### direction chosen by trial & error...
	good=(target>=min(yvals) & target<=max(yvals))
	lout=rep(NA,k)
	rout=lout
	lout[good]=parapolate(yvals,lbounds,xout=target[good],upward=TRUE)
	rout[good]=parapolate(yvals,rbounds,xout=target[good],upward=FALSE)
} else {
	lout=approx(yvals,lbounds,xout=target,rule=1)$y
	rout=approx(yvals,rbounds,xout=target,rule=1)$y
}
return(cbind(lout,rout))
}



