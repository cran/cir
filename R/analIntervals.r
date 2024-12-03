##' Backend utility to calculate analytical CIR/IR interval estimates, given the point estimates
#'
#'
#' For confidence intervals at design points (x values with obesrvations), this function calls \code{intfun} to do the work. In addition, CIs for any x value are calculated using linear interpolation between design points (note that for CIR, this differs from the interpolation of point estimates which is carried out between shrinkage points, as explained in \code{\link{quickIsotone}}). The interval estimation method is presented and discussed by Oron and Flournoy (2017).
#' 
#' @note All provided algorithms and formulae are for binary/Binomial data only. For other data, write your own \code{intfun}, returning a two-column matrix. 
#'
#' @note Interval coverage for extreme percentiles with adaptive designs may be lacking: use \code{adaptiveCurve=TRUE} whenever the \code{target} is outside (0.4, 0.6). This should work as far as the 10th or 90th percentile, but not for more extreme targets. 
#'
#' @seealso \code{\link{quickIsotone}},\code{\link{quickInverse}},\code{\link{morrisCI}},
#' @example inst/examples/fwdCiExamples.r
#'
##' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}
#' @export
#' @references Oron, A.P. and Flournoy, N., 2017. Centered Isotonic Regression: Point and Interval Estimation for Dose-Response Studies. Statistics in Biopharmaceutical Research 3, 258-267.

#' @return a data frame with two variables \code{ciLow, ciHigh} containing the estimated lower and upper confidence bounds, respectively.
#' 
#' @param isotPoint The output of an estimation function such as \code{\link{cirPAVA}}  with the option \code{full=TRUE}. Should be a list of 3 \code{\link{doseResponse}} objects named \code{input, output, shrinkage}.
#' @param outx vector of x values for which estimates will be made. If \code{NULL} (default), this will be set to the set of unique values in isotPoint$x argument (or equivalently in y$x).
#' @param conf numeric, the interval's confidence level as a fraction in (0,1). Default 0.9.
#' @param intfun the function to be used for interval estimation. Default \code{\link{morrisCI}} (see help on that function for additional options).
#' @param ... additional arguments passed on to \code{intfun}

isotInterval<-function(isotPoint,outx=isotPoint$output$x,conf=0.9,intfun=morrisCI,...)
{
## Validation
if(conf<=0 || conf>=1) stop("Confidence must be between 0 and 1.\n")
#if(!is.doseResponse(isotPoint)) stop("Point-estimate data must be in doseResponse format.\n")
if(min(outx)<min(isotPoint$output$x) || max(outx)>max(isotPoint$output$x)) stop("Cannot predict outside data boundaries.\n")

ycount=round(isotPoint$shrinkage$weight*isotPoint$shrinkage$y)
rawInt=intfun(phat=isotPoint$shrinkage$y,n=isotPoint$shrinkage$weight,y=ycount,conf=conf,...)

#if(all(outx %in% isotPoint$x)) return(designInt[match(outx,isotPoint$x),])

if(length(unique(rawInt[,1]))==1 || length(unique(rawInt[,2]))==1)
{ # degenerate case: only one y value
	lcl=rep(NA,length(outx))
	ucl=rep(NA,length(outx))
} else {
	n=isotPoint$shrinkage$weight
	lcl=approx(isotPoint$shrinkage$x[n>0],rawInt[n>0,1],xout=outx,rule=2)$y
	ucl=approx(isotPoint$shrinkage$x[n>0],rawInt[n>0,2],xout=outx,rule=2)$y
}
return(data.frame(ciLow=lcl,ciHigh=ucl))
}

####################### Inverse

#' Backend utility to calculate inverse (dose-finding) intervals, using local inversion and the Delta method
#'
#'
#' Calculate left-bound to right-bound intervals for the dose point estimates, using local slopes at design points (places where observations exist) to invert the forward lower-upper bounds.
#'
#'
#' This function is the "backend engine" for calculating confidence intervals for inverse (dose-finding) estimation. Methodologically this might be the most challenging task in the package. It is expected that most users will not interact with this function directly, but rather indirectly via the convenience wrapper \code{\link{quickInverse}}.
#' 
#' The Delta method in this application boils down to dividing the distance to the forward (vertical) bounds, by the slope, to get the left/right interval width. Both forward intervals and slopes are calculated across a standard set of \eqn{x} values, then interpolated at horizontal cross-sections determined by `target`. Slope estimates are performed by \code{\link{slope}}. 
#' 
#' Starting version 2.3.0, by default the slope estimate is different to the right and left of target. The intervals should now better accommodate the sharp slope changes that often happen with discrete dose-response datasets. Operationally, the intervals are first estimated via the single-slope approach described above. Then using a finer grid of \eqn{x} values, weighted-average slopes to the right and left of the point estimate separately are calculated over the first-stage's half-intervals. The weights are hard-coded as quadratic (Epanechnikov). 
#' 
#' An alternative and much simpler interval method (dubbed "global") is hard-coded into \code{\link{quickInverse}}, and can be chosen from there as an option. It is not recommended being far too conservative, and sometimes not existing. It is now also (since version 2.4.0) used in this function as a fallback upper bound on interval width.
#' 
#' 


#' @return two-column matrix with the left and right bounds, respectively

#' @param isotPoint The output of an estimation function such as \code{\link{cirPAVA},\link{doseFind}},  with the option \code{full=TRUE}. Should be at least a list of 3 \code{\link{doseResponse}} objects named \code{input, output, shrinkage}.
#' @param target A vector of target response rate(s), for which the interval is needed. Default (since version 2.3.0) is the 3 quartiles (`(1:3) / 4`). If changed to \code{NULL}, interval will be returned for the \eqn{y} values of `isotPoint$output`. 
#' @param intfun the function to be used for initial (forward) interval estimation. Default \code{\link{morrisCI}} (see help on that function for additional options).
#' @param conf numeric, the interval's confidence level as a fraction in (0,1). Default 0.9.
#' @param adaptiveCurve logical, should the CIs be expanded by using a parabolic curve between estimation points rather than straight interpolation? Default \code{FALSE}. Switch to `TRUE` recommended when adaptive design was used, and \code{target} is outside of (0.4, 0.6).
#' @param minslope minimum local slope (subsequently normalized by the dose-spacing unit) considered positive, passed on to \code{\link{slope}}. Needed to avoid unrealistically broad intervals. Default 0.01.
#' @param slopeRefinement **(new to 2.3.0)** logical: whether to allow refinement of the slope estimate, including different slopes to the left and right of target. Default `TRUE`. See Details.
#' @param globalCheck **(new to 2.4.0)** logical: whether to allow narrowing of the bound, in case the "global" bound *(obtained via inverting the forward interval, and generally more conservative)* is narrower. Default `TRUE`.
#' @param finegrid a numerical value used to guide how fine the grid of `x` values will be during slope estimation. Should be in (0,1) (preferably much less than 1). Default 0.05.
#' @param ... additional arguments passed on to \code{\link{quickIsotone}}

#' @seealso \code{\link{quickIsotone}},\code{\link{quickInverse}},\code{\link{isotInterval}},
#' \code{\link{slope}}; \code{\link{DRshrink}} for the shrinkage fix.
#' @example inst/examples/invCiExamples.r


#' @export

deltaInverse <- function(isotPoint, target=(1:3)/4, intfun = morrisCI, conf = 0.9,
	adaptiveCurve = FALSE, minslope = 0.01, slopeRefinement = TRUE, finegrid = 0.05, globalCheck = TRUE, ...)
{
k=length(target)
isotPoint$shrinkage$y=round(isotPoint$shrinkage$y,8) ### avoid rounding errors from PAVA
xvals=isotPoint$shrinkage$x
yvals=isotPoint$shrinkage$y
yval0=sort(unique(isotPoint$shrinkage$y))
n=isotPoint$shrinkage$weight
m = length(xvals)

# Normalizing minslope
minslope = minslope / mean(diff(isotPoint$input$x))


## degenerate case, completely flat or otherwise useless
if(sum(n>0)<2 || is.null(yval0) || length(yval0)<=1 || 
    var(yvals)<.Machine$double.eps*1e3) return(cbind(rep(NA,k),rep(NA,k))) 

### Using reference grid to calculate forward and then slope
# New 2.2.2! outx is not only at design/shrinkage points anymore
xgaps = diff(xvals)
gridx = sort(unique(c(xvals, seq(min(xvals), max(xvals), finegrid * diff(range(xvals)) / (m-1) ), max(xvals) ) ) )
xout = gridx
 xout = sort( unique( c(xvals, xvals[-m] + finegrid*xgaps, xvals[-1] - finegrid*xgaps) ) )
# return(xout)

festimate=approx(isotPoint$shrinkage$x,isotPoint$shrinkage$y, xout=xout)$y
cestimate=isotInterval(isotPoint, conf=conf, intfun=intfun, outx=xout, ...)
fslopes=slope(isotPoint$shrinkage$x,isotPoint$shrinkage$y, outx=xout, tol=minslope)

# inverse widths raw; note that lwidths are negative, rwidths positive
rwidths=(festimate-cestimate$ciLow)/fslopes
lwidths=(festimate-cestimate$ciHigh)/fslopes
# return(cbind(xout, lwidths, rwidths))


### New since 2.3.0: slope refinement to right/left slope, and generally wider CIs
if(slopeRefinement)
{
  # Now we actually need the point estimates and the slopes at them
#  testimate=approx(isotPoint$shrinkage$y,isotPoint$shrinkage$x, xout=target)$y
  
#  tslopes=slope(isotPoint$shrinkage$x,isotPoint$shrinkage$y, outx=testimate, tol=minslope)
  
  gridslopes=slope(isotPoint$shrinkage$x,isotPoint$shrinkage$y, outx=gridx, tol=minslope)
#    return(gridslopes)
    
  # "Long coding" this part for clarity?
  newslopes = data.frame(left=rep(NA, length(xout)), right=rep(NA, length(xout)) )
  for(a in seq_along(xout))
  {
    tmp = gridslopes[gridx >= xout[a]+lwidths[a] & gridx <= xout[a] ]
#        return(tmp)
    ntmp = length(tmp)
#  cat(a, ntmp,'\n')
    newslopes$left[a]  = ifelse(ntmp==0, fslopes[a], weighted.mean(tmp, w = (1:ntmp)^2) )
    
    tmp = gridslopes[gridx <= xout[a]+rwidths[a] & gridx >= xout[a] ]
    ntmp = length(tmp)
#  cat(a, ntmp,'\n')
    
    newslopes$right[a] = ifelse(ntmp==0, fslopes[a], weighted.mean(tmp, w = (ntmp:1)^2) )
  }
#  return(newslopes)

  rwidths = rwidths * fslopes / newslopes$right
  lwidths = lwidths * fslopes / newslopes$left

# print(cbind(xout, festimate, cestimate, fslopes, lwidths, rwidths)) ; stop()
  
}


# Adding the widths to the mean curve, self-consistently
rbounds=rev(cummin(rev(xout+rwidths)))
lbounds=cummax(xout+lwidths)


### Calculating at the requested targets
# No target specified: using design points
if (is.null(target)) target = isotPoint$output$y

# Note we use approx() with rule=1, forcing NAs when specified target is outside bounds
if(length(unique(lbounds))==1 || length(unique(rbounds))==1)
{ # degenerate case: only one y value. No interval can be calculated
	nout= length(target)
	return( cbind(rep(NA,nout), rep(NA,nout)) )
}	

if(adaptiveCurve) {
# Otherwise: target was specified
# First, curved case for adaptive design with target!=0.5
	lout=parapolate(unique(festimate),lbounds[!duplicated(festimate)],xout=target, upward=TRUE)
	rout=parapolate(unique(festimate),rbounds[!duplicated(festimate)],xout=target, upward=FALSE)
} else {
	lout=approx(festimate,lbounds,xout=target,rule=1,ties='ordered')$y
	rout=approx(festimate,rbounds,xout=target,rule=1,ties='ordered')$y
}

# New 2.4.0! Optional narrowing via global interval
if(globalCheck)
{
  
  # print(cbind(gridx,cestimate2)); stop() 
  rglob = approx(cestimate$ciLow, y=xout, xout = target, rule=1, ties='ordered')$y
  lglob = approx(cestimate$ciHigh, y=xout, xout = target, rule=1, ties='ordered')$y
#  return(cbind(lout,rout))
  
  #   rfine = approx(x=rbounds, y=xout, xout=festimate2, ties = 'ordered')$y
  #   lfine = approx(x=lbounds, y=xout, xout=festimate2, ties = 'ordered')$y
  #   return(cbind(festimate2, lfine, rfine))
  rout = ifelse(is.finite(rglob) &  rglob < rout, rglob, rout) 
  lout = ifelse(is.finite(lglob) &  lglob > lout, lglob, lout) 
}


dout = cbind(lout,rout)
rownames(dout) = target
return(dout)
}



