#' Inverse (dose-finding) point estimate (e.g., estimating a percentile)
#'
#'
#' Inverse ("dose-finding") point estimation of a dose (x) for a specified target y value (e.g., a response rate), 
#'  using a user-specified forward-estimation algorithm (default is CIR).
#'
#'
#' The function works by calling \code{estfun} for forward estimation of the x-y relationship, then using \code{\link{approx}} with the x and y roles reversed for inverse estimation. It is expected that most users will not interact with this function directly, but rather indirectly via the convenience wrapper \code{\link{quickInverse}}.
#' 
#' The \code{extrapolate} option sets the \code{rule} argument for this second call: 
#
#'  - \code{extrapolate=TRUE} translates to \code{rule=2}, which actually means that the x value on the edge of the estimated y range will be assigned.
#'  - \code{extrapolate=FALSE} (default) translates to \code{rule=1}, which means an \code{NA} will be returned for any target y value lying outside the estimated y range.
#' 
#' Note also that the function is set up to work with a vector of targets.
#'
#' If the data were obtained from an adaptive dose-finding design and you seek to estimate a dose other than the experiment's target, note that away from the target the estimates are likely biased (Flournoy and Oron, 2019). Use \code{adaptiveShrink=TRUE} to mitigate the bias. In addition, either provide the true target as \code{starget}, or a vector of values to \code{target}, with the first value being the true target.
#' 
#' Tie-breaking - the `tiemeth` argument passed on as the `ties` argument for `approx()` - provides yet another complication: as of 2.5.0, the default is `"decide"`, which means that the function chooses the most interior `x` value if `target` falls on the boundary of `y` estimates. Inside the boundaries the argument becomes `mean`, but with CIR this is generally ignored because there are no interior ties. Otherwise, if traditional isotonic regression (\code{\link{oldPAVA}}) is used, then the `"decide"` algorithm will pass `ties = "ordered"` on to `approx()`, respecting IR's flat stretches. A user-chosen value for `tiemeth` will override all of that; see `?approx` for options.


 
##' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}

#' @param y  can be either of the following: y values (response rates), a 2-column matrix with positive/negative response counts by dose, a \code{\link{DRtrace}} object or a \code{\link{doseResponse}} object. 
#' @param x dose levels (if not included in y). 
#' @param wt weights (if not included in y).
#' @param target A vector of target response rate(s), for which the percentile dose estimate is needed. See Note.
#' @param full logical, is a more complete output desired (relevant only for doseFind)? if \code{FALSE} (default), only a point estimate of the dose (x) for the provided target rate is returned.
#' @param extrapolate logical: should extrapolation beyond the range of estimated y values be allowed? Default \code{FALSE}.
#' @param dec (relevant only for doseFind) logical, is the true function is assumed to be monotone decreasing? Default \code{FALSE}.
#' @param estfun the name of the dose-response estimation function. Default \code{\link{cirPAVA}}.
#' @param errOnFlat logical: in case the forward estimate is completely flat making dose-finding infeasible, should an error be returned? Under default (\code{FALSE}), \code{NA}s are returned for the target estimate. 
#' @param adaptiveShrink logical, should the y-values be pre-shrunk towards an experiment's target? Recommended if data were obtained via an adaptive dose-finding design. See \code{\link{DRshrink}} and the Note.
#' @param starget The shrinkage target. Defaults to \code{target[1]}.
#' @param tiemeth The method to resolve ties. Default `"decide"`, meaning the function chooses based on context. See Details.
#' @param ...	Other arguments passed on to \code{\link{doseResponse}} and \code{estfun}.

#' @return under default, returns point estimate(s) of the dose (x) for the provided target rate(s). With \code{full=TRUE}, returns a list with
#'
#'  - `targest`: The said point estimate of x
#'  - `input`:    a \code{doseResponse} object summarizing the input data
#'  - `output`:  a \code{doseResponse} object with the forward estimate at design points
#'  - `shrinkage`:  a \code{doseResponse} object which is the \code{alg} output of the forward-estimation function
#' 

#' @seealso \code{\link{oldPAVA}},\code{\link{cirPAVA}}. If you'd like point and interval estimates together, use \code{\link{quickInverse}}.

#' @references Flournoy N and Oron AP, 2020. Bias Induced by Adaptive Dose-Finding Designs. Journal of Applied Statistics 47, 2431-2442.


#' @export

doseFind<-function(y,x=NULL,wt=NULL,estfun=cirPAVA,target=NULL,full=FALSE, dec=FALSE, extrapolate=FALSE,errOnFlat=FALSE, adaptiveShrink=FALSE, starget=target[1], tiemeth = "decide", ...) {


### converting to doseResponse object 
### Basically it's a numeric data frame with x,y,weight, and x increasing
if(is.null(target)) stop("Must provide target to dose-find for.\n")
if(!is.doseResponse(y)) dr=doseResponse(y=y,x=x,wt=wt,...) else dr=y
if (any(is.na(dr))) stop ("Missing values are not allowed.\n")
# Adaptive-design shrinkage fix
if(adaptiveShrink) dr=DRshrink(y=dr,target=starget,...)


# We start via forward estimation

pavout<-estfun(y=dr,full=TRUE,dec=dec,target=target,...)

newx<-pavout$shrinkage$x
newy<-pavout$shrinkage$y
# Degenerate case
if(min(newy)==max(newy))
{
	if(errOnFlat) stop("Cannot dose-find; curve completely flat.\n")
	if(full) return (list(targest=NA,input=dr,shrinkage=pavout$shrinkage,output=pavout$output))
	return(rep(NA,length(target)))
}
### New for 2.5.0! Context-sensitive ties method
if(tiemeth == 'decide')
{
  tiemeth = mean # safe-ish default
  # Old PAVA override
  if(identical(estfun, oldPAVA)) tiemeth = "ordered"
  # boundary behavior under both methods
  if(max(newy) == target) tiemeth = min
  if(min(newy) == target) tiemeth = max

}


### The estimate is generated by using 'approx' with x and y interchanged
#return(pavout)
tout=approx(x=newy, y=newx, xout=target, ties=tiemeth, rule=1)$y
# targets outside y range are provisionally NA'ed... now if so desired, extrapolate
if(any(is.na(tout)) && extrapolate)
{
	lowTargets=which(target<min(newy))
	if(length(lowTargets)>0) tout[lowTargets]=extrapol(c(newy[1],newx[1]),c(newy[2],newx[2]),xout=target[lowTargets])
	hiTargets=which(target>max(newy))
	m=length(newx)
	if(length(hiTargets)>0) tout[hiTargets]=extrapol(c(newy[m-1],newx[m-1]),c(newy[m],newx[m]),xout=target[hiTargets])
	tout[!is.finite(tout)]=NA
}
	
#### output

if (!full)  return(tout) 

return (list(targest=tout, input=dr, shrinkage=pavout$shrinkage, output=pavout$output))
}

#' Convenient point and Interval Inverse Estimation ("Dose-Finding"), using CIR or IR
#'
#'
#' Convenience wrapper for point and interval estimation of the "dose" that would generate a \code{target} "response" value, using CIR and IR.
#' 
#' 
#' The inverse point estimate is calculated in a straightforward manner from a forward estimate, using \code{\link{doseFind}}. For the inverse interval, the default option (\code{delta=TRUE}) calls \code{\link{deltaInverse}} for a "local" (Delta) inversion of the forward intervals. 

#' If \code{delta=FALSE}, a second call to \code{\link{quickIsotone}} generates a high-resolution grid outlining the forward intervals. Then the algorithm "draws a horizontal line" at \code{y=target} to find the right and left bounds on this grid. Note that the right (upper) dose-finding confidence bound is found on the lower forward confidence bound, and vice versa. This approach is not recommended, tending to produce CIs that are too wide.
#'
#' If the data were obtained from an adaptive dose-finding design and you seek to estimate a dose other than the experiment's target, note that away from the target the estimates are likely biased (Flournoy and Oron, 2019). Use \code{adaptiveShrink=TRUE} to mitigate the bias. In addition, either provide the true target as \code{starget}, or a vector of values to \code{target}, with the first value being the true target.


#' @param y  can be either of the following: y values (response rates), a 2-column matrix with positive/negative response counts by dose, a \code{\link{DRtrace}} object or a \code{\link{doseResponse}} object. 
#' @param x dose levels (if not included in y). 
#' @param wt weights (if not included in y).
#' @param target A vector of target response rate(s), for which the percentile dose estimate is needed. See Note.
#' @param estfun the function to be used for point estimation. Default \code{\link{cirPAVA}}.
#' @param intfun the function to be used for interval estimation. Default \code{\link{morrisCI}} (see help on that function for additional options).
#' @param delta logical: should intervals be calculated using the delta ("local") method (default, \code{TRUE}) or back-drawn directly from the forward bounds? See Details.
#' @param conf numeric, the interval's confidence level as a fraction in (0,1). Default 0.9.
#' @param resolution numeric: how fine should the grid for the inverse-interval approximation be? Default 100, which seems to be quite enough. See 'Details'.
#' @param extrapolate logical: should extrapolation beyond the range of estimated y values be allowed? Default \code{FALSE}. Note this affects only the point estimate; interval boundaries are not extrapolated.
#' @param adaptiveShrink logical, should the y-values be pre-shrunk towards an experiment's target? Recommended when the data were obtained via an adaptive dose-finding design. See \code{\link{DRshrink}} and the Note below.
#' @param starget The shrinkage target. Defaults to \code{target[1]}.
#' @param adaptiveCurve logical, should the CIs be expanded by using a parabolic curve between estimation points rather than straight interpolation (default \code{FALSE})? Recommended when adaptive design was used and \code{target} is not 0.5.
#' @param ...	Other arguments passed on to \code{\link{doseFind}} and \code{\link{quickIsotone}}, and onwards from there.

#' @return A data frame with 4 elements:
#' 
#'  - `target`: The user-provided target values of y, at which x is estimated
#'  - `point`: The point estimates of x
#'  - `lowerPPconf,upperPPconf`: the interval-boundary estimates for a 'PP'=\code{100*conf} confidence interval
#' 

#' @references Flournoy N and Oron AP, 2020. Bias Induced by Adaptive Dose-Finding Designs. Journal of Applied Statistics 47, 2431-2442.


#' @seealso \code{\link{quickIsotone}},\code{\link{doseFind}},\code{\link{deltaInverse}}
#' @example inst/examples/invExamples.r
#' @export

quickInverse<-function(y, x=NULL, wt=NULL, target, estfun=cirPAVA, intfun = morrisCI,
	delta=TRUE,conf = 0.9, resolution=100, extrapolate=FALSE,
	adaptiveShrink=FALSE, starget=target[1], adaptiveCurve = FALSE, ...)
{

#### Point estimate first
dr=doseResponse(y,x,wt)
# Adaptive-design shrinkage fix. We must do it here b/c used for both point and interval
if(adaptiveShrink) dr=DRshrink(y=dr,target=starget,...)
m=length(dr$x)
## adaptiveShrink set to FALSE to avoid double-shrinking
pestimate=doseFind(y=dr,estfun=estfun,target=target,full=TRUE,extrapolate=extrapolate,adaptiveShrink=FALSE,...)
dout=data.frame(target=target,point=pestimate$targest,low=-Inf,high=Inf)

# All point estimates are NA -> exit here
if(all(is.na(pestimate$targest))) {
	names(dout)[3:4]=paste(c("lower","upper"),round(100*conf),"conf",sep="")
	return(dout)
}

if(delta) { ## Default, delta-method ("local") intervals
	dout[,3:4]=deltaInverse(pestimate,target=target,intfun=intfun,conf=conf, adaptiveCurve = adaptiveCurve,...)
} else {
#### CI, using "global" interpolation; generally too conservative and not recommended
#	if(adaptiveShrink) dr=DRshrink(y=dr,target=starget,...)

	## Calculate forward CIs at high-rez grid
	higrid=seq(dr$x[1],dr$x[m],length.out=1+resolution*(m-1))
	fwdCIgrid=quickIsotone(dr,outx=higrid,conf=conf,intfun=intfun,...)
	fwdCIdesign=fwdCIgrid[match(dr$x,fwdCIgrid$x),]

	### CI by finding 'right points' on grid
	for (a in 1:length(dout$target))
	{	
		if(dout$target[a]<=max(fwdCIgrid[,3]))  # upper inverse interval taken from fwd LCL 
		{
			dout$high[a]=fwdCIgrid[min(which(fwdCIgrid[,3]>=dout$target[a])),1]
		} 
		if(dout$target[a]>=min(fwdCIgrid[,4]))  # and vice versa
		{
			dout$low[a]=fwdCIgrid[max(which(fwdCIgrid[,4]<=dout$target[a])),1]
		} 
	}
}
names(dout)[3:4]=paste(c("lower","upper"),round(100*conf),"conf",sep="")
return(dout)
}






