#### Some standard CI calculation functions
#' Standard unordered-Binomial confidence interval utilities.
#' 
#' Standard small-sample Binomial confidence interval utilities, using the methods of Wilson, Agresti-Coull and Jeffrys.
#' 
#' These functions implement the basic (uncorrected) three intervals which are seen by the consensus of literature as the "safest" off-the-shelf formulae. None of them account for ordering or monotonicity; therefore the \code{cir} package default is \code{\link{morrisCI}} which does account for that, with the 3 unordered formulae used for optional narrowing of the interval at individual points.

#' @aliases agcouCI jeffCI
#' @seealso \code{\link{isotInterval}} for more details about how forward CIs are calculated, \code{\link{quickInverse}} for inverse (dose-finding) intervals.
#' @example inst/examples/fwdCiExamples.r
#' @export

#' @return A two-column matrix with the same number of rows as \code{length(phat)}, containing the calculated lower and upper bounds, respectively.
#' 
#' @param phat numeric vector, point estimates for which an interval is sought
#' @param n integer vector of same length, of pointwise sample sizes
#' @param conf numeric in (0,1), the confidence level
#' @param w1,w2 numeric, weights used in \code{jeffCI} only
#' @param ... pass-through for compatibility with a variety of calling functions
#' 
wilsonCI<-function(phat,n,conf=0.9,...)
{
zalpha=qnorm(1-(1-conf)/2)
wid=zalpha*sqrt((phat*(1-phat)+zalpha^2/(4*n))/n)
return(cbind(pmax(0,(phat+zalpha^2/(2*n)-wid)/(1+zalpha^2/n)),
	pmin(1,(phat+zalpha^2/(2*n)+wid)/(1+zalpha^2/n))))
}

#' @rdname wilsonCI
#' @export
agcouCI<-function(phat,n,conf=0.9,...)
{
zalpha=qnorm(1-(1-conf)/2)
ptilde=(phat+zalpha^2/(2*n))/(1+zalpha^2/n)
ntilde=n+zalpha^2
wid=zalpha*sqrt(ptilde*(1-ptilde)/ntilde)
return(cbind(pmax(0,ptilde-wid),
	pmin(1,ptilde+wid)))
}

#' @rdname wilsonCI
#' @export
jeffCI<-function(phat,n,conf=0.9,w1=0.5,w2=w1,...)
{
x=n*phat
clow=ifelse(phat==0,0,qbeta((1-conf)/2,x+w1,n-x+w2))
chigh=ifelse(phat==1,1,qbeta(1-(1-conf)/2,x+w1,n-x+w2))

return(cbind(clow,chigh))
}


########### Interpolation, Extrapolation, Parapolation....

# Linearly interpolate/extrapolate from a single segment.

extrapol<-function(point1,point2,xout)
{
slope=(point2[2]-point1[2])/(point2[1]-point1[1])
return(point2[2]+(xout-point2[1])*slope)
}

# Locally monotone quadratic interpolation between points on a plane

parapolate<-function(x,y,xout,upward,full=FALSE)
{
## validations
if(any(diff(x)<=0)) stop("x must be monotone strictly increasing.\n")
m=length(x)
if(m<2) stop ("Need at least 2 points.\n")
#if(min(xout)<x[1] || max(xout)>x[m]) stop("This function does not extrapolate.\n")
if(length(y)!=m) stop("mismatched x,y lengths.\n")

yout=y[match(xout,x)]
if(!any(is.na(yout))) return(yout)  ### no interpolation to do

## These, to help make vectorized and interpretable formulae
x1=x[-m]
y1=y[-m]
x2=x[-1]
y2=y[-1]
denom=diff(x)^2
# A bit of boolean finesse is needed:
if(length(upward)==1) upward=rep(upward,m-1)
oneflat=xor(upward,diff(y)<0)
# Explanation: 'oneflat' determines whether parabola is 'flat' at point 1 or 2 of each segment

## parabola coefficients
a=(y1-y2)/denom
a[oneflat]=-a[oneflat]

b=ifelse(oneflat,-2*x1*a,-2*x2*a)
cee=y1-a*x1^2-b*x1  # per parabola definitions :)

intchoose=findInterval(xout,x)
intchoose[intchoose==0]=1
candy=a[intchoose]*xout^2+b[intchoose]*xout+cee[intchoose]
yout[!(xout %in% x)]=candy[!(xout %in% x)]
# Forcing NA outside boundaries
yout[ xout<x[1] | xout>x[m] ] = NA
if(!full) return(yout)
return(list(a=a,b=b,c=cee,outdat=data.frame(x=xout,y=yout)))
}

#####################
#' Piecewise-linear local slopes given a (non-strictly) monotone x-y sequence
#'
#'
#' Estimate monotone piecewise-linear slopes, with the default behavior forbidding zero slope. This behavior is due to the fact that the function is used to invert confidence intervals using the Delta method. The input interval has to be strictly increasing in \code{x}, and (non-strictly) monotone in \code{y} (increasing or decreasing).
#'
#'
#' At design points (i.e., the input \code{x} values), the function takes the average between the left and right slopes (on the edges the inside slope is technically replicated to the outside). If \code{allowZero=FALSE} (default), the algorithm gradually expands the x range over which slope is observed (by increments of one average \code{x} spacing), until a positive slope results. If the input is completely flat in \code{y} and \code{allowZero=FALSE}, the function returns \code{NA}s. 

#' @param x numeric or integer: input x values, must be strictly increasing
#' @param y numeric: input y values, must be monotone (can be non-strict) and in line with the direction specified by \code{decreasing}
#' @param outx numeric or integer: x values at which slopes are desired (default: same as input values)
#' @param allowZero logical: should zero be allowed in the output? Default \code{FALSE}
#' @param tol tolerance level: when \code{allowZero=FALSE}, slope below that value is considered zero. Default 1e-2. Might need to change if you use unusual units for x or y.
#' @param full logical: should a more detailed output be provided? Default \code{FALSE} (see details under 'Value').
#' @param decreasing logical: is input supposed to be monotone decreasing rather than increasing? Default \code{FALSE}

#' @return If \code{full=FALSE}, returns a vector of slopes at the points specified by \code{outx}. 
#' @return If \code{full=TRUE}, returns a list with slopes at the design point (\code{rawslopes}), the initial guess at output slopes (\code{initial}), and the official final ones (\code{final}). 

#' @seealso \code{\link{deltaInverse}}, which uses this function.
#' @export

slope <- function(x, y, outx=x, allowZero=FALSE, tol=1e-2, full=FALSE, decreasing=FALSE)
{
### Validation (might be mostly redundant if using doseResponse as input)
y=round(y,8)  # underflow error prevention
if (any(outx>max(x) | outx<min(x))) stop("No extrapolation allowed in 'slopes'.\n")
m=length(x)
if(length(y)!=m) stop("Mismatched lengths in 'slopes'.\n")
if(decreasing) y=-y
xdiffs=diff(x)
ydiffs=diff(y)
if (any(xdiffs<=0 | ydiffs<0)) stop("Monotonicity violation in slope().\n")
if (y[1]==y[m] && !allowZero)  return(rep(NA,length(outx))) # degenerate flat case; no solution
slopes=ydiffs/xdiffs
sslopes=c(slopes[1],slopes,slopes[m-1])  ### so that the edges get only the inward-side slope

interval=findInterval(outx,x)
## The trivial ones
candidate=slopes[interval]
## ones falling on design points get the average of the two segments
design=which(outx %in% x)
if (length(design)>0) {
	for(a in design) candidate[a]=(sslopes[interval[a]]+sslopes[interval[a]+1])/2
}
candidate0=candidate

## tougher nut: zero slope
if(!allowZero && any(candidate<tol))
{
  candidate[candidate<tol] = tol
	# xstep=mean(xdiffs)
	# for (a in which(candidate<tol))
	# {	
	# 	b=0
	# 	while(candidate[a]<tol)
	# 	{
	# 		b=b+1
	# 		xends=c(max(x[1],outx[a]-b*xstep),min(x[m],outx[a]+b*xstep))
	# 		yends=approx(x,y,xout=xends)$y
	# 		candidate[a]=diff(yends)/diff(xends)
	# 	}
	# }
}
if(decreasing) y=-y
if(!full) return (candidate)
return(list(rawslopes=slopes,initial=candidate0,final=candidate))		
}

#' Shrinkage fix to mitigate bias in observed rates, under adaptive dose-finding designs
#'
#' Adaptive dose-finding designs induce a bias on observed rates
#' away from the target dose. This is well-known in other adaptive-design fields,
#' but has been overlooked by the dose-finding research community.
#' Flournoy and Oron (2020) examine the bias in the dose-finding context,
#' and suggest a simple shrinkage fix that reduces both bias and variance.

#' The fix is analogous to the empirical-logit fix for zero counts in binary data, 
#' but instead of adding 0.5 to each cell, \code{target} is added to the 1's at each dose, 
#' and \code{1-target} to the 0's.
#' The shrinkage is applied to the raw observation, so CIR or IR are carried out
#' on the shrunk data.

#' @references Flournoy N and Oron AP, 2020. Bias Induced by Adaptive Dose-Finding Designs. Journal of Applied Statistics 47, 2431-2442.
#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}
#' 
#' @param y  can be either of the following: y values (response rates), a 2-column matrix with positive/negative response counts by dose, a \code{\link{DRtrace}} object or a \code{\link{doseResponse}} object. 
#' @param x dose levels (if not included in y). 
#' @param wt0 weights (if not included in y).
#' @param target the balance point (between 0 and 1) around which the design concentrates allocations.
#' @param swt the weight of the shrinkage. Default 1 (a single observation)
#' @param nmin the minimum n at each dose, for the shrinkage to be applied. Default 2.
#' @param ... parameters passed on to \code{doseResponse()} 
#' 
#' @example inst/examples/shrinkExamples.r
#' 
#' @export
DRshrink<-function(y, x=NULL, wt0=NULL, target, swt=1, nmin=2, ...) 
{
if(length(target)>1) stop('Shrinkage target must be a single constant.\n')
dr = doseResponse(y=y,x=x,wt=wt0,...)
dr$y=ifelse(dr$weight<nmin,dr$y,round((dr$y*dr$weight+target*swt)/(dr$weight+swt),8))
#print(dr)
return(dr)
}

