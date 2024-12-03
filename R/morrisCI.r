### Morris and Morris-style CI: functions and accessories
# Evaluation functions used by recursive algorithm (not exported)
# This is G(t_j,theta_j) in Morris (1988) equation (4.3), with typo correction
Gupper<-function(theta,y,n,j)
# conversion of Morris (4.3) to this function:
# s in Morris -> y here
# j index is the same
# theta is the same
# there's no need for Morris' t, it is a hypothetical variable
{
if(j==length(y)) return (pbinom(q=y[j],size=n[j],prob=theta))
pbinom(q=y[j]-1,size=n[j],prob=theta)+dbinom(x=y[j],size=n[j],prob=theta)*Gupper(theta=theta,y=y,n=n,j=j+1) ## The typo is the use of G_{j+1} rather than G_j throughout
}	

Glower<-function(theta,y,n,j)
{
if(j==1) return (pbinom(q=y[j]-1,size=n[j],prob=theta,lower.tail=FALSE))
pbinom(q=y[j],size=n[j],prob=theta,lower.tail=FALSE)+dbinom(x=y[j],size=n[j],prob=theta)*Glower(theta=theta,y=y,n=n,j=j-1)
}	


# Confidence intervals for ordered Binomial, based on Morris 1988: recursive algorithm (not exported)

 
morrisUCL<-function(y,n,halfa=0.05)
{
m=length(y)
if(length(n)!=m) stop("Mismatched lengths in Morris.\n")
# weird prep...
uout=rep(1,m)
a=m
### At uppermost doses as long as phat=1, no need for algorithms
while((y[a]==n[a] || n[a]==0) && a>=1) a=a-1
if(a<1) return(uout)

for(b in a:1) ### core algorithm
{
	uout[b]=uniroot(function(theta,h,d,alpha,...) h(theta=theta,...)-alpha,interval=c(0,1),
		alpha=halfa,j=b,h=Gupper,n=n,y=y)$root
}
if(any(n==0)) # Degenerate case with "phantom boundaries"
{
uout = approx((1:m)[n>0],uout[n>0],xout=1:m,rule=2)$y
}
return(uout)
}
#####

morrisLCL<-function(y,n,halfa=0.05)
{
m=length(y)
if(length(n)!=m) stop("Mismatched lengths in Morris.\n")
# weird prep...
uout=rep(0,m)
a=1
### At lowermost doses as long as phat=0, no need for algorithms
while((y[a]==0 || n[a]==0) && a<=m) a=a+1
if(a>m) return(uout)

for(b in a:m)  ### core algorithm
{
	uout[b]=uniroot(function(theta,h,d,alpha,...) h(theta=theta,...)-alpha,interval=c(0,1),
		alpha=halfa,j=b,h=Glower,n=n,y=y)$root
}
if(any(n==0)) # Degenerate case with "phantom boundaries"
{
uout = approx((1:m)[n>0],uout[n>0],xout=1:m,rule=2)$y
}
return(uout)
}
#####

##' Analytical ordered-binary-Y confidence intervals, using the Morris (1988) algorithm
##'
##'
##' Analytical confidence intervals for CIR and IR, using the recursive algorithm by Morris (1988), equation (4.3), for ordered-binary-Y point estimates. Optionally, the intervals are narrowed further using a backup (unordered) interval estimate at each individual x value.
##'
##' 
##' The default for backup is Wilson's (\code{wilconCI}). Also available are Jeffrys' (\code{jeffCI}) and Agresti-Coull (\code{agcouCI}).

#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}
#' @export
#' @seealso \code{\link{isotInterval}}
#' @example inst/examples/fwdCiExamples.r
#' 
#' @references Morris, M., 1988. Small-sample confidence limits for parameters under inequality constraints with application to quantal bioassay. Biometrics 44, 1083-1092.

#' @note This function found and corrected a typo in equation (4.3), namely the use of G_(j+1) in the recursion. The recursion cannot start in this way. Rather, it is the use of theta_(j+1) that delivers information from adjacent doses. Or perhaps in other words, there is only one G function rather than a different one for each dose. The correction has been verified by reproducing the numbers in the Morris (1988) example (Table 1), and also approved by the original author.
#' 
#' @return A two-column matrix with the same number of rows as \code{length(phat)}, containing the calculated lower and upper bounds, respectively.

#' @param y integer or numeric vector, the pointwise Binomial counts
#' @param n integer or numeric vector, the pointwise sample sizes
#' @param phat numeric vector, the point estimates. Defaults to \code{y/n}, but when called by \code{\link{isotInterval}} is overridden by the actual CIR/IR point estimate.
#' @param conf numeric, the interval's confidence level as a fraction in (0,1). Default 0.9.
#' @param narrower logical, if the \code{alternate}-produced interval is narrower at any point, should it replace the Morris result? Also, can we enforce straightforward monotonocity to narrow the bounds? Default \code{TRUE}.
#' @param alternate function to use for alternate pointwise interval. Default \code{wilconCI}.
#' @param ... parameters passed on to \code{alternate}.
morrisCI<-function(y,n,phat=y/n,conf=0.9,narrower=TRUE,alternate=wilsonCI,...)
{
if(conf<=0 || conf>=1) stop("Confidence must be between 0 and 1.\n")
tailp=(1-conf)/2
lcl=morrisLCL(y=y,n=n,halfa=tailp)
ucl=morrisUCL(y=y,n=n,halfa=tailp)

if(narrower) # Optional pointwise narrowing via alternate interval function
# See Oron and Flournoy (2017) for justification/performance
{
	relevants=which(n>0) # Avoiding n=0 boundary where alternate() produces NaNs
	altout=alternate(phat=phat[relevants],n=n[relevants],conf=conf,...)
# The cummax, cummin (added Dec. 2015) ensure monotonicity of boundaries.
# Monotonicity is not for the optics, but rather another way to pool information
# from where it is plentiful to where it might be lacking.

	lcl[relevants]=cummax(pmax(lcl[relevants],altout[,1]))
	ucl[relevants]=rev(cummin(rev(pmin(ucl[relevants],altout[,2]))))
}

return(cbind(lcl,ucl))
}

	
#uout[[a]]=wilsonCI(phat=x[a]/n[a],n=n[a],conf=1-2*halfa)[2]
#if(a<=1) return(uout)

	