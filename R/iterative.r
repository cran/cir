#' Centered-isotonic-regression (CIR): iterative version for better bias correction
#'
#'
#' EXPERIMENTAL: Nonparametric forward point estimation of a monotone response (y) as a function of dose (x), using an iterative version of the centered-isotonic-regression (CIR) algorithm. The code works, but delivers marginal improvement at greater computational cost (an issue if you simulate a large ensemble), and somewhat convoluted interpretation. Use at your own risk.
#' For explanation, see Oron and Flournoy (2017), Section 3.2.
#' 
#' @param y See \code{\link{cirPAVA}}
#' @param outx vector of x values for which predictions will be made. If \code{NULL} (default), this will be set to the set of unique values in the \code{x} argument (or equivalently in \code{y$x}).
#' @param tol The iteration's convergence tolerance level (default 1e-3)
#' @param maxit integer, maximum number of iterations (default 10)
#' @param full logical, is a more complete output desired? if \code{FALSE} (default), only a vector of point estimates for y at the provided dose levels is returned
#' @param ...	Other arguments passed on to \code{\link{cirPAVA}}

#' @return under default, returns a vector of y estimates at unique x values. With \code{full=TRUE}, returns a list of 3 \code{\link{doseResponse}} objects named \code{output,input,shrinkage} for the output data at dose levels, the input data, and the function as fit at algorithm-generated points, respectively.

#' @seealso \code{\link{cirPAVA}},\code{\link{quickIsotone}}
#' @export


iterCIR<-function(y,outx=NULL,tol=1e-3,maxit=10,full=FALSE,...){
  if(is.null(outx)) outx=y$x
  outx0=outx
  cand=cirPAVA(y,...)
  cand0=y$y
  y0=y
  it=0
  
  while(max(abs(cand-cand0))>tol) 
  {
    
 #   cat(cand,'\n')
    cand0=cand
    y0=y
    y0$weight=y$weight/(cand*(1-cand))
    y0$weight[!is.finite(y0$weight)]=y$weight[!is.finite(y0$weight)]
    cand=cirPAVA(y0,...)
    it=it+1
    if(it>maxit) return(cirPAVA(y,outx=outx,full=full,...))
  }
  tmp=cirPAVA(y0,full=TRUE,outx=outx0,...)
  if(!full) return(tmp$output$y[match(outx,tmp$output$x)])
  tmp$input=y
  tmp$output$weight=y$weight
  return(tmp)
}
