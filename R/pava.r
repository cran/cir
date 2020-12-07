##' Returns standard isotonic-regression estimate, with flexible dose-response input
#'
#'
#' Nonparametric forward point estimation of a monotone response (y), using the standard isotonic-regression pool-adjacent-violators algorithm (PAVA). Core code from Raubertas (1994) with many modifications.
#'
#'
#'  Compute the isotonic regression of a numeric vector 'y', with
#'  weights 'wt', with respect to simple order. The core algorithm is still the one
#' coded by R.F. Raubertas, dated 02 Sep 1994. However, the input and output modules have been
#' modified to allow more flexible formats in either direction.
#' note that unlike centered-isotonic-regression (CIR, see \code{\link{cirPAVA}}), this algorithm does not use the dose (x) values at all. For a discussion why CIR is preferred over "plain-vanilla" PAVA, see Oron and Flournoy (2017).
 
#  
##' @author C.R. Raubertas, Assaf P. Oron \code{<aoron.at.idmod.org>}
#' @references Oron, A.P. and Flournoy, N., 2017. Centered Isotonic Regression: Point and Interval Estimation for Dose-Response Studies. Statistics in Biopharmaceutical Research, In Press (author's public version available on arxiv.org).


#' @param y  can be either of the following: y values (response rates), a 2-column matrix with positive/negative response counts by dose, a \code{\link{DRtrace}} object or a \code{\link{doseResponse}} object. 
#' @param x dose levels (if not included in y). Note that the PAV algorithm doesn't really use them. 
#' @param wt weights (if not included in y).
#' @param outx vector of x values for which predictions will be made. If \code{NULL} (default), this will be set to the set of unique values in the x argument (or equivalently in y$x). Non-NULL inputs are relevant only if \code{full=TRUE}.
#' @param full logical, is a more complete output desired? if \code{FALSE} (default), only a vector of point estimates for y at the provided dose levels is returned
#' @param dec logical, is the true function is assumed to be monotone decreasing? Default \code{FALSE}.
#' @param adaptiveShrink logical, should the y-values be pre-shrunk towards an experimental target? May be relevant if data were obtain via an adaptive dose-finding design. See \code{\link{DRshrink}}.
#' @param ...	Other arguments passed on to the constructor functions that pre-process the input.

#' @return under default, returns a vector of y estimates at unique x values. With \code{full=TRUE}, returns a list of 3 \code{\link{doseResponse}} objects named \code{output,input,shrinkage} for the output data at dose levels, the input data, and the function as fit at algorithm-generated points, respectively. For this function, the first and third objects are identical.

#' @seealso \code{\link{cirPAVA}}
#' @export

oldPAVA<-function (y,x=NULL,wt=rep(1,length(x)),outx=NULL,full=FALSE,dec=FALSE,adaptiveShrink=FALSE,...) {

### converting to doseResponse object 
### Basically it's a numeric data frame with x,y,weight, and x increasing

# in case of old-style input of y only:
if(is.null(x)) { x=1:length(y); wt=rep(1,length(x))}

dr=doseResponse(y,x,wt,...)
if (any(is.na(dr))) stop ("Missing values are not allowed.\n")  
# Optional pre-shrinking of y for adaptive designs
if(adaptiveShrink) dr=DRshrink(y=dr,...)

### Predictions will be delivered for x=outx
if(is.null(outx)) outx=dr$x
if(min(outx)<min(dr$x) || max(outx)>max(dr$x)) stop("Cannot predict outside design boundaries.\n")

m <- dim(dr)[1]
if (m <= 1) {  ## degenerate case: only one dose level
if (!full) return (dr$y)
else return(list(output=dr,input=dr,shrinkage=dr))
}


dr0=dr ## clean copy of input data
### Decreasing monotone case: simple fix
if (dec) dr$y = -dr$y

#### Core algorithm

lvlsets <- (1:m)
 
repeat {
      viol <- (as.vector(diff(dr$y)) < 0)  # Find adjacent violators
      if (!(any(viol))) break
 
      i <- min( (1:(m-1))[viol])        # Pool first pair of violators
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i+1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      dr$y[ilvl] <- sum(dr$y[ilvl]*dr$weight[ilvl]) / sum(dr$weight[ilvl])
      lvlsets[ilvl] <- lvl1
}

if (dec) dr$y = -dr$y

outy=approx(dr$x,dr$y,outx,rule=2)$y

if (!full) {
	return(outy)
} else {
	
	if(all(outx %in% dr0$x)) {
		dr1=dr0
		dr1$y=outy
		dr1=dr1[match(outx,dr1$x),]
	} else dr1=doseResponse(y=outy,x=outx,wt=rep(0,length(outy)))
	return(list(output=dr1,input=dr0,shrinkage=dr))   }
}
