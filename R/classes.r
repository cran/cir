
##' @rdname DRtrace
#' @export

is.DRtrace<-function(dr)
{
if(!inherits(dr,"DRtrace")) return(FALSE)
if(!inherits(dr,"data.frame")) return(FALSE)
if(!all(c("x","y","weight") %in% names(dr))) return(FALSE)
if(!all(sapply(dr[,c("x","y","weight")],is.numeric))) return(FALSE)
if(any(dr$weight<0)) return(FALSE)
return(TRUE)
}
##' @rdname DRtrace
#' @export

is.doseResponse<-function(dr)
{
if(!inherits(dr,"doseResponse")) return(FALSE)
if(any(duplicated(dr$x))) return(FALSE)
if (any(dr$x!=sort(dr$x))) return(FALSE)
return(TRUE)
}

#########################

#' Constructor functions and class-checking functions for DRtrace and doseResponse classes
#'
#'
#' Functions to create and sanity-check objects of the \code{DRtrace} (dose-response experiment trace/trajectory) and \code{doseResponse} (dose-response raw summary) classes. Note that the latter inherits from the former, purely for programming-convenience reasons.
#'
#'
#' The input argument \code{y} can include the entire information, or as little as the \code{y} vector of responses (for a \code{DRtrace} object) or response rates (\code{doseResponse}). When including the entire information, it has to be a data frame with at least \code{y} (both y and x for \code{DRtrace}), or a two-column matrix with 'yes' and 'no' responses (assumed in this order, but can be the reverse with \code{noyes=TRUE}). In this case the doses \code{x} can be provided as a separate vector, or as the matrix row names. \code{doseResponse} will return an error if there are any duplicates in \code{x}.
#'
#' Even though both \code{DRtrace} and \code{doseResponse} accept two-column yes/no matrix input, the interpretation is different. For the former, this form of input is intended mostly to enable shorthand input when the experiment was run in cohorts. Each row represents a cohort's results, and rows must be in the order the experiment was run. For the latter, the yes-no table is a summary tabulation of responses and is treated accordingly, including rearrangement of rows to increasing \code{x}.

#'
#'
#' @aliases doseResponse is.doseResponse is.DRtrace
##' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}
##' 
#' @example inst/examples/classExamples.r
##' @seealso \code{\link{cirPAVA}}, \code{\link{plot.doseResponse}},\code{\link{plot.DRtrace}} 

#' @param y,x,wt  see Details.
#' @param noyes logical, in case of a 2-column input is the 1st column 'no'? Default \code{FALSE}, meaning the 1st column is 'yes'.
#' @param dr the object being checked
#' @param ... parameters passed on to \code{DRtrace()}, or ignored.
#' 
#' @return For constructor functions, the relevant object. For checking functions, a logical value indicating whether the object meets class definition.
#' @export

DRtrace<-function(y,x=NULL,wt=NULL,noyes=FALSE,...)
{
if(is.DRtrace(y)) return(y)

if(is.data.frame(y)) # data frame input, e.g, from read.csv
{
	if(!('x' %in% names(y)) || !('y' %in% names(y))) stop('Data frame must include variables named x and y.\n')
	x=y$x
# Code will accept either 'wt' or 'weight', with priority for the latter if both exist
	if('wt' %in% names(y)) wt=y$wt
	if('weight' %in% names(y)) wt=y$weight
	y=y$y  # warning: this overwrites y the data frame (inside the function)
}
ll=dim(y)

if (length(ll)>2 || (length(ll)==2 && ll[2]!=2)) stop ("y can be a vector, data frame, or yes-no table.\n")

if (length(ll)==2 && ll[2]==2) { # converting a yes-no table

	xvals<-x
	if(is.null(xvals)) if(is.null(rownames(y))) xvals=1:ll[1] else xvals=rownames(y)

	if(length(xvals) != ll[1]) stop("Mismatched lengths.\n")
	
	yntab=y
	if(noyes) yntab=yntab[,2:1]  # if the table is no-yes rather than yes-no, we reverse it here.
	y=unlist(apply(yntab,1,function(a,b)rep(b,a),b=1:0))
	x=rep(xvals,rowSums(yntab))
	
 #   y<-c(rep(1,sum(yntab[,1])),rep(0,sum(yntab[,2])))
 #   x<-c(rep(xvals,yntab[,1]),rep(xvals,yntab[,2]))
   wt<-rep(1,length(y))  ## in case of yes-no table we ignore incoming weights 
#	warning("Raw data is a yes-no table; therefore, observation order is arbitrary.\n")
}

n=length(y)
if(is.null(wt)) wt=rep(1,n)
if(length(x)!=n || length(wt)!=n) stop("Mismatched lengths.\n")

tout<-data.frame(x=x,y=y,weight=wt)
attr(tout,'class')<-c('DRtrace','data.frame')
return(tout)
}


#############
##' @rdname DRtrace
#' @export

doseResponse<-function(y,x=NULL,wt=rep(1,length(y)),noyes=FALSE,...)
{
if(is.doseResponse(y)) return(y)
if(!is.DRtrace(y) && is.data.frame(y)) # data frame input, e.g, from read.csv
{
	if(!('y' %in% names(y))) stop('Data frame must include variables named x and y.\n')
	if('x' %in% names(y)) x=y$x
# Code will accept either 'wt' or 'weight', with priority for the latter if both exist
	if('wt' %in% names(y)) wt=y$wt
	if('weight' %in% names(y)) wt=y$weight
	y=y$y  # warning: this overwrites y the data frame (inside the function)
}
ll=dim(y)
## When y is given as vector and no x, then x defaults to dose indices
if(is.null(ll) && is.null(x)) x=1:length(y)

## DRtrace conversion, or cases for doing DRtrace first
if(is.DRtrace(y) || any(duplicated(x)) || any(diff(x)<0))  
 {
	if(is.DRtrace(y)) {z=y} else z<-suppressWarnings(DRtrace(y=y,x=x,wt=wt,...))
	tout<-data.frame(x=sort(unique(z$x)),y=tapply(z$y,z$x,mean),weight=tapply(z$weight,z$x,sum))

} else if(length(ll)==2 && ll[2]==2) 
## Now, two-column yes-no matrix case
{
	if(noyes) y=y[,2:1]  # if the table is no-yes rather than yes-no, we reverse it here.

	if(is.null(x)) if(is.null(rownames(y))) x=1:ll[1] else x=rownames(y)
	tout=data.frame(x=x,y=(y[,1]/rowSums(y)),weight=rowSums(y))
	
} else if (length(x)==length(y) && (length(wt)==length(y) || is.null(wt)))  # straightforward x-y-wt input
{
	if(is.null(wt)) wt=rep(1,length(y))
	tout<-data.frame(x=x,y=y,weight=wt)
} else stop("Incompatible input data. Check the help.\n")

attr(tout,'class')<-c('doseResponse','data.frame')
# Final sanity check: duplicates in x?
if(any(duplicated(tout$x))) stop('X values must be unique.\n')
# Reordering rows if given in non-monotone order
tout=tout[order(tout$x),]
return(tout)
}








