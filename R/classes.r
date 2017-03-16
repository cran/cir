
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
if(!is.DRtrace(dr)) return(FALSE)
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
#' @aliases doseResponse is.doseResponse is.DRtrace
##' @author Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}
##' 
#' @example inst/examples/classExamples.r
##' @seealso \code{\link{cirPAVA}}, \code{\link{plot.doseResponse}},\code{\link{plot.DRtrace}} 

#' @param y,x,wt  see help to \code{\link{cirPAVA}}.
#' @param noyes logical, in case of a 2-column input is the 1st column 'no'? Default \code{FALSE}, meaning the 1st column is 'yes'.
#' @param dr the object being checked
#' @param ... (\code{doseResponse()} only) parameters passed on to \code{DRtrace()}
#' 
#' @return For constructor functions, the relevant object. For checking functions, a logical value indicating whether the object meets class definition.
#' @export

DRtrace<-function(y,x=NULL,wt=NULL,noyes=FALSE)
{
if(is.DRtrace(y)) return(y)

ll=dim(y)

if (length(ll)>2 || (length(ll)==2 && ll[2]!=2)) stop ("For conversion to dose-response, y can only be a vector or yes-no table.\n")

if (length(ll)==2 && ll[2]==2) { # converting a yes-no table

	xvals<-x
	if(is.null(xvals)){
		warning("No dose (x) values provided; using natural numbers.\n")
		xvals=1:ll[1]
	}
	if(length(xvals) != ll[1]) stop("Mismatched lengths.\n")
	
	yntab=y
	if(noyes) yntab=yntab[,2:1]  # if the table is no-yes rather than yes-no, we reverse it here.
    y<-c(rep(1,sum(yntab[,1])),rep(0,sum(yntab[,2])))
    x<-c(rep(xvals,yntab[,1]),rep(xvals,yntab[,2]))
    wt<-rep(1,length(y))  ## in case of yes-no table we ignore incoming weights 
	warning("Raw data is a yes-no table; therefore, observation order is arbitrary.\n")
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

doseResponse<-function(y,x=NULL,wt=rep(1,length(y)),...)
{
if(is.doseResponse(y)) return(y)

if(is.null(x)) x=1:length(y)


if(any(duplicated(x)) || any(diff(x)<0))  # Cases for doing DRtrace first
{
	z<-suppressWarnings(DRtrace(y=y,x=x,wt=wt,...))
	tout<-data.frame(x=sort(unique(z$x)),y=tapply(z$y,z$x,mean),weight=tapply(z$weight,z$x,sum))

} else if (length(x)==length(y) && (length(wt)==length(y) || is.null(wt)))  # straightforward x-y input
{
	if(is.null(wt)) wt=rep(1,length(y))
	tout<-data.frame(x=x,y=y,weight=wt)
} else stop("Incompatible input data. Check the help.\n")

attr(tout,'class')<-c('doseResponse','DRtrace','data.frame')
return(tout)
}








