\name{cir.pava}
\alias{cir.pava}
\alias{Centered Isotonic Regression (CIR)}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function for Centered Isotonic Regression (CIR)}
\description{
Performs a modified version of isotonic regression (IR), more appropriate
when the true function is assumed strictly monotone and smooth. No
parameteric assumptions or smoothing parameters are needed. Output is
piecewise-linear like IR, but avoids 'flat' stretches.
}
\usage{
cir.pava(y, x, wt = rep(1, length(x)), boundary = 2, full = FALSE, dec =
FALSE,wt.overwrite=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{y values (responses). Can be a vector or a two-column yes-no
    table (for binary responses). y is used as first argument, both for compatibility with 'pava' 
 and to enable parallel running via 'apply' type routines. Order of y
 must match that of x. }
  \item{x}{x values (treatments). Need to be pre-sorted in increasing order. }
  \item{wt}{Weights. Will be overwritten with observation counts (row
    sums) in case of a yes-no input format for y.}
  \item{boundary}{Action on boundaries. See 'Details' below.}
  \item{full}{If FALSE, only point estimates at x values are returned;
    otherwise, a more detailed list. See 'Value' below.}
  \item{dec}{Is the true function monotone decreasing (defaults to FALSE)?}
  \item{wt.overwrite}{Should the variable 'wt' be recalculated as the
    observation counts in each row? Defaults to TRUE. Applicable only for yes-no table input.}
}
\details{
Isotonic regression (IR, Barlow et al. 1972) replaces
monotonicity-violating sequences of observations with a 'flat' stretch
whose y value is the weighted average of the original observations. This
is the non-parametric MLE under order restrictions. IR is implemented as
\code{\link{pava}} in this package (PAVA stands for Pooled Adjacent
Violators Algorithms -- a fancy name for a very simple procedure); and also in a somewhat-crippled
version as \code{\link[stats]{isoreg}} in the stats package.

If it is known that the original function is strictly increasing and
reasonably smooth (i.e., has at least twice continuously
differentiable), then IR's performance can be improved by replacing the
'flat' stretches with a strictly-increasing estimate. CIR does precisely
this, in the simplest way: the weighted-average estimate is placed at a
point that is the weighted-average of corresponding x values, and
function values between points are estimated via linear
interpolation. When there are no monotonicity violations in the input
data, CIR provides an identical output to IR, which is simply to return
the original y values. More details are in Oron (2007), Chapter 3. 

Data can be provided as paired x-y values (x,y in two separate vectors) or,
for dose-response style applications, with y as a table summarizing 'yes'
and 'no' responses, with each dose summarized on one row ('yes' would be
column 1), and a matched x vector giving the doses. For the latter, it is okay to give
all-zero rows (i.e., rows with no observations); the function will get
rid of them and of the redundant x values.

If y is a yes-no table, the weights will be set as the observation
  counts in each row - UNLESS 'wt.overwrite' is set to FALSE, in which
  case it will use the given value of 'wt'.


If n is large (observations at many distinct values of x), then you
might do better using more sophisticated "smoothers" (kernel estimators,
splines, etc.). However, besides the extra complication in choosing parameter
values for the smoothing, as this goes to print (early 2008) I am not
aware of a reliable enough, plug-and-play R code for this task. Your best bet may
be Jim Ramsay's \code{\link[fda]{smooth.monotone}} function in the \texttt{fda} package (Ramsay developed a spline algorithm on a
transformed scale in a way that ensures monotonicity; see Ramsay
(1998)). But this requires some understanding of data structures specific to
that package.

In any case, CIR provides a nice, no-moving parts, nonparametric benchmark to compare
against any such "smoothers", even when you use them. Of course, if the
original function is expected to be staircase-like, then neither CIR nor
most "smoothers" are preferable over plain IR.

The most common potential caveat with CIR are boundaries. When there is a
monotonicity violation involving the largest or smallest x values, CIR
needs to be told how to extrapolate via the "boundary" parameter. The default option is
"boundary=2", which creates 'flat' intervals near the boundaries in such a
case (an output identical to IR; so you can do no worse). "boundary=1"
does linear interpolation (not recommended in general).

In case you have meaningful prior information you'd prefer to use on the
boundaries (e.g., cumulative menopause rate at age zero must be
zero), use the default boundary option but fix the situation outside the function: if
the information is deterministic (i.e., has infinite weight), use
"full=TRUE" to get the output.x and output.y, and augment them with the
information (say, x=0,y=0) before interpolating to get the final answer.
If you consider the information to have a finite
weight (e.g., equal to the sample size), add it to the *input* using the
appropriate weight.

If you need calibration/inverse estimation (guessing an x value), use
"full=TRUE" and \code{\link{approx}}, or in the dose-response case use \code{\link{cir.upndown}}
}
\value{

  If full==FALSE, you get only the new y estimates at the x values.

  If full==TRUE, you get a list:
  \item{output.y}{estimates of y at original x values. Same as the
  output in case full==FALSE; produced from alg.y via interpolation }
  \item{original.x}{original x values}
  \item{alg.x}{x values - *NOT* the original ones, but the ones calculated
  at the algorithm's final stage}
  \item{alg.y}{corresponding final estimates of y}
  \item{alg.wt}{corresponding final weights}
}
\references{Barlow R.E., Bartholomew D.J., Bremner J.M. and Brunk H.D.,
  Statistical Inference under Order Restriction. John Wiley & Sons
  1972. 
Oron A.P., Up-and-Down and the Percentile-Finding Problem. Doctoral
Dissertation, University of Washington. 2007

Ramsay, J. O. (1998) Estimating smooth monotone functions. Journal of
the Royal Statistical Society, Series B, 60, 365-375.
}
\author{ Assaf Oron (assaf@u.washington.edu,aoron@fhcrc.org}
\note{ Note that unlike \code{\link{pava}}, \code{\link{cir.pava}}
  requires the x values as input. 
}
\section{WARNING }{If you provide y as a yes-no table, do *NOT* set
  'wt.overwrite' to FALSE unless you really want to input different
  weights to the 'wt' variable. If you leave 'wt.overwrite' as is, the
  function will calculate the correct weights, i.e., obesrvation
  counts. Be aware that any weights other than observation counts for
  binary data will
  yield a non-standard solution, so tinker with them only if you know
  what you are doing.} 

\seealso{ Compare with \code{\link{pava}} for plain IR (a more limited
  IR version is available in \code{\link[stats]{isoreg}}), and with the
  sophisticated smoothing of \code{\link[fda]{smooth.monotone}}. For
  percentile (inverse) estimation in the dose-response case, see
  \code{\link{cir.upndown}}. }

\examples{
### In the 'stackloss' dataset, escape of ammonia through some plant's
### chimney appears driven mostly by plant operation rate with a clearly
### monotone dependence. Linearity is questionable, though, and there
### are monotonicity violations in the data.
### There are are 21 observations at 8 distinct rates, and the original dataset is not ordered.
### "pava" and "cir.pava" require unique and ordered x values.
### So this example also shows how to prepare such data for input to
### "pava" or "cir.pava" (not difficult):

data(stackloss)
attach(stackloss)

meanrate=sort(unique(Air.Flow))
meanloss=sapply(split(stack.loss,Air.Flow),mean)/10 ## according to stackloss documentation, this turns the data into percent loss
weights=sapply(split(stack.loss,Air.Flow),length) ### we don't want to lose the effect of multiple observations at certain points

### Raw data shows overall monotone pattern, linearity questionable, but
### perhaps not enough points for fancy smoothers 
plot(meanrate,meanloss,main="CIR Example (Stack Loss data)",xlab="Plant Operation Rate (Air Flow)",ylab="Mean Ammonia Loss Through Stack (percent)")

### PAVA gives a staircase solution in black
lines(meanrate,pava(meanloss,wt=weights))

### try CIR for a much more realistic curve in red
lines(meanrate,cir.pava(y=meanloss,x=meanrate,wt=weights),col=2)
 
### Compare with standard linear regression line in blue
abline(lsfit(meanrate,meanloss,wt=weights),col=4)


### This just to display what the "full=T" option provides:
cir.pava(y=meanloss,x=meanrate,wt=weights,full=TRUE)

######## yes-no table example #####
### Taken from Lacassie and Columb
### Anesth. Analg. 97, 1509-1513, 2003.

levo=cbind(c(0,2,2,4,2,1,0,1),c(3,3,5,3,2,1,1,0))

levo

### you should get this table:

###     [,1] [,2]
#[1,]    0    3
#[2,]    2    3
#[3,]    2    5
#[4,]    4    3
#[5,]    2    2
#[6,]    1    1
#[7,]    0    1
#[8,]    1    0

### Note that all doses except the lowest and highest are involved in
### some monotonicity violation (in terms of observed frequency of 'yes' responses)

pava(levo)

### Since the experiment's goal was to estimate the ED50 of the drug
### abbreviated here as 'levo', pava's solution is highly problematic as
### you can pick and choose your favorite ED50 from any of doses 4
### through 7!
###
### We call 'cir.pava' to our aid, meaning we need to specify x values
### for the doses:

levdoses=seq(0.25,0.425,0.025) ### values taken from the article

cir.pava(levo,x=levdoses)

### Now the ED50 will be unique (though hard to directly pinpoint from
### the default vector output of 'cir.pava'; try 'full=TRUE')
### see 'cir.upndown' for direct estimation of ED50 and its confidence
### interval on the same data using CIR.

### Also, play with 'wt.overwrite' to see how it affects the solutions

  }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nonparametric }
\keyword{ survival }
\keyword{smooth}
