\name{cir.upndown}
\alias{cir.upndown}

\title{ Nonparametric Percentile Estimation for dose-response
  experiments, using Centered Isotonic Regression (CIR)}
\description{
 Provides a unique point estimate for a specified percentile of a
 binary-response threshold CDF, and also an approximate interval
 estimate for the same percentile. This is done using linear
 interpolation and, in the event of monotonicity violations in observed
 response frequencies, centered isotonic regression (CIR).
}
\usage{
cir.upndown(yesno, xseq, target, xbounds = c(0, 1), ybounds = c(0, 1), full = FALSE, cioption = "poisson", plist = c(0.025, 0.975))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yesno}{ 2-column response tally table, with 'yes' counts in the left column and 'no' in the right }
  \item{xseq}{ vector with the corresponding treatments. Must be
    strictly increasing and with the same length as the number of rows
    in 'yesno'.}
  \item{target}{Target response frequency for the desired percentile,
    e.g. for median estimation set target=0.5}
  \item{xbounds,ybounds}{ parameters needed if estimate falls near boundary. See
  "Details".}
  \item{full}{Whether a full-detailed output is requested. If not, only
    point and interval estimates are returned.}
  \item{cioption}{ Type of interval estimation. See "Details." }
  \item{plist}{ Percent points for interval estimation, as fractions,
    e.g. for symmetric 95\% CI plist=c(0.025,0.975) (default).}
}
\details{
  This function uses the CIR function  \code{\link{cir.pava}}, and then
  performs inverse linear interpolation, that is: 'cir.pava's y-output
  becomes the x-coordinate when 'approx' is called to find the target percentile.

  Boundary options: any linear interpolation, CIR included, can run into
  trouble if the target to interpolate is outside the available range to
  interpolate upon. 'cir.upndown' solves that by forcing a reference x value at
  y=0 and 1. The user can set this reference according to prior
  knowledge (e.g., if there is zero risk at x=0 and 100\% risk at x=50,
  then set 'xbounds=c(0,50)'). In any case, the function then checks
  the input x values and expands the bounds to accommodate them if
  needed. The entire boundary-crossing response is only triggered if needed.

  Interval estimation options: interval estimation is done via
  linearly-smoothed quantile functions. The default option (a symmetrized
  Poisson) is adequate for sequential dose-response designs such as
  up-and-down, in which treatment allocation is random. If you use the
  function to analyze a nonrandom-assignment experiment, set
  "cioption='binomial'". Another ,less reommended but available option is 't', for
  t-distribution based quantiles. The utilities providing these
  functions are internal to the package, but some detail is available in
  Oron (2007), Section 3.3.
}
\value{
If full=FALSE, a vector with the first element being the point estimate,
and the other elements being interval estimate cutpoints (in the order given by 'plist')

If full=TRUE, a list with these objects:
  \item{raw}{The input yes-no table}
  \item{paved}{The output of the internal call to 'cir.pava', which is itself a
    list (see that function's documentation for details}
  \item{out}{The percentile point estimate}
  \item{ci}{The interval cutpoint estimates, in the same order as 'plist'}
}
\references{
  Barlow R.E., Bartholomew D.J., Bremner J.M. and Brunk H.D.,
  Statistical Inference under Order Restriction. John Wiley & Sons
  1972.

  Lacassie, H.J. and Columb, M.O. The Relative Motor Blocking Potencies
  of Bupivacaine and Levobupivacaine in Labor. Anesth. Analg. 97,
  1509-1513, 2003. (used in the example below)

  Oron A.P., Up-and-Down and the Percentile-Finding Problem. Doctoral
Dissertation, University of Washington. 2007
}
\author{ Assaf P. Oron}
\note{ 

}
\seealso{  \code{\link{cir.pava}}, \code{\link{pava}}}
\examples{

### Up-and-down example from Lacassie and Columb (2003):
### We are estimating the median (ED50) threshold for motor blocking
### potency of levobupivacaine in labor.
### Yes-no table is tallied from their Figure 1.

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

# The corresponding x values:

levdoses=seq(0.25,0.425,0.025)

### function call
levo.cir.outcome=cir.upndown(yesno=levo,xseq=levdoses,target=0.5,full=TRUE)

### This example is interesting, in that regular isotonic regression
###    would fit y-hat=0.5 for the entire interval between x=0.325 and x=0.4,
###    thereby not providing a unique point estimate for the ED50

### The authors' old-fashioned averaging estimator for ED50 yields
### 0.31 as the point estimate, and a probit regression (as reported in the article) 
### places it at 0.37

### The CIR point estimate (below) is around 0.345, almost exactly in the
### middle between these two

levo.cir.outcome$out

### The authors' 95\% confidence interval estimate is amazingly optimistic
### at (0.29,0.34). CIR's estimate is more conservative
### because it indirectly accounts for the gross monotonicity violations in the data
### It is more in line with the probit-regression estimate of (0.30,0.45)

levo.cir.outcome$ci

}
\keyword{ nonparametric }
\keyword{ survival }
\keyword{smooth}

