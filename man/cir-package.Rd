\name{cir-package}

\alias{cir-package}

\alias{Isotonic Regression and Centered Isotonic Regression}

\alias{cir}

\docType{package}

\title{Nonparametric estimation of monotone
functions via isotonic regression and centered isotonic regression }

\description{
The 'cir' package provides a
documented version of the well-known isotonic regression (IR)
algorithm (function 'pava'), and an improvement to IR in case the true
function is known to be smooth and strictly monotone. This improvement called Centered Isotonic Regression (CIR) is available via the function 'cir.pava'. Additionally, the function 'cir.upndown' provides percentile estimation for dose-response experiments (e.g., ED50 estimation of a medication) using CIR.
}

\details{

  \tabular{ll}{
Package: \tab cir\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-01-31\cr
License: \tab                          GPL (version 2 or later)\cr
}
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

'cir.upndown' provides an easy interface for direct percentile
estimation following a binary-response dose-finding experiment (such as
'up-and-down'). The function also includes an interval estimate for the
target percentile.

}


\author{ Assaf P. Oron

Maintainer: Assaf P. Oron <aoron@fhcrc.org> 

}

\references{
Barlow R.E., Bartholomew D.J., Bremner J.M. and Brunk H.D.,
  Statistical Inference under Order Restriction. John Wiley & Sons
  1972. 

Robertson T., Wright F.T. and Dykstra R.L. Order Restricted Statistical
Inference, Wiley, Chichester (1988).

Oron A.P., Up-and-Down and the Percentile-Finding Problem. Doctoral
Dissertation, University of Washington. 2007



}


\keyword{ package }
\keyword{ nonparametric }
\keyword{ survival }
\seealso{
 \code{\link[stats:]{isoreg}} 
}

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

### Play with 'wt.overwrite' to see how it affects the solutions

### With CIR, the ED50 will be unique (though hard to directly pinpoint from
### the default vector output of 'cir.pava')
### We can use 'cir.upndown' for direct estimation of ED50 and its confidence
### interval on the same data using CIR.

levo.cir.outcome=cir.upndown(yesno=levo,xseq=levdoses,target=0.5,full=TRUE)

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
