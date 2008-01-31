\name{pava}
\alias{pava}
\alias{Isotonic Regression}

\title{Univariate Isotonic Regression, a.k.a. Pool Adjacent Violators Algorithm (PAVA)}
\description{
  Isotonic regression is a standard nonparametric method to produce an
  estimate for a function known to be monotone (Barlow et al. (1972),
  Robertson et al. (1988)). This is a documentation and slight upgrade
  of a code that has been circulating since 1994, but to my knowledge
  has never made it into an R package. This version also accepts 'yes-no'
  tables from binary experiments as input.

}
\usage{
pava(y, wt = rep(1, length(y)),dec=FALSE,wt.overwrite=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{The data to be monotonized. Can be a vector or a two-column yes-no
    table (for binary response experiments).}
  \item{wt}{vector, weights paired with the y values. Defaults to
    rep(1,length(y)).}
  \item{dec}{Is the true function monotone decreasing (defaults to FALSE)?}
  \item{wt.overwrite}{Should the variable 'wt' be recalculated as the
    observation counts in each row? Defaults to TRUE. Applicable only for yes-no table input.}
}
\details{
Isotonic regression is a standard nonparametric method to produce point
  estimates for a function known to be monotone (Barlow et al. (1972),
  Robertson et al. (1988)). It replaces each stretch of monotonicity-violating
  observations with the weighted average of the original y values in that
  stretch. This creates the flat (i.e., piecewise-constant) intervals characteristic of IR's
  output.

  If y is a yes-no table, the weights will be set as the observation
  counts in each row - UNLESS 'wt.overwrite' is set to FALSE, in which
  case it will use the given value of 'wt'.

  If there are no monotonicity-violating stretches, the output equals y.

 
}
\value{
  A vector with the same length as y. 
 
}
\references{ Barlow R.E., Bartholomew D.J., Bremner J.M. and Brunk
H.D., Statistical Inference under Order Restriction. John Wiley & Sons
1972.

Robertson T., Wright F.T. and Dykstra R.L. Order Restricted Statistical
Inference, Wiley, Chichester (1988).

}

\author{Richard F. Raubertas (documented and extended by Assaf P. Oron)}

\note{ 1. Even though the term "isotonic" literally means
monotone-increasing, we have included an option to model
monotone-decreasing functions as well, by setting 'dec=TRUE'. }

\section{WARNING }{If you provide y as a yes-no table, do *NOT* set
  'wt.overwrite' to FALSE unless you really want to input different
  weights to the 'wt' variable. If you leave 'wt.overwrite' as is, the
  function will calculate the correct weights, i.e., obesrvation
  counts. Be aware that any weights other than observation counts for
  binary data will
  yield a non-standard solution, so tinker with them only if you know
  what you are doing.} 

\seealso{ A more limited version of IR, written as a recursive function, is available in
\code{\link[stats]{isoreg}}. Compare IR with \code{\link{cir.pava}}, a
simple modification that avoids the flat-stretch output; and with the
sophisticated smoothing of \code{\link[fda]{smooth.monotone}}. For
percentile (inverse) estimation in a dose-response setting, see
\code{\link{cir.upndown}}. }

\examples{ 
### In the 'stackloss' dataset,escape of ammonia through some plant's
### chimney appears driven mostly by plant operation rate with a clearly
### monotone dependence. Linearity is questionable, though, and there
### are monotonicity violations in the data.
### There are are 21 observations at 8 distinct rates, and the original
### dataset is not ordered.
### "pava" and "cir.pava" require unique and ordered x values
### (implicitly in "pava").
### So this example also shows how to prepare such data for input to
### "pava" or"cir.pava" (not difficult):

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
