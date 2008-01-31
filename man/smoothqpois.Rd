\name{smoothqpois}
\alias{smoothqpois}
\alias{smoothqbinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Linearly Smoothed Quantile Functions}
\description{
  Provide linearly-interpolated and somewhat modified quantile functions
  for the  binomial and Poisson (as an approximation to the binomial) distributions, for use in interval
  estimation via 'cir.upndown'. All functions are used to approximate the binomial.
}
\usage{
smoothqpois(p, size, prob, add = FALSE)
smoothqbinom(p, size, prob, add = TRUE,half=FALSE)
}
\arguments{
  \item{p}{ The inteval-cutoff quantile points (between 0 and 1). Could be a vector. }
  \item{size}{ Sample size (for the original binomial problem) }
  \item{prob}{ True probabilty (center of the interval)}
  \item{add}{ Whether to add a 'padding' to the output in order to
    account for randomness in sample size (defaults to FALSE, since it
    is already accounted for in using Poisson)}
  \item{half}{(only for 'smoothqbinom') Should the final outcome shifted
    by half a unit, for symmetry? Defaults to FALSE}
}


\details{
These functions are utilities for 'cir.upndown'. They are tailored to help provide realistic interval
estimation for percentile-finding binary experiments in which the
treatment allocation is *random*. They make little sense outside this
context. See example for how the functions look. For theoretical
details, see Oron (2007), Section 3.3.
}
\value{
A vector of the same length as p
}
\references{   Oron A.P., Up-and-Down and the Percentile-Finding Problem. Doctoral
  Dissertation, University of Washington. 2007
}
  \author{Assaf P. Oron}

\seealso{ \code{\link{cir.pava}}, \code{\link{cir.upndown}} }
\examples{

pvec=seq(0,1,.01)

n=10
##### This is for a median-targeting application
targ=0.5

### Plain vanilla binomial (does not account for random allocation, and
### is asymmetric in location of vertical jumps)

plot(100*pvec,qbinom(pvec,size=n,prob=targ),xlab="Cutoff point for
Confidence Interval (percent)",ylab="Confidence Quantile Estimate",type='l')

### Binomial, but linearly interpolated and 'padded'
lines(100*pvec,smoothqbinom(pvec,size=n,prob=targ),col=4)

### Poisson (the preferred approach since it has heavy tails)
lines(100*pvec,smoothqpois(pvec,size=n,prob=targ),col=2)

########## Now the whole thing when the target is 20th percentile

targ=0.2
plot(100*pvec,qbinom(pvec,size=n,prob=targ),xlab="Cutoff point for
Confidence Interval (percent)",ylab="Confidence Quantile Estimate",type='l',ylim=c(-2,10))

lines(100*pvec,smoothqbinom(pvec,size=n,prob=targ),col=4)
lines(100*pvec,smoothqpois(pvec,size=n,prob=targ),col=2)

### A negative value here may seem counter-intuitive, but it makes sense
###in a random-allocation setting
  }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution }
