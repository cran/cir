smoothqbinom=function(p,size,prob,add=TRUE,half=FALSE) {
q1=qbinom(p,size,prob)
p1=pbinom(q1,size,prob)

### We smooth out the binomial quantile function, because n is random

p2=pbinom(q1-1,size,prob)
out=q1-(p1-p)/(p1-p2)
if(add==TRUE) {
		out[p>0.5]=out[p>0.5]+2
		out[p<0.5]=out[p<0.5]-1

}
if (half==TRUE) out=out-0.5

#### This part to ensure we are still conservative compared with the traditional quantile function

out[p<0.5]=ifelse(out[p<0.5]<q1[p<0.5],out[p<0.5],q1[p<0.5]) 
out[p>0.5]=ifelse(out[p>0.5]>q1[p>0.5],out[p>0.5],q1[p>0.5]) 

out 
}
####

smoothqpois=function(p,size,prob,add=FALSE) { 

refp=c(0.5,ifelse(p<0.5,1-p,p))
q1=qpois(refp,lambda=size*prob)
p1=ppois(q1,lambda=size*prob)
p2=ppois(q1-1,lambda=size*prob)

### We smooth out the Poisson quantile function, because n is random

out=q1-(p1-refp)/(p1-p2)

extra=ifelse(add==TRUE,1,0)

ifelse(p>0.5,out[2:length(out)]+extra,2*out[1]-out[2:length(out)]-extra)+size*prob-out[1] 

}
####
