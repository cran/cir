`cir.pava` <-
function (y,x, wt=rep(1,length(x)),boundary=2,full=FALSE,dec=FALSE,wt.overwrite=TRUE) {
# Returns centered-isotonic-regression y values at original x points #
# Assaf Oron, 10/2007
# 
### ARGUMENTS:
# y: y values (responses). Can be a vector or a yes-no table (for binary responses)
#    it is Given as first argument, both for compatibility with 'pava' 
# and to enable parallel running via 'apply' type routines
# x: treatments. Need to be pre-sorted in increasing order, with order matching y's
# wt: weights. Will be overwritten in case of a yes-no input for y
# boundary: action on boundaries. Defaults to 2, 
# which is analogous to 'rule=2' on function 'approx', i.e.
# returned y values are constant outside the boundaries.
# boundary=1 does linear extrapolation.
# In addition, one can impose boundaries as inputs or
# augment the output with boundaries, as discussed in
# the dissertation text.
# full: if FALSE, only point estimates at x values are returned
# dec: Whether the true function is assumed to be
#    decreasing. Defaults to FALSE
# wt.overwrite: whether, in case of yes-no input, the weights should be recalculated as row
# observation counts

### adapting input in case of yes-no table, and some validation

ll=dim(y)

if (length(ll)>2) stop ("y values can only be a vector or yes-no table.")

if (length(ll)==2 && ll[2]>2) stop ("y values can only be a vector or yes-no table.")

if (length(ll)==2 && ll[2]==2) { # converting a yes-no table
    
    n.u<-y[,1]+y[,2]
    x<-x[n.u>0]
    y<-y[n.u>0,1]/n.u[n.u>0]
    if (wt.overwrite) wt<-n.u[n.u>0] 

}
### More validation stuff
n <- length(x)
if (n <= 1) {
if (!full) return (y)
else return(list(x=x,y=y,z=x))
}
if (any(is.na(x)) || any(is.na(y))) {       
    stop ("Missing values in 'x' or 'y' not allowed")    }
if (any(diff(x)<=0)) {stop ("x must be strictly increasing")}

if (dec) y = -y


z<-x  # Keeping a 'clean' copy of x for final stage

lvlsets <- (1:n)
repeat {
    viol <- (as.vector(diff(y)) <= 0) # Find adjacent violators
    if (!(any(viol))) break
    i <- min( (1:(n-1))[viol]) # Pool first pair of violators
    y[i] <- (y[i]*wt[i]+y[i+1]*wt[i+1]) / (wt[i]+wt[i+1])
    x[i] <- (x[i]*wt[i]+x[i+1]*wt[i+1]) / (wt[i]+wt[i+1])  # new x is calculated
    wt[i]<-wt[i]+wt[i+1]  # weights are combined

# Deleting the i-1-th element
    y<-y[-(i+1)]
    x<-x[-(i+1)]
    wt<-wt[-(i+1)]
    n <- length(y)
    if (n <= 1) break
  }

if (boundary==1) {

### Utilize this option if you wish to use linear extrapolation 
### outside the boundary
### (the 'approx' function does not have this option)
### In general, this is *not* recommended; 
### rather, impose boundary conditions whenever possible
### (as inputs or after output of this function)
### or use the default, constant boundary conditions 

    if (x[n]<max(z)) {
        x<-c(x,max(z))
        y<-c(y,y[n]+(y[n]-y[n-1])*(x[n+1]-x[n])/(x[n]-x[n-1])) }
    if (x[1]>min(z)) {
        x<-c(min(z),x)
        y<-c(y[1]-(y[2]-y[1])*(x[2]-x[1])/(x[3]-x[2]),y) }
}

# Now we re-interpolate to original x values, stored in z
# If we didn't set boundary=1, then this will give constant
# y values for x values falling outside new range of x

if (dec) y = -y

if (!full) return(approx(x,y,z,rule=2)$y)

else return(list(output.y=approx(x,y,z,rule=2)$y,original.x=z,alg.x=x,alg.y=y,alg.wt=wt))
}

