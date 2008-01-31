pava <- function (y, wt=rep(1,length(y)),dec=FALSE,wt.overwrite=TRUE)
#  Compute the isotonic regression of numeric vector 'y', with
#  weights 'wt', with respect to simple order.  The pool-adjacent-
#  violators algorithm is used.  Returns a vector of the same length
#  as 'y' containing the regression.
 
#  02 Sep 1994 / R.F. Raubertas
# Slightly modified, A.P. Oron, Jan. 2008
{
ll=dim(y)

if (length(ll)>2) stop ("y values can only be a vector or yes-no table.")

if (length(ll)==2 && ll[2]>2) stop ("y values can only be a vector or yes-no table.")

if (length(ll)==2 && ll[2]==2) { # converting a yes-no table
    
    n.u<-y[,1]+y[,2]
    y<-y[n.u>0,1]/n.u[n.u>0]
    if (wt.overwrite) wt<-n.u[n.u>0] 
}
   n <- length(y)
   if (n <= 1) return (y)
 
   if (any(is.na(y)) || any(is.na(wt))) {
      stop ("Missing values in 'y' or 'wt' not allowed")
   }
   if (dec) y = -y
   lvlsets <- (1:n)
 
   repeat {
      viol <- (as.vector(diff(y)) < 0)  # Find adjacent violators
      if (!(any(viol))) break
 
      i <- min( (1:(n-1))[viol])        # Pool first pair of violators
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i+1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      y[ilvl] <- sum(y[ilvl]*wt[ilvl]) / sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
   }

   if (dec) y = -y
   y
}
 
