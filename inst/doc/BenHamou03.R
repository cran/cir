## ----setup, include = FALSE---------------------------------------------------
#### Note to self: Use this whenever editing the vignette for resubmission/rebuild
#   (after the final CRAN-compatibility check)
# tools::buildVignettes(dir = '.',  tangle=TRUE)
#  (You can also replace the 2 lines below here with point/click directory creation and file copy/paste)
# dir.create("inst/doc")
# file.copy(dir("vignettes", full.names=TRUE), "inst/doc", overwrite=TRUE)
# (from https://community.rstudio.com/t/browsevignettes-mypackage-saying-no-vignettes-found/68656/6 )

options(rmarkdown.html_vignette.check_title = FALSE)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#"
)

## ----data---------------------------------------------------------------------
# For brevity, we initially use integers to denote the doses. 

xropi = c(11:9,10:8,9,10,9,10:7,8:11,10:12,11:7,8,7:10,9,8,9,8:10,9,10,9,10)
xlevo = c(11,10,11,10,11:9,10:7,8,7,8:5,6:8,7,8:6,7,6,7,6,7:5,6,7,6:12)

## ----DRtrace------------------------------------------------------------------
library(cir)
bhamou03ropi = DRtrace(x=xropi[-40]/100, y=(1-diff(xropi))/2)
bhamou03levo = DRtrace(x=xlevo[-40]/100, y=(1-diff(xlevo))/2)

## ----tracefig,fig.width=9,fig.height=4.5,out.height=400,out.width=800---------
par(mfrow=c(1,2), las=1, mar=c(4,4,4,1)) # image format parameters
doserange = c(5,12)/100

plot(bhamou03ropi, ylim=doserange, ylab="Concentration (%)", main='Ropivacaine Arm')
legend('bottomright',legend=c('Effective','Ineffective'),pch=c(19,1),bty='n')
plot(bhamou03levo, ylim=doserange, ylab="Concentration (%)", main='Levobupivacaine Arm')

## ----doseResponse-------------------------------------------------------------
bhamou03ropiRates = doseResponse(bhamou03ropi)
bhamou03levoRates = doseResponse(bhamou03levo)
knitr::kable(bhamou03ropiRates, row.names=FALSE,align='ccr',digits=c(2,4,0))
knitr::kable(bhamou03levoRates, row.names=FALSE,align='ccr',digits=c(2,4,0))

## ----drfig0,fig.width=10,fig.height=5,out.height=400,out.width=800------------
par(mfrow=c(1,2), las=1, mar=c(4,4,4,1)) # image format parameters
plot(bhamou03ropiRates, xlab="Concentration (%)", 
ylab="Proportion Effective", main='Ropivacaine Arm', ylim=0:1)
# Showing the target response rate
abline( h=0.5, col='purple', lty=2)

plot(bhamou03levoRates, xlab="Concentration (%)", 
ylab="Proportion Effective", main='Levobupivacaine Arm', ylim=0:1)

abline( h=0.5, col='purple', lty=2)

## ----fcir---------------------------------------------------------------------
quickIsotone(bhamou03ropiRates, target=0.5, adaptiveShrink=TRUE)
quickIsotone(bhamou03levoRates, target=0.5, adaptiveShrink=TRUE)

## ----cir----------------------------------------------------------------------
ropiTargetCIR = quickInverse(bhamou03ropiRates, target=0.5, adaptiveShrink=TRUE)
ropiTargetCIR 
levoTargetCIR = quickInverse(bhamou03levoRates, target=0.5, adaptiveShrink=TRUE)
levoTargetCIR

## ----ci83---------------------------------------------------------------------
ropi83 = quickInverse(bhamou03ropiRates, target=0.5, adaptiveShrink=TRUE, conf=0.83)
ropi83
levo83 = quickInverse(bhamou03levoRates, target=0.5, adaptiveShrink=TRUE, conf=0.83)
levo83

## ----pacesty------------------------------------------------------------------
PSropiEstimate = quickInverse(bhamou03ropiRates, target=0.5, estfun=oldPAVA, conf=0.83)
PSropiEstimate
PSlevoEstimate =quickInverse(bhamou03levoRates, target=0.5, estfun=oldPAVA, conf=0.83)
PSlevoEstimate

## ----forward------------------------------------------------------------------
ropiCurveIR = oldPAVA(bhamou03ropiRates, full=TRUE)
ropiCurveCIR = cirPAVA(bhamou03ropiRates, target=0.5, adaptiveShrink=TRUE, full=TRUE)
levoCurveIR = oldPAVA(bhamou03levoRates, full=TRUE)
levoCurveCIR = cirPAVA(bhamou03levoRates, target=0.5, adaptiveShrink=TRUE, full=TRUE)

## ----forward2-----------------------------------------------------------------
levoCurveCIR

## ----drfig,fig.width=6,fig.height=12,out.height=1000,out.width=500------------
par(mfrow=c(2,1), las=1, mar=c(4,4,4,1)) # image format parameters
plot(bhamou03ropiRates, xlab="Concentration (%)", 
  ylab="Proportion Effective", main='Ropivacaine Arm', 
  ylim=0:1, xlim=c(0.05, 0.12), dosevals = (5:12)/100)
# Adding IR and CIR lines with the same colors/lines as article’s Fig. 4
lines(y~x, data=ropiCurveIR$output, lty=2)
lines(y~x, data=ropiCurveCIR$shrinkage, col='blue',lwd=2)
# Showing the CIR estimate, and confidence interval as a horizontal line
points(target ~ point, data=ropiTargetCIR, col='purple', pch=19, cex=2)
lines(c(ropi83$lower83conf,ropi83$upper83conf), rep(0.5,2), col='purple', lwd=2)
# The estimate appearing in Pace and Stylianou 2007, nudged 0.01 units down:
points(I(target-0.01) ~ point, data=PSropiEstimate, cex=2)
lines(c(0.087, 0.097), rep(0.49,2))

# Adding legend:
legend('topleft', legend=c("Observed Proportions", 'Isotonic Regression',
                              'Centered Isotonic Regression', paste(c('CIR', 2007), 'Estimate +/- 83% CI')),
       bty='n',pch=c(4,NA,NA,16,1), col=c('black','black','blue','purple', 'black'), lty=c(0,2,1,1,1))

### Now, second plot for Levobupivacaine
plot(bhamou03levoRates, xlab="Concentration (%)", 
  ylab="Proportion Effective", main='Levobupivacaine Arm', 
  ylim=0:1, xlim=c(0.05, 0.12), dosevals = (5:12)/100)

lines(y~x, data=levoCurveIR$output, lty=2)
lines(y~x, data=levoCurveCIR$shrinkage, col='blue',lwd=2)
points(target ~ point, data=levoTargetCIR, col='purple', pch=19, cex=2)
lines(c(levo83$lower83conf,levo83$upper83conf), rep(0.5,2), col='purple', lwd=2)
points(I(target-0.01) ~ point, data=PSlevoEstimate, cex=2)
lines(c(0.059, 0.081), rep(0.49,2))


