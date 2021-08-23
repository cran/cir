## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#"
)

## ----data---------------------------------------------------------------------
# For brevity, we initially use integers to denote the doses. 
# We make use of R’s shorthand for consecutive sequences, 
# e.g., 1:3 is really 1,2,3
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

## ----cir----------------------------------------------------------------------
ropiTargetCIR=quickInverse(bhamou03ropiRates, target=0.5, adaptiveShrink=TRUE)
ropiTargetCIR 
levoTargetCIR=quickInverse(bhamou03levoRates, target=0.5, adaptiveShrink=TRUE)
levoTargetCIR

## ----ci83---------------------------------------------------------------------
quickInverse(bhamou03ropiRates, target=0.5, adaptiveShrink=TRUE, conf=0.83)
quickInverse(bhamou03levoRates, target=0.5, adaptiveShrink=TRUE, conf=0.83)

## ----pacesty------------------------------------------------------------------
quickInverse(bhamou03ropiRates, target=0.5, estfun=oldPAVA, conf=0.83)
quickInverse(bhamou03levoRates, target=0.5, estfun=oldPAVA, conf=0.83)

## ----forward------------------------------------------------------------------
ropiCurveIR = oldPAVA(bhamou03ropiRates, full=TRUE)
ropiCurveCIR = cirPAVA(bhamou03ropiRates, target=0.5, adaptiveShrink=TRUE, full=TRUE)
levoCurveIR = oldPAVA(bhamou03levoRates, full=TRUE)
levoCurveCIR = cirPAVA(bhamou03levoRates, target=0.5, adaptiveShrink=TRUE, full=TRUE)

## ----forward2-----------------------------------------------------------------
ropiCurveCIR

## ----drfig,fig.width=10,fig.height=5,out.height=400,out.width=800-------------
par(mfrow=c(1,2), las=1, mar=c(4,4,4,1)) # image format parameters
plot(bhamou03ropiRates, xlab="Concentration (%)", 
ylab="Proportion Effective", main='Ropivacaine Arm')
# Adding IR and CIR lines with the same colors/lines as article’s Fig. 4
lines(y~x, data=ropiCurveIR$output, lty=2)
lines(y~x, data=ropiCurveCIR$shrinkage, col='blue',lwd=2)
# Showing the CIR estimate, and confidence interval as a horizontal line
points(target~point, data=ropiTargetCIR, col='purple', pch=19, cex=2)
lines(c(ropiTargetCIR$lower90conf,ropiTargetCIR$upper90conf), rep(0.5,2), col='purple', lwd=2)

# Adding legend:
legend('bottomright', legend=c("Observed Proportions",'Isotonic Regression',
                              'Centered Isotonic Regression','Estimate +/- 90% CI'),
       bty='n',pch=c(4,rep(NA,2),16),col=c('black','black','blue','purple'),lty=c(0,2,1,1))

### Now, second plot for Levobupivacaine
plot(bhamou03levoRates, xlab="Concentration (%)", 
ylab="Proportion Effective", main='Levobupivacaine Arm', ylim=0:1)
lines(y~x, data=levoCurveIR$output, lty=2)
lines(y~x, data=levoCurveCIR$shrinkage, col='blue',lwd=2)
points(target~point, data=levoTargetCIR, col='purple', pch=19, cex=2)
lines(c(levoTargetCIR$lower90conf,levoTargetCIR$upper90conf), rep(0.5,2), col='purple', lwd=2)

