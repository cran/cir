## Summary of raw data from the notorious Neuenschwander et al. (Stat. Med., 2008) trial
## Note the use of the 'cohort' argument to specify the cohort order
neundatTrace = DRtrace(x = c(rep(1:4, c(3,4,5,4) ), 7, 7, rep(6,9)),
                       y = c(rep(0,16), 1,1, rep(c(0,0,1),2), 0,0,0), 
                       cohort = rep(1:8, c(3,4,5,4, 2, 3,3,3)) )
par(mar=c(3,3,3,1), mgp=c(2,.5,0), tcl=-0.25)
layout(t(1:2))
plot(neundatTrace ,main="N. et al. (2008) Trajectory", xlab = 'Cohort', 
          ylab="Ordinal Dose Level" ,cex.main=1.5)

## Same data, in 'doseResponse' format with actual doses rather than dose levels
neundatDose = doseResponse(x=c(1,2.5,5,10,20,25), y = c(rep(0,4),2/9,1), wt = c(3,4,5,4,9,2) )
plot(neundatDose ,main="N. et al. (2008) Final Dose-Toxicity", ylim=c(0,1),
	xlab="Dose (mg/sq.m./wk)", ylab="Toxicity Response Curve (F)", cex.main=1.5)
## We can also convert the DRtrace object to doseResponse...
neundatLevel = doseResponse(neundatTrace)

### Now plotting the former, vs. IR/CIR estimates
neunCIR0 = cirPAVA(neundatDose,full=TRUE, adaptiveShrink = TRUE, target = 0.3)
lines(neunCIR0$shrinkage$x, neunCIR0$shrinkage$y)
legend(1,1, pch=c(4,NA), lty = 0:1, legend=c('Observations', 'CIR w/bias corr.'), bty='n')
