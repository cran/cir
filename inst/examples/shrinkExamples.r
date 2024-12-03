## Summary of raw data from the notorious Neuenschwander et al. (Stat. Med., 2008) trial
## Note the use of the 'cohort' argument to specify the cohort order
neundatDose = doseResponse(x=c(1,2.5,5,10,20,25), y = c(rep(0,4),2/9,1), wt = c(3,4,5,4,9,2) )

neundatDose

# Compare to this:
DRshrink(neundatDose, target = 0.3)