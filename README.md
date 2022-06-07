# pdR
Package of R
    pdR 1.8 is updated and is on CRAN now (install.packages("pdR") works).  In the previous version, pdR contains several econometric routines, for example, panel threshold model and panel unit root (seasonal). Version 1.6 includes the "Hansen B. E. (2000) Sample Splitting and Threshold Estimation. Econometrica, 68, 575-603."  Hansen (2000) is an econometric version of regression tree, but has solid econometric foundation. Version 1.8 modifies Hansen(2000) to make its outputs reusable.

     This version has added three functions:

          SMPLSplit_het() tests the statistical significance of thresholds(sample splitting) by LM statistic and bootstrap p-value and

          SMPLSplit_est() estimates the model. Finally, 

          SMPLSplit_example() offers a learning example.

     I modified Hansen's original code and gives a learning function to use Sample Splitting model. SMPLSplit_example() is a function containing learning codes, which executes Hansen's example of Durlauf-Johnson growth data, but has a detailed explanation of  how sample is classified.

library(pdR)
data("dur_john.rda")
rep <- 500
trim_per <- 0.15
dep <- "gdpGrowth"
indep <- colnames(dur_john)[c(2,3,4,5)]
th1 <- "GDP60"
th2 <- "Literacy"

SMPLSplit_est(data=dur_john,dep,indep,th=th1,plot=1,h=1,nonpar=2) 

I also make Hansen's example as a function, you may also type "SMPLSplit_example" to view function in the colsole.
OUT=SMPLSplit_example(data=dur_john,dep,indep,th1,th2,trim_per,rep,plot=0)
OUT$TEST
OUT$Hypothesis
OUT$Threshold
stat=matrix(as.numeric(OUT$TEST),byrow = TRUE,8,2)
colnames(stat)=c("F-Stat","P-value")
rownames(stat)=OUT$Hypothesis
stat
