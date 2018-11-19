# pdR
Package of R
    pdR 1.6 is updated and is on CRAN now (install.packages("pdR") works).  In the previous version, pdR contains several econometric routines, for example, panel threshold model and panel unit root (seasonal). Version 1.6 includes the "Hansen B. E. (2000) Sample Splitting and Threshold Estimation. Econometrica, 68, 575-603."  Hansen (2000) is an econometric version of regression tree, but has solid econometric foundation.     

     This version has added three functions:

          SMPLSplit_het() tests the statistical significance of thresholds(sample splitting) by LM statistic and bootstrap p-value and

          SMPLSplit_est() estimates the model. Finally, 

          SMPLSplit_example() offers a learning example.



     I modified Hansen's original code and gives a learning function to use Sample Splitting model. SMPLSplit_example() is a function containing learning codes, which executes Hansen's example of Durlauf-Johnson growth data, but has a detailed explanation of  how sample is classified.



Step 1. Prepare the data and declare options

library(pdR)

data("dur_john.rda")

rep <- 200

trim_per <- 0.15

dep=1  # Column number of  dependent  variables

indep=c(2,3,4,5) # Column number of regime-dependent independent  variables

th1=6  #column of the first threshold variable

th2=7  #column of the second threshold variable


Step 2. We then execute the function SMPLSplit_example() below

OUT=SMPLSplit_example(data=dur_john, dep, indep, th1, th2, trim_per, rep, plot=1)

     I gave a output function to collect important information, which can be show below

OUT$TEST

OUT$Hypothesis

OUT$Threshold

stat=matrix(as.numeric(OUT$TEST), byrow = TRUE,8,2)

colnames(stat)=c("F-Stat","P-value")

rownames(stat)=OUT$Hypothesis

stat

In this example, it tests 8 hypothesis and estimates two models, the test can be summarized below:

> stat

                                                     F-Stat           P-value

Testing for a First Sample Split, Using GDP60        12.601835   0.085

Testing for a First Sample Split, Using Literacy       10.786273   0.186

Testing for a Second Sample Split, Using GDP60    11.009342   0.161

Testing for a Second Sample Split, Using Literacy  12.091351   0.073

Testing for a Third Sample Split, Using GDP60        7.423359   0.558

Testing for a Third Sample Split, Using Literacy       9.233321   0.186

Testing for a Third Sample Split, Using GDP60        9.415477   0.223

Testing for a Third Sample Split, Using Literacy       9.327215   0.176



     Moreover, the computational details are listed below, which details the implication of sample splitting; for individual need, users can retrieve the source code of SMPLSplit_example() to modify. This code has three levels and I add a conclusion note right below each stage in blue. 



Level 1: Testing for a First Sample Split



<1-1> Testing for a First Sample Split, Using GDP60

Test of Null of No Threshold Against Alternative of Threshold

Allowing Heteroskedastic Errors (White Corrected)

Number of Bootstrap Replications  1000

Trimming Percentage               0.15

Threshold Estimate                833

LM-test for no threshold          12.60184

Bootstrap P-Value                 0.085



<1-2> Testing for a First Sample Split, Using Literacy

Test of Null of No Threshold Against Alternative of Threshold

Allowing Heteroskedastic Errors (White Corrected)

Number of Bootstrap Replications  1000

Trimming Percentage               0.15

Threshold Estimate                10

LM-test for no threshold          10.78627

Bootstrap P-Value                 0.186



=== Estimate First Sample Split, Using GDP60 as Threshold ===



1. Global OLS Estimation

Dependent Variable:      gdpGrowth

Heteroskedasticity Correction Used

Variable       Estimate        St Error

----------------------------------------

Constant       2.8262496      0.72472989

logGDP60      -0.2815736      0.05335405

Inv_GDP       0.4918965      0.10434991

popGrowth      -0.5467001      0.23277830

School       0.2395413      0.06560082



Observations:                       96

Degrees of Freedom:                 91

Sum of Squared Errors:              9.622743

Residual Variance:                  0.1057444

R-squared:                          0.4805041

Heteroskedasticity Test (P-Value):  0.2008698

****************************************************



2. Threshold Estimation

Threshold Variable:                 GDP60

Threshold Estimate:                 863

0.95 Confidence Interval:           [594.0000, 1794.0000]

Sum of Squared Errors:              8.024881

Residual Variance:                  0.09331257

Joint R-squared:                    0.5667667

Heteroskedasticity Test (P-Value):  0.1565995

****************************************************



     Regime 1: GDP60<=863.000000



Parameter Estimates

Variable       Estimate        St Error

----------------------------------------

Constant       4.3120283      1.62679940

logGDP60      -0.6569710      0.21761579

Inv_GDP       0.2277417      0.07160391

popGrowth      -0.2948695      0.33677597

School       0.0180607      0.09685598

----------------------------------------

0.95 Confidence Regions for Parameters

Variable       Low               High

----------------------------------------

Constant       0.68755168       9.5624430

logGDP60      -1.25006934      -0.1464676

Inv_GDP       0.02470878       0.5740113

popGrowth      -1.51316271       0.9224831

School      -0.24700520       0.4397419



Observations:                       18

Degrees of Freedom:                 13

Sum of Squared Errors:              0.6742722

Residual Variance:                  0.05186709

R-squared:                          0.5165107

****************************************************



     Regime 2: GDP60>863.000000



Parameter Estimates

Variable       Estimate        St Error

----------------------------------------

Constant       3.6630685      0.71904747

logGDP60      -0.3233915      0.06144147

Inv_GDP       0.4957500      0.14497428

popGrowth      -0.4876940      0.25532245

School       0.3569407      0.08996972

----------------------------------------

0.95 Confidence Regions for Parameters

Variable       Low               High

----------------------------------------

Constant       1.84478944       5.79544463

logGDP60      -0.52300712      -0.18203224

Inv_GDP       0.18229032       0.95436145

popGrowth      -1.06852144       0.03368612

School      -0.08479933       0.54918938



Observations:                       78

Degrees of Freedom:                 73

Sum of Squared Errors:              7.350609

Residual Variance:                  0.1006933

R-squared:                          0.5492007





####################################################

####################################################

We check output above to determine which way to go.

Because regime ' GDP60 > 863 ' has more obs, Level 2 continues

#####################################################

#####################################################



Level 2

Sub-Sample GDP60  > 863



<2-1> Testing for a Second Sample Split, Using GDP60



Test of Null of No Threshold Against Alternative of Threshold

Allowing Heteroskedastic Errors (White Corrected)

Number of Bootstrap Replications  1000

Trimming Percentage               0.15

Threshold Estimate                1410

LM-test for no threshold          11.00934

Bootstrap P-Value                 0.161





<2-2> Testing for a Second Sample Split, Using Literacy



Test of Null of No Threshold Against Alternative of Threshold

Allowing Heteroskedastic Errors (White Corrected)

Number of Bootstrap Replications  1000

Trimming Percentage               0.15

Threshold Estimate                57

LM-test for no threshold          12.09135

Bootstrap P-Value                 0.073







=== Estimate Second Sample Split, Using Literacy as Threshold ===



1. Global OLS Estimation

Dependent Variable:      gdpGrowth

Heteroskedasticity Correction Used

Variable       Estimate        St Error

----------------------------------------

Constant       3.6630685      0.71904747

logGDP60      -0.3233915      0.06144147

Inv_GDP       0.4957500      0.14497428

popGrowth      -0.4876940      0.25532245

School       0.3569407      0.08996972



Observations:                       78

Degrees of Freedom:                 73

Sum of Squared Errors:              7.350609

Residual Variance:                  0.1006933

R-squared:                          0.5492007

Heteroskedasticity Test (P-Value):  0.362994

****************************************************



2. Threshold Estimation

Threshold Variable:                 Literacy

Threshold Estimate:                 45

0.95 Confidence Interval:           [19.0000, 57.0000]

Sum of Squared Errors:              6.198249

Residual Variance:                  0.09115072

Joint R-squared:                    0.6198728

Heteroskedasticity Test (P-Value):  0.551557

****************************************************



     Regime 1: Literacy<=45.000000



Parameter Estimates

Variable       Estimate        St Error

----------------------------------------

Constant       2.0922814      1.8701913

logGDP60      -0.1163168      0.1648712

Inv_GDP       0.1726423      0.2161463

popGrowth      -0.3902176      0.5161286

School       0.4525333      0.1168289

----------------------------------------

0.95 Confidence Regions for Parameters

Variable       Low               High

----------------------------------------

Constant      -1.8538485      6.1610409

logGDP60      -0.4708739      0.3447183

Inv_GDP      -0.3181809      0.6554835

popGrowth      -1.5932356      0.8683475

School       0.1953806      0.7459717



Observations:                       30

Degrees of Freedom:                 25

Sum of Squared Errors:              2.609994

Residual Variance:                  0.1043998

R-squared:                          0.5780366

****************************************************



     Regime 2: Literacy>45.000000



Parameter Estimates

Variable       Estimate        St Error

----------------------------------------

Constant       4.31048461      0.96515976

logGDP60      -0.39502686      0.06102754

Inv_GDP       0.83363986      0.13936990

popGrowth      -0.41800939      0.26956672

School       0.09457796      0.13489238

----------------------------------------

0.95 Confidence Regions for Parameters

Variable       Low               High

----------------------------------------

Constant       1.5818922       6.2291568

logGDP60      -0.5354340      -0.2495195

Inv_GDP       0.4271144       1.1318791

popGrowth      -1.1133443       0.1145814

School      -0.3443145       0.4020814



Observations:                       48

Degrees of Freedom:                 43

Sum of Squared Errors:              3.588255

Residual Variance:                  0.08344778

R-squared:                          0.5821175





####################################################

####################################################

We check outputs above to determine which way to go.

Because both sub-regimes by Literacy have similar obs, Level 3 continues

############################################################

############################################################



Level 3

3A: Given GDP60 > 863 , Sub-Sample Literacy <= 45



<3A-1> Testing for a Third Sample Split, Using GDP60

Test of Null of No Threshold Against Alternative of Threshold

Allowing Heteroskedastic Errors (White Corrected)

Number of Bootstrap Replications  1000

Trimming Percentage               0.15

Threshold Estimate                1618

LM-test for no threshold          7.423359

Bootstrap P-Value                 0.558





<3A-2> Testing for a Third Sample Split, Using Literacy

Test of Null of No Threshold Against Alternative of Threshold

Allowing Heteroskedastic Errors (White Corrected)

Number of Bootstrap Replications  1000

Trimming Percentage               0.15

Threshold Estimate                26

LM-test for no threshold          9.233321

Bootstrap P-Value                 0.186





3B: Given GDP60 > 863 , Sub-Sample Literacy > 45



<3B-1> Testing for a Third Sample Split, Using GDP60

Test of Null of No Threshold Against Alternative of Threshold

Allowing Heteroskedastic Errors (White Corrected)

Number of Bootstrap Replications  1000

Trimming Percentage               0.15

Threshold Estimate                3493

LM-test for no threshold          9.415477

Bootstrap P-Value                 0.223





<3B-2> Testing for a Third Sample Split, Using Literacy

Test of Null of No Threshold Against Alternative of Threshold

Allowing Heteroskedastic Errors (White Corrected)

Number of Bootstrap Replications  1000

Trimming Percentage               0.15

Threshold Estimate                83

LM-test for no threshold          9.327215

Bootstrap P-Value                 0.176



####################################################

####################################################

Because, by the Bootstrap P-Value, the LM tests of both sub-regimes do not have

statistical significance, we then stop here without estimating the third level estimation.

#####################################################################

#####################################################################
