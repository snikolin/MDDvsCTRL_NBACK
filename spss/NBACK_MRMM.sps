* Encoding: UTF-8.

FILTER BY outlier.

*** Study outcomes ***
* d-prime
* response time
* P2 amplitude
* P3 amplitude
* theta

*** Check MRMM covariance matrix ***
* replace COVTYPE(UN) with covariances below, one at a time, to determine matrix with best model fit *
* UN, UNR, VC, TPH, TP, ID, HF, FAH1, FA1, DIAG, CSR, CSH, CS, ARMA11, ARH1, AR1, AD1 *
* run for each outcome measure (d-prime, response time, etc...) *

MIXED dprime BY group load
  /FIXED = group load group*load | SSTYPE(3)
  /METHOD = ML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= load | SUBJECT(PID) COVTYPE(UN) 
  /EMMEANS = TABLES(group*load) COMPARE(group) ADJ(LSD).


*** Check MRMM covariates (hierarchical model) ***
* use information criterion to see whether covariate improves model fit *
* run for each outcome measure (d-prime, response time, etc...) *

MIXED dprime BY group load
  /FIXED = group load group*load | SSTYPE(3)
  /METHOD = ML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= load | SUBJECT(PID) COVTYPE(AD1) 
  /EMMEANS = TABLES(group*load) COMPARE(group) ADJ(LSD).

MIXED dprime BY group load sex
  /FIXED = group load group*load sex | SSTYPE(3)
  /METHOD = ML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= load | SUBJECT(PID) COVTYPE(AD1) 
  /EMMEANS = TABLES(group*load) COMPARE(group) ADJ(LSD).

MIXED dprime BY group load WITH age
  /FIXED = group load group*load age | SSTYPE(3)
  /METHOD = ML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= load | SUBJECT(PID) COVTYPE(AD1) 
  /EMMEANS = TABLES(group*load) COMPARE(group) ADJ(LSD).

MIXED dprime BY group load WITH edu
  /FIXED = group load group*load edu | SSTYPE(3)
  /METHOD = ML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= load | SUBJECT(PID) COVTYPE(AD1) 
  /EMMEANS = TABLES(group*load) COMPARE(group) ADJ(LSD).


*** Final MRMMs ***

MIXED dprime BY group load
  /FIXED = group load group*load | SSTYPE(3)
  /METHOD = REML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= load | SUBJECT(PID) COVTYPE(AD1) 
  /EMMEANS = TABLES(group*load) COMPARE(group) ADJ(LSD)
  /SAVE RESID (RESID_dprime).

MIXED rt BY group load
  /FIXED = group load group*load | SSTYPE(3)
  /METHOD = REML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= load | SUBJECT(PID) COVTYPE(ARH1) 
  /EMMEANS = TABLES(group*load) COMPARE(group) ADJ(LSD)
  /SAVE RESID (RESID_rt).

MIXED P2_amp BY group load
  /FIXED = group load group*load | SSTYPE(3)
  /METHOD = REML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= load | SUBJECT(PID) COVTYPE(AR1) 
  /EMMEANS = TABLES(group*load) COMPARE(group) ADJ(LSD)
  /SAVE RESID (RESID_P2).

MIXED P3_amp BY group load
  /FIXED = group load group*load | SSTYPE(3)
  /METHOD = REML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= load | SUBJECT(PID) COVTYPE(ARH1) 
  /EMMEANS = TABLES(group*load) COMPARE(group) ADJ(LSD)
  /SAVE RESID (RESID_P3).

MIXED theta BY group load
  /FIXED = group load group*load | SSTYPE(3)
  /METHOD = REML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= load | SUBJECT(PID) COVTYPE(CS) 
  /EMMEANS = TABLES(group*load) COMPARE(group) ADJ(LSD)
  /SAVE RESID (RESID_theta).


*** Examine MRMM residuals ***

EXAMINE VARIABLES = RESID_dprime BY group
  /ID=PID
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES = RESID_rt BY group
  /ID=PID
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES = RESID_P2 BY group
  /ID=PID
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES = RESID_P3 BY group
  /ID=PID
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES = RESID_theta BY group
  /ID=PID
  /PLOT BOXPLOT HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.


***** Check Normality of Residuals *****

EXAMINE VARIABLES= RESID_dprime RESID_rt RESID_P2 RESID_P3 RESID_theta BY group
/PLOT BOXPLOT NPPLOT
/COMPARE GROUPS   
/STATISTICS DESCRIPTIVES
/CINTERVAL 95
/MISSING PAIRWISE
/NOTOTAL.


***** USE WIDEFORM - COMPARE BY NBACK LEVEL *****

***** Check Normality of Outcome Measures per n-Back Level *****

EXAMINE VARIABLES= dprime_0 rt_0 dprime_1 rt_1 dprime_2 rt_2 
p2_amp_0back_group_all p2_amp_1back_group_all p2_amp_2back_group_all 
p3_amp_0back_group_all p3_amp_1back_group_all p3_amp_2back_group_all
theta_0back theta_1back theta_2back  BY group
/PLOT BOXPLOT NPPLOT HISTOGRAM
/COMPARE GROUPS   
/STATISTICS DESCRIPTIVES
/CINTERVAL 95
/MISSING PAIRWISE
/NOTOTAL.


* Nonparametric Tests: Independent Samples *

NPAR TESTS
  /M-W= dprime_0 dprime_1 dprime_2 BY group(0 1)
  /MISSING ANALYSIS.

NPAR TESTS
  /M-W= rt_0 rt_1 rt_2 BY group(0 1)
  /MISSING ANALYSIS.

NPAR TESTS
  /M-W= p2_amp_0back_group_all p2_amp_1back_group_all p2_amp_2back_group_all BY group(0 1)
  /MISSING ANALYSIS.

NPAR TESTS
  /M-W= p3_amp_0back_group_all p3_amp_1back_group_all p3_amp_2back_group_all BY group(0 1)
  /MISSING ANALYSIS.

NPAR TESTS
  /M-W= theta_0back theta_1back theta_2back BY group(0 1)
  /MISSING ANALYSIS.
