**============================================================================
* Paper:       Beyond the Screen: Disentangling Causal Effects and Familial Confounding in Adolescent Screen Time and Cognition
* Authors:     Alex Campbell - alexander.campbell1@student.unimelb.edu.au & Josue Linarte Teran - josue.teran.linarte@hhu.de
* Script:      07_ICEFALCON_OLS.do
* Purpose:     Baseline OLS

**----------------------------------------------------------------------------
**# Section 1.0: Setup & Globals
**----------------------------------------------------------------------------
version 17
frames reset
macro drop _all
set linesize 80

// Relative root path: works automatically wherever the project folder is placed
global route "`c(pwd)'"
cd "$route"

global datasets    "$route/Stata"
global master_data "$datasets/ScreentimeCognition_Master_vx.dta"
global results_dir "$datasets/Results"

use "$master_data", clear


// 1. Restrict to Monozygotic Twins only
keep if zyg0102 == 1

// 2. Generate the Individual Dimensions (Raw Hours)
capture drop TV_Hours Internet_Hours Gaming_Hours Pc_Games_Hours Computer_Hours Total_Hours
generate TV_Hours       = (med0900_w01 * 5 + med1000_w01 * 2) / 7
generate Internet_Hours = (med0300_w01 * 5 + med0400_w01 * 2) / 7
generate Gaming_Hours   = (med0700_w01 * 5 + med0800_w01 * 2) / 7
generate Pc_Games_Hours = (med0500_w01 * 5 + med0600_w01 * 2) / 7
generate Computer_Hours = (med0100_w01 * 5 + med0200_w01 * 2) / 7
generate Total_Hours    = medTotal_w01

// 3. Log-Transform all these new variables
capture drop log_TV log_Net log_Game log_PcG log_Computer log_Total
generate log_TV       = log(TV_Hours + 1)
generate log_Net      = log(Internet_Hours + 1)
generate log_Game     = log(Gaming_Hours + 1)
generate log_PcG      = log(Pc_Games_Hours + 1)
generate log_Computer = log(Computer_Hours + 1)
generate log_Total    = log_medTotal_w01

// 4. Define the Sample 
generate mysamp = 1
replace mysamp = 0 if missing(log_Total, log_TV, log_Net, log_Game, log_PcG, log_Computer)
replace mysamp = 0 if missing(Cog_Z_Cohort_w01, age0100_w01, sex_w01, svw0100)

// 6. Standardize 
capture drop z_*
egen z_Cog_Z_Cohort_w01 = std(Cog_Z_Cohort_w01)
foreach var in log_Total log_TV log_Net log_Game log_PcG log_Computer {
    egen z_`var' = std(`var')
}

// Save
tempfile _master_standardized
save `_master_standardized'


// Simple OLS with Clustered Standard Errors
reg Cog_Z_Cohort_w01 c.z_log_Total c.age0100_w01 i.sex_w01, vce(cluster fid)
reg z_log_Total c.Cog_Z_Cohort_w01 c.age0100_w01 i.sex_w01, vce(cluster fid)

reg Cog_Z_Cohort_w01 c.z_log_Total c.age0100_w01 i.sex_w01 if mysamp == 1, vce(cluster fid)
reg z_log_Total c.Cog_Z_Cohort_w01 c.age0100_w01 i.sex_w01 if mysamp == 1, vce(cluster fid)

reg Cog_Z_Cohort_w01 c.medTotal_w01 if mysamp==1, vce(cluster fid)



