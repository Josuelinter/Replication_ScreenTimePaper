**============================================================================
* Paper:       Beyond the Screen: Disentangling Causal Effects and Familial Confounding in Adolescent Screen Time and Cognition
* Authors:     Alex Campbell - alexander.campbell1@student.unimelb.edu.au & Josue Linarte Teran - josue.teran.linarte@hhu.de
* Script:      07_Knots.do
* Purpose:     Knot Composition Analysis
* Software:    Stata 17 (or higher)
**============================================================================

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

local min_knot = 1
local max_knot = 8

local num_models = `max_knot' - `min_knot' + 1
matrix bic_scores = J(`num_models', 3, .)
matrix colnames bic_scores = knots BIC AIC

local i = 1
forvalues k = `min_knot'/`max_knot' {
	
	capture drop spl*
	
	local n_pieces = `k' + 1
	
	mkspline spl `n_pieces' = Total_Hours
	
	quietly regress Cog_Z_Cohort_w01 spl* age0100_w01 sex_w01
	
	estat ic
	
	matrix temp_S = r(S)
	matrix bic_scores[`i', 1] = `k'
	matrix bic_scores[`i', 2] = temp_S[1, 6]  // BIC
	matrix bic_scores[`i', 3] = temp_S[1, 5]  // AIC
	
	local i = `i' + 1
}

matlist bic_scores

summarize Total_Hours, detail

capture drop spl5_*
mkspline spl5_1 5 spl5_2 = Total_Hours
regress Cog_Z_Cohort_w01 spl5_* age0100_w01 sex_w01
estat ic