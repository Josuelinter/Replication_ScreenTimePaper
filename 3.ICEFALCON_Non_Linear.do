**============================================================================
* Paper:       Beyond the Screen: Disentangling Causal Effects and Familial Confounding in Adolescent Screen Time and Cognition
* Authors:     Alex Campbell - alexander.campbell1@student.unimelb.edu.au & Josue Linarte Teran - josue.teran.linarte@hhu.de
* Script:      03_ICEFALCON_Non_Linear.do
* Purpose:     Main ICE FALCON estimation (MZ Twins, N = 1,942 / 971 pairs)
*              Estimates Models I, II, III via GEE, Wright's path tracing,
*              10,000 cluster-bootstrap resamples, and exports Table 2.
* Software:    Stata 17 (or higher)
* Note:        10,000 cluster-bootstrap iterations across 12 models require 
*              substantial computation time (~2-4 hours).
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
global results_dir "$datasets/Results"

// Specific filenames for the output
global results_csv "$datasets\Results\Split_Screentime_ICEFALCON.xlsx"
global output_word "$datasets\Results\Split_Screentime_ICEFALCON.docx"

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

// 4. Define the Sample (Listwise Deletion across ALL domains)
generate mysamp = 1
replace mysamp = 0 if missing(log_Total, log_TV, log_Net, log_Game, log_PcG, log_Computer)
replace mysamp = 0 if missing(Cog_Z_Cohort_w01, age0100_w01, sex_w01, svw0200)

// 5. Apply Complete-Pair Deletion
bysort fid: egen valid_twins = total(mysamp)
replace mysamp = 0 if valid_twins != 2
drop valid_twins
keep if mysamp == 1

// 6. Standardize
capture drop z_*
egen z_Cog_Z_Cohort_w01 = std(Cog_Z_Cohort_w01)
foreach var in log_Total log_TV log_Net log_Game log_PcG log_Computer {
    egen z_`var' = std(`var')
}


// <<<<<< THE 5-HOUR NON-LINEAR SAMPLE SPLIT >>>>>>
// Calculate the family's average raw screen time to preserve twin pairs
bysort fid: egen family_screen = mean(Total_Hours)

capture drop screen_split
generate screen_split = .
replace screen_split = 1 if family_screen <= 5
replace screen_split = 2 if family_screen > 5
drop family_screen

// Save
tempfile _master_standardized
save `_master_standardized'

**----------------------------------------------------------------------------
**# Section 3.0: ICE FALCON ANALYSIS SETUP
**----------------------------------------------------------------------------
mata: all_results = J(0,3, "")
local __i = 1

local screen_domains "log_Total log_TV log_Net log_Game log_PcG log_Computer"
local directions "Forward Reverse"

// ===========================================================================
// BEGIN DOMAIN LOOP
// ===========================================================================
foreach domain in `screen_domains' {

    // =======================================================================
    // BEGIN DIRECTION LOOP
    // =======================================================================
    foreach dir in `directions' {
		
		foreach split_group in 1 2 {
        
        if "`dir'" == "Forward" {
            local x_target "`domain'"
            local y_target "Cog_Z_Cohort_w01"
        }
        else if "`dir'" == "Reverse" {
            local x_target "Cog_Z_Cohort_w01"
            local y_target "`domain'"
        }
        
        local w_target "svw0200" 
		
        use `_master_standardized', clear
        tempfile _tempfile_for_boostrap

            keep if screen_split == `split_group'
            
            if `split_group' == 1 {
                local ses_tertile "`domain' (<= 5 Hrs)"
            }
            else if `split_group' == 2 {
                local ses_tertile "`domain' (> 5 Hrs)"
            }    
			
        local x_analysis "z_`x_target'"
        local y_analysis "z_`y_target'"

	**************************************************************************
	// Step 1: Generating co-family member value
	**************************************************************************
	sort fid ptyp                       
	
	bysort fid: gen `x_analysis'cot = `x_analysis'[_N-_n+1]
	bysort fid: gen `y_analysis'cot = `y_analysis'[_N-_n+1]
	save `_tempfile_for_boostrap'

	**************************************************************************
	// Step 2: Within-pair correlations & Background Stats
	**************************************************************************
	mata: all_results = all_results \  								///			
		st_local("x_target") , st_local("y_target"), st_local("ses_tertile") \ 	///
		st_local("x_target") , st_local("y_target"), st_local("ses_tertile")

	scalar i = `__i'
	scalar n = _N
	quietly: levelsof fid
		scalar pairs_n = r(r) 
	mata: background_stats = st_numscalar("i") , st_numscalar("n") , st_numscalar("pairs_n")

	corr `x_analysis' `x_analysis'cot
		mata: background_stats = background_stats , st_matrix("r(C)")[2,1]
		scalar z_xcorr = 1/2 * (log( ( 1+r(C)[2,1]) / (1-r(C)[2,1]) ) )

	corr `y_analysis' `y_analysis'cot
		mata: background_stats = background_stats , st_matrix("r(C)")[2,1]
		scalar z_ycorr = 1/2 * (log( (1+r(C)[2,1])/(1-r(C)[2,1]) ) )
		
	scalar corr_diff_fisher = (z_xcorr  - z_ycorr) / ( sqrt(1/ (pairs_n - 3)) + sqrt(1/ (pairs_n - 3)) )
	mata:  background_stats = background_stats , st_numscalar("corr_diff_fisher")

	**************************************************************************
	// Step 3: Performing 3 regressions for ICE FALCON
	**************************************************************************
	xtset fid

	* Model I
	xtgee `y_analysis' `x_analysis' age0100_w01 sex_w01 [pw = `w_target'] , corr(exchangeable)
	mata: model1Self = st_matrix("r(table)")[1, 1], st_matrix("r(table)")[5, 1] , st_matrix("r(table)")[6, 1]

	* Model II
	xtgee `y_analysis' `x_analysis'cot age0100_w01 sex_w01 [pw = `w_target'] , corr(exchangeable)
	mata: model2Cotw = st_matrix("r(table)")[1, 1], st_matrix("r(table)")[5, 1] , st_matrix("r(table)")[6, 1]

	* Model III
	xtgee  `y_analysis' `x_analysis' `x_analysis'cot age0100_w01 sex_w01 [pw = `w_target'] , corr(exchangeable)
	mata: model3Self = st_matrix("r(table)")[1, 1], st_matrix("r(table)")[5, 1] , st_matrix("r(table)")[6, 1] 
	mata: model3Cotw = st_matrix("r(table)")[1, 2], st_matrix("r(table)")[5, 2] , st_matrix("r(table)")[6, 2]

	mata: loop_results = ( (background_stats) \ (background_stats) ) , 	///
							(model1Self \ model2Cotw) , 				///
							(model3Self \ model3Cotw) 
							
	//////////////////////////////////////////////////////////////////////////
	**# Sub-program 1: Wright's Path Tracing
	//////////////////////////////////////////////////////////////////////////
	mata: BetaSelf 			= model1Self[1,1]
	mata: BetaSelf_LCI 		= model1Self[1,2] 
	mata: BetaSelf_UCI 		= model1Self[1,3] 
	mata: BetaCotwin 		= model2Cotw[1,1]
	mata: ChangeBetaSelf 	= model3Self[1,1] - model1Self[1,1]
	mata: ChangeBetaCotwin 	= model3Cotw[1,1] - model2Cotw[1,1]
	mata: CorrelationX 		= loop_results[1,4]

	mata: isValid = ((BetaSelf_LCI * BetaSelf_UCI > 0) & abs(BetaSelf) > 0.001 & abs(CorrelationX) > 0.001)

	mata: Pr 					= isValid ? (((ChangeBetaCotwin - (ChangeBetaSelf/BetaSelf) * BetaCotwin) / CorrelationX ) / BetaSelf) : .
	mata: CausalEffect 			= isValid ? (BetaSelf * Pr) : .
	mata: FamilialConfounding 	= isValid ? (1 - Pr) : .

	//////////////////////////////////////////////////////////////////////////
	**# Sub-program 2: Bootstrap re-sampling
	//////////////////////////////////////////////////////////////////////////
	set seed 1307
	local NumberOfBootstrapResamples = 10000

	forvalue __bootstrap_i = 1/`NumberOfBootstrapResamples' {
		quietly {
			use `_tempfile_for_boostrap', clear

			quietly: levelsof fid
			local number_of_pairs = r(r)
			bsample `number_of_pairs', cluster(fid) idcluster(fid2)

			xtset fid2

			xtgee `y_analysis' `x_analysis' age0100_w01 sex_w01 [pw = `w_target'] , corr(exchangeable)
			mata: model1SelfBoot = st_matrix("r(table)")[1, 1]

			xtgee `y_analysis' `x_analysis'cot age0100_w01 sex_w01 [pw = `w_target'] , corr(exchangeable)
			mata: model2CotwBoot = st_matrix("r(table)")[1, 1]

			xtgee  `y_analysis' `x_analysis' `x_analysis'cot age0100_w01 sex_w01 [pw = `w_target'] , corr(exchangeable)
			mata: model3SelfBoot = st_matrix("r(table)")[1, 1]
			mata: model3CotwBoot = st_matrix("r(table)")[1, 2]

			if `__bootstrap_i' == 1 {
				mata: total_boots = strtoreal(st_local("NumberOfBootstrapResamples"))
				mata: model1Self_bootstrap = J(total_boots, 1, .)
				mata: model2Cotw_bootstrap = J(total_boots, 1, .)
				mata: model3Self_bootstrap = J(total_boots, 1, .)
				mata: model3Cotw_bootstrap = J(total_boots, 1, .)
			}

			mata: current_row = strtoreal(st_local("__bootstrap_i"))
			mata: model1Self_bootstrap[current_row, 1] = model1SelfBoot
			mata: model2Cotw_bootstrap[current_row, 1] = model2CotwBoot
			mata: model3Self_bootstrap[current_row, 1] = model3SelfBoot
			mata: model3Cotw_bootstrap[current_row, 1] = model3CotwBoot

			mata: mata drop model1SelfBoot model2CotwBoot model3SelfBoot model3CotwBoot

			if mod(`__bootstrap_i', 500) == 0 {
				display "Completed bootstrap `__bootstrap_i' of `NumberOfBootstrapResamples'..."
			}
		}
	}
		
	// Bootstrap Math
	mata: self_boot_diffs = model3Self_bootstrap - model1Self_bootstrap
	mata: cotw_boot_diffs = model3Cotw_bootstrap - model2Cotw_bootstrap 

	mata: self_boot_mean_diff 	= mean(self_boot_diffs)
	mata: self_boot_se 			= sqrt(variance(self_boot_diffs)) 
	mata: self_boot_lci 		= self_boot_mean_diff - 1.96 * self_boot_se
	mata: self_boot_uci 		= self_boot_mean_diff + 1.96 * self_boot_se
	mata: self_boot_zval 		= self_boot_mean_diff / self_boot_se
	mata: self_boot_p 			= 2 * normal(-abs(self_boot_zval)) 

	mata: cotw_boot_mean_diff 	= mean(cotw_boot_diffs)
	mata: cotw_boot_se 			= sqrt(variance(cotw_boot_diffs)) 
	mata: cotw_boot_lci 		= cotw_boot_mean_diff - 1.96 * cotw_boot_se
	mata: cotw_boot_uci 		= cotw_boot_mean_diff + 1.96 * cotw_boot_se
	mata: cotw_boot_zval 		= cotw_boot_mean_diff / cotw_boot_se
	mata: cotw_boot_p 			= 2 * normal(-abs(cotw_boot_zval)) 

	mata: loop_results = ( loop_results , 									///
			( 	(self_boot_mean_diff , self_boot_lci , self_boot_uci, self_boot_p) \		/// 
				(cotw_boot_mean_diff , cotw_boot_lci , cotw_boot_uci, cotw_boot_p) ) ) 	

	mata: mata drop 									///
		self_boot_mean_diff self_boot_lci self_boot_uci ///
		self_boot_se self_boot_zval self_boot_p 		///
		cotw_boot_mean_diff cotw_boot_lci cotw_boot_uci ///
		cotw_boot_se cotw_boot_zval cotw_boot_p

	mata: loop_results = (loop_results , (				///
						 (Pr , FamilialConfounding ) \ 	///
						 (Pr , FamilialConfounding ) 	///
						 ) )
						 
	if `__i' == 1 {
		mata: all_loops = loop_results
	}
	if `__i' > 1 {
		mata: all_loops = all_loops \ loop_results
	}
	local ++__i			

} 
}
}
// ===========================================================================
// END DIRECTION LOOP
// ===========================================================================
**----------------------------------------------------------------------------
**# Section 4.0: Generating Results Dataset & Exporting
**----------------------------------------------------------------------------
frame change default
capture: frame drop results_frame 
frame create results_frame
frame change results_frame 
getmata (x_var y_var SES) = all_results

getmata (loop_num n_indiv n_pairs x_corr y_corr corr_diff_fisher		///
		 beta_m12 beta_lci_m12 beta_uci_m12 beta_m3 					///
		 beta_lci_m3 beta_uci_m3 										///
		 boot_mean_diff boot_lci boot_uci boot_p						///
		 causal_pr familial_pr 											///
		 ) 																///
		 = all_loops 

generate direction = ""
replace direction = x_var + " (Forward)" if y_var == "Cog_Z_Cohort_w01"
replace direction = y_var + " (Reverse)" if x_var == "Cog_Z_Cohort_w01"

label var direction "Causal Direction"
label var n_indiv		 "Individuals N"
label var n_pairs		 "Pairs N"
label var x_corr		 "X within-pair correlation"
label var y_corr		 "Y within-pair correlation"
label var beta_m12		 "Beta models 1 & 2"
label var beta_lci_m12	 "LCI models 1 & 2"
label var beta_uci_m12	 "UCI models 1 & 2"
label var beta_m3		 "Beta model 3"
label var beta_lci_m3	 "LCI model 3"
label var beta_uci_m3	 "UCI model 3"
label var boot_mean_diff "Bootstrap change"
label var boot_lci 		 "LCI bootstrap change"
label var boot_uci 		 "UCI bootstrap change"
label var boot_p		 "Bootstrap p-value"
label var causal_pr		 "Causal effect (%)"
label var familial_pr	 "Familial confounding (%)"

bysort loop_num: generate model = "Self" if _n == 1, after(loop_num)
bysort loop_num: replace model = "Co-twin" if _n == 2

export excel using "$results_csv", sheet("ICEFALCON_Results", replace) firstrow(variables)

// Word Formatting
local decimal_places = 2
foreach model in 12 3 {
	generate Model`model' = string(beta_m`model', "%4.`decimal_places'f") + " (" + string(beta_lci_m`model', "%4.`decimal_places'f") + ", " + string(beta_uci_m`model', "%4.`decimal_places'f") + ")"
	drop beta_m`model' beta_lci_m`model' beta_uci_m`model' 
}
generate bootstrap = string(boot_mean_diff, "%4.`decimal_places'f") + " (" + string(boot_lci, "%4.`decimal_places'f") + ", " + string(boot_uci, "%4.`decimal_places'f") + "), p=" + string(boot_p, "%4.3f")

drop boot_mean_diff boot_lci boot_uci

foreach var of varlist x_corr y_corr corr_diff_fisher causal_pr familial_pr {
	format `var' %3.`decimal_places'f
}

// Output Document
putdocx clear 
putdocx begin, pagesize(A4) margin(all, 1cm) landscape

putdocx table mytable = data(direction x_var y_var model n_pairs n_indiv x_corr y_corr corr_diff_fisher Model12 Model3 bootstrap causal_pr familial_pr), 	///
		varnames 													///
		border(insideV, nil) border(insideH, nil) 					///
		border(left, nil) border(right, nil) 						///
	    title("Table 2. ICE FALCON Forward & Reverse Results", bold )			 		///
		layout(autofitc)
		
putdocx table mytable(.,.), font(Times New Roman, 8, black)
putdocx table mytable(1,.), border(bottom, single)
putdocx table mytable(2,.), border(bottom, single)

putdocx save "$output_word", replace

**********************************************
***** END OF ICE FALCON ANALYSIS DO FILE *****
**********************************************