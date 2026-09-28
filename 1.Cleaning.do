**============================================================================
* Paper: Beyond the Screen: Disentangling Causal Effects and Familial Confounding in Adolescent Screen Time and Cognition
* Authors: Alex Campbell - alexander.campbell1@student.unimelb.edu.au & Josue Linarte Teran - josue.teran.linarte@hhu.de
* Script:      01_Cleaning.do
* Purpose: Data cleaning, harmonization across waves, and master file creation
* Software: Stata 17 (or higher)
* Inputs: ZA6701 TwinLife datasets (Wave 1, Wave 3, Wave 7, and Weights v9-0-0)
* Outputs: ScreentimeCognition_Master_vX.dta
**============================================================================

**----------------------------------------------------------------------------
**# Section 1.0: Setup & Working Directory
**----------------------------------------------------------------------------
version 17
frames reset
macro drop _all
set linesize 80
set maxvar 20000, permanently

// Use current working directory dynamically for full public replicability
global route "`c(pwd)'"
cd "$route"
global datasets "$route\Stata"

**----------------------------------------------------------------------------
**# Section 2.0: Variable Definitions
**----------------------------------------------------------------------------
* ID & constant demographic variables
local IdVars "wid fid cgr zyg0102"

* Time-varying demographic variables
local DemographicVars "age0100 age0102 age0200 age0201 sex nbi0115"

* Educational variables (math & German grades, prepared for auxiliary/positive control checks)
local EducationVars "cer2200 cer2201" 	

* Cognitive abilities (CFT 20-R and CFT 1-R items)
local CognitiveAbilities "igf0182 igf0282 igf0382 igf0482 igf0582 igf0682 igf0782"

* Screentime variables (Wave 1, Wave 3)
local ComputerWeekday 	"med0100" 
local ComputerWeekend 	"med0200" 
local InternetWeekday 	"med0300" 
local InternetWeekend 	"med0400" 
local PcWeekday 		"med0500" 
local PcWeekend 		"med0600" 
local ConsoleWeekday 	"med0700" 
local ConsoleWeekend 	"med0800" 
local TvWeekday 		"med0900" 
local TvWeekend 		"med1000" 

local ScreentimeVars 	`ComputerWeekday' `ComputerWeekend' ///
						`InternetWeekday' `InternetWeekend'	///
						`PcWeekday' `PcWeekend'				///
						`ConsoleWeekday' `ConsoleWeekend'	///
						`TvWeekday' `TvWeekend'

local FixedScreenRule 	"med2000 med2000t med2000u"

* Socioeconomic status variables (ISEI; prepared for background/stratification checks)
local SesVars 			"emp0505" 

local Pers_Limit        "med1202"


**----------------------------------------------------------------------------
**# Section 3.0: Generating Wave Datasets (Long Format)
**----------------------------------------------------------------------------
foreach WaveNumber in 1 3 7 {

    // Wave 3 does not contain cognitive test items
    if `WaveNumber' == 3 {
        local CurrentCogVars "" 
        local CurrentReshapeList sex age0100 age0101 `SesVars' `ScreentimeVars' `EducationVars'
    }
    else {
        local CurrentCogVars "`CognitiveAbilities'"
        local CurrentReshapeList sex age0100 age0101 `SesVars' `CognitiveAbilities' `ScreentimeVars' `EducationVars'
    }

    use `IdVars' *_t_`WaveNumber' *_u_`WaveNumber' *_m_`WaveNumber' *_f_`WaveNumber' ///
        using "$datasets\ZA6701_family_wide_wid`WaveNumber'_v9-0-0.dta", clear

    quietly mvdecode _all, mv(-99/-1)

    foreach stump in t u m f {
        capture rename *_`stump'_`WaveNumber' *_`stump'
    }

    // -----------------------------------------------------------------------
    // Note: Parental SES harmonization (Highest ISEI across parents)
    // -----------------------------------------------------------------------
    foreach var of local SesVars { 
        capture confirm variable `var'_m `var'_f 
        if !_rc {
            egen `var'_parent = rowmax(`var'_m `var'_f)
            replace `var'_t = `var'_parent 
            replace `var'_u = `var'_parent 
            drop `var'_parent
        }
    }

    foreach var of local SesVars {
        capture confirm variable `var'_t
        if !_rc {
            xtile `var'_tert_t = `var'_t, nq(3)
            generate `var'_tert_u = `var'_tert_t
            capture label define tert_lbl 1 "Low SES" 2 "Medium SES" 3 "High SES"
            label values `var'_tert_t tert_lbl
            label values `var'_tert_u tert_lbl
            local CurrentReshapeList "`CurrentReshapeList' `var'_tert"
        }
    }

    // -----------------------------------------------------------------------
    // Screentime parent-child information retrieval (if self-report was missing)
    // -----------------------------------------------------------------------
    foreach var of local ScreentimeVars {
        capture confirm variable `var't_m `var't_f
        if !_rc {
            egen `var't_parent = rowmax(`var't_m `var't_f)
            capture confirm variable `var'_t
            if !_rc {
                replace `var'_t = `var't_parent if missing(`var'_t)
            }
            else {
                generate `var'_t = `var't_parent
            }
            drop `var't_parent
        }

        capture confirm variable `var'u_m `var'u_f
        if !_rc {
            egen `var'u_parent = rowmax(`var'u_m `var'u_f)
            capture confirm variable `var'_u
            if !_rc {
                replace `var'_u = `var'u_parent if missing(`var'_u)
            }
            else {
                generate `var'_u = `var'u_parent
            }
            drop `var'u_parent
        }
    }

    // Reshaping data to long format
    reshape long `CurrentReshapeList', i(fid) j(ptyp_str) string

    generate ptyp = . 
    recode ptyp . = 1 if ptyp_str == "_t"
    recode ptyp . = 2 if ptyp_str == "_u"
    keep if inlist(ptyp, 1, 2)
    label define ptyp_lbl 1 "1: Twin 1" 2 "2: Twin 2"
    label values ptyp ptyp_lbl 
    drop ptyp_str

    rename wid_`WaveNumber' wid 
    bysort fid: generate pid = fid * 10 + _n 
    order pid ptyp, after(fid)

    // -----------------------------------------------------------------------
    // Section 3.1: Screentime Aggregation
    // -----------------------------------------------------------------------
    capture confirm variable `ScreentimeVars'
    if !_rc {
        foreach var of local ScreentimeVars {
            replace `var' = . if `var' > 25 & !missing(`var')
        }
    }

    // Weekday screentime (excluding internet/pc games to avoid double counting)
    egen n_miss_screen = rowmiss(med0100 med0700 med0900)
    egen medWeekdayTotal = rowtotal(med0100 med0700 med0900) if n_miss_screen < 3
    label var medWeekdayTotal "Total daily weekday screentime"
    
    generate log_medWeekday = log(medWeekdayTotal + 1)
    label var log_medWeekday "Log-transformed daily weekday screentime" 
    
    // Weekend screentime
    egen n_miss_weekend = rowmiss(med0200 med0800 med1000)
    egen medWeekendTotal = rowtotal(med0200 med0800 med1000) if n_miss_weekend < 3
    label var medWeekendTotal "Total daily weekend screentime"
    
    generate log_medWeekend = log(medWeekendTotal + 1)
    label var log_medWeekend "Log-transformed daily weekend screentime" 
    
    // Total weekly average (hours per day)
    generate medTotal = (medWeekdayTotal * 5 + medWeekendTotal * 2) / 7
    label var medTotal "Avg daily screentime (hours)"

    // Censor aggregate multitasking 
    replace medTotal = . if medTotal > 25

    generate log_medTotal = log(medTotal + 1)
    label variable log_medTotal "Log-transformed avg daily screentime"
        
    drop n_miss_screen n_miss_weekend

    // -----------------------------------------------------------------------
    // Section 3.2: Educational Achievement (School Math Grade)
    // Note: Prepared for auxiliary/positive control analyses
    // -----------------------------------------------------------------------
    capture confirm variable cer2200
    if !_rc {
        replace cer2200 = . if cer2200 < 1 | cer2200 > 6
    }

    // -----------------------------------------------------------------------
    // Section 3.3: Cognitive Ability Scoring (CFT 20-R and CFT 1-R)
    // Note: CFT 20-R applies to Cohorts 2-4 (aged 10+), which forms the primary
    // analytical sample. CFT 1-R applies only to Cohort 1 (aged 5).
    // -----------------------------------------------------------------------
    if `WaveNumber' == 1 | `WaveNumber' == 7 {
        capture confirm variable `CognitiveAbilities'
        if !_rc {
            // CFT 20-R (Cohorts 2-4, aged 10+)
            egen n_miss_cft20 = rowmiss(igf0182 igf0282 igf0382 igf0482)
            egen CogTotal_CFT20 = rowtotal(igf0182 igf0282 igf0382 igf0482) if n_miss_cft20 == 0
            
            // CFT 1-R (Cohort 1, aged 5)
            egen n_miss_cft1 = rowmiss(igf0582 igf0682 igf0782)
            egen CogTotal_CFT1 = rowtotal(igf0582 igf0682 igf0782) if n_miss_cft1 == 0
            
            // Cohort-standardized scores
            bysort cgr: egen Z_CFT20 = std(CogTotal_CFT20)
            bysort cgr: egen Z_CFT1  = std(CogTotal_CFT1)

            // Combined standardized score across cohorts
            generate Cog_Z_Cohort = Z_CFT20
            replace  Cog_Z_Cohort = Z_CFT1 if missing(Cog_Z_Cohort)
            label variable Cog_Z_Cohort "Cognitive Ability (Cohort/Test-Standardized Z-Score)"
            
            drop n_miss_cft20 n_miss_cft1 CogTotal_CFT20 CogTotal_CFT1 Z_CFT20 Z_CFT1
        }
    }

    // -----------------------------------------------------------------------
    // Section 3.4: Pruning Variables and Saving Wave Data
    // -----------------------------------------------------------------------
    local GeneratedVars "pid medTotal log_medTotal emp0505_tert Cog_Z_Cohort log_medWeekend log_medWeekday"
    local StructureVars "wid fid cgr ptyp zyg0102"
    local AllWantedVars "`StructureVars' `DemographicVars' `ScreentimeVars' `SesVars' `CognitiveAbilities' `GeneratedVars' `EducationVars'"
    
    local FinalKeepList ""
    foreach var of local AllWantedVars {
        capture confirm variable `var'
        if !_rc {
            local FinalKeepList "`FinalKeepList' `var'"
        }
    }
    keep `FinalKeepList'

    ds wid fid pid cgr ptyp zyg0102, not
    local TimeChangingVars `r(varlist)'
    foreach var of local TimeChangingVars {
        local var_label: variable label `var'
        rename `var' `var'_w0`WaveNumber'
        label var `var'_w0`WaveNumber' "`var_label' wave `WaveNumber'"
    }

    save "$datasets\ScreentimeCognition_w0`WaveNumber'.dta", replace
}


**----------------------------------------------------------------------------
**# Section 4.0: Merging Waves and Adding Weights
**----------------------------------------------------------------------------
use "$datasets\ScreentimeCognition_w01.dta", clear

merge 1:1 pid using "$datasets\ScreentimeCognition_w03.dta"
rename _merge followup_w03
label define _merge_w3 1 "W1 only" 2 "W3 only" 3 "W1 and W3"
label values followup_w03 _merge_w3

merge 1:1 pid using "$datasets\ScreentimeCognition_w07.dta"
rename _merge followup_w07
label define _merge_w7 1 "W1 only" 2 "W7 only" 3 "W1 and W7"
label values followup_w07 _merge_w7

// Merge baseline survey weights
merge m:1 fid using "$datasets\ZA6701_weights_v9-0-0.dta", ///
    keepusing(svw0100 svw0200 svw0301_3 svw0301_5 svw0301_7) ///
    keep(master match) ///
    nogenerate
    
generate w1_baseline_weight = svw0100 * svw0200
label variable w1_baseline_weight "Combined Design & Nonresponse Weight for W1"

// Recode sex to dummy 0/1 (0 = Male, 1 = Female)
recode sex_w01 (1 = 0) (2 = 1)
label define sex_lbl 0 "Male" 1 "Female", replace
label values sex_w01 sex_lbl
label variable sex_w01 "Sex (0 = Male, 1 = Female)"

// Save master analytical file
save "$datasets\ScreentimeCognition_Master_vx.dta", replace

**********************************************
****** END OF DO FILE FOR CLEANING DATA ******
**********************************************