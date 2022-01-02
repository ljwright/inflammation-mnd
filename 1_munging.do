cd "L:\David Batty\Data"

// clear all
// insheet using "Cgale_211117_data.csv", comma clear names
// save "df_raw.dta", replace
use "df_raw", clear

foreach var of varlist *{
	local lbl: variable label `var'
	rename `var' v`lbl'
}
gen id = vEID

* Baseline Diagnosis
gen baseline_diag = 0
foreach var of varlist v20002_0*{
	replace baseline_diag = 1 if `var' == 1259
}
gen date_baseline = date(v53_0_0, "YMD")

* Hospital Diagnosis
gen max_date_hosp = 0
gen event_hosp = 0
gen date_hosp = .
qui capture forval i = 0/74{	
	gen date_hosp_`i' = date(v41262_0_`i', "YMD")
	
	egen xx = max(date_hosp_`i')
	replace max_date_hosp = xx if xx > max_date_hosp
	drop xx	
	
	replace date_hosp = date_hosp_`i' if  ///
		v41202_0_`i' == "G122" & event_hosp == 0
	replace date_hosp = date_hosp_`i' if  ///
		v41202_0_`i' == "G122" & !missing(date_hosp) & date_hosp > date_hosp_`i'
		
	replace event_hosp = 1 if v41202_0_`i' == "G122"	
}

drop date_hosp_*

* Death
gen died = 0
gen died_date = .
gen max_date_dead = 0
gen event_dead = 0
gen date_dead = .
qui forval i = 0/1{
	gen date_dead_`i' = date(v40000_`i'_0, "YMD")
	
	egen xx = max(date_dead_`i')
	replace max_date_dead = xx if xx > max_date_dead
	drop xx
	
	replace died = 1 if !missing(date_dead_`i') | v40001_`i'_0 != ""
	replace died_date = date_dead_`i' if missing(died_date) & ///
		(!missing(date_dead_`i') | v40001_`i'_0 != "")
	
	replace date_dead = date_dead_`i' if  ///
		v40001_`i'_0 == "G122" & event_dead == 0
	replace date_dead = date_dead_`i' if  ///
		v40001_`i'_0 == "G122" & !missing(date_dead) & date_dead > date_dead_`i'
		
	replace event_dead = 1 if v40001_`i'_0 == "G122"	
}

drop date_dead_*

replace date_hosp = max_date_hosp if event_hosp == 0
replace date_hosp = died_date if died == 1 & event_hosp == 0
gen time_hosp = (date_hosp - date_baseline) / 365.25
replace date_dead = max_date_hosp if event_dead == 0
replace date_dead = died_date if died == 1 & event_dead == 0
gen time_dead = (date_dead - date_baseline) / 365.25

replace baseline_diag = 1 if time_hosp < 0 & event_hosp == 1

drop max_date_*

* Other Variables
gen townsend = v189_0_0

gen bmi = v21001_0_0

gen age = v21003_0_0

gen diabetes = v2443_0_0 if inrange(v2443_0_0, 0, 1)

gen vasculardis = 0 if v6150_0_0 == -7
replace vasculardis = 1 if inrange(v6150_0_0, 1, 4)

gen cancer = v2453_0_0 if inrange(v2453_0_0, 0, 1)

gen seen_psychiatrist = v2100_0_0 if inrange(v2100_0_0, 0, 1)

gen smoke_status = v20116_0_0 if inrange(v20116_0_0, 0, 2)
label define smoke_status 0 "Never" 1 "Former" 2 "Current"
label values smoke_status smoke_status
gen smoke_curr = smoke_status == 2 if !missing(smoke_status)

gen low_activity = 0 if inrange(v6164_0_0, 1, 5)
replace low_activity = 1 if v6164_0_0 == -7

gen crp = v30710_0_0
gen crp_fup = v30710_1_0 // change if below < 0.15

gen sex = v31_0_0 if inrange(v31_0_0, 0, 1) 
label define sex 0 "Female" 1 "Male"
label values sex sex
gen female = sex == 0 if inrange(sex, 0, 1)

gen ethnic_group = cond(v21000_0_0 < 1000, v21000_0_0, floor(v21000_0_0/1000))
replace ethnic_group = . if ethnic_group < 0
label define ethnic_group 1 "White" 2 "Mixed" 3 "Asian" 4 "Black" 5 "Chinese" 6 "Other"
label values ethnic_group ethnic_group
gen non_white = ethnic_group != 1 if !missing(ethnic_group)

gen height = v50_0_0

gen hba1c = v30750_0_0
    
gen hdl = v30760_0_0
    
gen fev_1 = max(v3063_0_0, v3063_0_1, v3063_0_2)
    
gen grip_strength = max(v46_0_0, v47_0_0)


* Choose Sample
keep id-grip_strength

local vlist time_dead event_dead time_hosp event_hosp ///
			crp age female vasculardis diabetes cancer ///
			fev_1 townsend bmi smoke_status non_white ///
			seen_psychiatrist low_activity
misstable sum `vlist'

drop if time_dead < 0
drop if baseline_diag == 1
drop if time_hosp < 3 & event_hosp == 1
drop if time_dead < 3 & event_dead == 1
drop if crp >= 10 & !missing(crp)

egen miss = rowmiss(`vlist')
drop if miss >= 1

tab female
save df_analysis.dta, replace