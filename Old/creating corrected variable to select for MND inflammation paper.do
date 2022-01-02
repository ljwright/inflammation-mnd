*creating variable to indicate has complete data on variables on interest & exclusions apply
gen corr_inMNDInflamm=0
replace corr_inMNDInflamm=1 if crpnew3~=. & age_assessment~=. & female~=. & vasculardis_diag~=. & diabetes_diag~=. & cancer_diag~=. & fev1~=. & townsend_score~=. & BMI~=. & currsmok~=. & nonwhite~=. & seen_psychiatrist~=. & rev_anyactivity~=. & MND_diag==0 & admitted_MND_beforebaseline==0 & MND_within3yrs==0
label var corr_inMNDInflamm "corrected variable - included in MND and inflammation analysis"
