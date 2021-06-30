*TABLE 1 (based on all those with data on crp)
ta crpnew3
su crpnew3 if crpnew3==0 detail
su crpnew3 if crpnew3==1, detail
su crpnew3 if crpnew3==2, detail

oneway fev1 crpnew3, ta  
oneway age_assessment crpnew3, ta  
oneway BMI crpnew3, ta  
oneway height crpnew3, ta  
ta disadvantaged crpnew3, col chi2

ta female crpnew3, col chi2
ta nonwhite crpnew3, col chi2
ta currsmok crpnew3, col chi2
ta rev_anyactivity crpnew3, col chi2
ta vasculardis_diag crpnew3, col chi2
ta diabetes_diag crpnew3, col chi2
ta cancer_diag crpnew3, col chi2
ta seen_psychiatrist crpnew3, col chi2
ta  crpnew3, col chi2

*TABLE 2 HRs for hospitalisation (based on those with complete data, no prior diagnosis and not admitted with/died of MND within 3 years)

*cases per crp group
ta update_MNDhospdiag rev_crpnew3 if inMNDinflamm==1

*declaring survival data
stset update_survivaltime_MNDhospdiag, failure(update_MNDhospdiag) id(eid)

*adjusted age and sex according to thirds of distribution and SD increase in logged CRP

xi:stcox i.rev_crpnew3 age_assessment female if inMNDinflamm==1
stcox rev_crpnew3 age_assessment female if  inMNDinflamm==1
stcox log_crpnew_sd age_assessment female if  inMNDinflamm==1

*multiple adjustment as above
xi:stcox i.rev_crp3 age_assessment female vasculardis_diag diabetes_diag cancer_diag fev1 townsend_score BMI currsmok nonwhite seen_psychiatrist rev_anyactivity if inMNDinflamm==1
stcox rev_crpnew3 age_assessment female vasculardis_diag diabetes_diag cancer_diag fev1 townsend_score BMI currsmok nonwhite seen_psychiatrist rev_anyactivity if inMNDinflamm==1
stcox log_crpnew_sd age_assessment female vasculardis_diag diabetes_diag cancer_diag fev1 townsend_score BMI currsmok nonwhite seen_psychiatrist rev_anyactivity if  inMNDinflamm==1

*TABLE 2 HRs for death (based on those with complete data, no prior diagnosis and not admitted with/died of MND within 3 years)

*cases per crp group
ta dead_MND rev_crpnew3 if inMNDinflamm==1


*declaring survival data
stset survivaltime_dead4oct2020, failure(dead_MND) id(eid)

*adjusted age and sex according to.. (as above)
xi:stcox i.rev_crpnew3 age_assessment female if inMNDinflamm==1
stcox rev_crpnew3 age_assessment female if inMNDinflamm==1
stcox log_crpnew_sd age_assessment female if  inMNDinflamm==1

*multiple adjustment as above
xi:stcox i.rev_crpnew3 age_assessment female vasculardis_diag diabetes_diag cancer_diag fev1 townsend_score BMI currsmok nonwhite seen_psychiatrist rev_anyactivity if inMNDinflamm==1
stcox rev_crpnew3 age_assessment female vasculardis_diag diabetes_diag cancer_diag fev1 townsend_score BMI currsmok nonwhite seen_psychiatrist rev_anyactivity if inMNDinflamm==1
stcox log_crpnew_sd age_assessment female vasculardis_diag diabetes_diag cancer_diag fev1 townsend_score BMI currsmok nonwhite seen_psychiatrist rev_anyactivity if inMNDinflamm==1

*FOR FIGURE - both outcomes according to 10ths, multiple adjustment

*hospitalisation

*cases per crp group
ta update_MNDhospdiag crp10 if inMNDinflamm==1

*declaring survival data
stset update_survivaltime_MNDhospdiag, failure(update_MNDhospdiag) id(eid)

*multiple adjustment
xi:stcox i.crp10 age_assessment female vasculardis_diag diabetes_diag cancer_diag fev1 townsend_score BMI currsmok nonwhite seen_psychiatrist rev_anyactivity if inMNDinflamm==1


*death

*cases per crp group
ta dead_MND crp10 if inMNDinflamm==1

*declaring survival data
stset survivaltime_dead4oct2020, failure(dead_MND) id(eid)

*multiple adjustment
xi:stcox i.crp10 age_assessment female vasculardis_diag diabetes_diag cancer_diag fev1 townsend_score BMI currsmok nonwhite seen_psychiatrist rev_anyactivity if inMNDinflamm==1
