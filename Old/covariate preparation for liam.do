
*creating number of types of physical activity in last 4 weeks

gen walking = A_61640_0
recode walking (1=1)(2/max=0)(-7=0)(-3=.)
label var walking "type of activity - walking"
gen heavydiy= 0 if A_61640_0==-7 | (A_61640_0<5 & A_61640_0>0)| (A_61640_1 <5 & A_61640_1>0) | (A_61640_2<5 & A_61640_2>0) | (A_61640_3<5 & A_61640_3>0) | (A_61640_4<5 & A_61640_4>0) 
replace heavydiy=1 if A_61640_0==5 | A_61640_1==5 | A_61640_2==5 | A_61640_3==5 | A_61640_4==5
label var heavydiy "type of activity - heavydiy"
gen lightdiy=0 if A_61640_0==-7 | A_61640_0~=.
replace lightdiy=1 if A_61640_0==4 | A_61640_1==4 | A_61640_2==4 | A_61640_3==4 | A_61640_4==4
label var lightdiy "type of activity - lightdiy"
gen strensport=0 if A_61640_0==-7 | A_61640_0<6
replace strensport=1 if A_61640_0==3 | A_61640_1==3 | A_61640_2==3 | A_61640_3==3 | A_61640_4==3
label var strensport "type of activity - strenuous sport"
gen othexercise=0 if A_61640_0==-7 | A_61640_0<6
replace othexercise=1 if A_61640_0==2 | A_61640_1==2 | A_61640_2==2 | A_61640_3==2 | A_61640_4==2
label var othexercise "type of activity - other exercise"
gen numbertypes_exercise=walking+heavydiy+lightdiy+strensport+othexercise
label var numbertypes_exercise "number of types of exercise taken"


*creating maximum grip strength
egen maxgrip=rowmax(A_460_0 A_470_0)
label var maxgrip "maximum grip strength"


*creating vascular/heart disease diagnosis (people could give multiple answers to this question, but as I only wanted to know whether individuals had any type of vascular/heart disease, I only needed to use the first variable to derive this)
gen vasculardis_diag=A_61500_0
recode vasculardis_diag (-7=0)(-3=.)(1/4=1)
label var vasculardis_diag "diagnosed with vascular disease or heart disease" 

*Steven's syntax for creating biomarkers for analysis

gen crpnew=crp
label var crpnew "crp adjusted according to Steven Bell's syntax"
replace crpnew=. if crpnew >=10
replace crpnew=0.075 if crpnew <0.15

*creating MND diagnosed at baseline using non-cancer illness codes
gen MND_diag=0
replace MND_diag=1 if A_200020_0==1259
replace MND_diag=1 if A_200020_1==1259
replace MND_diag=1 if A_200020_2==1259
replace MND_diag=1 if A_200020_3==1259
replace MND_diag=1 if A_200020_4==1259
replace MND_diag=1 if A_200020_5==1259
replace MND_diag=1 if A_200020_6==1259
label var MND_diag "MND diagnosed at baseline"
