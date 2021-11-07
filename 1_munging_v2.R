library(tidyverse)
library(labelled)
library(glue)
library(summarytools)

rm(list = ls())

options(scipen=999)

# 1. Load Data ----
df <- read_csv("Data/Cgale_211027_Data.csv", col_types = cols(.default = "d"))

miss <- map_dbl(df, ~ 100*sum(is.na(.x))/length(.x))

fields <- read_csv("Data/CGale_271022_Fieldlist.csv") %>%
  mutate(miss = miss[FID]) %>%
  mutate(lbl_lower = str_to_lower(FNAME)) %>%
  arrange(desc(miss)) %>%
  rename_with(str_to_lower)

write_csv(fields, file = "Data/missingness.csv")

# 2. Clean Dataset ----
fields %>%
  select(fid, fname, miss, lbl_lower) %>%
  mutate(fid = glue("v{fid}")) %>%
  filter(str_detect(lbl_lower, "activ")) %>% View()

df %>%
  rename_with(~ glue("v{.x}"), .cols = -EID) %>%
  # select(id = EID, matches("v6164_0"))  %>%
  # pivot_longer(-id) %>%
  # count(value)
  count(v6164_0_0)

df %>%
  rename_with(~ glue("v{.x}"), .cols = -EID) %>%
  mutate(ethnic_group = ifelse(v21000_0_0 < 0, NA, v21000_0_0) %>%
           as.character() %>% str_sub(1, 1) %>%
           factor(levels = 1:6, 
                  labels = c("White", "Mixed", "Asian", "Black", "Chinese", "Other"))) %>%
  count(ethnic)


#https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20116

clean_binary <- function(var){
  ifelse(var %in% 0:1, var, NA)
}

df %>%
  rename_with(~ glue("v{.x}"), .cols = -EID) %>%
  mutate(eid = EID, townsend = v189_0_0, 
         bmi = v21001_0_0, age = v21003_0_0,
         diabetes = clean_binary(v2443_0_0),
         vasculardis = case_when(v6150_0_0 == -7 ~ 0,
                                 between(v6150_0_0, 1, 4) ~ 1),
         cancer = clean_binary(v2453_0_0),
         seen_psychiatrist = clean_binary(v2100_0_0),
         smok_status = factor(v20116_0_0, levels = 0:2, labels = c("Never", "Former", "Current")),
         smok_curr = ifelse(smok_status == "Current", 1, 0),
         low_activity = case_when(v6164_0_0 == -7 ~ 1,
                                  between(v6164_0_0, 1, 5) ~ 0),
         crp = v30710_0_0, crp_fup = v30710_1_0,
         sex = factor(v31_0_0, 0:1, c("Female", "Male")),
         female = ifelse(sex == "Female", 1, 0),
         ethnic_group = ifelse(v21000_0_0 < 0, NA, v21000_0_0) %>%
           as.character() %>% 
           str_sub(1, 1) %>%
           factor(levels = 1:6, 
                  labels = c("White", "Mixed", "Asian", "Black", "Chinese", "Other")),
         non_white = ifelse(ethnic_group == "White", 0, 1),
         fev_1 = max(v3063_0_0, v3063_0_1), #  ONLY SEE TWO MEASUREMENTS
         height = v50_0_0, hba1c = v30750_0_0, hdl = v30760_0_0, # OUTLIERS FOR EACH?
         grip_strength = max(v46_0_0, v47_0_0),
         .before = 1) %>% count(low_activity, v6164_0_0)
set_variable_labels(eid = "ID",
                    townsend = "Townsend Score",
                    BMI = "Body Mass Index", age = "Age",
                    diabetes = "Diabetes", vasculardis = "Vascular Disease",
                    cancer = "Cancer", seen_psychiatrist = "Psychiatric Consulation",
                    currsmok = "Current Smoker", low_activity = "Low Physical Activity",
                    crp = "C-Reactive Protein", female = "Female",
                    time_dead = "Follow-Up Time (Death)",
                    event_dead = "Died of MND", nonwhite = "Non-White",
                    log_crp = "(Log) C-Reactive Protein",
                    event_hosp = "Hospital Diagnosis of MND", 
                    time_hosp = "Follow-Up Time (Hospitalisation)",
                    fev1 = "FEV1", inflam = "Included in Analysis",
                    height = "Height (cm)",
                    disadvantaged = "Disadvantaged Neighbourhood",
                    rev_crp3 = "CRP Tertiles", crp10 = "CRP Deciles",
                    crp_fup = "CRP at Follow-Up", rev_crp3_f = "CRP Tertiles",
                    log_crp_fup = "(Log) CRP at Follow-Up",
                    hba1c = "HbA1c", hdl = "HDL Cholesterol",
                    grip_strength = "Grip Strength",
                    .strict = FALSE) %>%
  select(eid:male) %>% lookfor()

names(df) 

df$`189_0_0`








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