library(tidyverse)
library(haven)
library(glue)
library(gallimaufr)
library(labelled)

rm(list = ls())

# 1. Load Data ----
df <- read_dta("Data/variables for CRP MND paper.dta") %>%
  rename_with(~ str_replace(.x, "crpnew", "crp")) %>%
  rename_with(~ str_replace(.x, "_diag$", "")) %>%
  rename(time_hosp = update_survivaltime_MNDhospdiag,
         event_hosp = update_MNDhospdiag,
         time_dead = survivaltime_dead4oct2020,
         event_dead = dead_MND,
         inflam = inMNDinflamm, 
         age = age_assessment, 
         low_activity = rev_anyactivity,
         townsend = townsend_score) %>%
  mutate(rev_crp3_f = factor(rev_crp3, 0:2),
         crp10 = factor(crp10),
         log_crp_fup = log(crp_fup) %>% wtd_scale(),
         log_crp = log(crp) %>% wtd_scale()) %>%
  select(-log_crp_sd, -admitted_MND_beforebaseline, -ethnic_group,
         -numbertypes_exercise, -smokestatus, -MND, -MND_within3yrs,
         -crp3) %>%
  zap_label() %>%
  zap_formats() %>%
  set_variable_labels(
    eid = "ID",
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
    log_crp_fup = "(Log) CRP at Follow-Up"
  ) %>%
  select(eid, matches("crp"), matches("(event|time)"), everything())

glimpse(df)
look_for(df)

save(df, file = "Data/df_analysis.Rdata")