library(tidyverse)
library(haven)
library(glue)
library(gallimaufr)
library(labelled)

rm(list = ls())

# 1. Load Data ----
df <- read_dta("Data/variables for CRP MND paper.dta") %>%
  full_join(read_dta("Data/liam_12aug21.dta"), by = "eid") %>%
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
  rename(hba1c = HbA1c, hdl = HDL, grip_strength = maxgrip) %>%
  mutate(hba1c = ifelse(hba1c >= 100, NA, hba1c)) %>%
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
    log_crp_fup = "(Log) CRP at Follow-Up",
    hba1c = "HbA1c", hdl = "HDL Cholesterol",
    grip_strength = "Grip Strength"
  ) %>%
  select(eid, matches("crp"), matches("(event|time)"), everything())

glimpse(df)
look_for(df)

# 2, Check CRP Data ----
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Distribution of CRP
ggplot(df) +
  aes(x = crp) +
  geom_density() +
  scale_x_log10()

df %>%
  count(crp) %>%
  ggplot() +
  aes(x = crp, y = n) +
  geom_col() +
  theme_minimal() +
  labs(x = "CRP", y = "Observations")
ggsave("Images/hist.png",
       width = 21, height = 9.9, units = "cm")


# Reliability: Check Linearity Assumption
get_linear <- function(low){
  df_fu <- df %>%
    filter(crp > !!low, 
           crp_fup > !!low)
  
  ggplot(df_fu) +
    aes(x = log_crp, y = log_crp_fup) +
    geom_jitter(alpha = 0.2, color = "grey60") +
    geom_abline(linetype = "dashed") +
    geom_smooth(method = "lm", color = cbPalette[7]) +
    geom_smooth(color = cbPalette[6]) +
    theme_minimal() +
    labs(x = "(Log) Baseline CRP", y = "(Log) CRP at follow-up")
}

get_linear(0.1)
ggsave("Images/measurement_scatter_01.png",
       width = 21, height = 9.9, units = "cm")


get_linear(-Inf)
ggsave("Images/measurement_scatter_inf.png",
       width = 21, height = 9.9, units = "cm")

# Reliability: Check Residuals
get_residual <- function(low){
  df_fu <- df %>%
    filter(crp > !!low, 
           crp_fup > !!low)
  
  lm(log_crp_fup ~ log_crp, df_fu) %>%
    augment() %>%
    ggplot() +
    aes(x = log_crp, y = .resid) +
    geom_hline(yintercept = 0) +
    geom_jitter(alpha = 0.2, color = "grey60") +
    geom_smooth(color = cbPalette[7]) +
    theme_minimal() +
    labs(x = "(Log) Baseline CRP", y = "Residual")
}

get_residual(0.1)
ggsave("Images/residual_scatter_01.png",
       width = 21, height = 9.9, units = "cm")

get_residual(-Inf)
ggsave("Images/residual_scatter_inf.png",
       width = 21, height = 9.9, units = "cm")


# 3. Create Datasets ----
df <- df %>%
  mutate(log_crp_lin = ifelse(log_crp <= 1.3, log_crp, NA))

save(df, file = "Data/df_analysis.Rdata")
