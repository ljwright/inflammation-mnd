library(tidyverse)
library(labelled)
library(glue)
library(summarytools)
library(gallimaufr)
library(broom)

rm(list = ls())

options(scipen = 999)

# 1. Load Data ----
df_raw <- read_csv("Data/Cgale_211117_data.csv", guess_max = 502421) %>%
  rename_with(~ glue("v{.x}"), .cols = -EID) %>%
  relocate(sort(current_vars()))

miss <- map_dbl(df_raw, ~ 100*sum(is.na(.x))/length(.x))

fields <- read_csv("Data/Fieldlist_211117.csv") %>%
  mutate(FID = glue("v{FID}"),
         miss = miss[FID]) %>%
  mutate(lbl_lower = str_to_lower(FNAME)) %>%
  arrange(desc(miss)) %>%
  rename_with(str_to_lower)

save(df_raw, fields, file = "Data/df_raw.Rdata")
write_csv(fields, file = "Data/missingness.csv")


# 2. Find Variables ----
load("Data/df_raw.Rdata")

# fields %>%
#   select(fid, fname, miss, lbl_lower) %>%
#   filter(str_detect(lbl_lower, "activ")) %>% View()

#https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20116

# 3. Clean Dataset ----
# Baseline Diagnosis
df_baseline <- df_raw %>%
  select(id = EID, matches("v20002_0")) %>%
  pivot_longer(-id) %>%
  filter(value == 1259) %>%
  distinct(id) %>%
  mutate(baseline_diag = 1) %>%
  full_join(df_raw %>%
              select(id = EID),
            by = "id") %>%
  mutate(baseline_diag = ifelse(is.na(baseline_diag), 0, 1)) %>%
  select(id, baseline_diag)

# Hospital Diagnosis of MND
df_hosp <- df_raw %>%
  select(id = EID, matches("v412(0|6)2_0")) %>%
  pivot_longer(-id,
               names_to = c(".value", "response"),
               names_pattern = "(.*)_0_(.*)") %>%
  rename(diag_hosp = v41202, date_hosp = v41262) %>%
  mutate(date_hosp = as.Date(date_hosp, format = "%Y-%m-%d"),
         date_max = max(date_hosp, na.rm = TRUE)) %>%
  filter(diag_hosp == "G122") %>%
  arrange(id, date_hosp) %>%
  group_by(id, date_max) %>%
  summarise(date_hosp = min(date_hosp),
            .groups = "drop") %>%
  mutate(event_hosp = 1) %>%
  full_join(df_raw %>%
              select(id = EID),
            by = "id") %>%
  mutate(event_hosp = ifelse(is.na(event_hosp), 0, 1),
         date_hosp = if_else(event_hosp == 1, 
                             date_hosp,
                             max(date_max, na.rm = TRUE))) %>%
  select(id, event_hosp, date_hosp)

# Death from MND
df_dead <- df_raw %>%
  select(id = EID, matches("v4000[0-1]")) %>%
  pivot_longer(-id,
               names_to = c(".value", "response"),
               names_pattern = "(.*)_(.*)_0") %>%
  rename(diag_dead = v40001, date_dead = v40000) %>%
  mutate(died = ifelse(!is.na(diag_dead) | !is.na(date_dead), 1, 0),
         event_dead = case_when(diag_dead == "G122" ~ 1,
                                TRUE ~ 0),
         date_dead = case_when(diag_dead == "G122" ~ date_dead,
                               !is.na(date_dead) ~ date_dead, 
                               TRUE ~ max(date_dead, na.rm = TRUE))) %>%
  pivot_wider(names_from = response, 
              values_from = -c(id, response)) %>%
  mutate(died = pmax(died_0, died_1),
         event_dead = pmax(event_dead_0, event_dead_1),
         date_dead = pmin(date_dead_0, date_dead_1, na.rm = TRUE)) %>%
  distinct(id, died, event_dead, date_dead)

# Other Variables
clean_binary <- function(var){
  ifelse(var %in% 0:1, var, NA)
}

df_cov <- df_raw %>%
  mutate(
    id = EID,
    date_baseline = v53_0_0,
    townsend = v189_0_0, 
    bmi = v21001_0_0, 
    age = v21003_0_0,
    diabetes = clean_binary(v2443_0_0),
    vasculardis = case_when(v6150_0_0 == -7 ~ 0,
                            between(v6150_0_0, 1, 4) ~ 1),
    cancer = clean_binary(v2453_0_0),
    seen_psychiatrist = clean_binary(v2100_0_0),
    smoke_status = factor(v20116_0_0, levels = 0:2, 
                          labels = c("Never", "Former", "Current")),
    smoke_curr = ifelse(smoke_status == "Current", 1, 0),
    low_activity = case_when(v6164_0_0 == -7 ~ 1,
                             between(v6164_0_0, 1, 5) ~ 0),
    crp = v30710_0_0, crp_fup = v30710_1_0, # change if below < 0.15
    sex = factor(v31_0_0, 0:1, c("Female", "Male")),
    female = ifelse(sex == "Female", 1, 0),
    ethnic_group = ifelse(v21000_0_0 < 0, NA, v21000_0_0) %>%
      as.character() %>% 
      str_sub(1, 1) %>%
      factor(levels = 1:6, 
             labels = c("White", "Mixed", "Asian", "Black", "Chinese", "Other")),
    non_white = ifelse(ethnic_group == "White", 0, 1),
    height = v50_0_0, 
    hba1c = v30750_0_0, 
    hdl = v30760_0_0,
    fev_1 = pmax(v3063_0_0, v3063_0_1, v3063_0_2, na.rm = TRUE),
    grip_strength = pmax(v46_0_0, v47_0_0, na.rm = TRUE),
    .before = 1
  ) %>%
  select(id:grip_strength)

# Combine
df_new <- df_baseline %>%
  full_join(df_hosp, by = "id") %>%
  full_join(df_dead, by = "id") %>%
  full_join(df_cov, by = "id") %>%
  mutate(
    time_dead = as.numeric(date_dead - date_baseline)/365.25,
    date_hosp = if_else(died == 1 & event_hosp == 0, date_dead, date_hosp),
    time_hosp = as.numeric(date_hosp - date_baseline)/365.25,
    baseline_diag = ifelse(time_hosp < 0, 1, baseline_diag)
  )

gc()

save(df_new, file = "Data/df_new.Rdata")


# 4. Exclude Participants ---- 
load("Data/df_new.Rdata")

covars <- lst("Age + Sex" = c("age", "female"),
              "Multiply-Adjusted" = c(`Age + Sex`,
                                      "vasculardis", "diabetes", "cancer",
                                      "fev_1", "townsend", "bmi",
                                      "smoke_status", "non_white",
                                      "seen_psychiatrist", "low_activity"),
              "Fully-Adjusted" = c(`Multiply-Adjusted`, 
                                   "hdl", "hba1c", "grip_strength"))
save(covars, file = "Data/model_covars.Rdata")

vars <- map(covars, ~ c(.x,c("time_dead", "event_dead",
                               "time_hosp", "event_hosp",
                               "crp")))

df_exclude <- df_new %>%
  add_count(name = "n_all") %>%
  filter(baseline_diag == 0) %>%
  add_count(name = "n_baseline") %>%
  mutate(within_3years = case_when(between(time_hosp, 0, 3) & event_hosp == 1 ~ 1,
                                   between(time_dead, 0, 3) & event_dead == 1 ~ 1,
                                   TRUE ~ 0)) %>%
  filter(within_3years == 0) %>%
  add_count(name = "n_3years") %>%
  filter(is.na(crp) | crp < 10) %>%
  add_count(name = "n_crp10") %>%
  filter(!is.na(crp)) %>%
  add_count(name = "n_crp_miss") %>%
  drop_na(all_of(vars[[2]])) %>%
  add_count(name = "n_miss") %>%
  select(id, matches("^n_"))

# Missingness
get_n <- function(df){
  df %>%
    distinct(across(matches("^n_"))) %>%
    pivot_longer(everything(), names_to = "sample", values_to = "n") %>%
    mutate(drop_n = lag(n) - n,
           drop_p = drop_n*100/lag(n),
           p = n*100/first(n),
           across(c(drop_p, p), round, 2))
}

get_n(df_exclude)

df_new %>%
  drop_na(all_of(vars[[3]])) %>%
  add_count(name = "n_miss_all") %>%
  select(id, n_miss_all) %>%
  inner_join(df_exclude, ., by = "id") %>%
  get_n()
  
df_new %>%
  drop_na(crp_fup) %>%
  count()

df_new %>%
  semi_join(df_exclude, by = "id") %>%
  drop_na(crp_fup) %>%
  count()

df_new %>%
  select(id, all_of(vars[[3]])) %>%
  mutate(across(-id, ~ ifelse(is.na(.x), 1, 0))) %>%
  pivot_longer(-id) %>%
  filter(value == 1) %>%
  count(name)

df_new %>%
  semi_join(df_exclude, by = "id") %>%
  count(female) %>%
  add_count(wt = n, name = "total")

gc()

# 5. Final Dataset ----
clean_lbls <- c(
  id = "ID", date_baseline = "Date at Baseline",
  event_hosp = "Hospital Diagnosis of MND", 
  time_hosp = "Follow-Up Time (Hospitalisation)",
  event_dead = "Died of MND",
  time_dead = "Follow-Up Time (Death)",
  crp = "C-Reactive Protein", log_crp = "Log CRP",
  crp_3 = "CRP Tertiles", crp_3f = "CRP Tertiles",
  crp_10f = "CRP Deciles",
  crp_fup = "CRP at Follow-Up",
  log_crp_fup = "(Log) CRP at Follow-Up",
  female = "Female", age = "Age",
  non_white = "Non-White", 
  townsend = "Townsend Score",
  disadvantaged = "Disadvantaged Neighbourhood",
  low_activity = "Low Physical Activity",
  smoke_status = "Smoking Status",
  diabetes = "Diabetes", vasculardis = "Vascular Disease",
  cancer = "Cancer", seen_psychiatrist = "Psychiatric Consulation",
  bmi = "Body Mass Index", height = "Height (cm)",
  hba1c = "HbA1c", hdl = "HDL Cholesterol",
  grip_strength = "Grip Strength", fev_1 = "FEV1"
)


df <- df_new %>%
  semi_join(df_exclude, by = "id") %>%
  mutate(
    crp_3f = cut_number(crp, 3),
    crp_3 = as.numeric(crp_3f) - 1,
    crp_10f = cut_number(crp, 10),
    log_crp = log(crp) %>% wtd_scale(),
    log_crp_fup = log(crp_fup) %>% wtd_scale(),
    disadvantaged = ifelse(percent_rank(townsend) >= 0.8, 1, 0)
  ) %>%
  select(all_of(names(clean_lbls))) %>%
  set_variable_labels(.labels = as.list(clean_lbls), .strict = FALSE)


# 6. Check CRP Data ----
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
