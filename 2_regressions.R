library(tidyverse)
library(haven)
library(glue)
library(broom)
library(survival)
library(gallimaufr)
library(tictoc)
library(furrr)
library(splines)

rm(list = ls())

# 1. Load Data ----
load("Data/df_analysis.Rdata")

covars <- list("Age + Sex" = c("age", "female"),
               "Multiple Adjusted" = c("age", "female", "vasculardis",
                                       "diabetes", "cancer", "fev1", 
                                       "townsend", "BMI", "currsmok",
                                       "nonwhite", "seen_psychiatrist",
                                       "low_activity"))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# 2. Check Data ----
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

df_cox <- df %>%
  filter(inflam == 1)

df_fu <- df %>%
  drop_na(crp, crp_fup)


# 4. Original Regressions ----
reg_grid <- expand_grid(outcome = c("hosp", "dead"),
                        crp = c("rev_crp3", "rev_crp3_f", "log_crp", "crp10", "log_crp_lin"),
                        mod = names(covars),
                        low = c(-Inf, 0.1))


get_cox <- function(outcome, crp, mod, low){
  df_reg <- df_cox %>%
    filter(crp > !!low)
  
  glue("Surv(time_{outcome}, event_{outcome}) ~ 
       {crp} + {glue_collapse(covars[[mod]], ' + ')}") %>%
    as.formula() %>%
    coxph(df_reg, method = "breslow") %>%
    tidy(conf.int = TRUE, conf.level = 0.95, 
         exponentiate = TRUE) %>%
    select(term, beta = estimate, p = p.value,
           lci = conf.low, uci = conf.high)
}

# Run Models
main_res <- reg_grid %>%
  mutate(res = pmap(list(outcome, crp, mod, low), get_cox)) %>%
  unnest(res) %>%
  filter(str_detect(term, crp))


# 3. Regression Dilution ----
# Model Function
get_boot <- function(boot, outcome, crp, mod, low, covars){
  set.seed(boot)
  
  boot_fu <- df_fu %>%
    filter(crp > !!low, crp_fup > !!low) %>%
    sample_frac(replace = TRUE)
  
  boot_df <- df_cox %>%
    filter(crp > !! low) %>%
    sample_frac(replace = TRUE)
  
  df_new <- glue("log_crp_fup ~ {crp} +
                 {glue_collapse(covars[[mod]], ' + ')}") %>% 
    as.formula() %>%
    lm(boot_fu) %>%
    broom::augment(newdata = boot_df)
  
  mod <- glue("Surv(time_{outcome}, event_{outcome}) ~ 
       .fitted + {glue_collapse(covars[[mod]], ' + ')}") %>%
    as.formula() %>%
    coxph(df_new, method = "breslow")
  
  coef(mod)[[".fitted"]]
}

get_ci <- function(estimates){
  quantile(estimates, c(.5, .025, .975)) %>%
    as_tibble_row() %>%
    rename(beta = 1, lci = 2, uci = 3)
}


# Run Models
tic()
plan(multisession, workers = 8)
dilut_res <- reg_grid %>% 
  filter(crp %in% c("log_crp", "log_crp_lin")) %>%
  distinct() %>%
  uncount(1000, .id = "boot") %>%
  sample_frac() %>%
  mutate(estimate = future_pmap_dbl(list(boot, outcome, crp, mod, low), 
                                    get_boot, covars,
                                    .progress = TRUE)) %>%
  arrange(outcome, crp, mod, low) %>%
  group_by(outcome, crp, mod, low) %>%
  summarise(get_ci(estimate),
            .groups = "drop")
future:::ClusterRegistry("stop")
toc()

save(main_res, dilut_res, covars,
     file = "Data/reg_results.Rdata")


# 3. Splines ----
df_cox <- df_cox %>%
  mutate(ns(crp, 3) %>%
           as_tibble() %>%
           mutate(across(everything(), as.numeric)) %>%
           rename_with(~ glue("ns_{.x}")))

df_s <- df_cox %>%
  select(crp, matches("^ns_")) %>%
  arrange(crp) %>%
  distinct() %>%
  pivot_longer(-crp)

get_cox_splines <- function(outcome, mod, low, covars, df_s){
  df_cox <- df_cox %>%
    filter(crp >= !!low)
  
  mod_form <- glue("Surv(time_{outcome}, event_{outcome}) ~
                   {glue_collapse(unique(df_s$name), ' + ')} +
              {glue_collapse(covars[[mod]], ' + ')}")
  
  get_coefs <- function(boot, df_cox, df_s, mod_form){
    set.seed(boot)
    
    df_reg <- sample_frac(df_cox, replace = TRUE)
    
    as.formula(mod_form) %>%
      coxph(df_reg) %>%
      coef() %>%
      enframe(value = "coef") %>%
      right_join(df_s, by = "name") %>%
      group_by(crp) %>%
      summarise(estimate = sum(value*coef),
                .groups = "drop")
  }
  
  map_dfr(1:1000, get_coefs, df_cox, df_s, mod_form, .id = "boot") %>%
    group_by(crp) %>%
    summarise(get_ci(estimate),
              .groups = "drop")
}

tic()
plan(multisession, workers = 8)
spline_res <- reg_grid %>%
  distinct(outcome, mod, low) %>%
  mutate(res = future_pmap(list(outcome, mod, low),
                           get_cox_splines, covars, df_s,
                           .progress = TRUE)) %>%
  unnest(res)
future:::ClusterRegistry("stop")
toc()

save(spline_res, file = "Data/spline_results.Rdata")
