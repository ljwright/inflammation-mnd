library(tidyverse)
library(haven)
library(glue)
library(broom)
library(survival)
library(gallimaufr)
library(tictoc)
library(furrr)
library(splines)
library(flextable)
library(officer)

rm(list = ls())

# 1. Load Data ----
df <- read_dta("Data/variables for CRP MND paper.dta") %>%
  rename_with(~ str_replace(.x, "crpnew", "crp")) %>%
  rename(time_hosp = update_survivaltime_MNDhospdiag,
         event_hosp = update_MNDhospdiag,
         time_dead = survivaltime_dead4oct2020,
         event_dead = dead_MND,
         inflam = inMNDinflamm) %>%
  mutate(rev_crp3_f = factor(rev_crp3, 0:2),
         crp10 = factor(crp10),
         log_crp_fup = log(crp_fup) %>% wtd_scale(),
         log_crp = log(crp) %>% wtd_scale())

reg_grid <- expand_grid(outcome = c("hosp", "dead"),
                        crp = c("rev_crp3", "rev_crp3_f", "log_crp", "crp10"),
                        covars = list(c("age_assessment", "female"),
                                      c("age_assessment", "female", "vasculardis_diag",
                                        "diabetes_diag", "cancer_diag", "fev1", 
                                        "townsend_score", "BMI", "currsmok",
                                        "nonwhite", "seen_psychiatrist",
                                        "rev_anyactivity")))

# 2. Original Regressions ----
df_cox <- df %>%
  filter(inflam == 1,
         crp > 0.1)

get_cox <- function(outcome, crp, covars){
  glue("Surv(time_{outcome}, event_{outcome}) ~ 
       {crp} + {glue_collapse(covars, ' + ')}") %>%
    as.formula() %>%
    coxph(df_cox, method = "breslow") %>%
    tidy(conf.int = TRUE, conf.level = 0.95) 
}

main_res <- reg_grid %>%
  mutate(res = pmap(list(outcome, crp, covars), get_cox))


# 3. Regression Dilution ----
# CHECK LINEARITY ASSUMPTION
df_fu <- df %>%
  filter(crp > 0.1, 
         crp_fup > 0.1)
rm(df)


# Check Linearity Assumption
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(df_fu) +
  aes(x = log_crp, y = log_crp_fup) +
  geom_jitter(alpha = 0.2, color = "grey60") +
  geom_abline(linetype = "dashed") +
  geom_smooth(method = "lm", color = cbPalette[7]) +
  geom_smooth(color = cbPalette[6]) +
  theme_minimal() +
  labs(x = "(Log) Baseline CRP", y = "(Log) CRP at follow-up")
ggsave("Images/measurement_scatter.png",
       width = 21, height = 9.9, units = "cm")

lm(log_crp_fup ~ log_crp, df_fu) %>%
  augment() %>%
  ggplot() +
  aes(x = log_crp, y = .resid) +
  geom_hline(yintercept = 0) +
  geom_jitter(alpha = 0.2, color = "grey60") +
  geom_smooth(color = cbPalette[7]) +
  theme_minimal() +
  labs(x = "(Log) Baseline CRP", y = "Residual")


df %>%
  count(crp) %>%
  ggplot() +
  aes(x = crp, y = n) +
  geom_col() +
  theme_minimal() +
  labs(x = "CRP", y = "Observations")
ggsave("Images/hist.png",
       width = 21, height = 9.9, units = "cm")

# Run Regressions
get_boot <- function(boot, outcome, covars, top_score){
  set.seed(boot)
  boot_fu <- df_fu %>%
    filter(log_crp_sd <= !!top_score) %>%
    sample_frac(replace = TRUE)
  boot_df <- df_cox %>%
    filter(log_crp_sd <= !!top_score) %>%
    sample_frac(replace = TRUE)
  
  df_new <- glue("log_crp_fup ~ log_crp_sd + {glue_collapse(covars, ' + ')}") %>% 
    as.formula() %>%
    lm(boot_fu) %>%
    broom::augment(newdata = boot_df)
  
  mod <- glue("Surv(time_{outcome}, event_{outcome}) ~ 
       .fitted + {glue_collapse(covars, ' + ')}") %>%
    as.formula() %>%
    coxph(df_new, method = "breslow")
  
  coef(mod)[[".fitted"]]
}


tic()
plan(multisession, workers = 4)
dilut_res <- reg_grid %>% 
  select(-crp) %>%
  distinct() %>%
  uncount(1000, .id = "boot") %>%
  expand_grid(top_score = c(1.3, Inf)) %>%
  sample_frac() %>%
  mutate(estimate = future_pmap_dbl(list(boot, outcome, covars, top_score), get_boot,
                               .progress = TRUE)) %>%
  arrange(outcome, covars, top_score) %>%
  group_by(outcome, covars, top_score) %>%
  summarise(beta = quantile(estimate, probs = .5),
            lci = quantile(estimate, probs = .025),
            uci = quantile(estimate, probs = .975))
future:::ClusterRegistry("stop")
toc()


# 3. Splines ----
df_all <- df_cox %>%
  mutate(ns(crp, 3) %>%
           as_tibble() %>%
           mutate(across(everything(), as.numeric)) %>%
           rename_with(~ glue("ns_{.x}")))

df_s <- df_all %>%
  select(crp, matches("^ns_")) %>%
  arrange(crp) %>%
  distinct() %>%
  nest(data = -crp)

get_cox_splines <- function(outcome, covars){
  ns_names <- str_subset(names(df_all), '^ns_') %>%
    c(covars) %>%
    glue_collapse(" + ")
  
  mod <- glue("Surv(time_{outcome}, event_{outcome}) ~ {ns_names}") %>%
    as.formula() %>%
    coxph(df_all, method = "breslow")
  
  get_spec <- function(row){
    data <- df_s$data[[row]]

    spec <- enframe(coef(mod)) %>%
      mutate(value = unlist(data)[name],
             value = ifelse(is.na(value), 0, value)) %>%
      deframe() %>%
      as_tibble_row() %>%
      as.matrix()

    multcomp::glht(mod, spec) %>%
      confint() %>%
      pluck("confint") %>%
      as_tibble() %>%
      rename(beta = 1, lci = 2, uci = 3)
  }

  map_dfr(1:nrow(df_s), get_spec) %>%
    mutate(crp = df_s$crp, .before = 1)
}

spline_res <- reg_grid %>%
  select(-crp) %>%
  distinct() %>%
  mutate(res = map2(outcome, covars, get_cox_splines)) %>%
  unnest(res)


save(main_res, dilut_res, spline_res,
     file = "Data/reg_results.Rdata")


# 4. Results ----
clean_res <- function(res_df){
  res_df %>%
    mutate(mod = ifelse(map_dbl(covars, length) == 2, "Age + Female", "Fully Adjusted"),
           dep_clean = ifelse(outcome == "hosp", "Hospitalisation", "Death")) %>%
    mutate(across(beta:uci, exp)) 
}

# Dilution
clean_res(dilut_res) %>%
  mutate(obs_clean = ifelse(top_score == 1.3, "Log CRP \U02264 1.3", "All Observations")) %>%
  ggplot() +
  aes(x = obs_clean, color = mod, shape = mod,
      y = beta, ymin = lci, ymax = uci) +
  facet_wrap(~ dep_clean, scales = "free") +
  geom_hline(yintercept = 1) +
  geom_pointrange(position = position_dodge(0.5)) +
  geom_ribbon(alpha = 0.2, color = NA) +
  scale_y_log10() +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(x = NULL, y = "Hazard Ratio (+ 95% CI)",
       color = "Model", fill = "Model", shape = "Model") +
  theme(legend.position = c(.1, .9))
ggsave("Images/dilute_plot.png",
       width = 21, height = 9.9, units = "cm")

dilut_tbl <- clean_res(dilut_res) %>%
  ungroup() %>%
  mutate(obs_clean = ifelse(top_score == 1.3, "Log CRP \U02264 1.3", "All Observations"),
         across(beta:uci, round, 2),
         string = glue("{beta} ({lci}, {uci})")) %>%
  select(dep_clean, obs_clean, mod, string) %>%
  arrange(dep_clean, obs_clean, mod) %>%
  make_flx(list(dep_clean = "Outcome", obs_clean = "Observations", 
                mod = "Model", string = "HR (95% CI)")) %>%
  merge_v(2) %>%
  fix_border_issues()
dilut_tbl
save_as_docx(dilut_tbl, path = "Tables/dilut_res.docx")

# Spline Plot
clean_res(spline_res) %>%
  filter(beta != 0) %>%
  ggplot() +
  aes(x = crp, y = beta, ymin = lci, ymax = uci,
      color = mod, fill = mod) +
  facet_wrap(~ dep_clean, scales = "free") +
  geom_hline(yintercept = 1) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_line() +
  scale_y_log10() +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(x = "Baseline CRP", y = "Hazard Ratio (+ 95% CI)",
       color = "Model", fill = "Model") +
  theme(legend.position = c(.1, .9))
ggsave("Images/splines_plot.png",
       width = 21, height = 9.9, units = "cm")


# Main Results
main_res %>%
  unnest(res) %>%
  select(outcome:term, beta = estimate,
         lci = conf.low, uci = conf.high) %>%
  clean_res() %>%
  filter(str_detect(term, crp),
         crp != "crp10",
         crp == "log_crp")


main_res %>%
  filter(crp == "rev_crp3_f") %>%
  unnest(res) %>%
  filter(str_detect(term, crp)) %>%
  mutate(across(c(estimate, conf.low, conf.high), exp))
