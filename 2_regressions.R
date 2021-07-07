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

df_cox <- df %>%
  filter(inflam == 1) %>%
  mutate(ns(crp, 3) %>%
           as_tibble() %>%
           mutate(across(everything(), as.numeric)) %>%
           rename_with(~ glue("ns_{.x}")))

df_fu <- df %>%
  drop_na(crp, crp_fup)

covars <- list("Age + Sex" = c("age", "female"),
               "Multiple Adjusted" = c("age", "female", "vasculardis",
                                       "diabetes", "cancer", "fev1", 
                                       "townsend", "BMI", "currsmok",
                                       "nonwhite", "seen_psychiatrist",
                                       "low_activity"))

reg_grid <- expand_grid(outcome = c("hosp", "dead"),
                        crp = c("rev_crp3", "rev_crp3_f", "log_crp", "crp10", "log_crp_lin"),
                        mod = names(covars),
                        low = c(-Inf, 0.1))

get_ci <- function(estimates){
  quantile(estimates, c(.5, .025, .975)) %>%
    as_tibble_row() %>%
    rename(beta = 1, lci = 2, uci = 3)
}

save(df_cox, df_fu, covars, reg_grid,
     file = "Data/df_regressions.Rdata")

rm(df)


# 2. Original Regressions ----
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
  
  converge <- 0L
  
  while (converge == 0L){
    out <- tryCatch(
      {
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
      },
      error = function(cond) {
        return(0L)
      },
      warning  =function(cond) {
        return(0L)
      }
    ) 
    
    if (!identical(out, 0L)){
      converge <- 1
    }
  }
  
  return(out)
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


# 4. Splines ----
df_s <- df_cox %>%
  select(crp, matches("^ns_")) %>%
  arrange(crp) %>%
  distinct() %>%
  pivot_longer(-crp)

# Bootstraps
get_splines <- function(outcome, mod, low){
  df_cox <- df_cox %>%
    filter(crp >= !!low)
  
  mod_form <- glue_collapse(covars[[mod]], ' + ')
  mod_form <- glue("Surv(time_{outcome}, event_{outcome}) ~
                   ns_1 + ns_2 + ns_3 + {mod_form}")
  
  get_coefs <- function(boot, df_cox, mod_form){
    set.seed(boot)
    
    converge <- 0L
    
    while (converge == 0L){
      out <- tryCatch(
        {
          df_reg <- sample_frac(df_cox, replace = TRUE)
          
          as.formula(mod_form) %>%
            coxph(df_reg) %>%
            coef() %>%
            enframe(value = "coef") %>%
            filter(str_detect(name, "^ns_"))
        },
        error = function(cond) {
          return(0L)
        },
        warning  =function(cond) {
          return(0L)
        }
      ) 
      
      if (!identical(out, 0L)){
        converge <- 1
      }
    }
    
   return(out)
  }
  
  map_dfr(1:1000, get_coefs,
          df_cox, mod_form,
          .id = "boot")
}

tic()
plan(multisession, workers = 8)
spline_res <- reg_grid %>%
  distinct(outcome, mod, low) %>%
  mutate(res = future_pmap(list(outcome, mod, low),
                           get_splines, 
                           .progress = TRUE,
                           .options = furrr_options(seed = NULL)))
future:::ClusterRegistry("stop")
toc()

get_splines_obs <- function(outcome, mod, low){
  df_cox <- df_cox %>%
    filter(crp >= !!low)
  
  mod_form <- glue_collapse(covars[[mod]], ' + ')
  mod_form <- glue("Surv(time_{outcome}, event_{outcome}) ~
                   ns_1 + ns_2 + ns_3 + {mod_form}")
  
  as.formula(mod_form) %>%
    coxph(df_cox) %>%
    coef() %>%
    enframe(value = "coef") %>%
    filter(str_detect(name, "^ns_"))
}

spline_obs <- reg_grid %>%
  distinct(outcome, mod, low) %>%
  mutate(res = pmap(list(outcome, mod, low),
                           get_splines_obs))

save(spline_res, spline_obs, df_s, get_ci,
     file = "Data/spline_results.Rdata")


# Delta Method
get_delta <- function(outcome, mod, low){
  df_reg <- df_cox %>%
    filter(crp >= !!low)
  
  mod_form <- glue("Surv(time_{outcome}, event_{outcome}) ~
                   ns_1 + ns_2 + ns_3 + {glue_collapse(covars[[mod]], ' + ')}")
  
  mod <- as.formula(mod_form) %>%
            coxph(df_reg)
  
  # Get Delta SE
  vars <- names(coef(mod))
  coefs <- coef(mod)[str_detect(vars, "^ns_")]
  vcovs <- vcov(mod)[str_detect(vars, "^ns_"),
                     str_detect(vars, "^ns_")]
  
  df_s %>%
    arrange(crp, name) %>%
    group_by(crp) %>%
    mutate(delta = glue("{value}*x{row_number()}") %>%
             glue_collapse(" + ") %>%
             paste0("~ exp(", ., ")")) %>%
    ungroup() %>%
    nest(lincom = c(name, value)) %>%
    mutate(lincom = map(lincom, deframe),
           delta = map(delta, as.formula)) %>%
    mutate(beta = map_dbl(lincom, ~ sum(.x * coefs) %>% exp()),
           se = deltamethod(delta, coefs, vcovs),
           lci = qnorm(.025, beta, se),
           uci = qnorm(.975, beta, se)) %>%
    select(-delta, -lincom)
}

delta_res <- reg_grid %>%
  distinct(outcome, mod, low) %>%
  mutate(res = pmap(list(outcome, mod, low), get_delta)) %>%
  unnest(res)

save(delta_res, df_s, file = "Data/delta_results.Rdata")
