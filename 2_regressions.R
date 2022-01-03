library(tidyverse)
library(haven)
library(glue)
library(broom)
library(survival)
library(gallimaufr)
library(tictoc)
library(furrr)
library(splines)
library(msm)

rm(list = ls())

# 1. Load Data ----
load("Data/df_analysis.Rdata")
load("Data/model_covars.Rdata")

df_cox <- df %>%
  mutate(ns(crp, 3) %>%
           as_tibble() %>%
           mutate(across(everything(), as.numeric)) %>%
           rename_with(~ glue("ns_{.x}")),
         age = age / 10)

df_fu <- df %>%
  drop_na(crp, crp_fup)

df_obs <- map(covars,
              ~ df %>%
                drop_na(all_of(.x)) %>%
                select(id))

reg_grid <- expand_grid(outcome = c("hosp", "dead"),
                        crp = c("crp_3", "crp_3f", "log_crp", "crp_10f", "log_crp_lin"),
                        mod = names(covars),
                        sample = names(covars)[2:3]) %>%
  filter(!(mod == names(covars)[[3]] & sample == names(covars)[[2]])) %>%
  mutate(spec_id = row_number(), .before = 1)

get_ci <- function(estimates){
  quantile(estimates, c(.5, .025, .975)) %>%
    as_tibble_row() %>%
    rename(beta = 1, lci = 2, uci = 3)
}

save(df_cox, df_fu, covars, reg_grid, df_obs,
     file = "Data/df_regressions.Rdata")

rm(df)
gc()


# 2. Original Regressions ----
load("Data/df_regressions.Rdata")

get_cox <- function(spec_id){
  spec <- slice(reg_grid, !!spec_id)
  
  df_reg <- df_cox %>%
    semi_join(df_obs[[spec$sample]], by = "id")
  
  glue("Surv(time_{spec$outcome}, event_{spec$outcome}) ~
       {spec$crp} + {glue_collapse(covars[[spec$mod]], ' + ')}") %>%
    as.formula() %>%
    coxph(df_reg, method = "breslow") %>%
    tidy(conf.int = TRUE, conf.level = 0.95,
         exponentiate = TRUE) %>%
    select(term, beta = estimate, p = p.value,
           lci = conf.low, uci = conf.high)
}

# Run Models
main_res <- reg_grid %>%
  mutate(res = map(spec_id, get_cox)) %>%
  unnest(res)

main_res_covars <- main_res %>%
  filter(!str_detect(term, crp))

main_res <- main_res %>%
  filter(str_detect(term, crp))


# 3. Regression Dilution ----
# Model Function
get_boot <- function(spec_id, boot){
  spec <- slice(reg_grid, !!spec_id)
  
  df_fu <- df_fu %>%
    semi_join(df_obs[[spec$sample]], by = "id")
  
  df_cox <- df_cox %>%
    semi_join(df_obs[[spec$sample]], by = "id")
  
  mod_covars <- glue_collapse(covars[[spec$mod]], ' + ')
  
  set.seed(boot)
  converge <- 0L
  while (converge == 0L){
    out <- tryCatch(
      {
        boot_fu <- df_fu %>%
          sample_frac(replace = TRUE)
        
        boot_df <- df_cox %>%
          sample_frac(replace = TRUE)
        
        df_new <- glue("log_crp_fup ~ {spec$crp} + {mod_covars}") %>% 
          as.formula() %>%
          lm(boot_fu) %>%
          broom::augment(newdata = boot_df)
        
        mod <- glue("Surv(time_{spec$outcome}, event_{spec$outcome}) ~ 
       .fitted + {mod_covars}") %>%
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
  gc()
  return(out)
}


# Run Models
tic()
plan(multisession, workers = 3)
dilut_res <- reg_grid %>% 
  filter(crp %in% c("log_crp", "log_crp_lin")) %>%
  select(spec_id) %>%
  uncount(1000, .id = "boot") %>%
  sample_frac() %>%
  # filter(boot %in% 1:2) %>%
  # mutate(estimate = map2_dbl(spec_id, boot, get_boot)) %>%
  mutate(estimate = future_map2_dbl(spec_id, boot, get_boot,
                                    .progress = TRUE)) %>%
  arrange(spec_id) %>%
  group_by(spec_id) %>%
  summarise(get_ci(estimate),
            .groups = "drop") %>%
  right_join(reg_grid, ., by = "spec_id")
future:::ClusterRegistry("stop")
toc()

save(main_res, main_res_covars,
     dilut_res, covars,
     file = "Data/reg_results.Rdata")


# 4. Splines ----
splines_grid <- reg_grid %>%
  distinct(outcome, mod, sample) %>%
  mutate(spec_id = row_number(), .before = 1)

df_s <- df_cox %>%
  select(crp, matches("^ns_")) %>%
  arrange(crp) %>%
  distinct() %>%
  pivot_longer(-crp)

get_objects <- function(spec_id){
  spec <- slice(splines_grid, !!spec_id)
  
  df_cox <- df_cox %>%
    semi_join(df_obs[[spec$sample]], by = "id")
  
  mod_form <- glue_collapse(covars[[spec$mod]], ' + ')
  mod_form <- glue("Surv(time_{spec$outcome}, event_{spec$outcome}) ~
                   ns_1 + ns_2 + ns_3 + {mod_form}")
  
  list(df_cox = df_cox, mod_form = mod_form)
}

# Bootstraps
get_splines <- function(spec_id){
  objs <- get_objects(spec_id)
  df_cox <- objs$df_cox
  mod_form <- objs$mod_form
  
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
        warning = function(cond) {
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
          .id = "boot") %>%
    mutate(boot = as.integer(boot))
}

tic()
plan(multisession, workers = 3)
spline_res <- splines_grid %>%
  mutate(res = future_map(spec_id, get_splines, 
                           .progress = TRUE,
                           .options = furrr_options(seed = NULL)))
future:::ClusterRegistry("stop")
toc()

# Main Result
get_splines_obs <- function(spec_id){
  objs <- get_objects(spec_id)
  df_cox <- objs$df_cox
  mod_form <- objs$mod_form
  
  as.formula(mod_form) %>%
    coxph(df_cox) %>%
    coef() %>%
    enframe(value = "coef") %>%
    filter(str_detect(name, "^ns_"))
}

spline_obs <- splines_grid %>%
  mutate(res = map(spec_id, get_splines_obs))

save(spline_res, spline_obs, df_s, get_ci,
     file = "Data/spline_results.Rdata")

# Delta Method
get_delta <- function(spec_id){
  objs <- get_objects(spec_id)
  df_cox <- objs$df_cox
  mod_form <- objs$mod_form
  
  mod <- as.formula(mod_form) %>%
    coxph(df_cox)
  
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

delta_res <- splines_grid %>%
  mutate(res = map(spec_id, get_delta)) %>%
  unnest(res)

save(delta_res, df_s, file = "Data/delta_results.Rdata")


# 5. Prospective Analysis of CRP
load("Data/df_new.Rdata")

get_pros <- function(max_years){
  df_mnd <- df_new %>% 
    filter(baseline_diag == 1 | (event_hosp == 1 & time_hosp <= !!max_years)) %>%
    mutate(log_crp = log(crp) %>% wtd_scale()) %>%
    select(id, died, event_dead, time_dead, 
           log_crp, age, female) %>%
    drop_na()
  
  mod <- coxph(Surv(time_dead, event_dead) ~ log_crp + age + female,
                    df_mnd, method = "breslow")
  
  tidy(mod, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(nobs = glance(mod)$n,
           events = glance(mod)$nevent)
}

pros_res <- tibble(max_years = 1:10) %>%
  mutate(res = map(max_years, get_pros)) %>%
  unnest(res)
save(pros_res, file = "Data/prospective_results.Rdata")
