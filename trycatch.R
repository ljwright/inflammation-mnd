get_cox_splines <- function(outcome, mod, low, covars, df_s){
}

outcome <- "hosp"
mod <- "Age + Sex"
low <- 0.1

df_cox <- df_cox %>%
  filter(crp >= !!low)

mod_form <- glue("Surv(time_{outcome}, event_{outcome}) ~
                   {glue_collapse(unique(df_s$name), ' + ')} +
              {glue_collapse(covars[[mod]], ' + ')}")

get_coefs <- function(boot, df_cox, df_s, mod_form){
  set.seed(boot)
  
  converg <- 0
  
  # while (converg == 0){
  df_reg <- sample_frac(df_cox, replace = TRUE)
  
  mod <- as.formula(mod_form) %>%
    coxph(df_reg) %>%
    coef() %>%
    enframe(value = "coef") %>%
    right_join(df_s, by = "name") %>%
    group_by(crp) %>%
    summarise(estimate = sum(value*coef),
              .groups = "drop")
  # }
  
  
}


catch_boot <- function(boot){
  set.seed(boot)
  
  converge <- 0L
  i <- 0
  
  while (converge == 0){
    i <- i + 1
    out <- tryCatch(
      {
        df_reg <- sample_frac(df_cox, replace = TRUE)
        
        as.formula(mod_form) %>%
          coxph(df_reg) %>% 
          coef()
      },
      error=function(cond) {
        return(0L)
      },
      warning=function(cond) {
        return(0L)
      }
    ) 
    
    if (!identical(out, 0L)){
      converge <- 1
    }
  }
   
  return(i)
}


run_cox <- function(boot){
  df_reg <- sample_frac(df_cox, size = 1, replace = TRUE)
  
  as.formula(mod_form) %>%
    coxph(df_reg) %>% 
    coef()
}

map(1:100, run_cox)

x <- safely(run_cox)

x$fail
run_cox2 <- function(boot){
  set.seed(boot)
  
  converg <- 0
  
  # while (converg == 0){
  df_reg <- sample_frac(df_cox, replace = TRUE)
  
  res <- as.formula(mod_form) %>%
    coxph(df_reg) %>%
    coef() %>%
    enframe(value = "coef") %>%
    right_join(df_s, by = "name") %>%
    group_by(crp) %>%
    summarise(estimate = sum(value*coef),
              .groups = "drop")
}


x <- map(100:500, run_cox)


x <- map_dfr(101:, run_cox2, .id = "boot") %>%
  group_by(crp) %>%
  summarise(get_ci(estimate),
            .groups = "drop") %>%
  mutate(across(beta:uci, exp)) %>%
  ggplot() +
  aes(x = crp, y = beta, ymin = lci, ymax = uci) +
  geom_ribbon()


get_coefs <- function(boot, df_cox, df_s, mod_form){
  set.seed(boot)
  
  converge <- 0
  
  while (converge == 0){
    df_reg <- sample_frac(df_cox, replace = TRUE)
    
    mod <- as.formula(mod_form) %>%
      coxph(df_reg)
    
    if (is.null(mod$fail)){
      converge <- 1
    }
  }
  
  coef(mod) %>%
    enframe(value = "coef") %>%
    right_join(df_s, by = "name") %>%
    group_by(crp) %>%
    summarise(estimate = sum(value*coef),
              .groups = "drop")
}


