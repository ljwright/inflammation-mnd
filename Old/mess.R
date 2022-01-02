# 2b. Stratified by Sex ----
get_cox_sex <- function(outcome, crp, female, low){
  df_reg <- df_cox %>%
    filter(crp > !!low, 
           female == !!female)
  
  glue("Surv(time_{outcome}, event_{outcome}) ~ {crp} + age") %>%
    as.formula() %>%
    coxph(df_reg, method = "breslow") %>%
    tidy(conf.int = TRUE, conf.level = 0.95, 
         exponentiate = TRUE) %>%
    select(term, beta = estimate, p = p.value,
           lci = conf.low, uci = conf.high)
}

sex_res <- reg_grid %>%
  select(-mod) %>%
  distinct() %>%
  expand_grid(female = c(0, 1)) %>%
  mutate(res = pmap(list(outcome, crp, female, low), get_cox_sex)) %>%
  unnest(res) %>%
  filter(str_detect(term, crp))

sex_res %>%
  filter(term == "log_crp", low == -Inf) %>%
  ggplot() +
  aes(x = female, y = beta, ymin = lci, ymax = uci) +
  facet_wrap(~ outcome, scales = "free_y") +
  geom_pointrange()