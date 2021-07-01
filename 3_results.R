library(tidyverse)
library(glue)
library(broom)
library(gallimaufr)
library(flextable)
library(officer)
library(infer)
library(magrittr)
library(labelled)

rm(list = ls())

# 1. Load Data ----
load("Data/df_analysis.Rdata")
load("Data/reg_results.Rdata")
load("Data/spline_results.Rdata")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

crp_clean <- df %>%
  filter(inflam == 1) %>%
  group_by(rev_crp3) %>%
  summarise(min = min(crp, na.rm = TRUE),
            max = max(crp, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(across(c(min, max), round, 2),
         crp3_clean = glue("Tertile {rev_crp3+1}\n({min}, {max})"),
         term = glue("rev_crp3_f{rev_crp3}")) %>%
  select(group_var = rev_crp3, term, crp3_clean)

# Questions: Why are descriptives and controls different?


# 2. Descriptives ----
# Get Columns
desc_vars <- c("rev_crp3", "age", "female", "fev1", "height", "BMI",
               "nonwhite", "disadvantaged", "currsmok", "low_activity",
               "vasculardis", "diabetes", "cancer", "seen_psychiatrist")

factor_vars <- c("female", "nonwhite", "disadvantaged", "currsmok",
                 "low_activity", "vasculardis", "diabetes", "cancer",
                 "seen_psychiatrist")

df_desc <- df %>%
  filter(inflam == 1) %>%
  select(all_of(desc_vars)) %>%
  mutate(across(all_of(factor_vars), as.factor))

pretty_lbls <- look_for(df) %>% 
  as_tibble() %>% 
  select(var = variable, label) %>%
  arrange(match(var, desc_vars)) %>%
  add_row(var = "n", label = "N",.before = 1)


# Create Table
desc <- list()


get_p <- function(var){
  df_desc$rev_crp3 <- factor(df_desc$rev_crp3)
  
  mod_form <- glue("{var} ~ rev_crp3") %>%
    as.formula()
  
  if (is.numeric(df_desc[[var]])){
    p <- aov(mod_form, df_desc) %>%
      tidy() %>%
      slice(1) %>%
      pull(p.value)
  } else{
    p <- chisq_test(df_desc, mod_form) %>%
      pull(p_value)
  }
  
  return(p)
}

desc$p <- tibble(var = names(df_desc)) %>%
  filter(var != "rev_crp3") %>%
  mutate(p = map_dbl(var, get_p),
         p = ifelse(p < 0.0001, "p < 0.0001", round(p, 2)))


desc$df <- get_desc(df_desc, group_var = "rev_crp3") %>%
  filter(var == cat | cat == "1") %>%
  left_join(crp_clean, by = "group_var") %>%
  select(crp3_clean, var, string) %>%
  pivot_wider(names_from = crp3_clean, values_from = string) %>%
  left_join(pretty_lbls, by = "var") %>%
  left_join(desc$p, by = "var") %>%
  mutate(Variable = factor(label, unique(pretty_lbls$label))) %>%
  select(Variable, 2:4, `p-value` = p)  %>%
  arrange(Variable)

desc$flx <- make_flx(desc$df) %>%
  align(j = 2:4, align = "center", part = "all")

# 3. Reliability ----
get_cor <- function(low, method){
  data <- df %>%
    filter(crp > !!low, 
           crp_fup > !!low) %>%
    select(crp, crp_fup)
  
  res <- cor.test(data$crp, data$crp_fup, method = method)
  
  tibble(rho = res$estimate, p = res$p.value,
         n = nrow(data))
}

expand_grid(low = c(-Inf, 0.1),
            method = c("pearson", "spearman")) %>%
  mutate(res = map2(low, method, get_cor)) %>%
  unnest(res)


# 4. Regression Plots ----
lin_res <- bind_rows(
  main = main_res,
  dilut = dilut_res %>%
    ungroup() %>%
    mutate(across(beta:uci, exp),
           term = crp),
  .id = "type"
) %>%
  mutate(dep_clean = ifelse(outcome == "hosp", "Hospitalisation", "Death"),
         term = ifelse(type == "dilut", 
                       str_replace(term, "log_crp", "dilut"),
                       term),
         term_clean = ifelse(type == "dilut", "Dilution Corrected", "Observed CRP") %>%
           fct_rev())

plot_point <- function(df, var){
  ggplot(df) +
    aes(x = {{ var }}, y = beta, ymin = lci, ymax = uci,
        color = mod) +
    facet_wrap(~ dep_clean, scales = "free_y") +
    geom_hline(yintercept = 1) +
    geom_pointrange(position = position_dodge(0.5)) +
    scale_y_log10() +
    scale_color_brewer(palette = "Set1") +
    theme_minimal() +
    theme(legend.position = c(0.1, 0.9)) +
    labs(x = NULL, y = "Hazard Rate Ratio(+ 95% CI)",
         color = "Model")
}

# Figure 1
plot_1 <- function(low){
  lin_res %>%
    filter(crp == "rev_crp3_f", 
           low == !!low) %>%
    uncount(ifelse(term == "rev_crp3_f1", 2, 1), .id = "id") %>%
    mutate(tertile = ifelse(id == 2, 0, as.numeric(str_sub(term, -1))),
           across(c(beta, lci, uci),  ~ ifelse(tertile == 0, 1, .x))) %>%
    left_join(crp_clean, by = c("tertile" = "group_var")) %>%
    plot_point(crp3_clean)
}

plot_1(-Inf)
plot_1(0.1)


# Figure 2
plot_2 <- function(low, crp){
  lin_res %>%
    filter(str_detect(crp, "log_crp"),
           low == !!low, 
           crp == !!crp) %>%
    plot_point(term_clean)
}

lin_res %>%
  filter(str_detect(crp, "log_crp")) %>%
  distinct(low, crp) %$%
  map2(low, crp, plot_2)

# Figure 3
## SPLINES!!
plot_3 <- function(low){
  spline_res %>%
    mutate(dep_clean = ifelse(outcome == "hosp", "Hospitalisation", "Death"),
           across(beta:uci, exp)) %>%
    filter(!(beta == 0 & lci == 0),
           low == !!low) %>%
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
}

plot_3(-Inf)
ggsave("Images/splines_plot_inf.png",
       width = 21, height = 9.9, units = "cm")

plot_3(0.1)
ggsave("Images/splines_plot_01.png",
       width = 21, height = 9.9, units = "cm")


# 5. Regression Tables
n_events <- df %>%
  filter(inflam == 1) %>%
  select(all_of(covars[[2]]), matches("(time|event)"), rev_crp3) %>%
  drop_na() %>%
  select(rev_crp3, matches("event")) %>%
  pivot_longer(-rev_crp3, names_to = "var") %>%
  count(rev_crp3, var, value) %>%
  add_count(rev_crp3, var, wt = n, name = "total") %>%
  filter(value == 1) %>%
  mutate(across(c(n, total), 
                ~ format(.x, big.mark = ",") %>%
                  trimws()),
         string = glue("{n} / {total}"),
         outcome = str_replace(var, "event_", ""),
         dep_clean = ifelse(outcome == "hosp", "Hospitalisation", "Death"),
         term = glue("rev_crp3_f{rev_crp3}"),
         mod = "# Events / N") %>%
  select(dep_clean, mod, term, string)

header_lbls <- crp_clean %>%
  select(2:3) %>% 
  deframe() %>% 
  as.list() %>%
  c(list(dep_clean = "Outcome",
         mod = "",
         rev_crp3 = "p-value for trend",
         log_crp = "Per 1 SD increase\n(Observed)",
         dilut = "Per 1 SD increase\n(Dilution Corrected)",
         log_crp_lin = "Per 1 SD increase\n(Observed)",
         dilut_lin = "Per 1 SD increase\n(Dilution Corrected)"))

get_tbl <- function(low, lin){
  if (lin){
    drop_vars <- c("log_crp", "dilut")
  } else{
    drop_vars <- c("log_crp_lin", "dilut_lin")
  }
  
  lin_res %>%
    filter(crp != "crp10", low == !!low) %>%
    uncount(ifelse(term == "rev_crp3_f1", 2, 1), .id = "id") %>%
    mutate(term = ifelse(id == 2, "rev_crp3_f0", term),
           across(c(beta, lci, uci), 
                  ~ ifelse(term == "rev_crp3_f0", 1, .x)),
           across(beta:uci, round, 2),
           string = ifelse(term == "rev_crp3", p, glue("{beta} ({lci}, {uci})"))) %>%
    select(dep_clean, mod, term, string) %>%
    bind_rows(n_events) %>%
    pivot_wider(names_from = term, values_from = string) %>%
    arrange(desc(dep_clean), mod) %>%
    select(dep_clean, mod, rev_crp3_f0, rev_crp3_f1, rev_crp3_f2,
           rev_crp3, log_crp, dilut, log_crp_lin, dilut_lin) %>%
    select(-all_of(drop_vars)) %>%
    make_flx(header_lbls)
}

expand_grid(low = c(0.1, -Inf),
            lin = c(TRUE, FALSE)) %$%
  map2(low, lin, get_tbl)
