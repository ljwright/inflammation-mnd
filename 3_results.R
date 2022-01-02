library(tidyverse)
library(glue)
library(broom)
library(gallimaufr)
library(flextable)
library(officer)
library(infer)
library(magrittr)
library(labelled)
library(summarytools)

rm(list = ls())

# 1. Load Data ----
load("Data/df_analysis.Rdata")
load("Data/reg_results.Rdata")
load("Data/spline_results.Rdata")
load("Data/delta_results.Rdata")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

crp_clean <- df %>%
  distinct(crp_3f) %>%
  arrange(crp_3f) %>%
  mutate(crp3_clean = str_replace(crp_3f, "\\[", "(") %>%
           str_replace("\\]", ")"),
         crp3_clean = glue("Tertile {row_number()}\n{crp3_clean}"),
         term = glue("crp_3f{crp_3f}")) %>%
  select(group_var = crp_3f, term, crp3_clean)

crp_dict <- crp_clean %>%
  select(-group_var) %>%
  deframe()


# 2. Descriptives ----
# Get Columns
desc_vars <- c("crp_3f", "crp", "age", "female", "fev_1", "bmi", "non_white",
               "townsend", "disadvantaged", "smoke_status", "low_activity",
               "vasculardis", "diabetes", "cancer", "seen_psychiatrist",
               "hba1c", "grip_strength", "hdl")

factor_vars <- c("female", "non_white", "disadvantaged",
                 "low_activity", "vasculardis", "diabetes", "cancer",
                 "seen_psychiatrist")

df_desc <- df %>%
  select(all_of(desc_vars)) %>%
  mutate(across(all_of(factor_vars), as.factor))

pretty_lbls <- look_for(df) %>% 
  as_tibble() %>% 
  select(var = variable, label) %>%
  arrange(match(var, desc_vars)) %>%
  add_row(var = "n", label = "N",.before = 1)

# Follow-Up Time
df %>%
  select(time_hosp, time_dead) %>%
  descr() %>%
  tb()
  summarise(mean_hosp = mean(time_hosp),
            mean_dead = mean(time_dead))

df %>%
  group_by(female) %>%
  summarise(sum_hosp = sum(event_hosp),
            sum_dead = sum(event_dead))

df %>%
  summarise(sum_hosp = sum(event_hosp),
            sum_dead = sum(event_dead))


# Create Table
desc <- list()

get_p <- function(var){
  mod_form <- glue("{var} ~ crp_3f") %>%
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
  filter(var != "crp_3f") %>%
  mutate(p = map_dbl(var, get_p),
         p = ifelse(p < 0.0001, "p < 0.0001", round(p, 2)))


desc$df <- get_desc(df_desc, group_var = "crp_3f") %>%
  filter(var == cat | cat == "1" | var == "smoke_status") %>%
  left_join(crp_clean, by = "group_var") %>%
  select(crp3_clean, var, string) %>%
  pivot_wider(names_from = crp3_clean, values_from = string) %>%
  left_join(pretty_lbls, by = "var") %>%
  left_join(desc$p, by = "var") %>%
  mutate(Variable = factor(label, unique(pretty_lbls$label)) )%>%
  select(Variable, 2:4, `p-value` = p)  %>%
  arrange(Variable) %>%
  unnest(2:4) %>%
  group_by(Variable) %>%
  mutate(Variable = ifelse(Variable == "Smoking Status", 
                           glue("{levels(df$smoke_status)} Smoker"),
                           as.character(Variable)),
         `p-value` = ifelse(row_number() == 1, `p-value`, NA))

desc$flx <- make_flx(desc$df) %>%
  align(j = 2:4, align = "center", part = "all")
desc$flx
save_as_docx(desc$flx, path = "Tables/descriptives.docx")

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

expand_grid(low = c(-Inf, 0.15),
            method = c("pearson", "spearman")) %>%
  mutate(res = map2(low, method, get_cor)) %>%
  unnest(res)


# 4. Regression Plots ----
save_gg <- function(file_name, p, height = 9.9, width = 21){
  glue("Images/{file_name}.png") %>%
    ggsave(plot = p, height = height, width = width,
           units = "cm")
}

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
           fct_rev(),
         mod = ifelse(mod == "Age + Sex", "Age and sex adjusted", "Multiply-adjusted"))

plot_point <- function(df, var){
  ggplot(df) +
    aes(x = {{ var }}, y = beta, ymin = lci, ymax = uci,
        color = mod, shape = mod) +
    facet_wrap(~ dep_clean) + #, scales = "free_y") +
    geom_hline(yintercept = 1) +
    geom_pointrange(position = position_dodge(0.5)) +
    scale_y_log10() +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold")) +
    labs(x = NULL, y = "Hazard Rate Ratio(+ 95% CI)",
         color = NULL, shape = NULL)
}

# Figure 1
plot_1 <- function(low, save_p = FALSE){
  p <- lin_res %>%
    filter(crp == "crp_3f", 
           low == !!low) %>%
    uncount(ifelse(term == glue("crp_3f{levels(df$crp_3f)[2]}"), 2, 1), .id = "id") %>%
    mutate(crp3_clean = ifelse(id == 2, crp_dict[1], crp_dict[term]),
           across(c(beta, lci, uci),  ~ ifelse(id == 2, 1, .x)),
           dep_clean = fct_rev(dep_clean)) %>%
    plot_point(crp3_clean)
  
  if (save_p == TRUE){
    ifelse(low == 0.15, "fig1_low", "fig1_all") %>%
      save_gg(p)
  }
  
  return(p)
}

plot_1(-Inf, TRUE)
plot_1(0.15, TRUE)


# Figure 2
plot_2 <- function(low, crp, save_p = FALSE){
  p <- lin_res %>%
    filter(str_detect(crp, "log_crp"),
           low == !!low, 
           crp == !!crp) %>%
    mutate(dep_clean = fct_rev(dep_clean)) %>%
    plot_point(term_clean)
  
  if (save_p == TRUE){
    case_when(low == 0.15 & crp == "log_crp" ~ "fig2_low_all",
              low == -Inf & crp == "log_crp" ~ "fig2_all_all",
              low == 0.15 & crp == "log_crp_lin" ~ "fig2_low_lin",
              low == -Inf & crp == "log_crp_lin" ~ "fig2_all_lin") %>%
      save_gg(p)
  }
  
  return(p)
}

lin_res %>%
  filter(str_detect(crp, "log_crp")) %>%
  distinct(low, crp) %$%
  map2(low, crp, plot_2, TRUE)

# Figure 3
clean_splines <- function(df_res){
  df_res %>%
    unnest(res) %>%
    left_join(df_s, by = "name") %>%
    group_by(across(-c(name, coef, value))) %>%
    summarise(estimate = sum(coef*value),
              .groups = "drop") %>%
    group_by(across(-c(crp, estimate))) %>%
    mutate(id = cur_group_id()) %>%
    ungroup() %>%
    mutate(dep_clean = ifelse(outcome == "hosp", "Hospitalisation", "Death"),
           estimate = exp(estimate))
}

spline_obs <- clean_splines(spline_obs)

spline_res <- clean_splines(spline_res)

spline_ci <- spline_res %>%
  group_by(dep_clean, mod, low, crp) %>%
  summarise(get_ci(estimate),
            .groups = "drop")

plot_3 <- function(low, save_p = FALSE){
  p <- spline_ci %>%
    filter(!(beta == 1 & lci == 1),
           low == !!low) %>%
    mutate(dep_clean = fct_rev(dep_clean),
           mod = ifelse(mod == "Age + Sex", 
                        "Age and sex adjusted", 
                        "Multiply-adjusted")) %>%
    ggplot() +
    aes(x = crp, y = beta, ymin = lci, ymax = uci,
        color = mod, fill = mod) +
    facet_wrap(~ dep_clean) +
    geom_hline(yintercept = 1) +
    geom_ribbon(alpha = 0.2, color = NA) +
    geom_line() +
    scale_y_log10() +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(x = "Baseline CRP", y = "Hazard Ratio (+ 95% CI)",
         color = NULL, fill = NULL) +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "bold"))
  
  if (save_p == TRUE){
    ifelse(low == 0.15, "fig3_low", "fig3_all") %>%
      save_gg(p)
  }
  
  return(p)
}

plot_3(-Inf, TRUE)
plot_3(0.15, TRUE)

# Figure 4
spline_x <- df_s %>%
  distinct(crp) %>%
  slice(which(row_number() %% 10 == 1)) %>%
  left_join(spline_res, by = "crp")

rm(spline_res, clean_splines)

plot_4 <- function(low, save_p = FALSE){
  p <- spline_x %>%
    filter(low == !!low) %>%
    ggplot() +
    aes(x = crp, y = estimate, color = mod, fill = mod, group = id) +
    facet_wrap(~ dep_clean, scales = "free") +
    geom_hline(yintercept = 1) +
    geom_line(alpha = 0.1, size = 0.3) +
    geom_line(data = filter(spline_obs, low == !!low),
              aes(linetype = mod), size = 1) +
    scale_y_log10() +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(x = "Baseline CRP", y = "Hazard Ratio (+ 95% CI)",
         color = "Model", fill = "Model", linetype = "Model") +
    theme(legend.position = c(.1, .9)) +
    guides(color = guide_legend(override.aes = list(alpha = 1)))
  
  if (save_p == TRUE){
    ifelse(low == 0.15, "fig4_low", "fig4_all") %>%
      save_gg(p)
  }
  
  return(p)
}

plot_4(-Inf, TRUE)
plot_4(0.15, TRUE)

rm(spline_x, spline_ci, spline_obs)

# 5. Regression Tables
n_events <- df %>%
  select(all_of(covars[[2]]), matches("(time|event)"), crp_3f) %>%
  drop_na() %>%
  select(crp_3f, matches("event")) %>%
  pivot_longer(-crp_3f, names_to = "var") %>%
  count(crp_3f, var, value) %>%
  add_count(crp_3f, var, wt = n, name = "total") %>%
  filter(value == 1) %>%
  mutate(across(c(n, total), 
                ~ format(.x, big.mark = ",") %>%
                  trimws()),
         string = glue("{n} / {total}"),
         outcome = str_replace(var, "event_", ""),
         dep_clean = ifelse(outcome == "hosp", "Hospitalisation", "Death"),
         term = glue("crp_3f{crp_3f}"),
         mod = "# Events / N") %>%
  select(dep_clean, mod, term, string)

header_lbls <- crp_clean %>%
  mutate(term = glue("crp_3f{row_number()}")) %>%
  select(term, crp3_clean) %>% 
  deframe() %>% 
  as.list() %>%
  c(list(dep_clean = "Outcome",
         mod = "",
         crp_3 = "p-value for trend",
         log_crp = "Per 1 SD increase\n(Observed)",
         dilut = "Per 1 SD increase\n(Dilution Corrected)",
         log_crp_lin = "Per 1 SD increase\n(Observed)",
         dilut_lin = "Per 1 SD increase\n(Dilution Corrected)"))

get_tbl <- function(low, lin, file_name = NULL){
  if (lin){
    drop_vars <- c("log_crp", "dilut")
  } else{
    drop_vars <- c("log_crp_lin", "dilut_lin")
  }
  
  tbl <- lin_res %>%
    filter(crp != "crp10", low == !!low) %>%
    uncount(ifelse(term == glue("crp_3f{levels(df$crp_3f)[2]}"), 2, 1), .id = "id") %>%
    mutate(term = ifelse(id == 2, crp_clean$term[1], term),
           across(c(beta, lci, uci),  ~ ifelse(id == 2, 1, .x)),
           across(beta:uci, round, 2),
           string = ifelse(term == "crp_3", p, glue("{beta} ({lci}, {uci})"))) %>%
    select(dep_clean, mod, term, string) %>%
    bind_rows(n_events) %>% 
    mutate(term = ifelse(str_detect(term, "^crp_3f"), 
                         glue("crp_3f{match(term, crp_clean$term)}"),
                         term)) %>%
    pivot_wider(names_from = term, values_from = string) %>%
    arrange(desc(dep_clean), mod) %>%
    select(dep_clean, mod, crp_3f1, crp_3f2, crp_3f3,
           crp_3, log_crp, dilut, log_crp_lin, dilut_lin) %>%
    select(-all_of(drop_vars)) %>%
    make_flx(header_lbls)
  
  if (!is.null(file_name)){
    glue("Tables/{file_name}.docx") %>%
      save_as_docx(tbl, path = .)
  }
  
  return(tbl)
}

expand_grid(low = c(0.15, -Inf),
            lin = c(TRUE, FALSE)) %>%
  mutate(file_name = case_when(low == 0.15 & lin == FALSE ~ "tbl2_low_all",
                               low == -Inf & lin == FALSE ~ "tbl2_all_all",
                               low == 0.15 & lin == TRUE ~ "tbl2_low_lin",
                               low == -Inf & lin == TRUE ~ "tbl2_all_lin")) %$%
  pmap(list(low, lin, file_name), get_tbl)
