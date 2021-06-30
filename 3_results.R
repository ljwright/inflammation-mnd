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

# 2. Descriptives ----
df %>%
  mutate(crp_int = cut_number(crp, 3)) %>%
  group_by(crp_int) %>%
  descr() %>%
  tb()

df %>%
  mutate(xx = cut_number(crp, 3)) %>%
  select(-matches("crp")) %>%
  drop_na()
  
  tbl_summary(by = x)


# 3. Reliability ----

# 4. Regression Results ----
