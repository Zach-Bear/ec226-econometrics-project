# =============================================================================
# EC226 ECONOMETRICS 1 — GROUP 17
# Does Homeownership Increase Life Satisfaction? Evidence from the UK, 2013-2022
# =============================================================================
#
# REPLICATION INSTRUCTIONS
#
# Quickstart (recommended):
#   1. Update the setwd() path below to the folder containing ec226-group017.RData
#   2. Run the load("ec226-group017.RData") line at the top of Section 1
#   3. Run from Section 1 onwards to replicate all results
#   The .RData file contains df_panel (fully cleaned and merged analytical
#   dataset) and all fitted model objects.
#
# Full replication from raw data:
#   The raw UKHLS wave files (SN 6614) must be downloaded from the UK Data
#   Service. Place all .dta files and Average-prices-2025-01.csv in the
#   folder referenced by setwd(), then run the script from the top.
#
# Structure:
#   Section 0    Raw data loading, cleaning, and HPI merge (skip if loading .RData)
#   Section 1    Panel structure declaration        <- Start here with .RData, line 290
#   Section 2    Summary statistics
#   Section 3    Regression models (OLS, FE, COVID interaction, age split)
#   Section 4    Diagnostic tests (pFtest, Hausman, Breusch-Godfrey)
#   Section 5    Instrumental variables (first stage and 2SLS)
#   Section 6    Output tables
#   Section 7    Robustness checks (GHQ outcome, tenure changers)
#
# Packages required:
#   haven, dplyr, ggplot2, plm, lmtest, sandwich, car, AER,
#   stargazer, readr, lubridate
#
# =============================================================================

rm(list = ls())

setwd("")  # <-- UPDATE THIS

library(haven)
library(dplyr)
library(ggplot2)
library(plm)
library(lmtest)
library(sandwich)
library(car)
library(AER)
library(stargazer)
library(readr)
library(lubridate)



  

# LOAD AND TRIM EACH WAVE
# Wave letters: e=2013, f=2014, g=2015, h=2016, i=2017, j=2018,
#               k=2019, l=2020, m=2021, n=2022

load_wave <- function(prefix) {
  
  # Load individual file and select variables, renaming to remove wave prefix
  indresp <- read_dta(paste0(prefix, "_indresp.dta"))
  indresp_small <- indresp %>% select(
    pidp,
    hidp        = paste0(prefix, "_hidp"),
    sclfsato    = paste0(prefix, "_sclfsato"),
    scghq1_dv   = paste0(prefix, "_scghq1_dv"),
    dvage       = paste0(prefix, "_dvage"),
    sex         = paste0(prefix, "_sex"),
    mastat_dv   = paste0(prefix, "_mastat_dv"),
    nchild_dv   = paste0(prefix, "_nchild_dv"),
    hiqual_dv   = paste0(prefix, "_hiqual_dv"),
    jbstat      = paste0(prefix, "_jbstat"),
    fimnnet_dv  = paste0(prefix, "_fimnnet_dv"),
    finnow      = paste0(prefix, "_finnow"),
    urban_dv    = paste0(prefix, "_urban_dv"),
    gor_dv      = paste0(prefix, "_gor_dv")
  ) %>%
    mutate(wave = prefix)  # Add wave label so we know which year each row is from
  
  # Load household file and select variables, renaming to remove wave prefix
  hhresp <- read_dta(paste0(prefix, "_hhresp.dta"))
  hhresp_small <- hhresp %>% select(
    hidp       = paste0(prefix, "_hidp"),
    tenure_dv  = paste0(prefix, "_tenure_dv"),
    hsrooms    = paste0(prefix, "_hsrooms"),
    xpmg       = paste0(prefix, "_xpmg"),
    rent       = paste0(prefix, "_rent")
  )
  
  # Merge individual and household data
  left_join(indresp_small, hhresp_small, by = "hidp")
}

# Load all waves and stack into one long dataframe
# Each person can appear up to 10 times (once per wave)
df_panel <- bind_rows(
  load_wave("e"),
  load_wave("f"),
  load_wave("g"),
  load_wave("h"),
  load_wave("i"),
  load_wave("j"),
  load_wave("k"),
  load_wave("l"),
  load_wave("m"),
  load_wave("n")
)

# CLEANING

# Recode negatives to NA (UKHLS uses negative codes for missing data)
df_panel <- df_panel %>%
  mutate(across(where(is.numeric), ~ifelse(. < 0, NA, .)))

# Collapse tenure into four meaningful categories
# Social rent combines local authority and housing association (standard in literature)
# Codes 5 (employer rent) and 8 (other) set to NA - too small to interpret
df_panel <- df_panel %>%
  mutate(
    tenure = case_when(
      tenure_dv == 1 ~ "Owned outright",
      tenure_dv == 2 ~ "Mortgage",
      tenure_dv %in% c(3, 4) ~ "Social rent",
      tenure_dv == 6 ~ "Private rent",
      TRUE ~ NA_character_
    ),
    # Numeric ownership dummy for regressions (1 = any ownership, 0 = renting)
    owner = ifelse(tenure_dv %in% c(1, 2), 1, 0),
    # Age group for heterogeneity analysis (under 40 vs 40+)
    age_group = ifelse(dvage < 40, "Under 40", "40 and over"),
    # COVID dummy (waves l and m cover 2020-2022)
    covid = ifelse(wave %in% c("l", "m"), 1, 0),
    # Approximate calendar year per wave (for merging house price data later)
    year = case_when(
      wave == "e" ~ 2013,
      wave == "f" ~ 2014,
      wave == "g" ~ 2015,
      wave == "h" ~ 2016,
      wave == "i" ~ 2017,
      wave == "j" ~ 2018,
      wave == "k" ~ 2019,
      wave == "l" ~ 2020,
      wave == "m" ~ 2021,
      wave == "n" ~ 2022
    )
  )

# SAVE

# save(df_panel, file = "housing_happiness_panel.RData")

# VERIFICATION CHECKS

nrow(df_panel)               # Total observations 
length(unique(df_panel$pidp)) # Number of unique individuals

table(df_panel$wave)          # Observations per wave
table(df_panel$tenure)        # Tenure categories across all waves
table(df_panel$sclfsato)      # Life satisfaction distribution

summary(df_panel[c("sclfsato", "tenure_dv", "dvage", "fimnnet_dv")])

# PLOT - Mean life satisfaction by tenure (across all waves)
# From project report

df_panel %>%
  filter(!is.na(tenure), !is.na(sclfsato)) %>%
  group_by(tenure) %>%
  summarise(mean_satisfaction = mean(sclfsato, na.rm = TRUE)) %>%
  ggplot(aes(x = reorder(tenure, mean_satisfaction), y = mean_satisfaction, fill = tenure)) +
  geom_col(show.legend = FALSE) +
  coord_cartesian(ylim = c(4, 6)) +
  labs(title = "Mean Life Satisfaction by Housing Tenure",
       subtitle = "Understanding Society, Waves e-n (2013-2022)",
       x = "Tenure", y = "Mean Life Satisfaction (1-7)") +
  theme_minimal()




# SECTION 0 cont.: PREPARE VARIABLES FOR REGRESSION

# Strip haven labels from key variables to avoid merge/type issues
df_panel$gor_dv  <- as.numeric(df_panel$gor_dv)
df_panel$year    <- as.numeric(df_panel$year)

# Set reference categories explicitly
# Private rent is first so this is the baseline in the factor variable
df_panel <- df_panel %>%
  mutate(
    tenure = factor(tenure, levels = c("Private rent", "Social rent",
                                       "Mortgage", "Owned outright")),
    
    # Binary ownership dummy for IV (1 = any owner, 0 = any renter)
    # The IV compares the change when people buy a house, so no need for 4 categories
    # Note to self on how this works for future dummies:
    # If value of tenure is in the vector 1,2, assign 1, else 0
    owner = ifelse(tenure_dv %in% c(1, 2), 1, 0),
    
    # Sex dummy (2 = female in UKHLS)
    female = ifelse(sex == 2, 1, 0),
    
    # Partnered dummy (1 = married, 2 = cohabiting)
    partnered = ifelse(mastat_dv %in% c(1, 2), 1, 0),
    
    # Employment dummies - reference category is employed
    unemployed = ifelse(jbstat == 3, 1, 0),
    retired    = ifelse(jbstat == 4, 1, 0),
    sick       = ifelse(jbstat == 8, 1, 0),
    
    # Log income - reduces influence of very high earners
    # Aus study squares and cubes income, Germany and Korea studies log it
    log_income = log(fimnnet_dv + 1),
    
    # Age squared - captures U-shape of wellbeing over life course
    # Explained in Korea study, used in all 3
    dvage_sq = dvage^2,
    
    # COVID dummy
    covid = ifelse(wave %in% c("l", "m"), 1, 0),
    
    # Age group for heterogeneity analysis
    age_group = ifelse(dvage < 40, "Under 40", "40 and over"),
    
    # Calendar year per wave for HPI merge
    year = case_when(
      wave == "e" ~ 2013,
      wave == "f" ~ 2014,
      wave == "g" ~ 2015,
      wave == "h" ~ 2016,
      wave == "i" ~ 2017,
      wave == "j" ~ 2018,
      wave == "k" ~ 2019,
      wave == "l" ~ 2020,
      wave == "m" ~ 2021,
      wave == "n" ~ 2022
    )
  )

# SECTION 0 cont.: MERGE HOUSE PRICE INDEX DATA

hpi_raw <- read_csv("Average-prices-2025-01.csv")

# Filter to January observation per year per region (annual snapshot)
hpi_clean <- hpi_raw %>%
  mutate(year = year(Date), month = month(Date)) %>%
  filter(
    month == 1,
    year >= 2013, year <= 2022,
    Region_Name %in% c(
      "North East", "North West", "Yorkshire and The Humber",
      "East Midlands", "West Midlands", "East of England",
      "London", "South East", "South West",
      "Wales", "Scotland", "Northern Ireland"
    )
  ) %>%
  select(year, Region_Name, Average_Price)

# Lookup table matching region names to gor_dv codes used in Understanding Society
region_lookup <- data.frame(
  gor_dv = 1:12,
  Region_Name = c(
    "North East", "North West", "Yorkshire and The Humber",
    "East Midlands", "West Midlands", "East of England",
    "London", "South East", "South West",
    "Wales", "Scotland", "Northern Ireland"
  )
)

hpi_clean <- left_join(hpi_clean, region_lookup, by = "Region_Name")

# Merge HPI into panel
df_panel <- left_join(df_panel, hpi_clean, by = c("gor_dv", "year"))

# Log transform house prices
# Same reason as income
df_panel$Average_Price <- as.numeric(df_panel$Average_Price)
df_panel$log_hpi <- log(df_panel$Average_Price)

# Save progress here
save(df_panel, file = "housing_happiness_merged_ready.RData")
load("housing_happiness_merged_ready.RData")  # Uncomment to restart from this point





# SECTION 1: DECLARE PANEL STRUCTURE
# Load submitted data file to skip Section 0:
load("ec226-group017.RData")


panel <- pdata.frame(df_panel, index = c("pidp", "wave"))

# How many people appear in multiple waves? Need within-person variation
obs_per_person <- df_panel %>%
  group_by(pidp) %>%
  summarise(n_waves = n_distinct(wave))
table(obs_per_person$n_waves)


# SECTION 2: Summary Statistics
# For appendices

# Basic summary of key variables
summary(df_panel[c("sclfsato", "dvage", "fimnnet_dv", "nchild_dv")])

# Tenure distribution
prop.table(table(df_panel$tenure, useNA = "no")) * 100

# Employment status counts
prop.table(table(df_panel$jbstat, useNA = "no")) * 100

# Sex
prop.table(table(df_panel$female, useNA = "no")) * 100

# Partnered
prop.table(table(df_panel$partnered, useNA = "no")) * 100

# Life satisfaction distribution (for the text description)
table(df_panel$sclfsato)

# Standard deviations separately since summary() doesn't show them
sd(df_panel$sclfsato, na.rm = TRUE)
sd(df_panel$dvage, na.rm = TRUE)
sd(df_panel$fimnnet_dv, na.rm = TRUE)
sd(df_panel$nchild_dv, na.rm = TRUE)

# Panel structure
length(unique(df_panel$pidp))   # Unique individuals
nrow(df_panel)                   # Total observations
table(df_panel$wave)             # Obs per wave

# SECTION 3: REGRESSION MODELS

# MODEL 1: Pooled OLS - baseline for later comparison

model_ols <- lm(
  sclfsato ~ tenure + dvage + dvage_sq + female + log_income +
    partnered + nchild_dv + unemployed + retired + sick +
    covid + factor(gor_dv) + factor(wave),
  data = df_panel,
  na.action = na.omit
)

summary(model_ols)

# MODEL 2: Fixed Effects - main result
# Removes all time-invariant unobservables
# Identified by people who changed tenure between waves

model_fe <- plm(
  sclfsato ~ tenure + dvage + dvage_sq + log_income +
    partnered + nchild_dv + unemployed + retired + sick + covid,
  data = panel,
  model = "within",
  effect = "twoways",
  na.action = na.omit
)

summary(model_fe)

# MODEL 3: Fixed Effects with COVID interaction
# Looking for interactions based on COVID hypothesis
# Did COVID ownership "buffer" life satisfaction loss?
# This did not end up being major in the paper due to word count and the null result
# Being less interesting than age heterogeneity findings
# Still an interesting finding though, so left in the paper

model_fe_covid <- plm(
  sclfsato ~ tenure * covid + dvage + dvage_sq + log_income +
    partnered + nchild_dv + unemployed + retired + sick,
  data = panel,
  model = "within",
  effect = "twoways",
  na.action = na.omit
)

summary(model_fe_covid)

# MODEL 4: Age heterogeneity - split sample under/over 40
# Under-40s are the generation "priced out of buying"
# Hypothesis expects smaller/no ownership effect for younger people

panel_u40 <- pdata.frame(filter(df_panel, age_group == "Under 40"),
                         index = c("pidp", "wave"))
panel_o40 <- pdata.frame(filter(df_panel, age_group == "40 and over"),
                         index = c("pidp", "wave"))

model_fe_u40 <- plm(
  sclfsato ~ tenure + dvage + dvage_sq + log_income +
    partnered + nchild_dv + unemployed + retired + sick + covid,
  data = panel_u40,
  model = "within",
  effect = "twoways",
  na.action = na.omit
)

model_fe_o40 <- plm(
  sclfsato ~ tenure + dvage + dvage_sq + log_income +
    partnered + nchild_dv + unemployed + retired + sick + covid,
  data = panel_o40,
  model = "within",
  effect = "twoways",
  na.action = na.omit
)

summary(model_fe_u40)
summary(model_fe_o40)


# SECTION 4: DIAGNOSTIC TESTS

# Test 1: Fixed effects vs pooled OLS
# Null: pooled OLS is fine. If p < 0.05, FE is justified
pFtest(model_fe, model_ols)

# Test 2: Fixed vs random effects (Hausman test)
# Null: random effects is consistent. If p < 0.05, use fixed effects
model_re <- plm(
  sclfsato ~ tenure + dvage + dvage_sq + log_income +
    partnered + nchild_dv + unemployed + retired + sick + covid,
  data = panel,
  model = "random",
  na.action = na.omit
)

phtest(model_fe, model_re)

# Test 3: Serial correlation
pbgtest(model_fe)

# SECTION 5: IV ESTIMATION
# Using log regional house price for instrument
# Turns out to be quite weak

# FIRST STAGE
# Does regional house price predict ownership?

first_stage <- plm(
  owner ~ log_hpi + dvage + dvage_sq + log_income +
    partnered + nchild_dv + unemployed + retired + sick + covid,
  data = panel,
  model = "within",
  effect = "twoways",
  na.action = na.omit
)

summary(first_stage)

# 2nd stage: IV regression/2SLS
# Syntax: outcome ~ endogenous + controls | instrument + controls
# Everything after | is the instrument set, exogenous controls appear on both sides

model_iv <- ivreg(
  sclfsato ~ owner + dvage + dvage_sq + log_income +
    partnered + nchild_dv + unemployed + retired + sick + covid +
    factor(gor_dv) + factor(wave)
  | log_hpi + dvage + dvage_sq + log_income +
    partnered + nchild_dv + unemployed + retired + sick + covid +
    factor(gor_dv) + factor(wave),
  data = df_panel,
  na.action = na.omit
)

summary(model_iv, diagnostics = TRUE)

# Clustered standard errors for IV ?
coeftest(model_iv, vcov = sandwich(model_iv))

# SECTION 6: COMPARISON TABLE
# OLS vs FE vs IV - the full story of bias correction in one table

# Simple binary owner dummy versions of OLS and FE for direct comparison with IV
model_ols_simple <- lm(
  sclfsato ~ owner + dvage + dvage_sq + female + log_income +
    partnered + nchild_dv + unemployed + retired + sick + covid +
    factor(gor_dv) + factor(wave),
  data = df_panel,
  na.action = na.omit
)

model_fe_simple <- plm(
  sclfsato ~ owner + dvage + dvage_sq + log_income +
    partnered + nchild_dv + unemployed + retired + sick + covid,
  data = panel,
  model = "within",
  effect = "twoways",
  na.action = na.omit
)

stargazer(model_ols_simple, model_fe_simple, model_iv,
          type = "text",
          title = "OLS vs FE vs IV: Effect of Homeownership on Life Satisfaction",
          column.labels = c("OLS", "Fixed Effects", "IV"),
          keep = c("owner", "log_hpi"),
          omit.stat = c("f", "ser"),
          notes = "Instrument: log regional house price index. Clustered SEs.")

# Full results table for the paper
stargazer(model_ols, model_fe, model_fe_covid,
          type = "text",
          title = "Housing Tenure and Life Satisfaction, UK 2013-2022",
          column.labels = c("OLS", "FE", "FE x COVID"),
          omit = c("gor_dv", "wave"),
          omit.stat = c("f", "ser"),
          notes = "")


# ROBUSTNESS 1: Alternative measure - GHQ mental health score
# Higher GHQ = more distress, -ve coefficient on ownership expected?


model_fe_ghq <- plm(
  scghq1_dv ~ tenure + dvage + dvage_sq + log_income +
    partnered + nchild_dv + unemployed + retired + sick + covid,
  data = panel,
  model = "within",
  effect = "twoways",
  na.action = na.omit
)

summary(model_fe_ghq)

# ROBUSTNESS 2: Tenure changers only
# Restricts sample to individuals who changed tenure at least once
# Tighter test of within-person identification

# Identify people who changed tenure
tenure_changers <- df_panel %>%
  filter(!is.na(tenure_dv)) %>%
  group_by(pidp) %>%
  summarise(n_tenure = n_distinct(tenure_dv, na.rm = TRUE)) %>%
  filter(n_tenure > 1) %>%
  pull(pidp)

df_changers <- df_panel %>% filter(pidp %in% tenure_changers)
panel_changers <- pdata.frame(df_changers, index = c("pidp", "wave"))

model_fe_changers <- plm(
  sclfsato ~ tenure + dvage + dvage_sq + log_income +
    partnered + nchild_dv + unemployed + retired + sick + covid,
  data = panel_changers,
  model = "within",
  effect = "twoways",
  na.action = na.omit
)

summary(model_fe_changers)

# How many changers?
length(tenure_changers)
nrow(df_changers)


save(df_panel, model_ols, model_fe, model_fe_covid, model_fe_u40, model_fe_o40,
     model_iv, model_fe_ghq, model_fe_changers,
     file = "ec226-group017.RData")
