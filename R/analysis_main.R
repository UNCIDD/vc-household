
# Option for whether to re-run stan models
# If FALSE, code will use pre-saved model results
# If TRUE, stan models will re-run, this will take a long time (recommed running on a computing cluster)
run_models <- FALSE
cores <- 4 # number of cores for running stan models (should be 4 if running in parallel or 1 if not running in parallel)

# Dependencies
library(tidyverse)

source("R/helper_scripts/cholera_utilities.R")
source("R/helper_scripts/figure_functions.R")

# Read in raw data
bl <- read.csv("data/dhaka_baseline_2024-11-26.csv")[,-1]
fu <- read.csv("data/dhaka_outcomes_2024-11-26.csv")[,-1] 

# Get model results
if(run_models) {
  source("R/helper_scripts/generate_stan_data.R") # get raw data into stan format
  source("R/run_models/nocov_main.R") # main analysis, no covariates
} else {
  nocov <- read.csv("model_results/pre_saved/res_nocov_320.csv") # baseline model
  cov <- read.csv("model_results/pre_saved/res_cov.csv") # covariates models
}

# Make figures

## Figure 1
### Note that the points for the two contact types overlap, I fix this manually in illustrator
make_fig1()

## Figure 2 is made manually in illustrator

## Figure 3
make_fig3(dat = nocov)

## Figure 4
make_fig4(dat = cov)

# Make tables

## Table 1 (all values are counts except for age which are mean and IQR as indicated)
make_table1(bl, fu)

## Table 2
make_table2(bl, fu)
