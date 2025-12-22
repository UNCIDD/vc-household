
# Option for whether to re-run stan models
# If FALSE, code will use pre-saved model results
# If TRUE, stan models will re-run, this will take a long time (recommed running on a computing cluster)
# If TRUE, assumes analysis_main.R has already been run with run_models = TRUE 
run_models <- FALSE
cores <- 4 # number of cores for running stan models (should be 4 if running in parallel or 1 if not running in parallel)

# Dependencies
library(tidyverse)

source("R/helper_scripts/cholera_utilities.R")
source("R/helper_scripts/figure_functions.R")

# Get model results
if(run_models) {
  
  nocov <- read.csv("model_results/new_runs/res_nocov_320.csv") # baseline model
  cov <- read.csv("model_results/new_runs/res_cov.csv") # covariates model
  ih_chains <- read.csv("model_results/new_runs/res_nocov_320_ihchains.csv") # chains for ih risks
  
} else {
  nocov <- read.csv("model_results/pre_saved/res_nocov_320.csv") # baseline model
  cov <- read.csv("model_results/pre_saved/res_cov.csv") # covariates models
  ih_chains <- read.csv("model_results/pre_saved/res_nocov_320_ihchains.csv") # chains for ih risks
}

# Read in raw data
bl <- read.csv("data/dhaka_baseline_2024-11-26.csv")[,-1]
fu <- read.csv("data/dhaka_outcomes_2024-11-26.csv")[,-1] 

# Make supplemental figures

## Figure S1
make_figS1(bl, fu)

## Figure S2
make_figS2(bl, fu)

## Figure S3
make_figS3(ih_chains)


