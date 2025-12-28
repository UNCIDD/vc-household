
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
  
  source("R/run_models/vib_sens.R")
  source("R/run_models/init_sens.R")
  source("R/run_models/phi_sens.R")
  source("R/run_models/simulations.R")
  
  nocov <- read.csv("model_results/new_runs/res_nocov_320.csv") # baseline model
  cov <- read.csv("model_results/new_runs/res_cov.csv") # covariates model
  ih_chains <- read.csv("model_results/new_runs/res_nocov_320_ihchains.csv") # chains for ih risks
  vib_sens <- read.csv("model_results/new_runs/res_vibsens.csv") # sensitivity to high vibriocidal titer cutoff
  init_sens <- read.csv("model_results/new_runs/res_initsens.csv") # sensitivity to proportion of index cases who start symptomatically infected
  phi_sens <- read.csv("model_results/new_runs/res_phisens.csv") # sensitivity to probability of observing symptoms in uninfected individuals
  sim_params <- read.csv("model_results/new_runs/sim_param_ests.csv") # parameter estimates from simulations
  
} else {
  nocov <- read.csv("model_results/pre_saved/res_nocov_320.csv") # baseline model
  cov <- read.csv("model_results/pre_saved/res_cov.csv") # covariates models
  ih_chains <- read.csv("model_results/pre_saved/res_nocov_320_ihchains.csv") # chains for ih risks
  vib_sens <- read.csv("model_results/pre_saved/res_vibsens.csv") # sensitivity to high vibriocidal titer cutoff
  init_sens <- read.csv("model_results/pre_saved/res_initsens.csv") # sensitivity to proportion of index cases who start symptomatically infected
  phi_sens <- read.csv("model_results/pre_saved/res_phisens.csv") # sensitivity to probability of observing symptoms in uninfected individuals
  sim_params <- read.csv("model_results/pre_saved/sim_param_ests.csv") # parameter estimates from simulations
  sim_truth <- read.csv("model_results/pre_saved/sim_truth.rds")
  sim_stateprobs <- read.csv("model_results/pre_saved/sim_stateprobs.rds")
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

## Figure S4
make_figS4(vib_sens, nocov)

## Figure S5
make_figS5(init_sens, nocov)

## Figure S6
make_figS6(phi_sens, nocov)

## Figure S7
make_figS7(sim_params)

## Figure S8
make_figS8(sim_stateprobs, sim_truth)
