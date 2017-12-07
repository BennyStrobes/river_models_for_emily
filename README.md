# RIVER MODELS FOR EMILY


## System Requirements
The code is written in R (`module load R/3.2.3`) and relies on several R packages that you will need to install:
1. library(glmnet)
2. library(methods)
3. library(stats)
4. library(utils)
5. library(Biobase)
6. library(pROC)
7. library(ggplot2)

## Running the code
Both versions of the code can be run with the command `sh driver_key.sh`.

There are two algorithms called from driver_key.sh:
1. `univariate_naive_bayes.R`. This is basically RIVER (with a few add-ons to make it more customizable/generalizable)
2. To-be-added (CRF-type model)