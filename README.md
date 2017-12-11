# RIVER MODELS FOR EMILY


## System Requirements
The code was developed on MARCC.

The code is written in R (`module load R/3.4.0`) and relies on several R packages that you will need to install:
1. library(glmnet)
2. library(methods)
3. library(stats)
4. library(utils)
5. library(Biobase)
6. library(pROC)
7. library(ggplot2)
8. library(cowplot)
9. library(sigmoid)
10. library(Rcpp)
11. library(RColorBrewer)
12. library(ggthemes)
13. library(lbfgs)
14. library(optimx)
15. library(numDeriv)

## Running the code
Both versions of the code can be run with the command `sh driver_key.sh`.
You will probably wish to change the variables `input_file` and `output_dir` in `driver_key.sh`.


There are two algorithms called from driver_key.sh:
1. `univariate_naive_bayes.R`. This is basically RIVER (with a few add-ons to make it more customizable/generalizable)
2. `watershed_ising_pseudo.R`. This is the WATERSHED model. I've added a valid_tissues variable that allows the user to control which tissues to use.
NOTE: `watershed_ising_pseudo.R` accomplishes most of its work by calling cpp functions in `crf_vi_updates.cpp`