#!/bin/sh
#SBATCH --time=10:00:00
input_file="/home-1/bstrobe1@jhu.edu/work/ben/rare_var/rare_splice/modeling/v2/data/compressed_features_correct_format_sorted.txt"

output_dir="/home-1/bstrobe1@jhu.edu/work/ben/rare_var/rare_splice/modeling/v_scripts_for_emily/output/"





Rscript univariate_naive_bayes.R $input_file $output_dir"univariate_naive_bayes_"



Rscript watershed_ising_pseudo.R $input_file $output_dir"crf_type_model_"



