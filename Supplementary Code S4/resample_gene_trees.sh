#!/bin/bash

#######################################
# Resample gene trees with replacement
# Tauana Cunha | https://github.com/tauanajc/Cunha_Reimer_Giribet_2021_SystBio
# August, 2021
# Usage: resample_gene_trees.sh

# Input is one file with all gene trees (same as input for astral), here named AllGeneTrees.trees
# -n : total number of gene trees, same as the number of lines in input file, here 1027
# -r : with replacement

# This will create as many files as the number of resamplings,
# each named genetrees_ followed by the id number of the resampling
#######################################

for i in {1..2000}; do # 2000 resamplings
shuf -r -n 1027 AllGeneTrees.trees > genetrees_$i
done
