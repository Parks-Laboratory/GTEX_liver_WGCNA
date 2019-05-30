#!/bin/bash

# untar your R installation
tar -xzf R.tar.gz

# make sure the script will use your R installation
export PATH=$(pwd)/R/bin:$PATH
export RHOME=$(pwd)/R

# get new outputs directary
mkdir outputs

# run R code
R CMD BATCH gtexR_GO.R gtex_GO_$1.Rout


# move all files in outputs to home directioary
mv outputs/* .
