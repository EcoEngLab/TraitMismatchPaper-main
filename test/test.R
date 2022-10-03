#__________________________________________________________________________________________________________
#   Using Bayesian approaches to fit TPCs for Glossina traits involved in Trypanosomiasis transmission      
#__________________________________________________________________________________________________________

rm(list=ls())
graphics.off()

# Load libraries
require(rjags)   # does the fitting
require(coda)    # makes diagnostic plots
require(R2jags)  # fitting
require(MCMCvis) # for trace plots
require(IDPmisc) # makes nice colored pairs plots to look at joint posteriors
require(plotrix)
require(tidyverse)
require(patchwork) # easy way of creating panels