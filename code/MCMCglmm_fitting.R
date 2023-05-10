#!/usr/bin/env Rscript

# This script reads the T_pk estimates for the four traits of the study 
# and fits a multivariate phylogenetic regression with MCMCglmm.
#
# It is used to estimate:
#
#	i)		the across-species average T_pk value,
#	ii)		the phylogenetic heritability of the T_pk of each trait,
#	iii)	the rate of evolution of the T_pk of each trait,
#	iv) 	the correlations between T_pk pairs.

library(doParallel)
library(MCMCglmm)

#####################
# F U N C T I O N S #
#####################

# This function checks that each parameter was adequately sampled from 
# the posterior, by calculating the Effective Sample Size.
check_ESS <- function(model_fits)
{
        all_ESS_vals <- c()

        for ( i in 1:length(model_fits) )
        {
                all_ESS_vals <- c(
                        all_ESS_vals,
                        effectiveSize(model_fits[[i]]$Sol),
                        effectiveSize(model_fits[[i]]$VCV),
                        effectiveSize(model_fits[[i]]$Deviance)
                )
        }

        all_ESS_vals <- all_ESS_vals[all_ESS_vals > 0]
        low_ESS <- which(all_ESS_vals < 1000)

        if ( length(low_ESS) > 0 )
        {
                stop("PROBLEM! The ESS is not big enough!")
        } else
        {
                return()
        }
}

# This function checks if the four MCMCglmm chains converged on statistically 
# indistinguishable posterior distributions, based on the Potential 
# Scale Reduction Factor diagnostic (Gelman & Rubin, Stat. Sci., 1992).
check_convergence <- function(model_fits)
{	
	Sols		<- 'mcmc.list('
	VCVs		<- 'mcmc.list('
	Deviances	<- 'mcmc.list('
	
	for ( i in 1:length(model_fits) )
	{
		Sols <- paste(
			Sols, 'model_fits[[', i, ']]$Sol, ', sep = ''
		)
		VCVs <- paste(
			VCVs, 'model_fits[[', i, ']]$VCV, ', sep = ''
		)
		Deviances <- paste(
			Deviances, 'model_fits[[', i, ']]$Deviance, ', sep = ''
		)
	}
	
	Sols		<- sub(', $', ')', Sols)
	VCVs		<- sub(', $', ')', VCVs)
	Deviances	<- sub(', $', ')', Deviances)
	
	psrf_vals <- c(
		gelman.diag(
			eval(parse(text = Sols)), multivariate = FALSE
		)$psrf[,1],
		gelman.diag(
			eval(parse(text = VCVs)), multivariate = FALSE
		)$psrf[,1],
		gelman.diag(
			eval(parse(text = Deviances)), multivariate = FALSE
		)$psrf[,1]
	)

	psrf_vals <- psrf_vals[!is.na(psrf_vals)]
	high_psrf <- which(psrf_vals >= 1.1)
	
	if ( length(high_psrf) > 0 )
	{
		stop("PROBLEM! The chains have not converged!")		
	} else
	{
		return()
	}
}

####################
# M A I N  C O D E #
####################

# Read the T_pk estimates for the four traits.
alpha_T_opts <- read.csv(
	'../data/alpha_Tpks_AllParams.csv', row.names = 1
)
alpha_T_opts <- alpha_T_opts[alpha_T_opts$param == 'topt',]

b_max_T_opts <- read.csv(
	'../data/bmax_Tpks_AllParams.csv', row.names = 1
)
b_max_T_opts <- b_max_T_opts[b_max_T_opts$param == 'topt',]

z_T_opts <- read.csv(
	'../data/z_Tpks_AllParams.csv', row.names = 1
)
z_T_opts <- z_T_opts[z_T_opts$param == 'topt',]

z_j_T_opts <- read.csv(
	'../data/zj_Tpks_AllParams.csv', row.names = 1
)
z_j_T_opts <- z_j_T_opts[z_j_T_opts$param == 'topt',]

# Get all the species that have at least one T_pk estimate.
all_sp <- unique(
	c(
		alpha_T_opts$species,	b_max_T_opts$species,
		z_T_opts$species,		z_j_T_opts$Species
	)
)

# Prepare a data frame to store all the data by species.
dataset <- data.frame(
	Species = all_sp,
	alpha_T_opt = rep(NA, length(all_sp)),
	alpha_T_opt_se = rep(0, length(all_sp)),
	b_max_T_opt = rep(NA, length(all_sp)),
	b_max_T_opt_se = rep(0, length(all_sp)),
	z_T_opt = rep(NA, length(all_sp)),
	z_T_opt_se = rep(0, length(all_sp)),
	z_j_T_opt = rep(NA, length(all_sp)),
	z_j_T_opt_se = rep(0, length(all_sp))
)

# Calculate the uncertainty around each estimate.
for ( i in 1:nrow(dataset) )
{
	if (
		length(
			alpha_T_opts$estimate[
				alpha_T_opts$species == dataset$Species[i]
			]
		) == 1 
	)
	{
		dataset$alpha_T_opt[i] <- alpha_T_opts$estimate[
			alpha_T_opts$species == dataset$Species[i]
		]
		
		dataset$alpha_T_opt_se[i] <- (
			alpha_T_opts$conf_upper[
				alpha_T_opts$species == dataset$Species[i]
			] - 
			alpha_T_opts$conf_lower[
				alpha_T_opts$species == dataset$Species[i]
			]
		) / 3.92
	}
	
	if (
		length(
			b_max_T_opts$estimate[
				b_max_T_opts$species == dataset$Species[i]
			]
		) == 1 
	)
	{
		dataset$b_max_T_opt[i] <- b_max_T_opts$estimate[
			b_max_T_opts$species == dataset$Species[i]
		]
		
		dataset$b_max_T_opt_se[i] <- (
			b_max_T_opts$conf_upper[
				b_max_T_opts$species == dataset$Species[i]
			] - 
			b_max_T_opts$conf_lower[
				b_max_T_opts$species == dataset$Species[i]
			]
		) / 3.92
	}
	
	if (
		length(
			z_T_opts$estimate[
				z_T_opts$species == dataset$Species[i]
			]
		) == 1 
	)
	{
		dataset$z_T_opt[i] <- z_T_opts$estimate[
			z_T_opts$species == dataset$Species[i]
		]
		
		dataset$z_T_opt_se[i] <- (
			z_T_opts$conf_upper[
				z_T_opts$species == dataset$Species[i]
			] - 
			z_T_opts$conf_lower[
				z_T_opts$species == dataset$Species[i]
			]
		) / 3.92
	}
	
	if (
		length(
			z_j_T_opts$estimate[
				z_j_T_opts$species == dataset$Species[i]
			]
		) == 1 
	)
	{
		dataset$z_j_T_opt[i] <- z_j_T_opts$estimate[
			z_j_T_opts$species == dataset$Species[i]
		]
		
		dataset$z_j_T_opt_se[i] <- (
			z_j_T_opts$conf_upper[
				z_j_T_opts$species == dataset$Species[i]
			] - 
			z_j_T_opts$conf_lower[
				z_j_T_opts$species == dataset$Species[i]
			]
		) / 3.92
	}
}

# Replace spaces with underscores in species' names.
dataset$Species <- gsub(' ', '_', dataset$Species)

# Prepare 3 threads for parallel execution.
cl <- makeCluster(3, outfile = '')
registerDoParallel(cl)

# Read each of the 100 alternative trees and fit the model with MCMCglmm 
# with 3 independent chains.
for ( i in 1:100 )
{
	cat('Now at run ', i, '...\n', sep = '')
	
	tree <- read.tree(
		paste(
			'../results/phylogeny/congruification/tree_', i, '.nwk', 
			sep = ''
		)
	)
	
	MCMCglmm_fits <- foreach(
		x = 1:3,
		.packages = c('MCMCglmm')
	) %dopar%
	{
		set.seed(x)
		fit <- MCMCglmm(
			cbind(alpha_T_opt, b_max_T_opt, z_T_opt, z_j_T_opt) ~ trait - 1,
			family = rep('gaussian', 4),
			random = ~us(trait):Species,
			ginverse = list(
				Species = inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv
			),
			prior = list(
				G = list(
					G1 = list(
						V = diag(4), nu = 1, alpha.mu = rep(0, 4),
						alpha.V = diag(rep(1000, 4))
					)
				),
				R = list(V = diag(4), nu = 0.002)
			),
			mev = c(
				dataset$alpha_T_opt_se,	dataset$b_max_T_opt_se,
				dataset$z_T_opt_se,		dataset$z_j_T_opt_se
			)^2,
			nitt = 200000000,
			burnin = 20000000,
			thin = 5000,
			rcov=~us(trait):units,
			data = dataset
		)
	}
	
	# Check that the 3 chains converged and that the posterior was 
	# adequately sampled.
	check_ESS(MCMCglmm_fits)
	check_convergence(MCMCglmm_fits)
	
	# Save the 3 chains to an .Rda file.
	save(
		MCMCglmm_fits, file = paste(
			'../results/MCMCglmm_runs/', i, '.Rda', sep = ''
		)
	)
}
