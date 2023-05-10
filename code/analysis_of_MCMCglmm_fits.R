#!/usr/bin/env Rscript

# This script reads the 100 MCMCglmm fits of the model, combines their 
# posterior estimates, and generates plots of:
#
#	i)		the across-species average T_pk value,
#	ii)		the phylogenetic heritability of the T_pk of each trait,
#	iii)	the rate of evolution of the T_pk of each trait,
#	iv) 	the correlations between T_pk pairs.

library(cowplot)
library(ggplot2)
library(MCMCglmm)
library(robustbase)

#####################
# F U N C T I O N S #
#####################

# This function combines the posterior samples from all chains from 
# the 100 runs of the model.
combine_samples <- function(run_files)
{
	
	# Prepare a data frame to store each estimate of interest.
	n_rows <- 36000 * 3 * length(run_files)
	
	sampled_vars <- data.frame(
		her_alpha = rep(NA, n_rows),
		her_beta = rep(NA, n_rows),
		her_z = rep(NA, n_rows),
		her_z_j = rep(NA, n_rows),
		var_alpha = rep(NA, n_rows),
		var_beta = rep(NA, n_rows),
		var_z = rep(NA, n_rows),
		var_z_j = rep(NA, n_rows),
		alpha = rep(NA, n_rows),
		beta = rep(NA, n_rows),
		z = rep(NA, n_rows),
		z_j = rep(NA, n_rows),
		cor_alpha_beta = rep(NA, n_rows),
		cor_alpha_z = rep(NA, n_rows),
		cor_alpha_z_j = rep(NA, n_rows),
		cor_beta_z = rep(NA, n_rows),
		cor_beta_z_j = rep(NA, n_rows),
		cor_z_z_j = rep(NA, n_rows)
	)
	
	current_row <- 1
	
	# Read each of the 100 model runs and store the estimates.
	for ( i in 1:length(run_files) )
	{
		load(
			paste('../results/MCMCglmm_runs/', run_files[i], sep = ''),
			.GlobalEnv
		)
	
		entries_to_replace <- current_row:(36000 * 3 + current_row - 1)
		
		current_hers <- get_phylogenetic_heritabilities()
		
		sampled_vars$her_alpha[entries_to_replace] <- current_hers$her_alpha
		sampled_vars$her_beta[entries_to_replace] <- current_hers$her_beta	
		sampled_vars$her_z[entries_to_replace] <- current_hers$her_z	
		sampled_vars$her_z_j[entries_to_replace] <- current_hers$her_z_j
		
		current_vars <- get_variances()
		
		sampled_vars$var_alpha[entries_to_replace] <- current_vars$var_alpha
		sampled_vars$var_beta[entries_to_replace] <- current_vars$var_beta
		sampled_vars$var_z[entries_to_replace] <- current_vars$var_z
		sampled_vars$var_z_j[entries_to_replace] <- current_vars$var_z_j
		
		current_means <- get_intercepts()
		
		sampled_vars$alpha[entries_to_replace] <- current_means$alpha
		sampled_vars$beta[entries_to_replace] <- current_means$beta
		sampled_vars$z[entries_to_replace] <- current_means$z
		sampled_vars$z_j[entries_to_replace] <- current_means$z_j
		
		current_correlations <- get_correlations()
		
		sampled_vars$cor_alpha_beta[entries_to_replace] <- current_correlations$cor_alpha_beta
		sampled_vars$cor_alpha_z[entries_to_replace] <- current_correlations$cor_alpha_z
		sampled_vars$cor_alpha_z_j[entries_to_replace] <- current_correlations$cor_alpha_z_j
		sampled_vars$cor_beta_z[entries_to_replace] <- current_correlations$cor_beta_z
		sampled_vars$cor_beta_z_j[entries_to_replace] <- current_correlations$cor_beta_z_j
		sampled_vars$cor_z_z_j[entries_to_replace] <- current_correlations$cor_z_z_j
	
		current_row <- current_row + 36000 * 3
	}
	
	return(sampled_vars)
}

# This function calculates the correlations between pairs of T_pks from 
# the variance-covariance matrices of the 3 chains per MCMCglmm fit.
get_correlations <- function()
{
	cor_alpha_beta <- c(
		(
			MCMCglmm_fits[[1]]$VCV[,'traitalpha_T_opt:traitb_max_T_opt.Species'] +
			MCMCglmm_fits[[1]]$VCV[,'traitalpha_T_opt:traitb_max_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[1]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species'] +
				MCMCglmm_fits[[1]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.units']
			) * (
				MCMCglmm_fits[[1]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species'] +
				MCMCglmm_fits[[1]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.units']
			)
		),
		(
			MCMCglmm_fits[[2]]$VCV[,'traitalpha_T_opt:traitb_max_T_opt.Species'] +
			MCMCglmm_fits[[2]]$VCV[,'traitalpha_T_opt:traitb_max_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[2]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species'] +
				MCMCglmm_fits[[2]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.units']
			) * (
				MCMCglmm_fits[[2]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species'] +
				MCMCglmm_fits[[2]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.units']
			)
		),
		(
			MCMCglmm_fits[[3]]$VCV[,'traitalpha_T_opt:traitb_max_T_opt.Species'] +
			MCMCglmm_fits[[3]]$VCV[,'traitalpha_T_opt:traitb_max_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[3]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species'] +
				MCMCglmm_fits[[3]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.units']
			) * (
				MCMCglmm_fits[[3]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species'] +
				MCMCglmm_fits[[3]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.units']
			)
		)
	)
	
	cor_alpha_z <- c(
		(
			MCMCglmm_fits[[1]]$VCV[,'traitalpha_T_opt:traitz_T_opt.Species'] +
			MCMCglmm_fits[[1]]$VCV[,'traitalpha_T_opt:traitz_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[1]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species'] +
				MCMCglmm_fits[[1]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.units']
			) * (
				MCMCglmm_fits[[1]]$VCV[,'traitz_T_opt:traitz_T_opt.Species'] +
				MCMCglmm_fits[[1]]$VCV[,'traitz_T_opt:traitz_T_opt.units']
			)
		),
		(
			MCMCglmm_fits[[2]]$VCV[,'traitalpha_T_opt:traitz_T_opt.Species'] +
			MCMCglmm_fits[[2]]$VCV[,'traitalpha_T_opt:traitz_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[2]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species'] +
				MCMCglmm_fits[[2]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.units']
			) * (
				MCMCglmm_fits[[2]]$VCV[,'traitz_T_opt:traitz_T_opt.Species'] +
				MCMCglmm_fits[[2]]$VCV[,'traitz_T_opt:traitz_T_opt.units']
			)
		),
		(
			MCMCglmm_fits[[3]]$VCV[,'traitalpha_T_opt:traitz_T_opt.Species'] +
			MCMCglmm_fits[[3]]$VCV[,'traitalpha_T_opt:traitz_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[3]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species'] +
				MCMCglmm_fits[[3]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.units']
			) * (
				MCMCglmm_fits[[3]]$VCV[,'traitz_T_opt:traitz_T_opt.Species'] +
				MCMCglmm_fits[[3]]$VCV[,'traitz_T_opt:traitz_T_opt.units']
			)
		)
	)
	
	cor_alpha_z_j <- c(
		(
			MCMCglmm_fits[[1]]$VCV[,'traitalpha_T_opt:traitz_j_T_opt.Species'] +
			MCMCglmm_fits[[1]]$VCV[,'traitalpha_T_opt:traitz_j_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[1]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species'] +
				MCMCglmm_fits[[1]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.units']
			) * (
				MCMCglmm_fits[[1]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species'] +
				MCMCglmm_fits[[1]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.units']
			)
		),
		(
			MCMCglmm_fits[[2]]$VCV[,'traitalpha_T_opt:traitz_j_T_opt.Species'] +
			MCMCglmm_fits[[2]]$VCV[,'traitalpha_T_opt:traitz_j_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[2]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species'] +
				MCMCglmm_fits[[2]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.units']
			) * (
				MCMCglmm_fits[[2]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species'] +
				MCMCglmm_fits[[2]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.units']
			)
		),
		(
			MCMCglmm_fits[[3]]$VCV[,'traitalpha_T_opt:traitz_j_T_opt.Species'] +
			MCMCglmm_fits[[3]]$VCV[,'traitalpha_T_opt:traitz_j_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[3]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species'] +
				MCMCglmm_fits[[3]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.units']
			) * (
				MCMCglmm_fits[[3]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species'] +
				MCMCglmm_fits[[3]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.units']
			)
		)
	)
		
	cor_beta_z <- c(
		(
			MCMCglmm_fits[[1]]$VCV[,'traitb_max_T_opt:traitz_T_opt.Species'] +
			MCMCglmm_fits[[1]]$VCV[,'traitb_max_T_opt:traitz_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[1]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species'] +
				MCMCglmm_fits[[1]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.units']
			) * (
				MCMCglmm_fits[[1]]$VCV[,'traitz_T_opt:traitz_T_opt.Species'] +
				MCMCglmm_fits[[1]]$VCV[,'traitz_T_opt:traitz_T_opt.units']
			)
		),
		(
			MCMCglmm_fits[[2]]$VCV[,'traitb_max_T_opt:traitz_T_opt.Species'] +
			MCMCglmm_fits[[2]]$VCV[,'traitb_max_T_opt:traitz_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[2]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species'] +
				MCMCglmm_fits[[2]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.units']
			) * (
				MCMCglmm_fits[[2]]$VCV[,'traitz_T_opt:traitz_T_opt.Species'] +
				MCMCglmm_fits[[2]]$VCV[,'traitz_T_opt:traitz_T_opt.units']
			)
		),
		(
			MCMCglmm_fits[[3]]$VCV[,'traitb_max_T_opt:traitz_T_opt.Species'] +
			MCMCglmm_fits[[3]]$VCV[,'traitb_max_T_opt:traitz_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[3]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species'] +
				MCMCglmm_fits[[3]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.units']
			) * (
				MCMCglmm_fits[[3]]$VCV[,'traitz_T_opt:traitz_T_opt.Species'] +
				MCMCglmm_fits[[3]]$VCV[,'traitz_T_opt:traitz_T_opt.units']
			)
		)
	)
	
	cor_beta_z_j <- c(
		(
			MCMCglmm_fits[[1]]$VCV[,'traitb_max_T_opt:traitz_j_T_opt.Species'] +
			MCMCglmm_fits[[1]]$VCV[,'traitb_max_T_opt:traitz_j_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[1]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species'] +
				MCMCglmm_fits[[1]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.units']
			) * (
				MCMCglmm_fits[[1]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species'] +
				MCMCglmm_fits[[1]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.units']
			)
		),
		(
			MCMCglmm_fits[[2]]$VCV[,'traitb_max_T_opt:traitz_j_T_opt.Species'] +
			MCMCglmm_fits[[2]]$VCV[,'traitb_max_T_opt:traitz_j_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[2]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species'] +
				MCMCglmm_fits[[2]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.units']
			) * (
				MCMCglmm_fits[[2]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species'] +
				MCMCglmm_fits[[2]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.units']
			)
		),
		(
			MCMCglmm_fits[[3]]$VCV[,'traitb_max_T_opt:traitz_j_T_opt.Species'] +
			MCMCglmm_fits[[3]]$VCV[,'traitb_max_T_opt:traitz_j_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[3]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species'] +
				MCMCglmm_fits[[3]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.units']
			) * (
				MCMCglmm_fits[[3]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species'] +
				MCMCglmm_fits[[3]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.units']
			)
		)
	)
	
	cor_z_z_j <- c(
		(
			MCMCglmm_fits[[1]]$VCV[,'traitz_T_opt:traitz_j_T_opt.Species'] +
			MCMCglmm_fits[[1]]$VCV[,'traitz_T_opt:traitz_j_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[1]]$VCV[,'traitz_T_opt:traitz_T_opt.Species'] +
				MCMCglmm_fits[[1]]$VCV[,'traitz_T_opt:traitz_T_opt.units']
			) * (
				MCMCglmm_fits[[1]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species'] +
				MCMCglmm_fits[[1]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.units']
			)
		),
		(
			MCMCglmm_fits[[2]]$VCV[,'traitz_T_opt:traitz_j_T_opt.Species'] +
			MCMCglmm_fits[[2]]$VCV[,'traitz_T_opt:traitz_j_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[2]]$VCV[,'traitz_T_opt:traitz_T_opt.Species'] +
				MCMCglmm_fits[[2]]$VCV[,'traitz_T_opt:traitz_T_opt.units']
			) * (
				MCMCglmm_fits[[2]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species'] +
				MCMCglmm_fits[[2]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.units']
			)
		),
		(
			MCMCglmm_fits[[3]]$VCV[,'traitz_T_opt:traitz_j_T_opt.Species'] +
			MCMCglmm_fits[[3]]$VCV[,'traitz_T_opt:traitz_j_T_opt.units']
		) / sqrt(
			(
				MCMCglmm_fits[[3]]$VCV[,'traitz_T_opt:traitz_T_opt.Species'] +
				MCMCglmm_fits[[3]]$VCV[,'traitz_T_opt:traitz_T_opt.units']
			) * (
				MCMCglmm_fits[[3]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species'] +
				MCMCglmm_fits[[3]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.units']
			)
		)
	)
	
	return(
		list(
			cor_alpha_beta = cor_alpha_beta,	cor_alpha_z = cor_alpha_z,
			cor_alpha_z_j = cor_alpha_z_j,		cor_beta_z = cor_beta_z,
			cor_beta_z_j = cor_beta_z_j,		cor_z_z_j = cor_z_z_j
		)
	)
}

# This function returns the estimates of the average T_pks across species 
# from the 3 chains per MCMCglmm fit.
get_intercepts <- function()
{
	alpha <- c(
		MCMCglmm_fits[[1]]$Sol[,'traitalpha_T_opt'],
		MCMCglmm_fits[[2]]$Sol[,'traitalpha_T_opt'],
		MCMCglmm_fits[[3]]$Sol[,'traitalpha_T_opt']
	)
	
	beta <- c(
		MCMCglmm_fits[[1]]$Sol[,'traitb_max_T_opt'],
		MCMCglmm_fits[[2]]$Sol[,'traitb_max_T_opt'],
		MCMCglmm_fits[[3]]$Sol[,'traitb_max_T_opt']
	)
	
	z <- c(
		MCMCglmm_fits[[1]]$Sol[,'traitz_T_opt'],
		MCMCglmm_fits[[2]]$Sol[,'traitz_T_opt'],
		MCMCglmm_fits[[3]]$Sol[,'traitz_T_opt']
	)
	
	z_j <- c(
		MCMCglmm_fits[[1]]$Sol[,'traitz_j_T_opt'],
		MCMCglmm_fits[[2]]$Sol[,'traitz_j_T_opt'],
		MCMCglmm_fits[[3]]$Sol[,'traitz_j_T_opt']
	)
	
	return(
		list(
			alpha = alpha,	beta = beta, 
			z = z,			z_j = z_j
		)
	)
}

# This function calculates the phylogenetic heritability for each T_pk, 
# based on the variance-covariance matrices of the 3 chains per MCMCglmm 
# fit.
get_phylogenetic_heritabilities <- function()
{
	her_alpha <- c(
		MCMCglmm_fits[[1]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species'] / (
			MCMCglmm_fits[[1]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species'] + 
			MCMCglmm_fits[[1]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.units']
		),
		MCMCglmm_fits[[2]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species'] / (
			MCMCglmm_fits[[2]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species'] + 
			MCMCglmm_fits[[2]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.units']
		),
		MCMCglmm_fits[[3]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species'] / (
			MCMCglmm_fits[[3]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species'] + 
			MCMCglmm_fits[[3]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.units']
		)
	)
	
	her_beta <- c(
		MCMCglmm_fits[[1]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species'] / (
			MCMCglmm_fits[[1]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species'] + 
			MCMCglmm_fits[[1]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.units']
		),
		MCMCglmm_fits[[2]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species'] / (
			MCMCglmm_fits[[2]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species'] + 
			MCMCglmm_fits[[2]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.units']
		),
		MCMCglmm_fits[[3]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species'] / (
			MCMCglmm_fits[[3]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species'] + 
			MCMCglmm_fits[[3]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.units']
		)
	)
	
	her_z <- c(
		MCMCglmm_fits[[1]]$VCV[,'traitz_T_opt:traitz_T_opt.Species'] / (
			MCMCglmm_fits[[1]]$VCV[,'traitz_T_opt:traitz_T_opt.Species'] + 
			MCMCglmm_fits[[1]]$VCV[,'traitz_T_opt:traitz_T_opt.units']
		),
		MCMCglmm_fits[[2]]$VCV[,'traitz_T_opt:traitz_T_opt.Species'] / (
			MCMCglmm_fits[[2]]$VCV[,'traitz_T_opt:traitz_T_opt.Species'] + 
			MCMCglmm_fits[[2]]$VCV[,'traitz_T_opt:traitz_T_opt.units']
		),
		MCMCglmm_fits[[3]]$VCV[,'traitz_T_opt:traitz_T_opt.Species'] / (
			MCMCglmm_fits[[3]]$VCV[,'traitz_T_opt:traitz_T_opt.Species'] + 
			MCMCglmm_fits[[3]]$VCV[,'traitz_T_opt:traitz_T_opt.units']
		)
	)
	
	her_z_j <- c(
		MCMCglmm_fits[[1]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species'] / (
			MCMCglmm_fits[[1]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species'] + 
			MCMCglmm_fits[[1]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.units']
		),
		MCMCglmm_fits[[2]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species'] / (
			MCMCglmm_fits[[2]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species'] + 
			MCMCglmm_fits[[2]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.units']
		),
		MCMCglmm_fits[[3]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species'] / (
			MCMCglmm_fits[[3]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species'] + 
			MCMCglmm_fits[[3]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.units']
		)
	)
	
	return(
		list(
			her_alpha = her_alpha,	her_beta = her_beta, 
			her_z = her_z,			her_z_j = her_z_j
		)
	)
}

# This function calculates the evolutionary rate for each T_pk from the 
# 3 chains per MCMCglmm fit.
get_variances <- function()
{
	var_alpha <- c(
		MCMCglmm_fits[[1]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species'],
		MCMCglmm_fits[[2]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species'],
		MCMCglmm_fits[[3]]$VCV[,'traitalpha_T_opt:traitalpha_T_opt.Species']
	)
	
	var_beta <- c(
		MCMCglmm_fits[[1]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species'],
		MCMCglmm_fits[[2]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species'],
		MCMCglmm_fits[[3]]$VCV[,'traitb_max_T_opt:traitb_max_T_opt.Species']
	)
	
	var_z <- c(
		MCMCglmm_fits[[1]]$VCV[,'traitz_T_opt:traitz_T_opt.Species'],
		MCMCglmm_fits[[2]]$VCV[,'traitz_T_opt:traitz_T_opt.Species'],
		MCMCglmm_fits[[3]]$VCV[,'traitz_T_opt:traitz_T_opt.Species']
	)
	
	var_z_j <- c(
		MCMCglmm_fits[[1]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species'],
		MCMCglmm_fits[[2]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species'],
		MCMCglmm_fits[[3]]$VCV[,'traitz_j_T_opt:traitz_j_T_opt.Species']
	)
	
	return(
		list(
			var_alpha = var_alpha,	var_beta = var_beta,
			var_z = var_z,			var_z_j = var_z_j
		)
	)
}

# This function combines all the plots into a single output file.
make_plots <- function(sampled_vars)
{
	theme_set(theme_cowplot())
	
	p_across_species_average <- plot_across_species_average(sampled_vars)
	
	p_herit <- plot_heritabilities(sampled_vars)
	
	p_rate <- plot_evolutionary_rates(sampled_vars)
	
	p_cor <- plot_correlations(sampled_vars)
	
	plots_1st_pass <- align_plots(
		p_across_species_average, p_cor, align = 'v'
	)
	top_row_plots <- plot_grid(plots_1st_pass[[1]], p_herit, p_rate, nrow = 1)
	
	ggsave(
		plot_grid(
			top_row_plots, plots_1st_pass[[2]], ncol = 1
		),
		file = '../results/phylo_summary_plots.pdf',
		width = 4.75, height = 2.8, units = 'in'
	)
}

# This function reads the estimates of the mean T_pk per trait across 
# species and generates a plot.
plot_across_species_average <- function(sampled_vars)
{
	data_to_plot <- data.frame(
		parameter = c('alpha', 'beta', 'z_j', 'z'),
		median_value = colMedians(
			as.matrix(
				sampled_vars[
					,c('alpha', 'beta', 'z_j', 'z')
				]
			)
		),
		lwr_value = HPDinterval(
			mcmc(
				sampled_vars[
					,c('alpha', 'beta', 'z_j', 'z')
				]
			)
		)[,1],
		upr_value = HPDinterval(
			mcmc(
				sampled_vars[
					,c('alpha', 'beta', 'z_j', 'z')
				]
			)
		)[,2]
	)
	
	data_to_plot_2 <- data.frame(
		parameter = c(
			rep('alpha', nrow(sampled_vars)),
			rep('beta', nrow(sampled_vars)),
			rep('z_j', nrow(sampled_vars)),
			rep('z', nrow(sampled_vars))
		),
		value = c(
			sampled_vars$alpha,	sampled_vars$beta,
			sampled_vars$z_j,	sampled_vars$z
		)
	)
	data_to_plot_2$parameter <- factor(
		data_to_plot_2$parameter, levels = c('alpha', 'beta', 'z_j', 'z')
	)
	
	# Here we only keep the values that lie within the 95% HPD interval.
	data_to_plot_2 <- data_to_plot_2[
		(
			data_to_plot_2$parameter == 'alpha' & 
			data_to_plot_2$value <= data_to_plot['alpha', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['alpha', 'lwr_value']
		) |
		(
			data_to_plot_2$parameter == 'beta' & 
			data_to_plot_2$value <= data_to_plot['beta', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['beta', 'lwr_value']
		) |
		(
			data_to_plot_2$parameter == 'z_j' & 
			data_to_plot_2$value <= data_to_plot['z_j', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['z_j', 'lwr_value']
		) |
		(
			data_to_plot_2$parameter == 'z' & 
			data_to_plot_2$value <= data_to_plot['z', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['z', 'lwr_value']
		),
	]
		
	p_across_species_average <- ggplot(data_to_plot_2, aes(x = parameter, y = value)) +
		geom_violin(
			position = 'dodge', fill = 'goldenrod1', width = 1,
			color = 'goldenrod1'
		) + 
		ggtitle(' A') +
		geom_point(
			data = data_to_plot, aes(x = parameter, y = median_value), 
			pch = 21, fill = 'white', size = 3
		) +
		scale_y_continuous(
			name = 'Across-species average (°C)',
			expand = c(0.005,0)
		) +
		scale_x_discrete(
			name = '',
			labels = c(
				expression(italic(T)[italic('pk')]^alpha),
				expression(italic(T)[italic('pk')]^'b'['max']),
				expression(italic(T)[italic('pk')]^italic('z')['j']),
				expression(italic(T)[italic('pk')]^italic('z'))
			),
			expand = c(0.19,0)
		) +
		theme(
			plot.title = element_text(margin = margin(b = -10), size = 8),
			axis.line = element_line(colour = NA),
			axis.title.y = element_text(size = 7, family = 'serif'),
			axis.text.x = element_text(size = 7, family = 'serif'),
			axis.text.y = element_text(size = 7, family = 'serif'),
			plot.margin = unit(c(0.2,0.2,-0.55,0.05), 'cm'),
			panel.border = element_rect(colour = "black", fill=NA, size=0.5)
		)	

	return(p_across_species_average)
}

# This function reads the estimates of the correlations between T_pk 
# pairs and generates a plot.
plot_correlations <- function(sampled_vars)
{
	data_to_plot <- data.frame(
		parameter = c(
			'cor_alpha_beta',	'cor_alpha_z',	'cor_alpha_z_j',
			'cor_beta_z',		'cor_beta_z_j',	'cor_z_z_j'
		),
		median_value = colMedians(
			as.matrix(
				sampled_vars[
					,c(
						'cor_alpha_beta',	'cor_alpha_z',	'cor_alpha_z_j',
						'cor_beta_z',		'cor_beta_z_j',	'cor_z_z_j'
					)
				]
			)
		),
		lwr_value = HPDinterval(
			mcmc(
				sampled_vars[
					,c(
						'cor_alpha_beta',	'cor_alpha_z',	'cor_alpha_z_j',
						'cor_beta_z',		'cor_beta_z_j',	'cor_z_z_j'
					)
				]
			)
		)[,1],
		upr_value = HPDinterval(
			mcmc(
				sampled_vars[
					,c(
						'cor_alpha_beta',	'cor_alpha_z',	'cor_alpha_z_j',
						'cor_beta_z',		'cor_beta_z_j',	'cor_z_z_j'
					)
				]
			)
		)[,2]
	)
	
	data_to_plot_2 <- data.frame(
		parameter = c(
			rep('cor_alpha_beta', nrow(sampled_vars)),
			rep('cor_alpha_z', nrow(sampled_vars)),
			rep('cor_alpha_z_j', nrow(sampled_vars)),
			rep('cor_beta_z', nrow(sampled_vars)),
			rep('cor_beta_z_j', nrow(sampled_vars)),
			rep('cor_z_z_j', nrow(sampled_vars))
		),
		value = c(
			sampled_vars$cor_alpha_beta,	sampled_vars$cor_alpha_z,
			sampled_vars$cor_alpha_z_j,		sampled_vars$cor_beta_z,
			sampled_vars$cor_beta_z_j,		sampled_vars$cor_z_z_j
		)
	)
	data_to_plot_2$parameter <- factor(
		data_to_plot_2$parameter, levels = c(
			'cor_alpha_beta', 'cor_alpha_z_j', 
			'cor_alpha_z', 'cor_beta_z_j', 
			'cor_beta_z', 'cor_z_z_j'
		)
	)
	
	# Here we only keep the values that lie within the 95% HPD interval.
	data_to_plot_2 <- data_to_plot_2[
		(
			data_to_plot_2$parameter == 'cor_alpha_beta' & 
			data_to_plot_2$value <= data_to_plot['cor_alpha_beta', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['cor_alpha_beta', 'lwr_value']
		) |
		(
			data_to_plot_2$parameter == 'cor_alpha_z' & 
			data_to_plot_2$value <= data_to_plot['cor_alpha_z', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['cor_alpha_z', 'lwr_value']
		) |
		(
			data_to_plot_2$parameter == 'cor_alpha_z_j' & 
			data_to_plot_2$value <= data_to_plot['cor_alpha_z_j', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['cor_alpha_z_j', 'lwr_value']
		) |
		(
			data_to_plot_2$parameter == 'cor_beta_z' & 
			data_to_plot_2$value <= data_to_plot['cor_beta_z', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['cor_beta_z', 'lwr_value']
		) |
		(
			data_to_plot_2$parameter == 'cor_beta_z_j' & 
			data_to_plot_2$value <= data_to_plot['cor_beta_z_j', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['cor_beta_z_j', 'lwr_value']
		) |
		(
			data_to_plot_2$parameter == 'cor_z_z_j' & 
			data_to_plot_2$value <= data_to_plot['cor_z_z_j', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['cor_z_z_j', 'lwr_value']
		),		
	]
		
	p_cor <- ggplot(data_to_plot_2, aes(x = parameter, y = value)) +
		geom_violin(
			position = 'dodge', fill = 'goldenrod1', width = 1,
			color = 'goldenrod1'
		) + 
		ggtitle(' D') +
		geom_point(
			data = data_to_plot, aes(x = parameter, y = median_value), 
			pch = 21, fill = 'white', size = 3
		) +
		scale_y_continuous(
			name = 'Correlation coefficient',
			expand = c(0.005,0),
			limits = c(-1,1),
		) +
		scale_x_discrete(
			name = '',
			labels = c(
				expression('corr('*italic(T)[italic('pk')]^alpha*', '*italic(T)[italic('pk')]^italic(b)['max']*')'),
				expression('corr('*italic(T)[italic('pk')]^alpha*', '*italic(T)[italic('pk')]^italic(z)['j']*')'),
				expression('corr('*italic(T)[italic('pk')]^alpha*', '*italic(T)[italic('pk')]^italic(z)*')'),
				expression('corr('*italic(T)[italic('pk')]^italic(b)['max']*', '*italic(T)[italic('pk')]^italic(z)['j']*')'),
				expression('corr('*italic(T)[italic('pk')]^italic(b)['max']*', '*italic(T)[italic('pk')]^italic(z)*')'),
				expression('corr('*italic(T)[italic('pk')]^italic(z)['j']*', '*italic(T)[italic('pk')]^italic(z)*')')
			),
			expand = c(0.105,0)
		) +
		geom_hline(yintercept = 0, linetype = 'dashed', color = 'dimgray', lwd = 0.4) +
		theme(
			plot.title = element_text(margin = margin(b = -10), size = 8),
			axis.line = element_line(colour = NA),
			axis.title.y = element_text(size = 7, family = 'serif'),
			axis.text.x = element_text(size = 7, family = 'serif'),
			axis.text.y = element_text(size = 7, family = 'serif'),
			plot.margin = unit(c(0.2,0.2,-0.55,0.05), 'cm'),
			panel.border = element_rect(colour = "black", fill=NA, size=0.5)
		)	

	return(p_cor)
}

# This function reads the estimates of the evolutionary rate of each 
# T_pk and generates a plot.
plot_evolutionary_rates <- function(sampled_vars)
{
	data_to_plot <- data.frame(
		parameter = c('alpha', 'beta', 'z_j', 'z'),
		median_value = colMedians(
			as.matrix(
				sampled_vars[
					,c('var_alpha', 'var_beta', 'var_z_j', 'var_z')
				]
			)
		),
		lwr_value = HPDinterval(
			mcmc(
				sampled_vars[
					,c('var_alpha', 'var_beta', 'var_z_j', 'var_z')
				]
			)
		)[,1],
		upr_value = HPDinterval(
			mcmc(
				sampled_vars[
					,c('var_alpha', 'var_beta', 'var_z_j', 'var_z')
				]
			)
		)[,2]
	)
	
	data_to_plot_2 <- data.frame(
		parameter = c(
			rep('alpha', nrow(sampled_vars)),
			rep('beta', nrow(sampled_vars)),
			rep('z_j', nrow(sampled_vars)),
			rep('z', nrow(sampled_vars))
		),
		value = c(
			sampled_vars$var_alpha,	sampled_vars$var_beta,
			sampled_vars$var_z_j,		sampled_vars$var_z
		)
	)
	data_to_plot_2$parameter <- factor(
		data_to_plot_2$parameter, levels = c('alpha', 'beta', 'z_j', 'z')
	)
	
	# Here we only keep the values that lie within the 95% HPD interval.
	data_to_plot_2 <- data_to_plot_2[
		(
			data_to_plot_2$parameter == 'alpha' & 
			data_to_plot_2$value <= data_to_plot['var_alpha', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['var_alpha', 'lwr_value']
		) |
		(
			data_to_plot_2$parameter == 'beta' & 
			data_to_plot_2$value <= data_to_plot['var_beta', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['var_beta', 'lwr_value']
		) |
		(
			data_to_plot_2$parameter == 'z_j' & 
			data_to_plot_2$value <= data_to_plot['var_z_j', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['var_z_j', 'lwr_value']
		) |
		(
			data_to_plot_2$parameter == 'z' & 
			data_to_plot_2$value <= data_to_plot['var_z', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['var_z', 'lwr_value']
		),		
	]
		
	p_rate <- ggplot(data_to_plot_2, aes(x = parameter, y = value)) +
		geom_violin(
			position = 'dodge', fill = 'goldenrod1', width = 1,
			color = 'goldenrod1'
		) + 
		ggtitle(' C') +
		geom_point(
			data = data_to_plot, aes(x = parameter, y = median_value), 
			pch = 21, fill = 'white', size = 3
		) +
		scale_y_continuous(
			name = bquote('Evolutionary rate (C'^{2} * ' · myr'^{-1}*')'),
			expand = c(0.005,0), breaks = seq(0, 140, 35), limits = c(0, 140)
		) +
		scale_x_discrete(
			name = '',
			labels = c(
				expression(italic(T)[italic('pk')]^alpha),
				expression(italic(T)[italic('pk')]^'b'['max']),
				expression(italic(T)[italic('pk')]^italic('z')['j']),
				expression(italic(T)[italic('pk')]^italic('z'))
			),
			expand = c(0.175,0)
		) +
		theme(
			plot.title = element_text(margin = margin(b = -10), size = 8),
			axis.line = element_line(colour = NA),
			axis.title.y = element_text(size = 7, family = 'serif'),
			axis.text.x = element_text(size = 7, family = 'serif'),
			axis.text.y = element_text(size = 7, family = 'serif'),
			plot.margin = unit(c(0.2,0.2,-0.55,0.05), 'cm'),
			panel.border = element_rect(colour = "black", fill=NA, size=0.5)
		)	

	return(p_rate)
}

# This function reads the estimates of the phylogenetic heritability of 
# each T_pk and generates a plot.
plot_heritabilities <- function(sampled_vars)
{
	data_to_plot <- data.frame(
		parameter = c('alpha', 'beta', 'z_j', 'z'),
		median_value = colMedians(
			as.matrix(
				sampled_vars[
					,c('her_alpha', 'her_beta', 'her_z_j', 'her_z')
				]
			)
		),
		lwr_value = HPDinterval(
			mcmc(
				sampled_vars[
					,c('her_alpha', 'her_beta', 'her_z_j', 'her_z')
				]
			)
		)[,1],
		upr_value = HPDinterval(
			mcmc(
				sampled_vars[
					,c('her_alpha', 'her_beta', 'her_z_j', 'her_z')
				]
			)
		)[,2]
	)
	
	data_to_plot_2 <- data.frame(
		parameter = c(
			rep('alpha', nrow(sampled_vars)),
			rep('beta', nrow(sampled_vars)),
			rep('z_j', nrow(sampled_vars)),
			rep('z', nrow(sampled_vars))
		),
		value = c(
			sampled_vars$her_alpha,	sampled_vars$her_beta,
			sampled_vars$her_z_j,	sampled_vars$her_z
		)
	)
	data_to_plot_2$parameter <- factor(
		data_to_plot_2$parameter, levels = c('alpha', 'beta', 'z_j', 'z')
	)
	
	# Here we only keep the values that lie within the 95% HPD interval.
	data_to_plot_2 <- data_to_plot_2[
		(
			data_to_plot_2$parameter == 'alpha' & 
			data_to_plot_2$value <= data_to_plot['her_alpha', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['her_alpha', 'lwr_value']
		) |
		(
			data_to_plot_2$parameter == 'beta' & 
			data_to_plot_2$value <= data_to_plot['her_beta', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['her_beta', 'lwr_value']
		) |
		(
			data_to_plot_2$parameter == 'z_j' & 
			data_to_plot_2$value <= data_to_plot['her_z_j', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['her_z_j', 'lwr_value']
		) |
		(
			data_to_plot_2$parameter == 'z' & 
			data_to_plot_2$value <= data_to_plot['her_z', 'upr_value'] &
			data_to_plot_2$value >= data_to_plot['her_z', 'lwr_value']
		),		
	]
		
	p_herit <- ggplot(data_to_plot_2, aes(x = parameter, y = value)) +
		geom_violin(
			position = 'dodge', fill = 'goldenrod1', width = 1,
			color = 'goldenrod1'
		) + 
		ggtitle(' B') +
		geom_point(
			data = data_to_plot, aes(x = parameter, y = median_value), 
			pch = 21, fill = 'white', size = 3
		) +
		scale_y_continuous(
			name = 'Phylogenetic heritability', limits = c(0,1),
			expand = c(0.005,0), breaks = seq(0, 1, 0.2)
		) +
		scale_x_discrete(
			name = '',
			labels = c(
				expression(italic(T)[italic('pk')]^alpha),
				expression(italic(T)[italic('pk')]^'b'['max']),
				expression(italic(T)[italic('pk')]^italic('z')['j']),
				expression(italic(T)[italic('pk')]^italic('z'))
			),
			expand = c(0,0)
		) +
		theme(
			plot.title = element_text(margin = margin(b = -10), size = 8),
			axis.line = element_line(colour = NA),
			axis.title.y = element_text(size = 7, family = 'serif'),
			axis.text.x = element_text(size = 7, family = 'serif'),
			axis.text.y = element_text(size = 7, family = 'serif'),
			plot.margin = unit(c(0.2,0.2,-0.55,0.05), 'cm'),
			panel.border = element_rect(colour = "black", fill=NA, size=0.5)
		)	

	return(p_herit)
}

####################
# M A I N  C O D E #
####################

# Read all the MCMCglmm fits, combine their posterior samples, and 
# generate plots.
make_plots(
	combine_samples(
		list.files(
			path = '../results/MCMCglmm_runs/', pattern = '\\d+.Rda'
		)
	)
)
