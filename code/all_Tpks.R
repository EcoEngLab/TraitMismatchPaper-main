# load packages

require('tidyverse')
require('nls.multstart')
require('broom')
require('rTPC')
require('data.table')
require('car')
require('patchwork')
require('minpack.lm')
require('boot')
require('doMC')
require('foreach')
require('ggpubr')
require('ggtext')
require('cowplot')

###Fitting TPC for all species for all traits

rm(list=ls())
graphics.off()

######Alpha######

df <- as_tibble(read.csv('mismatch_dat.csv'))

df<- df %>%
  rename(temp = ambienttemp,species = interactor1)

dv <- df %>% filter(standardisedtraitname == "Juvenile Development Rate")

dv <- dv %>% select(curve_ID, species, temp, standardisedtraitname, standardisedtraitvalue, stage)

dv <- dv %>% rename(rate = standardisedtraitvalue)

dvzero      <- dv %>% filter(rate == 0)
dvzero$rate <- 10^-6

dv          <- dv %>% filter(rate != 0)

dv <- rbind(dv,dvzero)

dv <- dv %>% arrange(dv, curve_ID)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# fit TPC model for each species

start_vals <- get_start_vals(dv$temp, dv$rate, model_name = 'pawar_2018')
low_lims <- get_lower_lims(dv$temp, dv$rate, model_name = 'pawar_2018')
upper_lims <- get_upper_lims(dv$temp, dv$rate, model_name = 'pawar_2018')


dv_fits <- nest(dv, data = c(temp, rate)) %>%
  mutate(fit = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                        data = .x,
                                        iter = c(3,3,3,3),
                                        start_lower = start_vals - 10,
                                        start_upper = start_vals + 10,
                                        lower = low_lims,
                                        upper=upper_lims,
                                        supp_errors = 'Y',
                                        convergence_count = FALSE)))


dv_preds <- mutate(dv_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100)))) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(fit)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(curve_ID, species, standardisedtraitname, stage, preds) %>%
  # unlist the preds list column
  unnest(preds)

glimpse(dv_preds)

ggplot(dv_preds) +
  geom_line(aes(temp, .fitted, col = curve_ID)) +
  geom_point(aes(temp, rate), dv) +
  facet_wrap(~species, scales = 'free_y', ncol = 6) +
  theme_bw() +
  theme(legend.position = 'none')

#get good names with line break
dv$GoodName <- str_replace(pattern = " ",replacement =  "\n", dv$species)

#fit models
FitModel <- function(ID){
  # browser()
  ID <<- ID
  df <- dv %>% filter(curve_ID == ID)
  df <<-df
  
  Model <- tryCatch({minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                       data = df,
                                       start = coef(dv_fits$fit[[which(dv_fits$curve_ID==ID)]]),
                                       lower = get_lower_lims(df$temp, df$rate, model_name = 'pawar_2018'),
                                       upper = get_upper_lims(df$temp, df$rate, model_name = 'pawar_2018'),
                                       weights = rep(1, times = nrow(df)))
  }, error = function(error){
    print(paste("error in:", ID, ",", error))
    minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                      data = df,
                      start = coef(dv_fits$fit[[which(dv_fits$curve_ID==ID)]]),
                      weights = rep(1, times = nrow(df)))
    
  })
  extra_params <- calc_params(Model) %>%
    pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
  
  ci_extra_params <- Boot(Model, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(Model)), R = 200, method = 'residual')%>%
    confint(., method = 'bca')%>%
    as.data.frame() %>%
    rename(conf_lower = 1, conf_upper = 2) %>%
    rownames_to_column(., var = 'param') %>%
    mutate(method = 'residual bootstrap')
  
  ci_extra_params <- left_join(ci_extra_params, extra_params)
  
  ##Get other parameters
  params <- broom::tidy(Model) %>% select(param = term, estimate)
  BootOut <- Boot(Model, method = 'residual')
  ## Get the param Names that has multiple values:
  paramName <- colnames(BootOut[[2]])[which(apply(BootOut[[2]], 2, function(x) length(unique(x))>1))]
  
  
  params_cis <- BootOut %>%
    confint(.,parm = paramName, method = 'bca') %>%
    as.data.frame() %>%
    rename(conf_lower = 1, conf_upper = 2) %>%
    rownames_to_column(., var = 'param') %>%
    mutate(method = 'residual bootstrap')
  
  params_cis <- bind_rows(params_cis) %>%
    left_join(., params)%>%
    filter(param!="topt")
  
  topt <- as_tibble(rbind(ci_extra_params,params_cis))
  topt$species <- as.character(df$species[1])
  
  #Plot fit
  Boot_conf <- BootOut$t %>%
    as.data.frame() %>%
    drop_na() %>%
    mutate(iter = 1:n()) %>%
    group_by_all() %>%
    do(data.frame(temp = seq(min(df$temp), max(df$temp), length.out = 100))) %>%
    ungroup() %>%
    mutate(pred =pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15))
  
  # calculate bootstrapped confidence intervals
  boot_conf_preds <- group_by(Boot_conf, temp) %>%
    summarise(conf_lower = quantile(pred, 0.025),
              conf_upper = quantile(pred, 0.975)) %>%
    ungroup()
  
  
  ##plot
  counter <- which(dv_fits$curve_ID==ID)
  print(counter)
  
  plotData <- filter(dv_preds, dv_preds$curve_ID==ID)
  plot <- ggplot(data=df, aes(x=temp, y=rate))+
    geom_point()+
    geom_line(data=plotData, mapping = aes(x=temp, y=.fitted)) +
    geom_ribbon(aes(temp, ymin=conf_lower, ymax=conf_upper), boot_conf_preds, fill="#e66101",alpha=0.3,
                inherit.aes = F)+
    theme_bw()+
    theme(text = element_text(size=10))+
    labs(title=paste(df$GoodName[1], sep=""),
         y="Development rate", x=expression(plain(paste(" Temperature, ",degree,"C"))))
  
  ggsave(plot,file=paste("../results/TPC/Alpha_",ID,".pdf",sep=""), 
         height=10,width=15,units="cm")
  
  AlphaPlots[[counter]] <<- plot
  
  return(topt)
}

##Run everything:
#run using apply:
AlphaPlots <- vector(mode="list", length=length(unique(dv$curve_ID)))

ModelOut <- sapply(unique(dv$curve_ID), FitModel)
# ModelOutDFList <- apply(ModelOut, 2, function(x) as.data.frame(do.call(cbind, x)))
ModelOutDF <- do.call(rbind, ModelOut)
ModelOutDF$trait <- "juvenile development rate"

#run in parallel:
# doMC::registerDoMC(cores = 4)
# ModelOutList <- foreach(ID = unique(dv$curve_ID)) %dopar%{ FitModel(ID)}
# 
# ModelOutDF <- do.call(rbind, ModelOutList)
# ModelOutDF$trait <- "juvenile development rate"

##Get plots:
alphaPlot <- plot_grid(plotlist = AlphaPlots)

save_plot(alphaPlot, file="../results/AlphaFits.pdf", 
          base_height=28, base_asp = 1,units="cm")




##### Peak fecundity #####
rm(list=ls())
graphics.off() 

# setwd("~/Dropbox/ph_thesis/Topt_paper/data")

#read in the trait data 

df <- as_tibble(read.csv('mismatch_dat.csv'))

df<- df %>%
  rename(temp = ambienttemp,species = interactor1)

dv <- df %>% filter(standardisedtraitname == 'Fecundity Rate')

dv <- dv %>% select(curve_ID, species, temp, standardisedtraitname, standardisedtraitvalue, stage)

dv <- dv %>% rename(rate = standardisedtraitvalue)

dv$rate <- dv$rate + 0.001

dv <- dv %>% arrange(dv, curve_ID)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# fit TPC model for each species

start_vals <- get_start_vals(dv$temp, dv$rate, model_name = 'pawar_2018')
low_lims <- get_lower_lims(dv$temp, dv$rate, model_name = 'pawar_2018')
upper_lims <- get_upper_lims(dv$temp, dv$rate, model_name = 'pawar_2018')


dv_fits <- nest(dv, data = c(temp, rate)) %>%
  mutate(fit = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                        data = .x,
                                        iter = c(3,3,3,3),
                                        start_lower = start_vals - 10,
                                        start_upper = start_vals + 10,
                                        lower = low_lims,
                                        upper=upper_lims,
                                        supp_errors = 'Y',
                                        convergence_count = FALSE)))


dv_preds <- mutate(dv_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100)))) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(fit)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(curve_ID, species, standardisedtraitname, stage, preds) %>%
  # unlist the preds list column
  unnest(preds)

glimpse(dv_preds)

ggplot(dv_preds) +
  geom_line(aes(temp, .fitted, col = curve_ID)) +
  geom_point(aes(temp, rate), dv) +
  facet_wrap(~species, scales = 'free_y', ncol = 6) +
  theme_bw() +
  theme(legend.position = 'none')

#name with line break
dv$GoodName <- str_replace(pattern = " ",replacement =  "\n", dv$species)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

######vectorized version########
FitModel <- function(ID){
  # browser()
  ID <<- ID
  df <- dv %>% filter(curve_ID == ID)
  df <<-df
  
  Model <- tryCatch({minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                       data = df,
                                       start = coef(dv_fits$fit[[which(dv_fits$curve_ID==ID)]]),
                                       lower = get_lower_lims(df$temp, df$rate, model_name = 'pawar_2018'),
                                       upper = get_upper_lims(df$temp, df$rate, model_name = 'pawar_2018'),
                                       weights = rep(1, times = nrow(df)))
  }, error = function(error){
    print(paste("error in:", ID, ",", error))
    minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                      data = df,
                      start = coef(dv_fits$fit[[which(dv_fits$curve_ID==ID)]]),
                      weights = rep(1, times = nrow(df)))
    
  })
  
  extra_params <- calc_params(Model) %>%
    pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
  
  ci_extra_params <- Boot(Model, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(Model)), R = 200, method = 'residual')%>%
    confint(., method = 'bca')%>%
    as.data.frame() %>%
    rename(conf_lower = 1, conf_upper = 2) %>%
    rownames_to_column(., var = 'param') %>%
    mutate(method = 'residual bootstrap')
  
  ci_extra_params <- left_join(ci_extra_params, extra_params)
  
  ##Get other parameters
  params <- broom::tidy(Model) %>% select(param = term, estimate)
  BootOut <- Boot(Model, method = 'residual')
  ## Get the param Names that has multiple values:
  paramName <- colnames(BootOut[[2]])[which(apply(BootOut[[2]], 2, function(x) length(unique(x))>1))]
  
  
  params_cis <- BootOut %>%
    confint(.,parm = paramName, method = 'bca') %>%
    as.data.frame() %>%
    rename(conf_lower = 1, conf_upper = 2) %>%
    rownames_to_column(., var = 'param') %>%
    mutate(method = 'residual bootstrap')
  
  params_cis <- bind_rows(params_cis) %>%
    left_join(., params)%>%
    filter(param!="topt")
  
  
  topt <- as_tibble(rbind(ci_extra_params,params_cis))
  topt$species <- as.character(df$species[1])
  
  ##Plot fit##
  Boot_conf <- BootOut$t %>%
    as.data.frame() %>%
    drop_na() %>%
    mutate(iter = 1:n()) %>%
    group_by_all() %>%
    do(data.frame(temp = seq(min(df$temp), max(df$temp), length.out = 100))) %>%
    ungroup() %>%
    mutate(pred =pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15))
  
  # calculate bootstrapped confidence intervals
  boot_conf_preds <- group_by(Boot_conf, temp) %>%
    summarise(conf_lower = quantile(pred, 0.025),
              conf_upper = quantile(pred, 0.975)) %>%
    ungroup()
  
  
  plotData <- filter(dv_preds, dv_preds$curve_ID==ID)
  plot <- ggplot(data=df, aes(x=temp, y=rate))+
    geom_point()+
    geom_line(data=plotData, mapping = aes(x=temp, y=.fitted)) +
    geom_ribbon(aes(temp, ymin=conf_lower, ymax=conf_upper), boot_conf_preds, fill="#fdb863",alpha=0.3,
                inherit.aes = F)+
    theme_bw()+
    theme(text = element_text(size=10))+
    labs(title=paste(df$GoodName[1], sep=""),         
         y="Fecundity",x=expression(plain(paste(" Temperature, ",degree,"C"))))
  
  counter <- which(dv_fits$curve_ID==ID)
  print(counter)
  BetaPlots[[counter]] <<- plot
  
  
  ggsave(plot,file=paste("../results/TPC/Bopt_",ID,".pdf",sep=""), 
         height=10,width=15,units="cm")
  
  return(topt)
}

##Run everything:
# doMC::registerDoMC(cores = 4)  
# ModelOut <- foreach(ID = unique(dv$curve_ID)) %dopar%{ FitModel(ID)}
BetaPlots <- vector(mode="list", length=length(unique(dv$curve_ID)))

ModelOut <- sapply(unique(dv$curve_ID), FitModel)

##Get plots:
BetaPlot <- plot_grid(plotlist = BetaPlots, ncol=5)

save_plot(BetaPlot, file="../results/BetaFits.pdf", 
          base_height=14, base_width = 28,units="cm")



# ModelOutDFList <- apply(ModelOut, 2, function(x) as.data.frame(do.call(cbind, x)))
ModelOutDF <- do.call(rbind, ModelOut)
ModelOutDF$trait <- "fecundity"
write.csv(ModelOutDF, "../data/bpk_Tpks_AllParams.csv")

##### Adult Mortality #####
rm(list=ls())
graphics.off()

# setwd("~/Dropbox/ph_thesis/Topt_paper/data")

#read in and prepare the trait data 

df <- as_tibble(read.csv('mismatch_dat.csv'))

df<- df %>%
  rename(temp = ambienttemp,species = interactor1)

dv <- df %>% filter(standardisedtraitname == "Female Mortality Rate")

dv <- dv %>% select(curve_ID, species, temp, standardisedtraitname, standardisedtraitvalue, stage)

dv <- dv %>% rename(rate = standardisedtraitvalue)

dv$rate <- 1/dv$rate

dv <- dv %>% arrange(dv, curve_ID)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# fit TPC model for each species

start_vals <- get_start_vals(dv$temp, dv$rate, model_name = 'pawar_2018')
low_lims <- get_lower_lims(dv$temp, dv$rate, model_name = 'pawar_2018')
upper_lims <- get_upper_lims(dv$temp, dv$rate, model_name = 'pawar_2018')


dv_fits <- nest(dv, data = c(temp, rate)) %>%
  mutate(fit = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                        data = .x,
                                        iter = c(3,3,3,3),
                                        start_lower = start_vals - 10,
                                        start_upper = start_vals + 10,
                                        lower = low_lims,
                                        upper=upper_lims,
                                        supp_errors = 'Y',
                                        convergence_count = FALSE)))


dv_preds <- mutate(dv_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100)))) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(fit)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(curve_ID, species, standardisedtraitname, stage, preds) %>%
  # unlist the preds list column
  unnest(preds)

glimpse(dv_preds)

ggplot(dv_preds) +
  geom_line(aes(temp, .fitted, col = curve_ID)) +
  geom_point(aes(temp, rate), dv) +
  facet_wrap(~species, scales = 'free_y', ncol = 6) +
  theme_bw() +
  theme(legend.position = 'none')


#name with line break
dv$GoodName <- str_replace(pattern = " ",replacement =  "\n", dv$species)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

######vectorized version########
FitModel <- function(ID){
  # browser()
  ID <<- ID
  df <- dv %>% filter(curve_ID == ID)
  df <<-df
  
  Model <- tryCatch({minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                       data = df,
                                       start = coef(dv_fits$fit[[which(dv_fits$curve_ID==ID)]]),
                                       lower = get_lower_lims(df$temp, df$rate, model_name = 'pawar_2018'),
                                       upper = get_upper_lims(df$temp, df$rate, model_name = 'pawar_2018'),
                                       weights = rep(1, times = nrow(df)))
  }, error = function(error){
    print(paste("error in:", ID, ",", error))
    minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                      data = df,
                      start = coef(dv_fits$fit[[which(dv_fits$curve_ID==ID)]]),
                      weights = rep(1, times = nrow(df)))
    
  })
  
  extra_params <- calc_params(Model) %>%
    pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
  
  ci_extra_params <- Boot(Model, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(Model)), R = 200, method = 'residual')%>%
    confint(., method = 'bca')%>%
    as.data.frame() %>%
    rename(conf_lower = 1, conf_upper = 2) %>%
    rownames_to_column(., var = 'param') %>%
    mutate(method = 'residual bootstrap')
  
  ci_extra_params <- left_join(ci_extra_params, extra_params)
  
  ##Get other parameters
  params <- broom::tidy(Model) %>% select(param = term, estimate)
  BootOut <- Boot(Model, method = 'residual')
  ## Get the param Names that has multiple values:
  paramName <- colnames(BootOut[[2]])[which(apply(BootOut[[2]], 2, function(x) length(unique(x))>1))]
  
  
  params_cis <- BootOut %>%
    confint(.,parm = paramName, method = 'bca') %>%
    as.data.frame() %>%
    rename(conf_lower = 1, conf_upper = 2) %>%
    rownames_to_column(., var = 'param') %>%
    mutate(method = 'residual bootstrap')
  
  params_cis <- bind_rows(params_cis) %>%
    left_join(., params)%>%
    filter(param!="topt")
  
  
  topt <- as_tibble(rbind(ci_extra_params,params_cis))
  topt$species <- as.character(df$species[1])
  
  ##Plot fit##
  Boot_conf <- BootOut$t %>%
    as.data.frame() %>%
    drop_na() %>%
    mutate(iter = 1:n()) %>%
    group_by_all() %>%
    do(data.frame(temp = seq(min(df$temp), max(df$temp), length.out = 100))) %>%
    ungroup() %>%
    mutate(pred =pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15))
  
  # calculate bootstrapped confidence intervals
  boot_conf_preds <- group_by(Boot_conf, temp) %>%
    summarise(conf_lower = quantile(pred, 0.025),
              conf_upper = quantile(pred, 0.975)) %>%
    ungroup()
  
  
  plotData <- filter(dv_preds, dv_preds$curve_ID==ID)
  plot <- ggplot(data=df, aes(x=temp, y=rate))+
    geom_point()+
    geom_line(data=plotData, mapping = aes(x=temp, y=.fitted)) +
    geom_ribbon(aes(temp, ymin=conf_lower, ymax=conf_upper), boot_conf_preds, fill="#a6cee3",alpha=0.3,
                inherit.aes = F)+
    theme_bw()+
    theme(text = element_text(size=10))+
    labs(title=paste(df$GoodName[1], sep=""),         
         y="1/Adult Mortality",x=expression(plain(paste(" Temperature, ",degree,"C"))))
  
  ggsave(plot,file=paste("../results/TPC/z_",ID,".pdf",sep=""), 
         height=10,width=15,units="cm")
  
  
  ##invert
  plotIN <- ggplot(data=df, aes(x=temp, y=1/rate))+
    geom_point()+
    geom_line(data=plotData, mapping = aes(x=temp, y=1/.fitted)) +
    geom_ribbon(aes(temp, ymin=1/conf_lower, ymax=1/conf_upper), boot_conf_preds, fill="#a6cee3",alpha=0.3,
                inherit.aes = F)+
    theme_bw()+
    labs(title=paste(df$GoodName[1],", e =", round(topt$estimate[which(topt$param=="e")],digits = 2),
                     ", Topt= ", round(topt$estimate[which(topt$param=="topt")],digits = 2), sep=""),         
         y="Adult Mortality",x=expression(plain(paste(" Temperature, ",degree,"C"))))
  
  
  counter <- which(dv_fits$curve_ID==ID)
  print(counter)
  ZetaPlots[[counter]] <<- plot
  
  
  ggsave(plotIN,file=paste("../results/TPC/IN_z_",ID,".pdf",sep=""), 
         height=10,width=15,units="cm")
  
  return(topt)
}

##Run everything:
ZetaPlots <- vector(mode="list", length=length(unique(dv$curve_ID)))

ModelOut <- sapply(unique(dv$curve_ID), FitModel)

##Get plots:
ZetaPlot <- plot_grid(plotlist = ZetaPlots, ncol=5)

save_plot(ZetaPlot, file="../results/ZetaFits.pdf", 
          base_height=21, base_width = 28,units="cm")



# ModelOut <- sapply(unique(dv$curve_ID), FitModel)
# ModelOutDFList <- apply(ModelOut, 2, function(x) as.data.frame(do.call(cbind, x)))
# doMC::registerDoMC(cores = 4)  
# ModelOutDFList <- foreach(ID = unique(dv$curve_ID)) %dopar%{ FitModel(ID)}

ModelOutDF <- do.call(rbind, ModelOutDFList)
ModelOutDF$trait <- "adult mortality rate"
write.csv(ModelOutDF, "../data/z_Tpks_AllParams.csv")

##### Juvenile Mortality #####
rm(list=ls())
graphics.off()

# setwd("~/Dropbox/ph_thesis/Topt_paper/data")

#read in the trait data 

df <- as_tibble(read.csv('mismatch_dat.csv'))

df<- df %>%
  rename(temp = ambienttemp,species = interactor1)

dv <- df %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

dv <- dv %>% select(curve_ID, species, temp, standardisedtraitname, standardisedtraitvalue, stage)

dv <- dv %>% rename(rate = standardisedtraitvalue)

dv$rate <- 1/dv$rate

dv <- dv %>% arrange(dv, curve_ID)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# fit TPC model for each species

start_vals <- get_start_vals(dv$temp, dv$rate, model_name = 'pawar_2018')
low_lims <- get_lower_lims(dv$temp, dv$rate, model_name = 'pawar_2018')
upper_lims <- get_upper_lims(dv$temp, dv$rate, model_name = 'pawar_2018')

dv_fits <- nest(dv, data = c(temp, rate)) %>%
  mutate(fit = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                        data = .x,
                                        iter = c(3,3,3,3),
                                        start_lower = start_vals - 10,
                                        start_upper = start_vals + 10,
                                        lower = low_lims,
                                        upper=upper_lims,
                                        supp_errors = 'Y',
                                        convergence_count = FALSE)))


dv_preds <- mutate(dv_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100)))) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(fit)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(curve_ID, species, standardisedtraitname, stage, preds) %>%
  # unlist the preds list column
  unnest(preds)

glimpse(dv_preds)

ggplot(dv_preds) +
  geom_line(aes(temp, .fitted, col = curve_ID)) +
  geom_point(aes(temp, rate), dv) +
  facet_wrap(~species, scales = 'free_y', ncol = 6) +
  theme_bw() +
  theme(legend.position = 'none')


#name with line break
dv$GoodName <- str_replace(pattern = " ",replacement =  "\n", dv$species)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

######vectorized version########
FitModel <- function(ID){
  # browser()
  ID <<- ID
  df <- dv %>% filter(curve_ID == ID)
  df <<-df
  
  Model <- tryCatch({minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                       data = df,
                                       start = coef(dv_fits$fit[[which(dv_fits$curve_ID==ID)]]),
                                       lower = get_lower_lims(df$temp, df$rate, model_name = 'pawar_2018'),
                                       upper = get_upper_lims(df$temp, df$rate, model_name = 'pawar_2018'),
                                       weights = rep(1, times = nrow(df)))
  }, error = function(error){
    print(paste("error in:", ID, ",", error))
    minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                      data = df,
                      start = coef(dv_fits$fit[[which(dv_fits$curve_ID==ID)]]),
                      weights = rep(1, times = nrow(df)))
    
  })
  
  extra_params <- calc_params(Model) %>%
    pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
  
  ci_extra_params <- Boot(Model, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(Model)), R = 200, method = 'residual')%>%
    confint(., method = 'bca')%>%
    as.data.frame() %>%
    rename(conf_lower = 1, conf_upper = 2) %>%
    rownames_to_column(., var = 'param') %>%
    mutate(method = 'residual bootstrap')
  
  ci_extra_params <- left_join(ci_extra_params, extra_params)
  
  ##Get other parameters
  params <- broom::tidy(Model) %>% select(param = term, estimate)
  BootOut <- Boot(Model, method = 'residual')
  ## Get the param Names that has multiple values:
  paramName <- colnames(BootOut[[2]])[which(apply(BootOut[[2]], 2, function(x) length(unique(x))>1))]
  
  
  params_cis <- BootOut %>%
    confint(.,parm = paramName, method = 'bca') %>%
    as.data.frame() %>%
    rename(conf_lower = 1, conf_upper = 2) %>%
    rownames_to_column(., var = 'param') %>%
    mutate(method = 'residual bootstrap')
  
  params_cis <- bind_rows(params_cis) %>%
    left_join(., params)%>%
    filter(param!="topt")
  
  
  topt <- as_tibble(rbind(ci_extra_params,params_cis))
  topt$species <- as.character(df$species[1])
  
  ##Plot fit##
  Boot_conf <- BootOut$t %>%
    as.data.frame() %>%
    drop_na() %>%
    mutate(iter = 1:n()) %>%
    group_by_all() %>%
    do(data.frame(temp = seq(min(df$temp), max(df$temp), length.out = 100))) %>%
    ungroup() %>%
    mutate(pred =pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15))
  
  # calculate bootstrapped confidence intervals
  boot_conf_preds <- group_by(Boot_conf, temp) %>%
    summarise(conf_lower = quantile(pred, 0.025),
              conf_upper = quantile(pred, 0.975)) %>%
    ungroup()
  
  
  plotData <- filter(dv_preds, dv_preds$curve_ID==ID)
  plot <- ggplot(data=df, aes(x=temp, y=rate))+
    geom_point()+
    geom_line(data=plotData, mapping = aes(x=temp, y=.fitted)) +
    geom_ribbon(aes(temp, ymin=conf_lower, ymax=conf_upper), boot_conf_preds, fill="#1f78b4",alpha=0.3,
                inherit.aes = F)+
    theme_bw()+
    labs(title=paste(df$GoodName[1],sep=""),         
         y="1/Juvenlie Mortality",x=expression(plain(paste(" Temperature, ",degree,"C"))))
  
  ggsave(plot,file=paste("../results/TPC/zj_",ID,".pdf",sep=""), 
         height=10,width=15,units="cm")
  
  ##invert
  plotIN <- ggplot(data=df, aes(x=temp, y=1/rate))+
    geom_point()+
    geom_line(data=plotData, mapping = aes(x=temp, y=1/.fitted)) +
    geom_ribbon(aes(temp, ymin=1/conf_lower, ymax=1/conf_upper), boot_conf_preds, fill="#a6cee3",alpha=0.3,
                inherit.aes = F)+
    theme_bw()+
    theme(text = element_text(size=10))+
    labs(title=paste(df$GoodName[1],sep=""),         
         y="Adult Mortality",x=expression(plain(paste(" Temperature, ",degree,"C"))))
  
  ggsave(plotIN,file=paste("../results/TPC/IN_zj_",ID,".pdf",sep=""), 
         height=10,width=15,units="cm")
  
  counter <- which(dv_fits$curve_ID==ID)
  print(counter)
  ZetaJPlots[[counter]] <<- plot
  
  return(topt)
}

##Run everything:
ZetaJPlots <- vector(mode="list", length=length(unique(dv$curve_ID)))

ModelOut <- sapply(unique(dv$curve_ID), FitModel)

##Get plots:
ZetaJPlot <- plot_grid(plotlist = ZetaJPlots, ncol=5)

save_plot(ZetaJPlot, file="../results/ZetaJFits.pdf", 
          base_height=21, base_width = 28,units="cm")

# ModelOut <- sapply(unique(dv$curve_ID), FitModel)
# ModelOutDFList <- apply(ModelOut, 2, function(x) as.data.frame(do.call(cbind, x)))
# doMC::registerDoMC(cores = 4)  
# ModelOutDFList <- foreach(ID = unique(dv$curve_ID)) %dopar%{ FitModel(ID)}

ModelOutDF <- do.call(rbind, ModelOutDFList)
ModelOutDF$trait <- "juvenile mortality rate"
write.csv(ModelOutDF, "../data/zj_Tpks_AllParams.csv")


