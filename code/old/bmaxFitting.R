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

###### 2. bmax (daily fecundity rate) ######

df <- as_tibble(read.csv('TraitData.csv')) 

dv <- df %>% rename(temp = interactor1temp, species = interactor1, rate = standardisedtraitvalue) %>%
  select(species, temp, standardisedtraitname, rate) %>%
  filter(standardisedtraitname == "bmax", rate != 'NA') %>% 
  mutate(curve_ID = case_when(species == 'Aedes albopictus' ~ '1',
                              species == 'Aedes aegypti' ~ '2',
                              species == 'Anthonomus grandis' ~ '3',
                              species == 'Paracoccus marginatus' ~ '4',
                              species == 'Acyrthosiphon pisum' ~ '5',
                              species == 'Aphis gossypii' ~ '6',
                              species == 'Harmonia axyridis' ~ '7',
                              species == 'Tribolium castaneum' ~ '8',
                              species == 'Aedes krombeini' ~ '9',
                              species == 'Bemisia tabaci' ~ '10',
                              species == 'Tetraneura nigriabdominalis' ~ '11',
                              species == 'Stethorus punctillum' ~ '12',
                              species == 'Tetranychus mcdanieli' ~ '13',
                              species == 'Tetranychus urticae' ~ '14',
                              species == 'Clavigralla tomentosicollis' ~ '15',
                              species == 'Planococcus citri' ~ '16',
                              species == 'Muscidifurax zaraptor' ~ '17',
                              species == 'Aphis nasturtii' ~ '18',
                              species == 'Rhopalosiphum maidis' ~ '19')) %>%
  arrange(curve_ID) %>% mutate(temp = as.numeric(temp))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# fit TPC model for each species

start_vals <- get_start_vals(dv$temp, dv$rate, model_name = 'pawar_2018')
low_lims <-   get_lower_lims(dv$temp, dv$rate, model_name = 'pawar_2018')
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


dv_preds <- mutate(dv_fits, new_data = map(data, ~tibble(temp = seq(min(0), max(45), length.out = 1000)))) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(fit)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(curve_ID, species, standardisedtraitname, preds) %>%
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
    do(data.frame(temp = seq(min(0), max(45), length.out = 1000))) %>%
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
    geom_point(size=0.5, alpha=0.4)+
    geom_line(data=plotData, mapping = aes(x=temp, y=.fitted)) +
    geom_ribbon(aes(temp, ymin=conf_lower, ymax=conf_upper), boot_conf_preds, fill="#e66101",alpha=0.3,
                inherit.aes = F)+
    theme_bw()+
    theme(text = element_text(size=10))+
    labs(title=paste(df$GoodName[1], sep=""),
         y="bmax", x=expression(plain(paste(" Temperature, ",degree,"C"))))
  
  ggsave(plot,file=paste("../results/TPC/Bopt_",ID,".pdf",sep=""), 
         height=10,width=15,units="cm")
  
  BetaPlots[[counter]] <<- plot
  
  return(topt)
}

##Run everything:
#run using apply:
BetaPlots <- vector(mode="list", length=length(unique(dv$curve_ID)))

ModelOut       <- sapply(unique(dv$curve_ID), FitModel)
ModelOutDFList <- apply(ModelOut, 2, function(x) as.data.frame(do.call(cbind, x)))
ModelOutDF     <- do.call(rbind, ModelOut)
ModelOutDF$trait <- "fecundity rate"

#run in parallel:
doMC::registerDoMC(cores = 4)
ModelOutList <- foreach(ID = unique(dv$curve_ID)) %dopar%{ FitModel(ID)}
ModelOutDF <- do.call(rbind, ModelOutList)
ModelOutDF$trait <- "fecundity rate"
write.csv(ModelOutDF, "../data/bmax_Tpks_AllParams.csv")

#$$$$ Compile dataset containing fits with confidence bounds
#$$ 
AalBetaCIs  <- BetaPlots[[1]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AalBetaFits <- dv_preds %>% filter(species == 'Aedes albopictus') %>% bind_cols(AalBetaCIs)
#$$ 
BtaBetaCIs  <- BetaPlots[[2]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
BtaBetaFits <- dv_preds %>% filter(species == 'Bemisia tabaci')   %>% bind_cols(BtaBetaCIs)
#$$ 
TniBetaCIs  <- BetaPlots[[3]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TniBetaFits <- dv_preds %>% filter(species == 'Tetraneura nigriabdominalis')   %>% bind_cols(TniBetaCIs)
#$$
SpuBetaCIs  <- BetaPlots[[4]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
SpuBetaFits <- dv_preds %>% filter(species == 'Stethorus punctillum')   %>% bind_cols(SpuBetaCIs)
#$$
TemBetaCIs  <- BetaPlots[[5]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TemBetaFits <- dv_preds %>% filter(species == 'Tetranychus mcdanieli')   %>% bind_cols(TemBetaCIs)
#$$
TeuBetaCIs  <- BetaPlots[[6]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TeuBetaFits <- dv_preds %>% filter(species == 'Tetranychus urticae')   %>% bind_cols(TeuBetaCIs)
#$$
CtoBetaCIs  <- BetaPlots[[7]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
CtoBetaFits <- dv_preds %>% filter(species == 'Clavigralla tomentosicollis') %>% bind_cols(CtoBetaCIs)
#$$
PciBetaCIs  <- BetaPlots[[8]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PciBetaFits <- dv_preds %>% filter(species == 'Planococcus citri')   %>% bind_cols(PciBetaCIs)
#$$
MzaBetaCIs  <- BetaPlots[[9]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
MzaBetaFits <- dv_preds %>% filter(species == 'Muscidifurax zaraptor')   %>% bind_cols(MzaBetaCIs)
#$$
AnaBetaCIs  <- BetaPlots[[10]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AnaBetaFits <- dv_preds %>% filter(species == 'Aphis nasturtii')   %>% bind_cols(AnaBetaCIs)
#$$
RmaBetaCIs  <- BetaPlots[[11]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
RmaBetaFits <- dv_preds %>% filter(species == 'Rhopalosiphum maidis')   %>% bind_cols(RmaBetaCIs)
#$$
AaeBetaCIs  <- BetaPlots[[12]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AaeBetaFits <- dv_preds %>% filter(species == 'Aedes aegypti')   %>% bind_cols(AaeBetaCIs)
#$$
AgrBetaCIs  <- BetaPlots[[13]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgrBetaFits <- dv_preds %>% filter(species == 'Anthonomus grandis')   %>% bind_cols(AgrBetaCIs)
#$$
PmaBetaCIs  <- BetaPlots[[14]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PmaBetaFits <- dv_preds %>% filter(species == 'Paracoccus marginatus')   %>% bind_cols(PmaBetaCIs)
#$$
ApiBetaCIs  <- BetaPlots[[15]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
ApiBetaFits <- dv_preds %>% filter(species == 'Acyrthosiphon pisum')   %>% bind_cols(ApiBetaCIs)
#$$
AgoBetaCIs  <- BetaPlots[[16]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgoBetaFits <- dv_preds %>% filter(species == 'Aphis gossypii')   %>% bind_cols(AgoBetaCIs)
#$$
HaxBetaCIs  <- BetaPlots[[17]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
HaxBetaFits <- dv_preds %>% filter(species == 'Harmonia axyridis')   %>% bind_cols(HaxBetaCIs)
#$$
TcaBetaCIs  <- BetaPlots[[18]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TcaBetaFits <- dv_preds %>% filter(species == 'Tribolium castaneum')   %>% bind_cols(TcaBetaCIs)
#$$
AkrBetaCIs  <- BetaPlots[[19]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AkrBetaFits <- dv_preds %>% filter(species == 'Aedes krombeini')   %>% bind_cols(AkrBetaCIs)


BetaPredictions <- bind_rows(AalBetaFits,BtaBetaFits,TniBetaFits,SpuBetaFits, TemBetaFits,
                              TeuBetaFits,CtoBetaFits,PciBetaFits, MzaBetaFits, AnaBetaFits,
                              RmaBetaFits, AaeBetaFits,AgrBetaFits, PmaBetaFits,
                              ApiBetaFits, AgoBetaFits, HaxBetaFits,TcaBetaFits, AkrBetaFits) %>%
  rename(bmax = .fitted, bmaxLwr = conf_lower, bmaxUpr = conf_upper)

write_csv(BetaPredictions, 'BetaPredictions.csv')

# plot all TPCs

#$$ load in raw data
df <- as_tibble(read.csv('TraitData.csv')) 
dv <- df %>% select(interactor1, interactor1temp, standardisedtraitname, standardisedtraitvalue) %>%
      rename(species = interactor1, temp = interactor1temp, bmax = standardisedtraitvalue) %>%
      filter(standardisedtraitname == 'bmax', bmax != 'NA') %>%
      mutate(temp = as.numeric(temp))


BetaPlot <- ggplot(BetaPredictions) +
  geom_line(aes(temp, bmax)) +
  geom_point(aes(temp, bmax), dv, size = 0.75, alpha =0.3) +
  facet_wrap(~species, ncol = 4) +
  scale_y_continuous(limits=c(-0.001,0.2))+
  theme_bw() +
  geom_ribbon(aes(temp, ymin=bmaxLwr, ymax=bmaxUpr), BetaPredictions, fill="#841B2D",alpha=0.3,
              inherit.aes = T)+
  theme(text = element_text(size=10))+
  labs(y=expression(italic(b)[max]), x=expression(plain(paste(" Temperature, ",degree,"C"))))+
  theme(legend.position = 'none'); BetaPlot


save_plot(BetaPlot, file="../results/BetaFits.pdf", 
          base_height=15, base_width = 17.5, base_asp = 1, units="cm")





