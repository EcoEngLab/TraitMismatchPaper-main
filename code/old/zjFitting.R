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

###### 2. z (adult mortality rate) ######

df <- as_tibble(read.csv('TraitData.csv')) 

dv <- df %>% rename(temp = interactor1temp, species = interactor1, rate = standardisedtraitvalue) %>%
  select(species, temp, standardisedtraitname, rate) %>%
  filter(standardisedtraitname == 'zj', rate != 'NA', 
         species == 'Aedes albopictus'   |   species == 'Anopheles gambiae'| 
           species == 'Aedes aegypti'    |   species == 'Harmonia axyridis'|
           species == 'Tribolium castaneum'| species == 'Bemisia tabaci'   |   
           species == 'Muscidifurax zaraptor'|
           species == 'Aphis nasturtii'  |   species == 'Rhopalosiphum maidis' |
           species == 'Tetranychus mcdanieli'  |  species == 'Tetranychus urticae' |
           species == 'Stethorus punctillum' | species == 'Tetraneura nigriabdominalis' |
           species == 'Anthonomus grandis'  | species == 'Aphis gossypii' | 
           species == 'Acyrthosiphon pisum' | species == 'Planococcus citri' |
           species == 'Paracoccus marginatus' | species == 'Clavigralla tomentosicollis') %>%
  mutate(curve_ID = case_when(species == 'Aedes albopictus' ~ '1',
                              species == 'Aedes aegypti' ~ '2',
                              species == 'Anopheles gambiae' ~ '3',
                              species == 'Harmonia axyridis' ~ '5',
                              species == 'Tribolium castaneum' ~ '6',
                              species == 'Bemisia tabaci' ~ '7',
                              species == 'Muscidifurax zaraptor' ~ '8',
                              species == 'Aphis nasturtii' ~ '9',
                              species == 'Rhopalosiphum maidis' ~ '10',
                              species == 'Tetranychus mcdanieli' ~ '11',
                              species == 'Tetranychus urticae' ~ '12',
                              species == 'Stethorus punctillum' ~ '13',
                              species == 'Tetraneura nigriabdominalis' ~ '14',
                              species == 'Anthonomus grandis' ~ '15',
                              species == 'Aphis gossypii' ~ '16',
                              species == 'Acyrthosiphon pisum' ~ '17',
                              species == 'Planococcus citri' ~ '18',
                              species == 'Paracoccus marginatus' ~ '19',
                              species == 'Clavigralla tomentosicollis' ~ '20')) %>% 
  arrange(curve_ID) %>% 
  mutate(temp = as.numeric(temp), rate = 1/rate)


# fit TPC model for each species

start_vals <- get_start_vals(dv$temp, dv$rate, model_name = 'pawar_2018')
low_lims <-   get_lower_lims(dv$temp, dv$rate, model_name = 'pawar_2018')
upper_lims <- get_upper_lims(dv$temp, dv$rate, model_name = 'pawar_2018')


dv_fits <- nest(dv, data = c(temp, rate)) %>%
  mutate(fit = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                        data = .x,
                                        iter = c(3,3,3,3),
                                        start_lower = start_vals - 1,
                                        start_upper = start_vals + 1,
                                        supp_errors = 'Y',
                                        convergence_count = FALSE)))


dv_preds <- mutate(dv_fits, new_data = map(data, ~tibble(temp = seq(min(0), max(45), length.out = 1000)))) %>%
  select(., -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(fit)) %>%
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  select(curve_ID, species, standardisedtraitname, preds) %>%
  unnest(preds)

glimpse(dv_preds)

ggplot(dv_preds) +
  geom_line(aes(temp, 1/.fitted, col = curve_ID)) +
  geom_point(aes(temp, 1/rate), dv) +
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
  plot <- ggplot(data=df, aes(x=temp, y=1/rate))+
    geom_point(size=0.65, alpha=0.4)+
    geom_line(data=plotData, mapping = aes(x=temp, y=1/.fitted)) +
    geom_ribbon(aes(temp, ymin=1/conf_lower, ymax=1/conf_upper), boot_conf_preds, fill="#e66101",alpha=0.3,
                inherit.aes = F)+
    theme_bw()+
    theme(text = element_text(size=10))+
    labs(title=paste(df$GoodName[1], sep=""),
         y="zj", x=expression(plain(paste(" Temperature, ",degree,"C"))))
  
  ggsave(plot,file=paste("../results/TPC/zj_",ID,".pdf",sep=""), 
         height=10,width=15,units="cm")
  
  ZetaJPlots[[counter]] <<- plot
  
  return(topt)
}

## Run everything:
#  run using apply:

ZetaJPlots <- vector(mode="list", length=length(unique(dv$curve_ID)))

ModelOut       <- sapply(unique(dv$curve_ID), FitModel)
ModelOutDFList <- apply(ModelOut, 2, function(x) as.data.frame(do.call(cbind, x)))
ModelOutDF     <- do.call(rbind, ModelOut)
ModelOutDF$trait <- "juvenile mortality rate"

#run in parallel:
doMC::registerDoMC(cores = 4)
ModelOutList <- foreach(ID = unique(dv$curve_ID)) %dopar%{ FitModel(ID)}
ModelOutDF   <- do.call(rbind, ModelOutList)
ModelOutDF$trait <- "juvenile mortality rate"
write.csv(ModelOutDF, "../data/zj_Tpks_AllParams.csv")


#$$$$ Compile dataset containing fits for all other species
#$$ 
AalZetaJCIs  <- ZetaJPlots[[1]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AalZetaJFits <- dv_preds %>% filter(species == 'Aedes albopictus') %>% bind_cols(AalZetaJCIs)
#$$ 
RmaZetaJCIs  <- ZetaJPlots[[2]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
RmaZetaJFits <- dv_preds %>% filter(species == 'Rhopalosiphum maidis') %>% bind_cols(RmaZetaJCIs)
#$$
TemZetaJCIs  <- ZetaJPlots[[3]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TemZetaJFits <- dv_preds %>% filter(species == 'Tetranychus mcdanieli')   %>% bind_cols(TemZetaJCIs)
#$$
TeuZetaJCIs  <- ZetaJPlots[[4]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TeuZetaJFits <- dv_preds %>% filter(species == 'Tetranychus urticae')   %>% bind_cols(TeuZetaJCIs)
#$$
SpuZetaJCIs  <- ZetaJPlots[[5]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
SpuZetaJFits <- dv_preds %>% filter(species == 'Stethorus punctillum')   %>% bind_cols(SpuZetaJCIs)
#$$ 
TniZetaJCIs  <- ZetaJPlots[[6]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TniZetaJFits <- dv_preds %>% filter(species == 'Tetraneura nigriabdominalis')   %>% bind_cols(TniZetaJCIs)
#$$
AgrZetaJCIs  <- ZetaJPlots[[7]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgrZetaJFits <- dv_preds %>% filter(species == 'Anthonomus grandis')   %>% bind_cols(AgrZetaJCIs)
#$$
AgoZetaJCIs  <- ZetaJPlots[[8]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgoZetaJFits <- dv_preds %>% filter(species == 'Aphis gossypii')   %>% bind_cols(AgoZetaJCIs)
#$$
ApiZetaJCIs  <- ZetaJPlots[[9]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
ApiZetaJFits <- dv_preds %>% filter(species == 'Acyrthosiphon pisum')   %>% bind_cols(ApiZetaJCIs)
#$$
PciZetaJCIs  <- ZetaJPlots[[10]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PciZetaJFits <- dv_preds %>% filter(species == 'Planococcus citri')   %>% bind_cols(PciZetaJCIs)
#$$
PmaZetaJCIs  <- ZetaJPlots[[11]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PmaZetaJFits <- dv_preds %>% filter(species == 'Paracoccus marginatus')  %>% bind_cols(PmaZetaJCIs)
#$$
AaeZetaJCIs  <- ZetaJPlots[[12]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AaeZetaJFits <- dv_preds %>% filter(species == 'Aedes aegypti')   %>% bind_cols(AaeZetaJCIs)
#$$
CtoZetaJCIs  <- ZetaJPlots[[13]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
CtoZetaJFits <- dv_preds %>% filter(species == 'Clavigralla tomentosicollis')   %>% bind_cols(CtoZetaJCIs)
#$$
AgaZetaJCIs  <- ZetaJPlots[[14]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgaZetaJFits <- dv_preds %>% filter(species == 'Anopheles gambiae')   %>% bind_cols(AgaZetaJCIs)
#$$
HaxZetaJCIs  <- ZetaJPlots[[15]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
HaxZetaJFits <- dv_preds %>% filter(species == 'Harmonia axyridis')   %>% bind_cols(HaxZetaJCIs)
#$$
TcaZetaJCIs  <- ZetaJPlots[[16]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TcaZetaJFits <- dv_preds %>% filter(species == 'Tribolium castaneum')   %>% bind_cols(TcaZetaJCIs)
#$$
BtaZetaJCIs  <- ZetaJPlots[[17]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
BtaZetaJFits <- dv_preds %>% filter(species == 'Bemisia tabaci')   %>% bind_cols(BtaZetaJCIs)
#$$
MzaZetaJCIs  <- ZetaJPlots[[18]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
MzaZetaJFits <- dv_preds %>% filter(species == 'Muscidifurax zaraptor')   %>% bind_cols(MzaZetaJCIs)
#$$
AnaZetaJCIs  <- ZetaJPlots[[19]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AnaZetaJFits <- dv_preds %>% filter(species == 'Aphis nasturtii')   %>% bind_cols(AnaZetaJCIs)


ZetaJPredictions <- bind_rows(AalZetaJFits,RmaZetaJFits,TemZetaJFits,TeuZetaJFits,
                             SpuZetaJFits,TniZetaJFits,AgrZetaJFits,
                             AgoZetaJFits,ApiZetaJFits,PciZetaJFits, PmaZetaJFits,
                             AaeZetaJFits,CtoZetaJFits, AgaZetaJFits, HaxZetaJFits,
                             TcaZetaJFits, BtaZetaJFits, MzaZetaJFits, AnaZetaJFits) %>%
  rename(zj = .fitted, zjLwr = conf_lower, zjUpr = conf_upper) %>%
  mutate(zj = 1/zj, zjLwr = 1/zjLwr, zjUpr = 1/zjUpr) 

write_csv(ZetaJPredictions, 'ZetaJPredictions.csv')

