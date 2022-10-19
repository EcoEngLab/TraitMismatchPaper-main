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
  filter(standardisedtraitname == 'z', rate != 'NA', 
           species == 'Aedes albopictus' |   species == 'Anopheles gambiae'| 
           species == 'Aedes aegypti'    |   species == 'Harmonia axyridis'|
           species == 'Aedes krombeini'  |   species == 'Tribolium castaneum'|
           species == 'Bemisia tabaci'   |   species == 'Muscidifurax zaraptor'|
           species == 'Aphis nasturtii'  |   species == 'Rhopalosiphum maidis' |
           species == 'Tetranychus mcdanieli'  |  species == 'Tetranychus urticae' |
           species == 'Stethorus punctillum' | species == 'Tetraneura nigriabdominalis' |
           species == 'Anthonomus grandis'  | species == 'Aphis gossypii' | 
           species == 'Acyrthosiphon pisum' | species == 'Planococcus citri' |
           species == 'Paracoccus marginatus') %>%
  mutate(curve_ID = case_when(species == 'Aedes albopictus' ~ '1',
                              species == 'Aedes aegypti' ~ '2',
                              species == 'Anopheles gambiae' ~ '3',
                              species == 'Aedes krombeini' ~ '4',
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
                              species == 'Paracoccus marginatus' ~ '19')) %>% 
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
         y="z", x=expression(plain(paste(" Temperature, ",degree,"C"))))
  
  ggsave(plot,file=paste("../results/TPC/z_",ID,".pdf",sep=""), 
         height=10,width=15,units="cm")
  
  ZetaPlots[[counter]] <<- plot
  
  return(topt)
}

##Run everything:
#run using apply:
ZetaPlots <- vector(mode="list", length=length(unique(dv$curve_ID)))

ModelOut       <- sapply(unique(dv$curve_ID), FitModel)
ModelOutDFList <- apply(ModelOut, 2, function(x) as.data.frame(do.call(cbind, x)))
ModelOutDF     <- do.call(rbind, ModelOut)
ModelOutDF$trait <- "adult mortality rate"

#run in parallel:
doMC::registerDoMC(cores = 4)
ModelOutList <- foreach(ID = unique(dv$curve_ID)) %dopar%{ FitModel(ID)}
ModelOutDF   <- do.call(rbind, ModelOutList)
ModelOutDF$trait <- "adult mortality rate"
write.csv(ModelOutDF, "../data/z_Tpks_AllParams.csv")


##Get plots:
ZetaPlot <- plot_grid(plotlist = ZetaPlots)

save_plot(ZetaPlot, file="../results/ZetaFits.pdf", 
          base_height=28, base_asp = 1,units="cm")

#______________________________________
#  Clavigralla tomentosicollis

###### 2. z (adult mortality rate) ######

Cto <- as_tibble(read.csv('TraitData.csv')) 

Cto <- Cto %>% rename(temp = interactor1temp, species = interactor1, rate = standardisedtraitvalue) %>%
  select(species, temp, standardisedtraitname, rate) %>%
  filter(standardisedtraitname == 'z', rate != 'NA', 
         species == 'Clavigralla tomentosicollis') %>%
  mutate(curve_ID = case_when(species == 'Clavigralla tomentosicollis' ~ '20')) %>% 
  mutate(temp = as.numeric(temp), rate = 1/rate)


# fit chosen model formulation in rTPC
start_vals <- get_start_vals(Cto$temp, Cto$rate, model_name = 'pawar_2018')

b_fits <- nest(Cto, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 1,
                                          start_upper = start_vals + 1,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________
b_preds <- mutate(b_fits, new_data = map(data, ~tibble(temp = seq(min(0), max(45), length.out = 1000))))%>% 
  select(., -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  select(model_name, preds) %>%
  unnest(preds)

# plot panel
ggplot(b_preds) +
  geom_line(aes(temp, 1/.fitted, col = model_name)) +
  geom_point(aes(temp, 1/rate),size=0.2,alpha=0.5, Cto)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________
fit_nlsLM3 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = Cto,
                                start = coef(b_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(Cto)))

# bootstrap using residual resampling
boot3 <- Boot(fit_nlsLM3, method = 'residual')

boot3_preds <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(0), max(45), length.out = 1000))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot3_conf_preds <- group_by(boot3_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

# plot bootstrapped CIs
ggplot() +
  geom_line(aes(temp, 1/.fitted), b_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = 1/conf_lower, ymax = 1/conf_upper), 
              boot3_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, 1/rate), Cto, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/z')

# calculate params with CIs
extra_params3 <- calc_params(fit_nlsLM3) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')%>%
  mutate(trait="z")

ci_extra_params3 <- Boot(fit_nlsLM3, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM3)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params3 <- left_join(ci_extra_params3, extra_params3)

##Get other parameters
params <- broom::tidy(fit_nlsLM3) %>% select(param = term, estimate)
BootOut <- Boot(fit_nlsLM3, method = 'residual')

## Get the param Names that has multiple values:
paramName <- colnames(BootOut[[2]])[which(apply(BootOut[[2]], 2, function(x) length(unique(x))>1))]

params_cis <- BootOut %>%
  confint(.,parm = paramName, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

params_cis <- bind_rows(params_cis) %>%
  left_join(., params) %>% mutate(trait = 'z') %>% filter(param != 'topt')

CtoParams   <- as_tibble(rbind(ci_extra_params3,params_cis)) %>%
  mutate(species = 'Clavigralla tomentosicollis')

write_csv(CtoParams, 'CtoParams.csv')

CtoZeta <- boot3_conf_preds %>% select(conf_lower,conf_upper)
CtoZetaFits <- b_preds %>% select(temp,.fitted) %>% bind_cols(CtoZeta) %>%
  mutate(standardisedtraitname = 'z', species = 'Clavigralla tomentosicollis', curve_ID = '20') %>%
  select(curve_ID, species, standardisedtraitname, temp, .fitted, conf_lower, conf_upper)
  
#$$$$ Compile dataset containing fits for all other species
#$$ 
AalZetaCIs  <- ZetaPlots[[1]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AalZetaFits <- dv_preds %>% filter(species == 'Aedes albopictus') %>% bind_cols(AalZetaCIs)
#$$ 
RmaZetaCIs  <- ZetaPlots[[2]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
RmaZetaFits <- dv_preds %>% filter(species == 'Rhopalosiphum maidis') %>% bind_cols(RmaZetaCIs)
#$$
TemZetaCIs  <- ZetaPlots[[3]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TemZetaFits <- dv_preds %>% filter(species == 'Tetranychus mcdanieli')   %>% bind_cols(TemZetaCIs)
#$$
TeuZetaCIs  <- ZetaPlots[[4]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TeuZetaFits <- dv_preds %>% filter(species == 'Tetranychus urticae')   %>% bind_cols(TeuZetaCIs)
#$$
SpuZetaCIs  <- ZetaPlots[[5]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
SpuZetaFits <- dv_preds %>% filter(species == 'Stethorus punctillum')   %>% bind_cols(SpuZetaCIs)
#$$ 
TniZetaCIs  <- ZetaPlots[[6]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TniZetaFits <- dv_preds %>% filter(species == 'Tetraneura nigriabdominalis')   %>% bind_cols(TniZetaCIs)
#$$
AgrZetaCIs  <- ZetaPlots[[7]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgrZetaFits <- dv_preds %>% filter(species == 'Anthonomus grandis')   %>% bind_cols(AgrZetaCIs)
#$$
AgoZetaCIs  <- ZetaPlots[[8]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgoZetaFits <- dv_preds %>% filter(species == 'Aphis gossypii')   %>% bind_cols(AgoZetaCIs)
#$$
ApiZetaCIs  <- ZetaPlots[[9]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
ApiZetaFits <- dv_preds %>% filter(species == 'Acyrthosiphon pisum')   %>% bind_cols(ApiZetaCIs)
#$$
PciZetaCIs  <- ZetaPlots[[10]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PciZetaFits <- dv_preds %>% filter(species == 'Planococcus citri')   %>% bind_cols(PciZetaCIs)
#$$
PmaZetaCIs  <- ZetaPlots[[11]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PmaZetaFits <- dv_preds %>% filter(species == 'Paracoccus marginatus')  %>% bind_cols(PmaZetaCIs)
#$$
AaeZetaCIs  <- ZetaPlots[[12]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AaeZetaFits <- dv_preds %>% filter(species == 'Aedes aegypti')   %>% bind_cols(AaeZetaCIs)
#$$
AgaZetaCIs  <- ZetaPlots[[13]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgaZetaFits <- dv_preds %>% filter(species == 'Anopheles gambiae')   %>% bind_cols(AgaZetaCIs)
#$$
AkrZetaCIs  <- ZetaPlots[[14]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AkrZetaFits <- dv_preds %>% filter(species == 'Aedes krombeini')   %>% bind_cols(AkrZetaCIs)
#$$
HaxZetaCIs  <- ZetaPlots[[15]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
HaxZetaFits <- dv_preds %>% filter(species == 'Harmonia axyridis')   %>% bind_cols(HaxZetaCIs)
#$$
TcaZetaCIs  <- ZetaPlots[[16]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TcaZetaFits <- dv_preds %>% filter(species == 'Tribolium castaneum')   %>% bind_cols(TcaZetaCIs)
#$$
BtaZetaCIs  <- ZetaPlots[[17]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
BtaZetaFits <- dv_preds %>% filter(species == 'Bemisia tabaci')   %>% bind_cols(BtaZetaCIs)
#$$
MzaZetaCIs  <- ZetaPlots[[18]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
MzaZetaFits <- dv_preds %>% filter(species == 'Muscidifurax zaraptor')   %>% bind_cols(MzaZetaCIs)
#$$
AnaZetaCIs  <- ZetaPlots[[19]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AnaZetaFits <- dv_preds %>% filter(species == 'Aphis nasturtii')   %>% bind_cols(AnaZetaCIs)


ZetaPredictions <- bind_rows(AalZetaFits,RmaZetaFits,TemZetaFits,TeuZetaFits,
                             SpuZetaFits,TniZetaFits,AgrZetaFits,
                             AgoZetaFits,ApiZetaFits,PciZetaFits, PmaZetaFits,
                             AaeZetaFits,AgaZetaFits,AkrZetaFits, HaxZetaFits,
                             TcaZetaFits, BtaZetaFits, MzaZetaFits,AnaZetaFits, CtoZetaFits) %>%
  rename(z = .fitted, zLwr = conf_lower, zUpr = conf_upper) %>%
  mutate(z = 1/z, zLwr = 1/zLwr, zUpr = 1/zUpr)

write_csv(ZetaPredictions, 'ZetaPredictions.csv')

# plot all TPCs

#$$ raw data
df <- as_tibble(read.csv('TraitData.csv')) 
dv <- df %>% rename(temp = interactor1temp, species = interactor1, z = standardisedtraitvalue) %>%
  select(species, temp, standardisedtraitname, z) %>% 
  filter(standardisedtraitname == 'z', z != 'NA') %>%
  mutate(temp = as.numeric(temp))

#$$ truncate for plotting 
ZetaPredictions <- ZetaPredictions %>%
  mutate_at(vars(c(z)), 
            ~ifelse(z > 0.2, 0.2, .)) %>%
  mutate_at(vars(c(zLwr)), 
            ~ifelse(zLwr > 0.2, 0.2, .)) %>%
  mutate_at(vars(c(zUpr)), 
            ~ifelse(zUpr > 0.2, 0.2, .)) 

ZetaPredictions <- ZetaPredictions %>% filter(z < 0.2)
  
ZetaPlot <- ggplot(ZetaPredictions) +
  geom_line(aes(temp, z)) +
  geom_point(aes(temp, z), dv, size = 0.75, alpha =0.3) +
  facet_wrap(~species, ncol = 4)+
  scale_y_continuous(limits=c(-0.001,0.2))+
  theme_bw() +
  geom_ribbon(aes(temp, ymin=zLwr, ymax=zUpr), ZetaPredictions, fill="#004225",alpha=0.3,
              inherit.aes = T)+
  theme(text = element_text(size=10))+
  labs(y=expression(italic(z)), x=expression(plain(paste(" Temperature, ",degree,"C"))))+
  theme(legend.position = 'none'); ZetaPlot

save_plot(ZetaPlot, file="../results/ZetaFits.pdf", 
          base_height=15,base_width = 20, base_asp = 1,units="cm")





