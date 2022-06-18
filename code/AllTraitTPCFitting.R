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

######1. Juvenile development rate (1/alpha) ######

df <- as_tibble(read.csv('../data/TraitData.csv')) 
      
dv <- df %>% rename(temp = interactor1temp, species = interactor1, rate = standardisedtraitvalue) %>%
             select(species, temp, standardisedtraitname, rate) %>%
             filter(standardisedtraitname == "1/alpha", rate != 'NA') %>% 
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
                                         species == 'Rhopalosiphum maidis' ~ '19',
                                         species == 'Anopheles gambiae' ~ '20',
                                         species == 'Anoplophora glabripennis' ~ '21',
                                         species == 'Amblyseius womersleyi' ~ '22',
                                         species == 'Trichogramma sp. nr. Lutea' ~ '23',
                                         species == 'Trichogramma bruni' ~ '24',
                                         species == 'Culex annulirostris' ~ '25',
                                         species == 'Macrocentrus iridescens' ~ '26',
                                         species == 'Otiorhynchus sulcatus' ~ '27',
                                         species == 'Drosophila suzukii' ~ '28',
                                         species == 'Gastrolina depressa' ~ '29',
                                         species == 'Laricobius nigrinus' ~ '30',
                                         species == 'Aubeonymus mariaefranciscae' ~ '31',
                                         species == 'Iphiseius degenerans' ~ '32',
                                         species == 'Amblyseius swirskii' ~ '33',
                                         species == 'Macrosiphum euphorbia' ~ '34',
                                         species == 'Myzus persicae' ~ '35',
                                         species == 'Tetranychus evansi' ~ '36',
                                         species == 'Helicoverpa armigera' ~ '37',
                                         species == 'Antestiopsis thunbergii' ~ '38',
                                         species == 'Monochamus leuconotus' ~ '39',
                                         species == 'Kampimodromus aberrans' ~ '40',
                                         species == 'Phenacoccus solenopsis' ~ '41',
                                         species == 'Leptinotarsa decemlineata' ~ '42',
                                         species == 'Halyomorpha halys' ~ '43',
                                         species == 'Muscidifurax raptorellus' ~ '44',
                                         species == 'Thrips hawaiiensis' ~ '45')) %>%
             arrange(curve_ID) %>% mutate(temp = as.numeric(temp))
             
dv %>% distinct(species, curve_ID) %>% print(n=60)

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
                                        upper= upper_lims,
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

##Run everything using apply:

AlphaPlots <- vector(mode="list", length=length(unique(dv$curve_ID)))

ModelOut       <- sapply(unique(dv$curve_ID), FitModel)
ModelOutDFList <- apply(ModelOut, 2, function(x) as.data.frame(do.call(cbind, x)))
ModelOutDF     <- do.call(rbind, ModelOut)
ModelOutDF$trait <- "juvenile development rate"

#run in parallel:
doMC::registerDoMC(cores = 4)
ModelOutList <- foreach(ID = unique(dv$curve_ID)) %dopar%{ FitModel(ID)}
ModelOutDF <- do.call(rbind, ModelOutList)
ModelOutDF$trait <- "juvenile development rate"
write.csv(ModelOutDF, "../data/alpha_Tpks_AllParams.csv")


#compile dataset containing fits with confidence bounds
 
AalAlphaCIs  <- AlphaPlots[[1]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AalAlphaFits <- dv_preds %>% filter(species == 'Aedes albopictus') %>% bind_cols(AalAlphaCIs)
#$$ 
BtaAlphaCIs  <- AlphaPlots[[2]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
BtaAlphaFits <- dv_preds %>% filter(species == 'Bemisia tabaci')   %>% bind_cols(BtaAlphaCIs)
#$$ 
TniAlphaCIs  <- AlphaPlots[[3]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TniAlphaFits <- dv_preds %>% filter(species == 'Tetraneura nigriabdominalis') %>% bind_cols(TniAlphaCIs)
#$$
SpuAlphaCIs  <- AlphaPlots[[4]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
SpuAlphaFits <- dv_preds %>% filter(species == 'Stethorus punctillum')   %>% bind_cols(SpuAlphaCIs)
#$$
TmcAlphaCIs  <- AlphaPlots[[5]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TmcAlphaFits <- dv_preds %>% filter(species == 'Tetranychus mcdanieli')   %>% bind_cols(TmcAlphaCIs)
#$$
TeuAlphaCIs  <- AlphaPlots[[6]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TeuAlphaFits <- dv_preds %>% filter(species == 'Tetranychus urticae')   %>% bind_cols(TeuAlphaCIs)
#$$
CtoAlphaCIs  <- AlphaPlots[[7]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
CtoAlphaFits <- dv_preds %>% filter(species == 'Clavigralla tomentosicollis') %>% bind_cols(CtoAlphaCIs)
#$$
PciAlphaCIs  <- AlphaPlots[[8]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PciAlphaFits <- dv_preds %>% filter(species == 'Planococcus citri')   %>% bind_cols(PciAlphaCIs)
#$$
MzaAlphaCIs  <- AlphaPlots[[9]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
MzaAlphaFits <- dv_preds %>% filter(species == 'Muscidifurax zaraptor')   %>% bind_cols(MzaAlphaCIs)
#$$
AnaAlphaCIs  <- AlphaPlots[[10]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AnaAlphaFits <- dv_preds %>% filter(species == 'Aphis nasturtii')   %>% bind_cols(AnaAlphaCIs)
#$$
RmaAlphaCIs  <- AlphaPlots[[11]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
RmaAlphaFits <- dv_preds %>% filter(species == 'Rhopalosiphum maidis')   %>% bind_cols(RmaAlphaCIs)
#$$
AaeAlphaCIs  <- AlphaPlots[[12]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AaeAlphaFits <- dv_preds %>% filter(species == 'Aedes aegypti')   %>% bind_cols(AaeAlphaCIs)
#$$
AgaAlphaCIs  <- AlphaPlots[[13]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgaAlphaFits <- dv_preds %>% filter(species == 'Anopheles gambiae')   %>% bind_cols(AgaAlphaCIs)
#$$
AglAlphaCIs  <- AlphaPlots[[14]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AglAlphaFits <- dv_preds %>% filter(species == 'Anoplophora glabripennis')   %>% bind_cols(AglAlphaCIs)
#$$
AwoAlphaCIs  <- AlphaPlots[[15]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AwoAlphaFits <- dv_preds %>% filter(species == 'Amblyseius womersleyi')   %>% bind_cols(AwoAlphaCIs)
#$$
TluAlphaCIs  <- AlphaPlots[[16]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TluAlphaFits <- dv_preds %>% filter(species == 'Trichogramma sp. nr. Lutea')   %>% bind_cols(TluAlphaCIs)
#$$
TbrAlphaCIs  <- AlphaPlots[[17]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TbrAlphaFits <- dv_preds %>% filter(species == 'Trichogramma bruni')   %>% bind_cols(TbrAlphaCIs)
#$$
CanAlphaCIs  <- AlphaPlots[[18]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
CanAlphaFits <- dv_preds %>% filter(species == 'Culex annulirostris')   %>% bind_cols(CanAlphaCIs)
#$$
MirAlphaCIs  <- AlphaPlots[[19]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
MirAlphaFits <- dv_preds %>% filter(species == 'Macrocentrus iridescens')   %>% bind_cols(MirAlphaCIs)
#$$
OsuAlphaCIs  <- AlphaPlots[[20]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
OsuAlphaFits <- dv_preds %>% filter(species == 'Otiorhynchus sulcatus')   %>% bind_cols(OsuAlphaCIs)
#$$
DsuAlphaCIs  <- AlphaPlots[[21]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
DsuAlphaFits <- dv_preds %>% filter(species == 'Drosophila suzukii')   %>% bind_cols(DsuAlphaCIs)
#$$
GdeAlphaCIs  <- AlphaPlots[[22]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
GdeAlphaFits <- dv_preds %>% filter(species == 'Gastrolina depressa')   %>% bind_cols(GdeAlphaCIs)
#$$
AgrAlphaCIs  <- AlphaPlots[[23]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgrAlphaFits <- dv_preds %>% filter(species == 'Anthonomus grandis')   %>% bind_cols(AgrAlphaCIs)
#$$
LniAlphaCIs  <- AlphaPlots[[24]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
LniAlphaFits <- dv_preds %>% filter(species == 'Laricobius nigrinus')   %>% bind_cols(LniAlphaCIs)
#$$
AmaAlphaCIs  <- AlphaPlots[[25]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AmaAlphaFits <- dv_preds %>% filter(species == 'Aubeonymus mariaefranciscae')   %>% bind_cols(AmaAlphaCIs)
#$$
IdeAlphaCIs  <- AlphaPlots[[26]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
IdeAlphaFits <- dv_preds %>% filter(species == 'Iphiseius degenerans')   %>% bind_cols(IdeAlphaCIs)
#$$
AswAlphaCIs  <- AlphaPlots[[27]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AswAlphaFits <- dv_preds %>% filter(species == 'Amblyseius swirskii')   %>% bind_cols(AswAlphaCIs)
#$$
MeuAlphaCIs  <- AlphaPlots[[28]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
MeuAlphaFits <- dv_preds %>% filter(species == 'Macrosiphum euphorbia')   %>% bind_cols(MeuAlphaCIs)
#$$
MpeAlphaCIs  <- AlphaPlots[[29]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
MpeAlphaFits <- dv_preds %>% filter(species == 'Myzus persicae')   %>% bind_cols(MpeAlphaCIs)
#$$
TevAlphaCIs  <- AlphaPlots[[30]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TevAlphaFits <- dv_preds %>% filter(species == 'Tetranychus evansi')   %>% bind_cols(TevAlphaCIs)
#$$
HarAlphaCIs  <- AlphaPlots[[31]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
HarAlphaFits <- dv_preds %>% filter(species == 'Helicoverpa armigera')   %>% bind_cols(HarAlphaCIs)
#$$
AthAlphaCIs  <- AlphaPlots[[32]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AthAlphaFits <- dv_preds %>% filter(species == 'Antestiopsis thunbergii')   %>% bind_cols(AthAlphaCIs)
#$$
MleAlphaCIs  <- AlphaPlots[[33]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
MleAlphaFits <- dv_preds %>% filter(species == 'Monochamus leuconotus')   %>% bind_cols(MleAlphaCIs)
#$$
PmaAlphaCIs  <- AlphaPlots[[34]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PmaAlphaFits <- dv_preds %>% filter(species == 'Paracoccus marginatus')   %>% bind_cols(PmaAlphaCIs)
#$$
KabAlphaCIs  <- AlphaPlots[[35]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
KabAlphaFits <- dv_preds %>% filter(species == 'Kampimodromus aberrans')   %>% bind_cols(KabAlphaCIs)
#$$
PsoAlphaCIs  <- AlphaPlots[[36]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PsoAlphaFits <- dv_preds %>% filter(species == 'Phenacoccus solenopsis')   %>% bind_cols(PsoAlphaCIs)
#$$
LdeAlphaCIs  <- AlphaPlots[[37]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
LdeAlphaFits <- dv_preds %>% filter(species == 'Leptinotarsa decemlineata')   %>% bind_cols(LdeAlphaCIs)
#$$
HhaAlphaCIs  <- AlphaPlots[[38]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
HhaAlphaFits <- dv_preds %>% filter(species == 'Halyomorpha halys')   %>% bind_cols(HhaAlphaCIs)
#$$
MraAlphaCIs  <- AlphaPlots[[39]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
MraAlphaFits <- dv_preds %>% filter(species == 'Muscidifurax raptorellus')   %>% bind_cols(MraAlphaCIs)
#$$
ThaAlphaCIs  <- AlphaPlots[[40]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
ThaAlphaFits <- dv_preds %>% filter(species == 'Thrips hawaiiensis')   %>% bind_cols(ThaAlphaCIs)
#$$
ApiAlphaCIs  <- AlphaPlots[[41]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
ApiAlphaFits <- dv_preds %>% filter(species == 'Acyrthosiphon pisum')   %>% bind_cols(ApiAlphaCIs)
#$$
AgoAlphaCIs  <- AlphaPlots[[42]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgoAlphaFits <- dv_preds %>% filter(species == 'Aphis gossypii')   %>% bind_cols(AgoAlphaCIs)
#$$
HaxAlphaCIs  <- AlphaPlots[[43]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
HaxAlphaFits <- dv_preds %>% filter(species == 'Harmonia axyridis')   %>% bind_cols(HaxAlphaCIs)
#$$
TcaAlphaCIs  <- AlphaPlots[[44]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TcaAlphaFits <- dv_preds %>% filter(species == 'Tribolium castaneum')   %>% bind_cols(TcaAlphaCIs)
#$$
AkrAlphaCIs  <- AlphaPlots[[45]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AkrAlphaFits <- dv_preds %>% filter(species == 'Aedes krombeini')   %>% bind_cols(AkrAlphaCIs)
#$$

AlphaPredictions <- bind_rows(AalAlphaFits, 
                              BtaAlphaFits, 
                              TniAlphaFits, 
                              SpuAlphaFits, 
                              TmcAlphaFits,
                              TeuAlphaFits, 
                              CtoAlphaFits, 
                              PciAlphaFits, 
                              MzaAlphaFits, 
                              AnaAlphaFits,
                              RmaAlphaFits, 
                              AaeAlphaFits, 
                              AgaAlphaFits, 
                              AglAlphaFits,
                              AwoAlphaFits,
                              TluAlphaFits,
                              TbrAlphaFits,
                              CanAlphaFits,
                              MirAlphaFits,
                              OsuAlphaFits,
                              DsuAlphaFits,
                              GdeAlphaFits,
                              AgrAlphaFits, 
                              LniAlphaFits,
                              AmaAlphaFits,
                              IdeAlphaFits,
                              AswAlphaFits,
                              MeuAlphaFits,
                              MpeAlphaFits,
                              TevAlphaFits,
                              HarAlphaFits,
                              AthAlphaFits,
                              MleAlphaFits,
                              PmaAlphaFits,
                              KabAlphaFits,
                              PsoAlphaFits,
                              LdeAlphaFits,
                              HhaAlphaFits,
                              MraAlphaFits,
                              ThaAlphaFits,
                              ApiAlphaFits, 
                              AgoAlphaFits, 
                              HaxAlphaFits, 
                              TcaAlphaFits, 
                              AkrAlphaFits) %>%
                    rename(alpha = .fitted, alphaLwr = conf_lower, alphaUpr = conf_upper) %>%
                    mutate(alpha = 1/alpha, alphaLwr = 1/alphaLwr, alphaUpr = 1/alphaUpr)

write_csv(AlphaPredictions, '../data/AlphaPredictions.csv')

# plot all TPCs

#$$ load in raw data
df <- as_tibble(read.csv('../data/TraitData.csv')) 
dv <- df %>% select(interactor1, interactor1temp, standardisedtraitname, standardisedtraitvalue) %>%
  rename(species = interactor1, temp = interactor1temp, alpha = standardisedtraitvalue) %>%
  filter(standardisedtraitname == '1/alpha', alpha != 'NA') %>%
  mutate(temp = as.numeric(temp))

#$$ load in predictions and invert alpha (1/alpha)

AlphaPredictions <- as_tibble(read_csv('AlphaPredictions.csv')) %>% 
  mutate(alpha = 1/alpha, alphaLwr = 1/alphaLwr, alphaUpr = 1/alphaUpr)

AlphaPlot <- ggplot(AlphaPredictions) +
  geom_line(aes(temp, alpha)) +
  geom_point(aes(temp, alpha), dv, size = 0.75, alpha =0.3) +
  facet_wrap(~species, scales = 'free_y', ncol = 4) +
  theme_bw() +
  geom_ribbon(aes(temp, ymin=alphaLwr, ymax=alphaUpr), AlphaPredictions, fill="#e66101",alpha=0.3,
              inherit.aes = T)+
  theme(text = element_text(size=6))+theme(strip.text = element_text(face = "italic"))+
  labs(y=expression(italic(1/alpha)), x=expression(plain(paste(" Temperature, ",degree,"C"))))+
  theme(legend.position = 'none'); AlphaPlot


save_plot(AlphaPlot, file="../results/AlphaFits.pdf", 
          base_height=18,base_width = 14, base_asp = 0.75,units="cm")


#±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±

rm(list=ls())
graphics.off()

# 2. zj (juvenile mortality rate) ######

df <- as_tibble(read.csv('../data/TraitData.csv')) 

dv <- df %>% rename(temp = interactor1temp, species = interactor1, rate = standardisedtraitvalue) %>%
  select(species, temp, standardisedtraitname, rate) %>%
  filter(standardisedtraitname == 'zj', rate != 'NA') %>%
  mutate(curve_ID = case_when(species == 'Thrips hawaiiensis' ~ '1',
                              species == 'Aedes aegypti' ~ '2',
                              species == 'Anthonomus grandis' ~ '3',
                              species == 'Paracoccus marginatus' ~ '4',
                              species == 'Acyrthosiphon pisum' ~ '5',
                              species == 'Aphis gossypii' ~ '6',
                              species == 'Harmonia axyridis' ~ '7',
                              species == 'Tribolium castaneum' ~ '8',
                              species == 'Tetranychus mcdanieli' ~ '9',
                              species == 'Bemisia tabaci' ~ '10',
                              species == 'Tetraneura nigriabdominalis' ~ '11',
                              species == 'Stethorus punctillum' ~ '12',
                              species == 'Tetranychus urticae' ~ '14',
                              species == 'Planococcus citri' ~ '16',
                              species == 'Muscidifurax zaraptor' ~ '17',
                              species == 'Aphis nasturtii' ~ '18',
                              species == 'Rhopalosiphum maidis' ~ '19',
                              species == 'Anopheles gambiae' ~ '20',
                              species == 'Anoplophora glabripennis' ~ '21',
                              species == 'Culex annulirostris' ~ '25',
                              species == 'Laricobius nigrinus' ~ '30',
                              species == 'Aubeonymus mariaefranciscae' ~ '31',
                              species == 'Helicoverpa armigera' ~ '37',
                              species == 'Halyomorpha halys' ~ '43',
                              species == 'Muscidifurax raptorellus' ~ '44',
                              species == 'Aedes albopictus' ~ '45')) %>%
  arrange(curve_ID) %>% filter(curve_ID != 'NA') %>%
  mutate(temp = as.numeric(temp), rate = 1/rate)

dv %>% distinct(species, curve_ID) %>% print(n=50)

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



#$$$$ Compile dataset of fits for all species
ThaZetaJCIs  <- ZetaJPlots[[1]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
ThaZetaJFits <- dv_preds %>% filter(species == 'Thrips hawaiiensis')   %>% bind_cols(ThaZetaJCIs)
#$$
BtaZetaJCIs  <- ZetaJPlots[[2]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
BtaZetaJFits <- dv_preds %>% filter(species == 'Bemisia tabaci')   %>% bind_cols(BtaZetaJCIs)
#$$ 
TniZetaJCIs  <- ZetaJPlots[[3]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TniZetaJFits <- dv_preds %>% filter(species == 'Tetraneura nigriabdominalis') %>% bind_cols(TniZetaJCIs)
#$$ 
SpuZetaJCIs  <- ZetaJPlots[[4]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
SpuZetaJFits <- dv_preds %>% filter(species == 'Stethorus punctillum')   %>% bind_cols(SpuZetaJCIs)
#$$
TeuZetaJCIs  <- ZetaJPlots[[5]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TeuZetaJFits <- dv_preds %>% filter(species == 'Tetranychus urticae')   %>% bind_cols(TeuZetaJCIs)
#$$
PciZetaJCIs  <- ZetaJPlots[[6]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PciZetaJFits <- dv_preds %>% filter(species == 'Planococcus citri')   %>% bind_cols(PciZetaJCIs)
#$$
MzaZetaJCIs  <- ZetaJPlots[[7]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
MzaZetaJFits <- dv_preds %>% filter(species == 'Muscidifurax zaraptor')   %>% bind_cols(MzaZetaJCIs)
#$$
AnaZetaJCIs  <- ZetaJPlots[[8]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AnaZetaJFits <- dv_preds %>% filter(species == 'Aphis nasturtii')   %>% bind_cols(AnaZetaJCIs)
#$$ 
RmaZetaJCIs  <- ZetaJPlots[[9]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
RmaZetaJFits <- dv_preds %>% filter(species == 'Rhopalosiphum maidis') %>% bind_cols(RmaZetaJCIs)
#$$
AaeZetaJCIs  <- ZetaJPlots[[10]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AaeZetaJFits <- dv_preds %>% filter(species == 'Aedes aegypti')   %>% bind_cols(AaeZetaJCIs)
#$$
AgaZetaJCIs  <- ZetaJPlots[[11]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgaZetaJFits <- dv_preds %>% filter(species == 'Anopheles gambiae')   %>% bind_cols(AgaZetaJCIs)
#$$
AglZetaJCIs  <- ZetaJPlots[[12]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AglZetaJFits <- dv_preds %>% filter(species == 'Anoplophora glabripennis')   %>% bind_cols(AglZetaJCIs)
#$$
CanZetaJCIs  <- ZetaJPlots[[13]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
CanZetaJFits <- dv_preds %>% filter(species == 'Culex annulirostris')   %>% bind_cols(CanZetaJCIs)
#$$
AgrZetaJCIs  <- ZetaJPlots[[14]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgrZetaJFits <- dv_preds %>% filter(species == 'Anthonomus grandis')   %>% bind_cols(AgrZetaJCIs)
#$$
LniZetaJCIs  <- ZetaJPlots[[15]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
LniZetaJFits <- dv_preds %>% filter(species == 'Laricobius nigrinus')  %>% bind_cols(LniZetaJCIs)
#$$
AmaZetaJCIs  <- ZetaJPlots[[16]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AmaZetaJFits <- dv_preds %>% filter(species == 'Aubeonymus mariaefranciscae')  %>% bind_cols(AmaZetaJCIs)
#$$
HarZetaJCIs  <- ZetaJPlots[[17]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
HarZetaJFits <- dv_preds %>% filter(species == 'Helicoverpa armigera')  %>% bind_cols(HarZetaJCIs)
#$$
PmaZetaJCIs  <- ZetaJPlots[[18]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PmaZetaJFits <- dv_preds %>% filter(species == 'Paracoccus marginatus')  %>% bind_cols(PmaZetaJCIs)
#$$
HhaZetaJCIs  <- ZetaJPlots[[19]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
HhaZetaJFits <- dv_preds %>% filter(species == 'Halyomorpha halys')  %>% bind_cols(HhaZetaJCIs)
#$$
MraZetaJCIs  <- ZetaJPlots[[20]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
MraZetaJFits <- dv_preds %>% filter(species == 'Muscidifurax raptorellus')  %>% bind_cols(MraZetaJCIs)
#$$ 
AalZetaJCIs  <- ZetaJPlots[[21]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AalZetaJFits <- dv_preds %>% filter(species == 'Aedes albopictus') %>% bind_cols(AalZetaJCIs)
#$$
ApiZetaJCIs  <- ZetaJPlots[[22]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
ApiZetaJFits <- dv_preds %>% filter(species == 'Acyrthosiphon pisum')   %>% bind_cols(ApiZetaJCIs)
#$$
AgoZetaJCIs  <- ZetaJPlots[[23]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgoZetaJFits <- dv_preds %>% filter(species == 'Aphis gossypii')   %>% bind_cols(AgoZetaJCIs)
#$$
HaxZetaJCIs  <- ZetaJPlots[[24]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
HaxZetaJFits <- dv_preds %>% filter(species == 'Harmonia axyridis')   %>% bind_cols(HaxZetaJCIs)
#$$
TcaZetaJCIs  <- ZetaJPlots[[25]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TcaZetaJFits <- dv_preds %>% filter(species == 'Tribolium castaneum')   %>% bind_cols(TcaZetaJCIs)
#$$
TmcZetaJCIs  <- ZetaJPlots[[26]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TmcZetaJFits <- dv_preds %>% filter(species == 'Tetranychus mcdanieli')   %>% bind_cols(TmcZetaJCIs)

ZetaJPredictions <- bind_rows(ThaZetaJFits,
                              BtaZetaJFits, 
                              TniZetaJFits,
                              SpuZetaJFits, 
                              TeuZetaJFits,
                              PciZetaJFits,
                              MzaZetaJFits, 
                              AnaZetaJFits,
                              RmaZetaJFits,
                              AaeZetaJFits,
                              AgaZetaJFits,
                              AglZetaJFits,
                              CanZetaJFits,
                              AgrZetaJFits,
                              LniZetaJFits,
                              AmaZetaJFits,
                              HarZetaJFits,
                              PmaZetaJFits,
                              HhaZetaJFits,
                              MraZetaJFits,
                              AalZetaJFits,
                              ApiZetaJFits,
                              AgoZetaJFits, 
                              HaxZetaJFits,
                              TcaZetaJFits,
                              TmcZetaJFits) %>%
                              rename(zj = .fitted, zjLwr = conf_lower, zjUpr = conf_upper) %>%
                              mutate(zj = 1/zj, zjLwr = 1/zjLwr, zjUpr = 1/zjUpr) 

write_csv(ZetaJPredictions, 'ZetaJPredictions.csv')

#$$$$$$$$$$$$$$$$

# plot fits

# raw data
df <- as_tibble(read.csv('../data/TraitData.csv')) 
dv <- df %>% rename(temp = interactor1temp, species = interactor1, zj = standardisedtraitvalue) %>%
  select(species, temp, standardisedtraitname, zj) %>% 
  filter(standardisedtraitname == 'zj', zj != 'NA') %>%
  filter(species != 'Amblyseius womersleyi' & species != 'Clavigralla tomentosicollis') %>%
  mutate(temp = as.numeric(temp))

# load in predictions
ZetaJPredictions <- as_tibble(read_csv('ZetaJPredictions.csv'))

# truncate fits for plotting 
ZetaJPredictions <- ZetaJPredictions %>%
  mutate_at(vars(c(zj)), 
            ~ifelse(zj > 0.2, 0.2, .)) %>%
  mutate_at(vars(c(zjLwr)), 
            ~ifelse(zjLwr > 0.2, 0.2, .)) %>%
  mutate_at(vars(c(zjUpr)), 
            ~ifelse(zjUpr > 0.2, 0.2, .)) 

ZetaJPredictions <- ZetaJPredictions %>% filter(zj < 0.2)

ZetaJPlot <- ggplot(ZetaJPredictions) +
  geom_line(aes(temp, zj)) +
  geom_point(aes(temp, zj), dv, size = 0.75, alpha =0.3) +
  scale_y_continuous(limits=c(-0.001,0.2))+
  facet_wrap(~species, ncol = 4) +
  theme_bw() +
  geom_ribbon(aes(temp, ymin=zjLwr, ymax=zjUpr), ZetaJPredictions, fill="#1f78b4",alpha=0.3,
              inherit.aes = F)+
  theme(text = element_text(size=6))+theme(strip.text = element_text(face = "italic"))+
  labs(y=expression(italic(z[J])), x=expression(plain(paste(" Temperature, ",degree,"C"))))+
  theme(legend.position = 'none'); ZetaJPlot

save_plot(ZetaJPlot, file="../results/ZetaJFits.pdf", 
          base_height=18,base_width = 14, base_asp = 0.75,units="cm")

#±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±

# 3. z (adult mortality rate)

rm(list=ls())
graphics.off()

df <- as_tibble(read.csv('../data/TraitData.csv')) 

dv <- df %>% rename(temp = interactor1temp, species = interactor1, rate = standardisedtraitvalue) %>%
  select(species, temp, standardisedtraitname, rate) %>%
  filter(standardisedtraitname == 'z', rate != 'NA') %>%
  mutate(curve_ID = case_when(species == 'Culex pipiens' ~ '1',
                              species == 'Plutella xylostella' ~ '2',
                              species == 'Thrips hawaiiensis' ~ '3',
                              species == 'Phenacoccus solenopsis' ~ '4',
                              species == 'Antestiopsis thunbergii' ~ '5',
                              species == 'Culex quinquefasciatus' ~ '6',
                              species == 'Halyomorpha halys' ~ '7',
                              species == 'Monochamus leuconotus' ~ '8',
                              species == 'Anthonomus grandis' ~ '9',
                              species == 'Paracoccus marginatus' ~ '10',
                              species == 'Acyrthosiphon pisum' ~ '11',
                              species == 'Aphis gossypii' ~ '12',
                              species == 'Tribolium castaneum' ~ '14',
                              species == 'Tetranychus mcdanieli' ~ '15',
                              species == 'Bemisia tabaci' ~ '16',
                              species == 'Tetraneura nigriabdominalis' ~ '17',
                              species == 'Stethorus punctillum' ~ '18',
                              species == 'Aedes krombeini' ~ '19',
                              species == 'Muscidifurax zaraptor' ~ '21',
                              species == 'Aphis nasturtii' ~ '22',
                              species == 'Rhopalosiphum maidis' ~ '23',
                              species == 'Anopheles gambiae' ~ '24',
                              species == 'Anoplophora glabripennis' ~ '25',
                              species == 'Helicoverpa armigera' ~ '26',
                              species == 'Aedes albopictus' ~ '27',
                              species == 'Trichogramma bruni' ~ '28',
                              species == 'Trichogramma sp. nr. Lutea' ~ '29',
                              species == 'Aedes aegypti' ~ '30')) %>%
  arrange(curve_ID) %>% filter(curve_ID != 'NA') %>%
  mutate(temp = as.numeric(temp), rate = 1/rate)

dv %>% distinct(species,curve_ID) %>% print(n=50)


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

# Run everything:

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


#$$$$ Compile dataset containing fits for all other species

CpiZetaCIs  <- ZetaPlots[[1]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
CpiZetaFits <- dv_preds %>% filter(species == 'Culex pipiens')   %>% bind_cols(CpiZetaCIs)
#$$
PmaZetaCIs  <- ZetaPlots[[2]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PmaZetaFits <- dv_preds %>% filter(species == 'Paracoccus marginatus')  %>% bind_cols(PmaZetaCIs)
#$$
ApiZetaCIs  <- ZetaPlots[[3]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
ApiZetaFits <- dv_preds %>% filter(species == 'Acyrthosiphon pisum')   %>% bind_cols(ApiZetaCIs)
#$$
AgoZetaCIs  <- ZetaPlots[[4]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgoZetaFits <- dv_preds %>% filter(species == 'Aphis gossypii')   %>% bind_cols(AgoZetaCIs)
#$$
TcaZetaCIs  <- ZetaPlots[[5]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TcaZetaFits <- dv_preds %>% filter(species == 'Tribolium castaneum')   %>% bind_cols(TcaZetaCIs)
#$$
TmcZetaCIs  <- ZetaPlots[[6]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TmcZetaFits <- dv_preds %>% filter(species == 'Tetranychus mcdanieli')   %>% bind_cols(TmcZetaCIs)
#$$
BtaZetaCIs  <- ZetaPlots[[7]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
BtaZetaFits <- dv_preds %>% filter(species == 'Bemisia tabaci')   %>% bind_cols(BtaZetaCIs)
#$$ 
TniZetaCIs  <- ZetaPlots[[8]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TniZetaFits <- dv_preds %>% filter(species == 'Tetraneura nigriabdominalis')   %>% bind_cols(TniZetaCIs)
#$$
SpuZetaCIs  <- ZetaPlots[[9]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
SpuZetaFits <- dv_preds %>% filter(species == 'Stethorus punctillum')   %>% bind_cols(SpuZetaCIs)
#$$
AkrZetaCIs  <- ZetaPlots[[10]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AkrZetaFits <- dv_preds %>% filter(species == 'Aedes krombeini')   %>% bind_cols(AkrZetaCIs)
#$$
PxyZetaCIs  <- ZetaPlots[[11]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PxyZetaFits <- dv_preds %>% filter(species == 'Plutella xylostella')   %>% bind_cols(PxyZetaCIs)
#$$
MzaZetaCIs  <- ZetaPlots[[12]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
MzaZetaFits <- dv_preds %>% filter(species == 'Muscidifurax zaraptor')   %>% bind_cols(MzaZetaCIs)
#$$
AnaZetaCIs  <- ZetaPlots[[13]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AnaZetaFits <- dv_preds %>% filter(species == 'Aphis nasturtii')   %>% bind_cols(AnaZetaCIs)
#$$
RmaZetaCIs  <- ZetaPlots[[14]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
RmaZetaFits <- dv_preds %>% filter(species == 'Rhopalosiphum maidis')   %>% bind_cols(RmaZetaCIs)
#$$
AgaZetaCIs  <- ZetaPlots[[15]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgaZetaFits <- dv_preds %>% filter(species == 'Anopheles gambiae')   %>% bind_cols(AgaZetaCIs)
#$$
AglZetaCIs  <- ZetaPlots[[16]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AglZetaFits <- dv_preds %>% filter(species == 'Anoplophora glabripennis')   %>% bind_cols(AglZetaCIs)
#$$
HarZetaCIs  <- ZetaPlots[[17]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
HarZetaFits <- dv_preds %>% filter(species == 'Helicoverpa armigera')   %>% bind_cols(HarZetaCIs)
#$$ 
AalZetaCIs  <- ZetaPlots[[18]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AalZetaFits <- dv_preds %>% filter(species == 'Aedes albopictus') %>% bind_cols(AalZetaCIs)
#$$
TbrZetaCIs  <- ZetaPlots[[19]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TbrZetaFits <- dv_preds %>% filter(species == 'Trichogramma bruni')   %>% bind_cols(TbrZetaCIs)
#$$
TluZetaCIs  <- ZetaPlots[[20]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TluZetaFits <- dv_preds %>% filter(species == 'Trichogramma sp. nr. Lutea')   %>% bind_cols(TluZetaCIs)
#$$
ThaZetaCIs  <- ZetaPlots[[21]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
ThaZetaFits <- dv_preds %>% filter(species == 'Thrips hawaiiensis')   %>% bind_cols(ThaZetaCIs)
#$$
AaeZetaCIs  <- ZetaPlots[[22]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AaeZetaFits <- dv_preds %>% filter(species == 'Aedes aegypti')   %>% bind_cols(AaeZetaCIs)
#$$
PsoZetaCIs  <- ZetaPlots[[23]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PsoZetaFits <- dv_preds %>% filter(species == 'Phenacoccus solenopsis')   %>% bind_cols(PsoZetaCIs)
#$$
AthZetaCIs  <- ZetaPlots[[24]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AthZetaFits <- dv_preds %>% filter(species == 'Antestiopsis thunbergii')   %>% bind_cols(AthZetaCIs)
#$$
CquZetaCIs  <- ZetaPlots[[25]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
CquZetaFits <- dv_preds %>% filter(species == 'Culex quinquefasciatus')   %>% bind_cols(CquZetaCIs)
#$$
HhaZetaCIs  <- ZetaPlots[[26]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
HhaZetaFits <- dv_preds %>% filter(species == 'Halyomorpha halys')   %>% bind_cols(HhaZetaCIs)
#$$
MleZetaCIs  <- ZetaPlots[[27]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
MleZetaFits <- dv_preds %>% filter(species == 'Monochamus leuconotus')   %>% bind_cols(MleZetaCIs)
#$$
AgrZetaCIs  <- ZetaPlots[[28]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgrZetaFits <- dv_preds %>% filter(species == 'Anthonomus grandis')   %>% bind_cols(AgrZetaCIs)

ZetaPredictions <- bind_rows(CpiZetaFits,
                             PmaZetaFits,
                             ApiZetaFits,
                             AgoZetaFits,
                             TcaZetaFits,
                             TmcZetaFits,
                             BtaZetaFits,
                             TniZetaFits,
                             SpuZetaFits,
                             AkrZetaFits,
                             PxyZetaFits,
                             CpiZetaFits,
                             MzaZetaFits,
                             AnaZetaFits,
                             RmaZetaFits,
                             AgaZetaFits,
                             AglZetaFits,
                             HarZetaFits,
                             AalZetaFits,
                             TbrZetaFits,
                             TluZetaFits,
                             ThaZetaFits,
                             AaeZetaFits,
                             PsoZetaFits,
                             AthZetaFits,
                             CquZetaFits,
                             HhaZetaFits,
                             MleZetaFits,
                             AgrZetaFits) %>%
  rename(z = .fitted, zLwr = conf_lower, zUpr = conf_upper) %>%
  mutate(z = 1/z, zLwr = 1/zLwr, zUpr = 1/zUpr)

write_csv(ZetaPredictions, 'ZetaPredictions.csv')

# plot all TPCs

# raw data
df <- as_tibble(read.csv('../data/TraitData.csv')) 
dv <- df %>% rename(temp = interactor1temp, species = interactor1, z = standardisedtraitvalue) %>%
  select(species, temp, standardisedtraitname, z) %>% 
  filter(standardisedtraitname == 'z', z != 'NA') %>%
  mutate(temp = as.numeric(temp))

#load in predictions
ZetaPredictions <- as_tibble(read_csv('ZetaPredictions.csv'))

# truncate for plotting 
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
  geom_ribbon(aes(temp, ymin=zLwr, ymax=zUpr), ZetaPredictions, fill="#a6cee3",alpha=0.5,
              inherit.aes = T)+
  theme(text = element_text(size=6))+theme(strip.text = element_text(face = "italic"))+
  labs(y=expression(italic(z)), x=expression(plain(paste(" Temperature, ",degree,"C"))))+
  theme(legend.position = 'none'); ZetaPlot

save_plot(ZetaPlot, file="../results/ZetaFits.pdf", 
          base_height=18,base_width = 14, base_asp = 0.75,units="cm")

#±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±

# 4. bmax (daily fecundity rate) ######

rm(list=ls())
graphics.off()

df <- as_tibble(read.csv('../data/TraitData.csv')) 

dv <- df %>% rename(temp = interactor1temp, species = interactor1, rate = standardisedtraitvalue) %>%
  select(species, temp, standardisedtraitname, rate) %>%
  filter(standardisedtraitname == "bmax", rate != 'NA') %>% 
  mutate(curve_ID = case_when(species == 'Aedes aegypti' ~ '2',
                              species == 'Anthonomus grandis' ~ '3',
                              species == 'Paracoccus marginatus' ~ '4',
                              species == 'Acyrthosiphon pisum' ~ '5',
                              species == 'Aphis gossypii' ~ '6',
                              species == 'Aedes krombeini' ~ '9',
                              species == 'Bemisia tabaci' ~ '10',
                              species == 'Tetraneura nigriabdominalis' ~ '11',
                              species == 'Stethorus punctillum' ~ '12',
                              species == 'Tetranychus mcdanieli' ~ '13',
                              species == 'Clavigralla tomentosicollis' ~ '15',
                              species == 'Planococcus citri' ~ '16',
                              species == 'Muscidifurax zaraptor' ~ '17',
                              species == 'Aphis nasturtii' ~ '18',
                              species == 'Rhopalosiphum maidis' ~ '19',
                              species == 'Anoplophora glabripennis' ~ '21',
                              species == 'Amblyseius womersleyi' ~ '22',
                              species == 'Trichogramma sp. nr. Lutea' ~ '23',
                              species == 'Trichogramma bruni' ~ '24',
                              species == 'Drosophila suzukii' ~ '28',
                              species == 'Iphiseius degenerans' ~ '32',
                              species == 'Tetranychus evansi' ~ '36',
                              species == 'Helicoverpa armigera' ~ '37',
                              species == 'Antestiopsis thunbergii' ~ '38',
                              species == 'Monochamus leuconotus' ~ '39',
                              species == 'Kampimodromus aberrans' ~ '40',
                              species == 'Phenacoccus solenopsis' ~ '41',
                              species == 'Halyomorpha halys' ~ '43',
                              species == 'Muscidifurax raptorellus' ~ '44',
                              species == 'Thrips hawaiiensis' ~ '45',
                              species == 'Hylobius transversovittatus' ~ '46',
                              species == 'Callosobruchus maculatus' ~ '47',   
                              species == 'Callosobruchus chinensis' ~ '48',   
                              species == 'Callosobruchus analis' ~ '49',      
                              species == 'Callosobruchus rhodesianus' ~ '50',
                              species == 'Sepedon spinipes' ~ '51',
                              species == 'Plutella xylostella' ~ '52')) %>%
                              arrange(curve_ID) %>% 
                              mutate(temp = as.numeric(temp)) %>%
                              filter(curve_ID != 'NA')
   
#dv %>% distinct(species) %>% print(n=40)


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
BtaBetaCIs  <- BetaPlots[[1]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
BtaBetaFits <- dv_preds %>% filter(species == 'Bemisia tabaci')   %>% bind_cols(BtaBetaCIs)
#$$ 
TniBetaCIs  <- BetaPlots[[2]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TniBetaFits <- dv_preds %>% filter(species == 'Tetraneura nigriabdominalis')   %>% bind_cols(TniBetaCIs)
#$$
SpuBetaCIs  <- BetaPlots[[3]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
SpuBetaFits <- dv_preds %>% filter(species == 'Stethorus punctillum')   %>% bind_cols(SpuBetaCIs)
#$$
TmcBetaCIs  <- BetaPlots[[4]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TmcBetaFits <- dv_preds %>% filter(species == 'Tetranychus mcdanieli')   %>% bind_cols(TmcBetaCIs)
#$$
CtoBetaCIs  <- BetaPlots[[5]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
CtoBetaFits <- dv_preds %>% filter(species == 'Clavigralla tomentosicollis') %>% bind_cols(CtoBetaCIs)
#$$
PciBetaCIs  <- BetaPlots[[6]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PciBetaFits <- dv_preds %>% filter(species == 'Planococcus citri')   %>% bind_cols(PciBetaCIs)
#$$
MzaBetaCIs  <- BetaPlots[[7]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
MzaBetaFits <- dv_preds %>% filter(species == 'Muscidifurax zaraptor')   %>% bind_cols(MzaBetaCIs)
#$$
AnaBetaCIs  <- BetaPlots[[8]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AnaBetaFits <- dv_preds %>% filter(species == 'Aphis nasturtii')   %>% bind_cols(AnaBetaCIs)
#$$
RmaBetaCIs  <- BetaPlots[[9]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
RmaBetaFits <- dv_preds %>% filter(species == 'Rhopalosiphum maidis')   %>% bind_cols(RmaBetaCIs)
#$$
AaeBetaCIs  <- BetaPlots[[10]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AaeBetaFits <- dv_preds %>% filter(species == 'Aedes aegypti')   %>% bind_cols(AaeBetaCIs)
#$$ 
AglBetaCIs  <- BetaPlots[[11]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AglBetaFits <- dv_preds %>% filter(species == 'Anoplophora glabripennis') %>% bind_cols(AglBetaCIs)
#$$ 
AwoBetaCIs  <- BetaPlots[[12]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AwoBetaFits <- dv_preds %>% filter(species == 'Amblyseius womersleyi') %>% bind_cols(AwoBetaCIs)
#$$
TluBetaCIs  <- BetaPlots[[13]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TluBetaFits <- dv_preds %>% filter(species == 'Trichogramma sp. nr. Lutea')   %>% bind_cols(TluBetaCIs)
#$$
TbrBetaCIs  <- BetaPlots[[14]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TbrBetaFits <- dv_preds %>% filter(species == 'Trichogramma bruni')   %>% bind_cols(TbrBetaCIs)
#$$
DsuBetaCIs  <- BetaPlots[[15]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
DsuBetaFits <- dv_preds %>% filter(species == 'Drosophila suzukii')   %>% bind_cols(DsuBetaCIs)
#$$
AgrBetaCIs  <- BetaPlots[[16]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgrBetaFits <- dv_preds %>% filter(species == 'Anthonomus grandis')   %>% bind_cols(AgrBetaCIs)
#$$
IdeBetaCIs  <- BetaPlots[[17]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
IdeBetaFits <- dv_preds %>% filter(species == 'Iphiseius degenerans')   %>% bind_cols(IdeBetaCIs)
#$$
TevBetaCIs  <- BetaPlots[[18]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
TevBetaFits <- dv_preds %>% filter(species == 'Tetranychus evansi')   %>% bind_cols(TevBetaCIs)
#$$
HarBetaCIs  <- BetaPlots[[19]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
HarBetaFits <- dv_preds %>% filter(species == 'Helicoverpa armigera')   %>% bind_cols(HarBetaCIs)
#$$
AthBetaCIs  <- BetaPlots[[20]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AthBetaFits <- dv_preds %>% filter(species == 'Antestiopsis thunbergii')   %>% bind_cols(AthBetaCIs)
#$$
MleBetaCIs  <- BetaPlots[[21]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
MleBetaFits <- dv_preds %>% filter(species == 'Monochamus leuconotus')   %>% bind_cols(MleBetaCIs)
#$$
PmaBetaCIs  <- BetaPlots[[22]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PmaBetaFits <- dv_preds %>% filter(species == 'Paracoccus marginatus')   %>% bind_cols(PmaBetaCIs)
#$$
KabBetaCIs  <- BetaPlots[[23]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
KabBetaFits <- dv_preds %>% filter(species == 'Kampimodromus aberrans')   %>% bind_cols(KabBetaCIs)
#$$
PsoBetaCIs  <- BetaPlots[[24]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PsoBetaFits <- dv_preds %>% filter(species == 'Phenacoccus solenopsis')   %>% bind_cols(PsoBetaCIs)
#$$
HhaBetaCIs  <- BetaPlots[[25]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
HhaBetaFits <- dv_preds %>% filter(species == 'Halyomorpha halys')   %>% bind_cols(HhaBetaCIs)
#$$
MraBetaCIs  <- BetaPlots[[26]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
MraBetaFits <- dv_preds %>% filter(species == 'Muscidifurax raptorellus')   %>% bind_cols(MraBetaCIs)
#$$
ThaBetaCIs  <- BetaPlots[[27]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
ThaBetaFits <- dv_preds %>% filter(species == 'Thrips hawaiiensis')   %>% bind_cols(ThaBetaCIs)
#$$
HtrBetaCIs  <- BetaPlots[[28]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
HtrBetaFits <- dv_preds %>% filter(species == 'Hylobius transversovittatus')   %>% bind_cols(HtrBetaCIs)
#$$
CmaBetaCIs  <- BetaPlots[[29]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
CmaBetaFits <- dv_preds %>% filter(species == 'Callosobruchus maculatus')   %>% bind_cols(CmaBetaCIs)
#$$
CchBetaCIs  <- BetaPlots[[30]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
CchBetaFits <- dv_preds %>% filter(species == 'Callosobruchus chinensis')   %>% bind_cols(CchBetaCIs)
#$$
CanBetaCIs  <- BetaPlots[[31]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
CanBetaFits <- dv_preds %>% filter(species == 'Callosobruchus analis')   %>% bind_cols(CanBetaCIs)
#$$
ApiBetaCIs  <- BetaPlots[[32]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
ApiBetaFits <- dv_preds %>% filter(species == 'Acyrthosiphon pisum')   %>% bind_cols(ApiBetaCIs)
#$$
CrhBetaCIs  <- BetaPlots[[33]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
CrhBetaFits <- dv_preds %>% filter(species == 'Callosobruchus rhodesianus')   %>% bind_cols(CrhBetaCIs)
#$$
SspBetaCIs  <- BetaPlots[[34]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
SspBetaFits <- dv_preds %>% filter(species == 'Sepedon spinipes')   %>% bind_cols(SspBetaCIs)
#$$
PxyBetaCIs  <- BetaPlots[[35]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
PxyBetaFits <- dv_preds %>% filter(species == 'Plutella xylostella')   %>% bind_cols(PxyBetaCIs)
#$$
AgoBetaCIs  <- BetaPlots[[36]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AgoBetaFits <- dv_preds %>% filter(species == 'Aphis gossypii')   %>% bind_cols(AgoBetaCIs)
#$$
AkrBetaCIs  <- BetaPlots[[37]][["plot_env"]][["boot_conf_preds"]] %>% select(conf_lower, conf_upper)
AkrBetaFits <- dv_preds %>% filter(species == 'Aedes krombeini')   %>% bind_cols(AkrBetaCIs)

BetaPredictions <- bind_rows(BtaBetaFits,
                             TniBetaFits,
                             SpuBetaFits,
                             TmcBetaFits,
                             CtoBetaFits,
                             PciBetaFits,
                             MzaBetaFits,
                             AnaBetaFits,
                             RmaBetaFits,
                             AaeBetaFits,
                             AglBetaFits,
                             AwoBetaFits,
                             TluBetaFits,
                             TbrBetaFits,
                             DsuBetaFits,
                             AgrBetaFits,
                             IdeBetaFits,
                             TevBetaFits,
                             HarBetaFits,
                             AthBetaFits,
                             MleBetaFits,
                             PmaBetaFits,
                             KabBetaFits,
                             PsoBetaFits,
                             HhaBetaFits,
                             MraBetaFits,
                             ThaBetaFits,
                             HtrBetaFits,
                             CmaBetaFits,
                             CchBetaFits,
                             CanBetaFits,
                             ApiBetaFits,
                             CrhBetaFits,
                             SspBetaFits,
                             PxyBetaFits,
                             AgoBetaFits,
                             AkrBetaFits) %>%
  rename(bmax = .fitted, bmaxLwr = conf_lower, bmaxUpr = conf_upper)

write_csv(BetaPredictions, 'BetaPredictions.csv')

# plot all TPCs

#$$ load in raw data
df <- as_tibble(read.csv('../data/TraitData.csv')) 
dv <- df %>% select(interactor1, interactor1temp, standardisedtraitname, standardisedtraitvalue) %>%
  rename(species = interactor1, temp = interactor1temp, bmax = standardisedtraitvalue) %>%
  filter(standardisedtraitname == 'bmax', bmax != 'NA') %>%
  mutate(temp = as.numeric(temp))

# load in predictions
BetaPredictions <- as_tibble(read_csv('BetaPredictions.csv'))

BetaPlot <- ggplot(BetaPredictions) +
  geom_line(aes(temp, bmax)) +
  geom_point(aes(temp, bmax), dv, size = 0.75, alpha =0.3) +
  facet_wrap(~species, scales = 'free_y', ncol = 4) +
  theme_bw() +
  geom_ribbon(aes(temp, ymin=bmaxLwr, ymax=bmaxUpr), BetaPredictions, fill="#fdb863",alpha=0.3,
              inherit.aes = T)+
  theme(text = element_text(size=6))+theme(strip.text = element_text(face = "italic"))+
  labs(y=expression(italic(b)[max]), x=expression(plain(paste(" Temperature, ",degree,"C"))))+
  theme(legend.position = 'none'); BetaPlot


save_plot(BetaPlot, file="../results/BetaFits.pdf", 
          base_height=18,base_width = 14, base_asp = 0.75,units="cm")

