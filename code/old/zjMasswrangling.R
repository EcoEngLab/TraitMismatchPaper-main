
# Load libraries
require(tidyverse)
require(patchwork) # easy way of creating panels
require(car)
require(forcats)

setwd("~/TraitTesting/data")

bodyMass <- as_tibble(read.csv('sizeMeans.csv')) %>% 
  rename(species = interactor1) %>% 
  mutate(curve_ID = case_when(species == 'Aedes albopictus' ~ '1',
                              species == 'Aedes aegypti' ~ '2',
                              species == 'Anthonomus grandis' ~ '3',
                              species == 'Paracoccus marginatus' ~ '4',
                              species == 'Acyrthosiphon pisum' ~ '5',
                              species == 'Aphis gossypii' ~ '6',
                              species == 'Harmonia axyridis' ~ '7',
                              species == 'Tribolium castaneum' ~ '8',
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
                              species == 'Amblyseius womersleyi' ~ '22',
                              species == 'Culex annulirostris' ~ '25',
                              species == 'Laricobius nigrinus' ~ '30',
                              species == 'Aubeonymus mariaefranciscae' ~ '31')) %>%
                              arrange(curve_ID) %>% 
                              rename(massspecies = species, masscurve_ID = curve_ID) %>% 
                              filter(masscurve_ID != 'NA')


zjPk <- as_tibble(read.csv('zj_Tpks_AllParams.csv')) %>%
  mutate(curve_ID = case_when(species == 'Aedes albopictus' ~ '1',
                              species == 'Aedes aegypti' ~ '2',
                              species == 'Anthonomus grandis' ~ '3',
                              species == 'Paracoccus marginatus' ~ '4',
                              species == 'Acyrthosiphon pisum' ~ '5',
                              species == 'Aphis gossypii' ~ '6',
                              species == 'Harmonia axyridis' ~ '7',
                              species == 'Tribolium castaneum' ~ '8',
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
                              species == 'Amblyseius womersleyi' ~ '22',
                              species == 'Culex annulirostris' ~ '25',
                              species == 'Laricobius nigrinus' ~ '30',
                              species == 'Aubeonymus mariaefranciscae' ~ '31')) %>%
                              arrange(curve_ID) 
                              
  
zj_pk <- zjPk %>% filter(param == 'rmax') %>% 
  rename(zjpk = estimate, zjpkLwr = conf_lower, zjpkUpr = conf_upper) %>% 
  select(zjpk, zjpkLwr, zjpkUpr, species, curve_ID) %>% 
  mutate(zjpk = 1/zjpk, zjpkLwr = 1/zjpkLwr, zjpkUpr = 1/zjpkUpr)

T_pk <- zjPk %>% filter(param == 'topt') %>% 
  rename(T_pk = estimate, T_pkLwr = conf_lower,T_pkUpr = conf_upper) %>% 
  select(T_pk, T_pkLwr, T_pkUpr)

bodyMass <- bodyMass %>% select(avg, masscurve_ID)

zj_data <- bind_cols(zj_pk, T_pk, bodyMass)

zj_data <- zj_data %>% rename(mass = avg) %>% select(-masscurve_ID)

write_csv(zj_data, 'zj_pksT_pksMass.csv')

zj_model <- lm(log(zjpk) ~ T_pk, data = zj_data); summary(zj_model)

zj_data %>%
  ggplot(aes(x=T_pk,y = log(zjpk)))+
  geom_point()+
  geom_smooth(method = 'lm')


#±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±

# 2. Juvenile mortality rate (zj)

zj_data <- as_tibble(read.csv('../data/zj_pksT_pksMass.csv')) # load in data 
zj_data$kT <- 1/(8.617333262145 * 10^-5 * (zj_data$T_pk+273.15))
zj_data

#plot uncorrected data in log-log scale
zj_data %>%
  ggplot(aes(x=log(mass), y = log(zjpk)))+
  geom_point()+
  geom_smooth(method = 'lm')

#plot zj_pk vs Tpk
zj_data %>%
  ggplot(aes(x = T_pk, y = zjpk)) +
  geom_point()+
  geom_smooth(method = 'lm')

# linear model (note the allometry is linear in log-log scale)
zj_model <- lm(log(zjpk) ~ log(mass) + kT, data = zj_data)
summary(zj_model)
coef(zj_model)

cf <-  confint(zj_model,level = .95)
anova(zj_model)

#plot a_pk in 1/kT, correcting for mass
MassCorrectedzjpkTpk <- 
  zj_data %>%
  ggplot(aes(x = T_pk, y = log(zjpk/mass^coef(zj_model)[2]))) +
  geom_smooth(method = 'lm', colour = '#636363', size=0.15)+
  geom_linerange(aes(y=log(zjpk/mass^coef(zj_model)[2]), xmin=T_pkLwr, xmax=T_pkUpr), 
                 size=0.1,
                 col="#1f78b4")+
  geom_linerange(aes(x=T_pk, ymin=log(zjpkLwr/mass^cf[2,2]), ymax=log(zjpkUpr/mass^cf[2])),
                 size=0.1, col="#1f78b4")+scale_y_continuous(expression(plain(paste("ln(", italic(z[J][pk])~"/",italic(M^-0.201),")"))))+
  scale_x_continuous(expression(plain(paste(italic("T"[pk])))))+
  theme_bw()+
  geom_point(size = 1, col="#000000",stroke=0.1, shape = 21, fill="#1f78b4")


save_plot(MassCorrectedzjpkTpk, file="../results/MassCorrectedzjpkTpk.pdf", 
          base_height=5,base_width = 6, base_asp = 0.75,units="cm")







