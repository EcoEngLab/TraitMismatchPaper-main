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
                              species == 'Tribolium castaneum' ~ '8',
                              species == 'Aedes krombeini' ~ '9',
                              species == 'Bemisia tabaci' ~ '10',
                              species == 'Tetraneura nigriabdominalis' ~ '11',
                              species == 'Stethorus punctillum' ~ '12',
                              species == 'Tetranychus mcdanieli' ~ '13',
                              species == 'Muscidifurax zaraptor' ~ '17',
                              species == 'Aphis nasturtii' ~ '18',
                              species == 'Rhopalosiphum maidis' ~ '19',
                              species == 'Anopheles gambiae' ~ '20',
                              species == 'Trichogramma sp. nr. Lutea' ~ '23',
                              species == 'Trichogramma bruni' ~ '24')) %>%
  arrange(curve_ID) %>% 
  rename(massspecies = species, masscurve_ID = curve_ID) %>% 
  filter(masscurve_ID != 'NA')


zPk <- as_tibble(read.csv('z_Tpks_AllParams.csv')) %>%
  mutate(curve_ID = case_when(species == 'Aedes albopictus' ~ '1',
                              species == 'Aedes aegypti' ~ '2',
                              species == 'Anthonomus grandis' ~ '3',
                              species == 'Paracoccus marginatus' ~ '4',
                              species == 'Acyrthosiphon pisum' ~ '5',
                              species == 'Aphis gossypii' ~ '6',
                              species == 'Tribolium castaneum' ~ '8',
                              species == 'Aedes krombeini' ~ '9',
                              species == 'Bemisia tabaci' ~ '10',
                              species == 'Tetraneura nigriabdominalis' ~ '11',
                              species == 'Stethorus punctillum' ~ '12',
                              species == 'Tetranychus mcdanieli' ~ '13',
                              species == 'Muscidifurax zaraptor' ~ '17',
                              species == 'Aphis nasturtii' ~ '18',
                              species == 'Rhopalosiphum maidis' ~ '19',
                              species == 'Anopheles gambiae' ~ '20',
                              species == 'Trichogramma sp. nr. Lutea' ~ '23',
                              species == 'Trichogramma bruni' ~ '24')) %>%
  arrange(curve_ID) 

z_pk <- zPk %>% filter(param == 'rmax') %>% 
  rename(zpk = estimate, zpkLwr = conf_lower, zpkUpr = conf_upper) %>% 
  select(zpk, zpkLwr, zpkUpr, species, curve_ID) %>% 
  mutate(zpk = 1/zpk, zpkLwr = 1/zpkLwr, zpkUpr = 1/zpkUpr)

T_pk <- zPk %>% filter(param == 'topt') %>% 
  rename(T_pk = estimate, T_pkLwr = conf_lower,T_pkUpr = conf_upper) %>% 
  select(T_pk, T_pkLwr, T_pkUpr)

bodyMass <- bodyMass %>% select(avg, masscurve_ID)

z_data <- bind_cols(z_pk, T_pk, bodyMass)

z_data <- z_data %>% rename(mass = avg) %>% select(-masscurve_ID)

write_csv(z_data, 'z_pksT_pksMass.csv')

z_model <- lm(log(zpk) ~ T_pk, data = z_data); summary(z_model)

z_data %>%
  ggplot(aes(x=T_pk,y = log(zpk)))+
  geom_point()+
  geom_smooth(method = 'lm')


#±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±

z_data <- as_tibble(read.csv('../data/z_pksT_pksMass.csv')) # load in data 
z_data$kT <- 1/(8.617333262145 * 10^-5 * (z_data$T_pk+273.15))
z_data

#plot uncorrected data in log-log scale
z_data %>%
  ggplot(aes(x=log(mass), y = log(zpk)))+
  geom_point()+
  geom_smooth(method = 'lm')

#plot z_pk vs Tpk
z_data %>%
  ggplot(aes(x = T_pk, y = zpk)) +
  geom_point()+
  geom_smooth(method = 'lm')

# linear model (note the allometry is linear in log-log scale)
z_model <- lm(log(zpk) ~ log(mass) + kT, data = z_data)
summary(z_model)
coef(z_model)

cf <-  confint(z_model,level = .95)
anova(z_model)

#plot a_pk in 1/kT, correcting for mass
MassCorrectedzpkTpk <- 
  z_data %>%
  ggplot(aes(x = T_pk, y = log(zpk/mass^coef(z_model)[2]))) +
  geom_smooth(method = 'lm', colour = '#636363', size=0.15)+
  geom_linerange(aes(y=log(zpk/mass^coef(z_model)[2]), xmin=T_pkLwr, xmax=T_pkUpr), 
                 size=0.1,
                 col="#a6cee3")+
  geom_linerange(aes(x=T_pk, ymin=log(zpkLwr/mass^cf[2,2]), ymax=log(zpkUpr/mass^cf[2])),
                 size=0.1, col="#a6cee3")+scale_y_continuous(expression(plain(paste("ln(", italic(z[pk])~"/",italic(M^-0.204),")"))))+
  scale_x_continuous(expression(plain(paste(italic("T"[pk])))))+
  theme_bw()+
  geom_point(size = 1, col="#000000",stroke=0.1, shape = 22, fill="#a6cee3")


save_plot(MassCorrectedzpkTpk, file="../results/MassCorrectedzpkTpk.pdf", 
          base_height=5,base_width = 6, base_asp = 0.75,units="cm")

