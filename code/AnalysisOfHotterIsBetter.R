
#=====================================================================
#             Analysis of the 'hotter-is-better' pattern             #
#=====================================================================

# Load libraries
require(tidyverse)
require(patchwork) # easy way of creating panels
require(car)
require(forcats)
require(cowplot)

rm(list=ls())
graphics.off()

#±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±

# 1. Development time (a)

bodyMass <- as_tibble(read.csv('../data/sizeMeans.csv')) %>% 
  rename(species = interactor1) %>%
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
  arrange(curve_ID) %>% 
  rename(massspecies = species, masscurve_ID = curve_ID) %>% 
  filter(masscurve_ID != 'NA')

Tc <- as_tibble(read.csv('../data/alpha_Tpks_AllParams.csv')) %>%
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
  
  arrange(curve_ID) 
  
a_pk <- Tc %>% filter(param == 'rmax') %>% 
        rename(a_pk = estimate, a_pkLwr = conf_lower,a_pkUpr = conf_upper) %>% 
        select(a_pk, a_pkLwr, a_pkUpr, species, curve_ID)

T_pk <- Tc %>% filter(param == 'topt') %>% 
  rename(T_pk = estimate, T_pkLwr = conf_lower,T_pkUpr = conf_upper) %>% 
  select(T_pk, T_pkLwr, T_pkUpr)
                                                                         
bodyMass <- bodyMass %>% select(avg, masscurve_ID)

Tc <- bind_cols(a_pk, T_pk, bodyMass)

Tc <- Tc %>% rename(mass = avg) %>% select(-masscurve_ID)

write_csv(Tc, '../data/a_pksT_pksMass.csv')

# linear model (note that alpha_pk is linear in log-log scale)
TcLm <- lm(log(a_pk) ~ T_pk, data = Tc); summary(TcLm)

Tc %>%
  ggplot(aes(x=T_pk,y = log(a_pk)))+
  geom_point()+
  geom_smooth(method = 'lm')

#±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±

a_data <- as_tibble(read.csv('../data/a_pksT_pksMass.csv')) # load in data 
a_data$kT <- 1/(8.617333262145 * 10^-5 * (a_data$T_pk+273.15))
a_data

#plot uncorrected data in log-log scale
a_data %>%
  ggplot(aes(x=log(mass), y = log(a_pk)))+
  geom_point()+
  geom_smooth(method = 'lm')

#plot a_pk vs Tpk
a_data %>%
  ggplot(aes(x = T_pk, y = a_pk)) +
  geom_point()+
  geom_smooth(method = 'lm')

# linear model (note the allometry is linear in log-log scale)
a_model <- lm(log(a_pk) ~ log(mass) + kT, data = a_data)
summary(a_model)
coef(a_model)

a_cf <-  confint(a_model,level = .95)
anova(a_model)

#plot a_pk in 1/kT, correcting for mass
MassCorrectedApkTpk <- a_data %>%
  ggplot(aes(x = T_pk, y = log(a_pk/mass^coef(a_model)[2]))) +
  geom_smooth(method = 'lm', colour = '#636363', size=0.15)+
  scale_y_continuous(expression(plain(paste("ln((", italic(1/alpha[pk])~")/",italic(M^-0.25),")"))))+
  scale_x_continuous(expression(plain(paste(italic("T"[pk])))),
                     limits =c(16.5,39),
                     expand = c(0, 0),
                     breaks=seq(20,35, by=5))+
  geom_linerange(aes(y=log(a_pk/mass^coef(a_model)[2]), xmin=T_pkLwr, xmax=T_pkUpr), 
                 size=0.1,
                 col="#e66101")+
  geom_linerange(aes(x=T_pk, ymin=log(a_pkLwr/mass^a_cf[2]), 
                     ymax=log(a_pkUpr/mass^a_cf[2,2])),
                 size=0.1, col="#e66101")+
  theme_bw()+
  geom_point(size = 1, col="#000000",stroke=0.1, shape = 24, fill="#e66101")
  
#save_plot(MassCorrectedApkTpk, file="../results/MassCorrectedApkTpk.pdf", 
#          base_height=5,base_width = 6, base_asp = 0.75,units="cm")


#±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±

# 2. Juvenile mortality rate (zj)

bodyMass <- as_tibble(read.csv('../data/sizeMeans.csv')) %>% 
  rename(species = interactor1) %>% 
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
  arrange(curve_ID) %>% 
  rename(massspecies = species, masscurve_ID = curve_ID) %>% 
  filter(masscurve_ID != 'NA')


zjPk <- as_tibble(read.csv('zj_Tpks_AllParams.csv')) %>%
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
  arrange(curve_ID) 


zj_pk <- zjPk %>% filter(param == 'rmax') %>% 
  rename(zjpk = estimate, zjpkLwr = conf_lower, zjpkUpr = conf_upper) %>% 
  select(zjpk, zjpkLwr, zjpkUpr, species, curve_ID) %>% mutate(zjpk = as.numeric(zjpk)) %>%
  mutate(zjpk = 1/zjpk, zjpkLwr = 1/zjpkLwr, zjpkUpr = 1/zjpkUpr)

T_pk <- zjPk %>% filter(param == 'topt') %>% 
  rename(T_pk = estimate, T_pkLwr = conf_lower,T_pkUpr = conf_upper) %>% 
  select(T_pk, T_pkLwr, T_pkUpr)

bodyMass <- bodyMass %>% select(avg, masscurve_ID)

zj_data <- bind_cols(zj_pk, T_pk, bodyMass)

zj_data <- zj_data %>% rename(mass = avg) %>% select(-masscurve_ID)

write_csv(zj_data, '../data/zj_pksT_pksMass.csv')

zj_model <- lm(log(zjpk) ~ T_pk, data = zj_data); summary(zj_model)

zj_data %>%
  ggplot(aes(x=T_pk,y = log(zjpk)))+
  geom_point()+
  geom_smooth(method = 'lm')


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

b_cf <-  confint(zj_model,level = .95)
anova(zj_model)

#plot a_pk in 1/kT, correcting for mass
MassCorrectedzjpkTpk <- 
  zj_data %>%
  ggplot(aes(x = T_pk, y = log(zjpk/mass^coef(zj_model)[2])))+
  geom_smooth(method = 'lm', colour = '#636363', size=0.15)+
  geom_linerange(aes(y=log(zjpk/mass^coef(zj_model)[2]), 
                 xmin=T_pkLwr, xmax=T_pkUpr), 
                 size=0.1,
                 col="#1f78b4")+
  geom_linerange(aes(x=T_pk, 
                     ymin=log(zjpkLwr/mass^b_cf[2,2]), 
                     ymax=log(zjpkUpr/mass^b_cf[2])),
                     size=0.1, col="#1f78b4")+
  scale_y_continuous(expression(plain(paste("ln(", italic(z[J][pk])~"/",italic(M^-0.193),")"))))+
  scale_x_continuous(expression(plain(paste(italic("T"[pk])))))+
  theme_bw()+
  geom_point(size = 1, col="#000000",stroke=0.1, shape = 21, fill="#1f78b4")


#save_plot(MassCorrectedzjpkTpk, file="../results/MassCorrectedzjpkTpk.pdf", 
#          base_height=5,base_width = 6, base_asp = 0.75,units="cm")


#==================================================================

# 3. Adult mortality rate (z)

bodyMass <- as_tibble(read.csv('../data/sizeMeans.csv')) %>% 
  rename(species = interactor1) %>% 
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
  arrange(curve_ID) %>% 
  rename(massspecies = species, masscurve_ID = curve_ID) %>% 
  filter(masscurve_ID != 'NA')


zPk <- as_tibble(read.csv('../data/z_Tpks_AllParams.csv')) %>%
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

write_csv(z_data, '../data/z_pksT_pksMass.csv')

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

z_cf <-  confint(z_model,level = .95)
anova(z_model)

#plot a_pk in 1/kT, correcting for mass
MassCorrectedzpkTpk <- 
  z_data %>%
  ggplot(aes(x = T_pk, y = log(zpk/mass^coef(z_model)[2]))) +
  geom_smooth(method = 'lm', colour = '#636363', size=0.15)+
  geom_linerange(aes(y=log(zpk/mass^coef(z_model)[2]), xmin=T_pkLwr, xmax=T_pkUpr), 
                 size=0.1,
                 col="#a6cee3")+
  geom_linerange(aes(x=T_pk, ymin=log(zpkLwr/mass^z_cf[2,2]), 
                     ymax=log(zpkUpr/mass^z_cf[2])),
                 size=0.1, col="#a6cee3")+scale_y_continuous(expression(plain(paste("ln(", italic(z[pk])~"/",italic(M^-0.124),")"))))+
  scale_x_continuous(expression(plain(paste(italic("T"[pk])))))+
  theme_bw()+
  geom_point(size = 1, col="#000000",stroke=0.1, shape = 22, fill="#a6cee3")


#save_plot(MassCorrectedzpkTpk, file="../results/MassCorrectedzpkTpk.pdf", 
#          base_height=5,base_width = 6, base_asp = 0.75,units="cm")


#===============================================================
# 4. Peak fecundity (bmax)

bodyMass <- as_tibble(read.csv('../data/sizeMeans.csv')) %>% 
  rename(species = interactor1) %>% 
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
  arrange(curve_ID) %>% filter(curve_ID != 'NA') %>%
  rename(massspecies = species, masscurve_ID = curve_ID) 


bPk <- as_tibble(read.csv('../data/bmax_Tpks_AllParams.csv')) %>% 
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
  arrange(curve_ID) %>% filter(curve_ID != 'NA')

b_pk <- bPk %>% filter(param == 'rmax') %>% 
  rename(bpk = estimate, bpkLwr = conf_lower, bpkUpr = conf_upper) %>% 
  select(bpk, bpkLwr, bpkUpr, species, curve_ID)

T_pk <- bPk %>% filter(param == 'topt') %>% 
  rename(T_pk = estimate, T_pkLwr = conf_lower,T_pkUpr = conf_upper) %>% 
  select(T_pk, T_pkLwr, T_pkUpr)

bodyMass <- bodyMass %>% select(avg, masscurve_ID)

bpk_data <- bind_cols(b_pk, T_pk, bodyMass)

bpk_data <- bpk_data %>% rename(mass = avg)  %>% select(-masscurve_ID)

write_csv(bpk_data, '../data/bmaxT_pksMass.csv')

#bpk_model <- lm(log(bpk) ~ T_pk, data = bpk_data); summary(bpk_model)

bpk_data %>%
  ggplot(aes(x=T_pk,y = log(bpk)))+
  geom_point()+
  geom_smooth(method = 'lm')

#±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±

bpk_data <- as_tibble(read.csv('../data/bmaxT_pksMass.csv')) # load in data 
bpk_data$kT <- 1/(8.617333262145 * 10^-5 * (bpk_data$T_pk+273.15))
bpk_data

#plot uncorrected data in log-log scale
bpk_data %>%
  ggplot(aes(x=log(mass), y = log(bpk)))+
  geom_point()+
  geom_smooth(method = 'lm')

#plot b_pk vs T
bpk_data %>%
  ggplot(aes(x = T_pk, y = bpk)) +
  geom_point()+
  geom_smooth(method = 'lm')

# linear model (note the allometry is linear in log-log scale)
bpk_model <- lm(log(bpk) ~ log(mass) + kT, data = bpk_data)
summary(bpk_model)
coef(bpk_model)

bpk_cf <-  confint(bpk_model,level = .95)
anova(bpk_model)

#plot a_pk in 1/kT, correcting for mass
MassCorrectedBmaxTpk <- 
  bpk_data %>%
  ggplot(aes(x = T_pk, y = log(bpk/mass^coef(bpk_model)[2]))) +
  geom_smooth(method = 'lm', colour = '#636363', size=0.15)+
  scale_y_continuous(expression(plain(paste("ln(", italic(b[max])~"/",italic(M^0.08),")"))))+
  scale_x_continuous(expression(plain(paste(italic("T"[pk])))))+
  geom_linerange(aes(y=log(bpk/mass^coef(bpk_model)[2]), xmin=T_pkLwr, xmax=T_pkUpr), 
                 size=0.1,
                 col="#fdb863")+
  geom_linerange(aes(x=T_pk, ymin=log(bpkLwr/mass^bpk_cf[2]), ymax=log(bpkUpr/mass^bpk_cf[2,2])),
                 size=0.1, col="#fdb863")+
  theme_bw()+
  geom_point(size = 1, col="#000000",stroke=0.1, shape = 23, fill="#fdb863")


#save_plot(MassCorrectedBmaxTpk, file="../results/MassCorrectedBmaxTpk.pdf", 
#          base_height=5,base_width = 6, base_asp = 0.75,units="cm")


#========= plot hotter-is-better panel 

#p1 <- (MassCorrectedApkTpk+MassCorrectedBmaxTpk)/(MassCorrectedzpkTpk+MassCorrectedzjpkTpk); p1

#save_plot(p1, file="../results/hotterIsbetter.pdf", 
#          base_height=10,base_width = 12, base_asp = 0.75,units="cm")




