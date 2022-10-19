# Load libraries
require(tidyverse)
require(patchwork) # easy way of creating panels
require(car)
require(forcats)
require(cowplot)

setwd("~/Dropbox/TraitTesting/data")

bodyMass <- as_tibble(read.csv('sizeMeans.csv')) %>% 
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
    species == 'Culex annulirostris' ~ '25',
    species == 'Drosophila suzukii' ~ '28',
    species == 'Iphiseius degenerans' ~ '32',
    species == 'Hylobius transversovittatus' ~ '33',
    species == 'Callosobruchus maculatus' ~ '34',
    species == 'Callosobruchus chinensis' ~ '35',
    species == 'Callosobruchus analis' ~ '36',
    species == 'Callosobruchus rhodesianus' ~ '37',
    species == 'Sepedon spinipes' ~ '38')) %>%
  arrange(curve_ID) %>% filter(curve_ID != 'NA') %>%
  rename(massspecies = species, masscurve_ID = curve_ID) 
  

bPk <- as_tibble(read.csv('bmax_Tpks_AllParams.csv')) %>% 
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
                              species == 'Culex annulirostris' ~ '25',
                              species == 'Drosophila suzukii' ~ '28',
                              species == 'Iphiseius degenerans' ~ '32',
                              species == 'Hylobius transversovittatus' ~ '33',
                              species == 'Callosobruchus maculatus' ~ '34',
                              species == 'Callosobruchus chinensis' ~ '35',
                              species == 'Callosobruchus analis' ~ '36',
                              species == 'Callosobruchus rhodesianus' ~ '37',
                              species == 'Sepedon spinipes' ~ '38')) %>%
  arrange(curve_ID)
  


b_pk <- bPk %>% filter(param == 'rmax') %>% 
  rename(bpk = estimate, bpkLwr = conf_lower, bpkUpr = conf_upper) %>% 
  select(bpk, bpkLwr, bpkUpr, species, curve_ID)

T_pk <- bPk %>% filter(param == 'topt') %>% 
  rename(T_pk = estimate, T_pkLwr = conf_lower,T_pkUpr = conf_upper) %>% 
  select(T_pk, T_pkLwr, T_pkUpr)

bodyMass <- bodyMass %>% select(avg, masscurve_ID)

bpk_data <- bind_cols(b_pk, T_pk, bodyMass)

bpk_data <- bpk_data %>% rename(mass = avg)  %>% select(-masscurve_ID)

write_csv(bpk_data, 'bmaxT_pksMass.csv')

bpk_model <- lm(log(bpk) ~ T_pk, data = bpk_data); summary(bpk_model)

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

#plot a_pk vs T
bpk_data %>%
  ggplot(aes(x = T_pk, y = bpk)) +
  geom_point()+
  geom_smooth(method = 'lm')

# linear model (note the allometry is linear in log-log scale)
bpk_model <- lm(log(bpk) ~ log(mass) + kT, data = bpk_data)
summary(bpk_model)
coef(bpk_model)

cf <-  confint(bpk_model,level = .95)
anova(bpk_model)

#plot a_pk in 1/kT, correcting for mass
MassCorrectedBmaxTpk <- 
  bpk_data %>%
  ggplot(aes(x = T_pk, y = log(bpk/mass^coef(bpk_model)[2]))) +
  geom_smooth(method = 'lm', colour = '#636363', size=0.15)+
  scale_y_continuous(expression(plain(paste("ln(", italic(b[max])~"/",italic(M^0.102),")"))))+
  scale_x_continuous(expression(plain(paste(italic("T"[pk])))))+
  geom_linerange(aes(y=log(bpk/mass^coef(bpk_model)[2]), xmin=T_pkLwr, xmax=T_pkUpr), 
                 size=0.1,
                 col="#fdb863")+
  geom_linerange(aes(x=T_pk, ymin=log(bpkLwr/mass^cf[2]), ymax=log(bpkUpr/mass^cf[2,2])),
                 size=0.1, col="#fdb863")+
  theme_bw()+
  geom_point(size = 1, col="#000000",stroke=0.1, shape = 23, fill="#fdb863")


save_plot(MassCorrectedBmaxTpk, file="../results/MassCorrectedBmaxTpk.pdf", 
          base_height=5,base_width = 6, base_asp = 0.75,units="cm")


