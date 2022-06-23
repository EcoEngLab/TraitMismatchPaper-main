
#===========================================================
# Figures for Supplementary Information 

require('tidyverse')
require('data.table')
require('rTPC')
require('car')
require('grDevices')
require('Cairo')
require('ggplot2')
require('ggpubr')
require('patchwork')
require('cowplot')
require('ggtext')

rm(list=ls())
graphics.off()

#======================================

# 1. Variance 

AllTpks <- as_tibble(read_csv('../data/AllTpkParams.csv'))
rm_data <- as_tibble(read_csv('../data/rm_optSizeScaling.csv'))

#===================================================
# mass-corrected r_m opt vs r_m Tpks
rm_data$kT <- 1/(8.617333262145 * 10^-5 * (rm_data$rmTpk+273.15))
rm_data

#plot uncorrected data in log-log scale
rm_data %>%
  ggplot(aes(x=log(mass), y = log(rm_opt)))+
  geom_point()+
  geom_smooth(method = 'lm')

#plot rm_opt vs T
rm_data %>%
  ggplot(aes(x = rmTpk, y = rm_opt)) +
  geom_point()+
  geom_smooth(method = 'lm')

# linear model (note the allometry is linear in log-log scale)
rm_model <- lm(log(rm_opt) ~ log(mass) + kT, data = rm_data)
summary(rm_model)
coef(rm_model)

cf <-  confint(rm_model,level = .95)
anova(rm_model)

#prepare data for plotting
rm_data <- rm_data %>% 
           mutate(rm_massCor = log(rm_opt/mass^coef(rm_model)[2]))
         
#Calculate variance and Sum

OptVar <- AllTpks %>% group_by(species) %>%
  summarise(variance = var(estimate))
OptSum <- AllTpks %>% group_by(species) %>%
  summarise(sum = sum(estimate))

VarPlot <- left_join(rm_data, OptVar, by="species")
VarLM <-   lm(VarPlot$rm_massCor~VarPlot$variance)

summary(VarLM)

Var_rm <- 
  ggplot(VarPlot, aes(x=variance, y=rm_massCor, colour = species, fill = species)) +
  theme_bw(base_size = 12.5) +
  scale_x_continuous(expression(plain(paste("Variance of ", italic(T)[pk]))),
                     limits =c(34.5,101.5),
                     expand = c(0, 0),
                     breaks=seq(40,100, by=10))+
  scale_y_continuous(expression(plain(paste("Log mass-corrected ",italic(r[m])," at its optimal temperature (ln(",italic(r[paste(m,",", opt)])," / ",
                                            italic(M^-0.16),"))"))),
                     limits=c(-3.8,-0.5),
                     expand = c(0.01, 0),
                     breaks=seq(-3,-1, by=1))+
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="#004225")+
  geom_point(aes(shape=species, 
                 fill=species),
             size=2,
             stroke=0.25)+
  scale_shape_manual(values = c(21,21,21,21,
                                22,22,22,22,
                                23,23,23,23,
                                24,24,24,24),
                     name=expression(bold("")),
                     guide = guide_legend(nrow=8,ncol=2,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  scale_fill_manual(values = c("#EEEEEE","#CCCCCC","#999999","#333333",
                               "#EEEEEE","#CCCCCC","#999999","#333333",
                               "#EEEEEE","#CCCCCC","#999999","#333333",
                               "#EEEEEE","#CCCCCC","#999999","#333333"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow=8,ncol=2,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(values = c("#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=8,ncol=2,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  theme_bw()+
  theme(text=element_text(family="Times"))+
  theme(text = element_text(size=5.25),
        legend.position = 'right',
        legend.text = element_text(size = 4, face = 'italic'),
        legend.background = element_rect(colour = "black", size = 0.125), 
        legend.margin=margin(t = 0.05, b = 0.1, r=0.1,l=0.1, unit='cm'),
        legend.key.size = unit(0.3, 'cm'))
  
  
save_plot(Var_rm, file="../results/SI/rmVariance.pdf", 
          base_height=6.5,base_width = 11, base_asp = 0.75,units="cm")


#===================================================
#    2. Activation Energy histogram 

rm(list=ls())
graphics.off()

alpha   <- as_tibble(read.csv('../data/alpha_Tpks_AllParams.csv', header = TRUE))
zj      <- as_tibble(read.csv('../data/zj_Tpks_AllParams.csv', header = TRUE))
z       <- as_tibble(read.csv('../data/z_Tpks_AllParams.csv', header = TRUE))
bmax    <- as_tibble(read.csv('../data/bmax_Tpks_AllParams.csv', header = TRUE))

topt  <- rbind(alpha,zj,z,bmax)

E <-  filter(topt, trait == 'juvenile development rate' | trait == 'fecundity rate', param =="e")  
eh <- filter(topt, trait == 'juvenile mortality rate' | trait == 'adult mortality rate', param =="eh") %>%
  mutate(param=replace(param, param=='eh', 'e')) %>%
  filter(estimate > -8 & estimate < 50)

eh[8,3] <- -7.36; eh[23,3] <- -4.97; eh[41,3] <- -4.55; eh[22,3] <- -4.97
eh[16,4] <- 36.393048; eh[3,4] <- 30.027643

activation_e <- bind_rows(E, eh) %>% 
  mutate(trait = case_when(trait =="adult mortality rate" ~ "Adult Mortality Rate",
                                         trait == "fecundity rate" ~ "Peak Fecundity",
                                         trait == "juvenile development rate" ~ "Development Rate",
                                         trait == "juvenile mortality rate" ~ "Juvenile Mortality Rate"))


eplot <- ggplot(activation_e, aes(x=estimate))+
  theme_bw()+
  scale_x_continuous(expression(plain(paste("Activation Energy (",italic(E),")"))),
                     limits =c(-12,42),
                     expand = c(0, 0),
                     breaks=seq(-10,40, by=10))+
  labs(y="Density")+
  geom_histogram(position="identity", binwidth=2, col='#333333', fill='#004225',alpha = 0.3, size=0.01)+
  geom_density(col='red', size=0.1)+
  facet_grid(.~trait)+
  theme(text = element_text(family = 'Times', size = 10))


ggsave("../results/SI/eHist.pdf",eplot, width = 16, height = 4, units = "cm")

#===============================================================

# Figure SI3*: Tpks + activation energy plot 

# Tpks for all species

alpha   <- as_tibble(read.csv('../data/alpha_Tpks_AllParams.csv', header = TRUE))
zj      <- as_tibble(read.csv('../data/zj_Tpks_AllParams.csv', header = TRUE))
z       <- as_tibble(read.csv('../data/z_Tpks_AllParams.csv', header = TRUE))
bmax    <- as_tibble(read.csv('../data/bmax_Tpks_AllParams.csv', header = TRUE))

topt  <- rbind(alpha,zj,z,bmax)
toptz <- subset(topt, topt$param=="topt")

toptz$species <- as.factor(toptz$species)
toptz$trait   <- as.factor(toptz$trait)

levels(toptz$species)

toptz$trait <- fct_relevel(toptz$trait, "juvenile mortality rate")

toptz <- toptz %>% filter(param!= "rmax")

SpCount <- table(toptz$species)
RmSp <- names(SpCount)[which(SpCount==1)]

toptz <- filter(toptz, !(species %in% RmSp))

# change name to add asterisk to match next fig
toptz$species <- as.character(toptz$species)
toptz$species[which(toptz$species=="Anthonomus grandis")] <- "**Anthonomus grandis***"
toptz$species[which(toptz$species=="Paracoccus marginatus")] <- "**Paracoccus marginatus***"
toptz$species[which(toptz$species=="Aphis nasturtii")] <- "**Aphis nasturtii***"
toptz$species[which(toptz$species=="Tetraneura nigriabdominalis")] <- "**Tetraneura nigriabdominalis***"
toptz$species[which(toptz$species=="Muscidifurax zaraptor")] <- "**Muscidifurax zaraptor***"
toptz$species[which(toptz$species=="Rhopalosiphum maidis")] <- "**Rhopalosiphum maidis***"
toptz$species[which(toptz$species=="Aedes aegypti")] <- "**Aedes aegypti***"
toptz$species[which(toptz$species=="Helicoverpa armigera")] <- "**Helicoverpa armigera***"
toptz$species[which(toptz$species=="Thrips hawaiiensis")] <- "**Thrips hawaiiensis***"
toptz$species[which(toptz$species=="Tetranychus mcdanieli")] <- "**Tetranychus mcdanieli***"
toptz$species[which(toptz$species=="Halyomorpha halys")] <- "**Halyomorpha halys***"
toptz$species[which(toptz$species=="Bemisia tabaci")] <- "**Bemisia tabaci***"
toptz$species[which(toptz$species=="Aphis gossypii")] <- "**Aphis gossypii***"
toptz$species[which(toptz$species=="Anoplophora glabripennis")] <- "**Anoplophora glabripennis***"
toptz$species[which(toptz$species=="Stethorus punctillum")] <- "**Stethorus punctillum***"
toptz$species[which(toptz$species=="Acyrthosiphon pisum")] <- "**Acyrthosiphon pisum***"

toptz <- toptz %>% filter(species != 'Plutella xylostella')

#order by development rate 
alp <- subset(toptz, toptz$trait=="juvenile development rate")
SPorder <- alp$species[order(alp$estimate)]
toptz$species <- factor(toptz$species, levels=SPorder)

# All traits

fig3a <- ggplot(toptz, aes(estimate, species, shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(xmin = conf_lower, xmax = conf_upper),width=0.35,size=0.35) +
  geom_point(size = 2.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("", italic(T)[pk]))),
                     limits =c(6,40),
                     expand = c(0, 0),
                     breaks=seq(12,36, by=4))+
  scale_fill_manual(labels = c(expression(plain(paste("Juvenile Mortality Rate (",italic(z[J]),")"))),
                               expression(plain(paste("Adult Mortality Rate (",italic(z),")"))),
                               expression(plain(paste("Fecundity (",italic(b[max]),")"))),
                               expression(plain(paste("Juvenile Development Rate (",italic(alpha),")")))),
                    values = c("#1f78b4","#a6cee3","#fdb863","#e66101"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow = 1,ncol =4,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(labels = c(expression(plain(paste("Juvenile Mortality Rate (",italic(z[J]),")"))),
                                 expression(plain(paste("Adult Mortality Rate (",italic(z),")"))),
                                 expression(plain(paste("Fecundity (",italic(b[max]),")"))),
                                 expression(plain(paste("Juvenile Development Rate (",italic(alpha),")")))),
                      values = c("#1f78b4","#a6cee3","#fdb863","#e66101"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=1,ncol=4,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  scale_shape_manual(labels = c(expression(plain(paste("Juvenile Mortality Rate (",italic(z[J]),")"))),
                                expression(plain(paste("Adult Mortality Rate (",italic(z),")"))),
                                expression(plain(paste("Fecundity (",italic(b[max]),")"))),
                                expression(plain(paste("Juvenile Development Rate (",italic(alpha),")")))),
                     values = c(21,22,23,24),
                     name=expression(bold("")),
                     guide = guide_legend(nrow = 1,ncol =4,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  theme(legend.position = c(0.35,-0.1),legend.text = element_text(size = 10),
        axis.text.y = element_markdown(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  theme(legend.margin=margin(t = -0.4, unit='cm'))+
  # coord_fixed(ratio = 1.5)
  theme(aspect.ratio=2)+
  ggtitle("A")+
  theme(plot.title=element_text(face="bold", size = 15, vjust = -1))

fig3a

#===============================================
# Fig 3b: Activation energies 

E <-  filter(topt, trait == 'juvenile development rate' | trait == 'fecundity rate', param =="e")  
eh <- filter(topt, trait == 'juvenile mortality rate' | trait == 'adult mortality rate', param =="eh") %>%
  mutate(param=replace(param, param=='eh', 'e')) %>%
  filter(estimate > -8 & estimate < 50)

eh[8,3] <- -7.36; eh[23,3] <- -4.97; eh[41,3] <- -4.55; eh[22,3] <- -4.97
eh[16,4] <- 36.393048; eh[3,4] <- 30.027643

activation_e <- bind_rows(E, eh) 

# change name to add asterisk to match next fig
activation_e$species <- as.character(activation_e$species)
activation_e$species[which(activation_e$species=="Anthonomus grandis")] <- "**Anthonomus grandis***"
activation_e$species[which(activation_e$species=="Paracoccus marginatus")] <- "**Paracoccus marginatus***"
activation_e$species[which(activation_e$species=="Aphis nasturtii")] <- "**Aphis nasturtii***"
activation_e$species[which(activation_e$species=="Tetraneura nigriabdominalis")] <- "**Tetraneura nigriabdominalis***"
activation_e$species[which(activation_e$species=="Muscidifurax zaraptor")] <- "**Muscidifurax zaraptor***"
activation_e$species[which(activation_e$species=="Rhopalosiphum maidis")] <- "**Rhopalosiphum maidis***"
activation_e$species[which(activation_e$species=="Aedes aegypti")] <- "**Aedes aegypti***"
activation_e$species[which(activation_e$species=="Helicoverpa armigera")] <- "**Helicoverpa armigera***"
activation_e$species[which(activation_e$species=="Thrips hawaiiensis")] <- "**Thrips hawaiiensis***"
activation_e$species[which(activation_e$species=="Tetranychus mcdanieli")] <- "**Tetranychus mcdanieli***"
activation_e$species[which(activation_e$species=="Halyomorpha halys")] <- "**Halyomorpha halys***"
activation_e$species[which(activation_e$species=="Bemisia tabaci")] <- "**Bemisia tabaci***"
activation_e$species[which(activation_e$species=="Aphis gossypii")] <- "**Aphis gossypii***"
activation_e$species[which(activation_e$species=="Anoplophora glabripennis")] <- "**Anoplophora glabripennis***"
activation_e$species[which(activation_e$species=="Stethorus punctillum")] <- "**Stethorus punctillum***"
activation_e$species[which(activation_e$species=="Acyrthosiphon pisum")] <- "**Acyrthosiphon pisum***"


activation_e <- mutate(activation_e,
                       trait = case_when(trait =="adult mortality rate" ~ "Adult Mortality Rate",
                                         trait == "fecundity rate" ~ "Peak Fecundity",
                                         trait == "juvenile development rate" ~ "Development Rate",
                                         trait == "juvenile mortality rate" ~ "Juvenile Mortality Rate"))



TraitLevel <- c("Juvenile Mortality Rate", "Adult Mortality Rate","Peak Fecundity","Development Rate" )
activation_e$trait <- as.character(activation_e$trait)
activation_e$trait <- factor(activation_e$trait, levels = TraitLevel)
activation_e$species <- as.character(activation_e$species)
activation_e$species <- factor(activation_e$species, levels=SPorder)
activation_e <- na.omit(activation_e)


fig3b <- ggplot(activation_e, aes(estimate, species, shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(xmin = conf_lower, xmax = conf_upper),width=0.35,size=0.35) +
  geom_point(size = 2.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("Activation Energy (",italic(E),")"))),
                     limits =c(-10,38),
                     expand = c(0, 0),
                     breaks=seq(-8,32, by=8))+
  scale_fill_manual(values = c("#1f78b4","#a6cee3","#fdb863","#e66101"),
                    name=expression(bold("")))+
  scale_colour_manual(values = c("#1f78b4","#a6cee3","#fdb863","#e66101"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=1,ncol=4,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  scale_shape_manual(values = c(21,22,23,24),
                     name=expression(bold("")),
                     guide = guide_legend(nrow = 1,ncol =4,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  theme(legend.position = c(0.35,-0.1),legend.text = element_text(size = 8.5),
        axis.text.y = element_markdown(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  theme(legend.margin=margin(t = -0.4, unit='cm'))+
  theme(axis.text.y=element_blank())+
  theme(plot.margin=margin(l=-4,unit="cm"))+
  theme(aspect.ratio=2)+
  ggtitle("B")+
  theme(plot.title=element_text(face="bold", size = 15, vjust = -1))+
  geom_vline(xintercept = 0.65, lty = 'dashed')

fig3b


plotMain <- fig3a+theme(legend.position="none")+fig3b+theme(legend.position="none")
# plotMain
legend <- get_legend(fig3a+theme(legend.position = "bottom"))


# plotMain <- plot_grid(plotMain, NULL, rel_widths = c(1,0))

fig3 <- plot_grid(plotMain, legend,ncol = 1, rel_heights = c(1, 0.05))

# 
save_plot(fig3, file="../results/SI/Fig3_SI.pdf",
          base_height=15,base_asp=1.5, units="cm")


#=======================================================================

# Relationship between r_m opt and latitude

alphaLat <- as_tibble(read_csv('../data/TraitData.csv')) %>% 
            select(interactor1, standardisedtraitname,latitude) %>%
            filter(standardisedtraitname == '1/alpha' & latitude != 'NA') %>%
            rename(species = interactor1) %>%
            distinct(species, latitude) %>%
            mutate(Latcurve_ID = case_when(species == 'Aedes albopictus' ~ '1',
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
  arrange(Latcurve_ID) %>% rename(latspecies = species)
            
alphaLat %>% print(n=60)



alphaMass <- as_tibble(read_csv('../data/a_pksT_pksMass.csv')) %>% 
  select(species, a_pk, mass) %>% 
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
  arrange(curve_ID) %>% filter(curve_ID != 'NA')
  
  
alphaLat_data <- bind_cols(alphaMass,alphaLat) %>% 
                 select(-latspecies, -Latcurve_ID) %>% 
                 mutate(alphamassCor = a_pk/mass^-0.265)

# linear model (note the allometry is linear in log-log scale)

alphalatmodel <- lm(log(alphamassCor) ~ latitude, data = alphaLat_data)
summary(alphalatmodel)
anova(alphalatmodel)

#plot a_pk vs latitude

alphalat_plot <-
  alphaLat_data %>%
  ggplot(aes(x= latitude, y = log(alphamassCor)))+
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="#004225")+
  scale_y_continuous(expression(plain(paste("ln(", italic(1/alpha[pk])~")/",
                                            italic(M^-0.265),")"))),
                     limits=c(-4.5,-0.8),
                     expand = c(0.01, 0),
                     breaks=seq(-4,-1, by=1))+
  scale_x_continuous(expression(plain(paste("latitude"))))+
  geom_point(size = 1,stroke=0.2, col = '#000000', shape=24, fill ='#e66101')+
  theme_bw()+
  theme(text=element_text(family="Times"))+
  theme(text = element_text(size=6),
        legend.position = 'none',
        legend.text = element_text(size = 4.5, face = 'italic'),
        legend.background = element_rect(colour = "white", size = 0.125), 
        legend.margin=margin(t = 0.01, b = 0.1, r=0.1,l=0.1, unit='cm'),
        legend.key.size = unit(0.3, 'cm'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  

save_plot(alphalat_plot, file="../results/SI/alphaLat.pdf", 
         base_height=3,base_width = 4, base_asp = 0.75,units="cm")


#==========================================

# Trait level TPCs for all species 

# 1. juvenile development rate 

rm(list=ls())
graphics.off()

# load in raw data
df <- as_tibble(read.csv('../data/TraitData.csv')) 
dv <- df %>% select(interactor1, interactor1temp, standardisedtraitname, standardisedtraitvalue) %>%
  rename(species = interactor1, temp = interactor1temp, alpha = standardisedtraitvalue) %>%
  filter(standardisedtraitname == '1/alpha', alpha != 'NA') %>%
  mutate(temp = as.numeric(temp))

#$$ load in predictions and invert alpha (1/alpha)

AlphaPredictions <- as_tibble(read_csv('../data/AlphaPredictions.csv')) %>% 
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


save_plot(AlphaPlot, file="../results/SI/AlphaFits.pdf", 
          base_height=18,base_width = 14, base_asp = 0.75,units="cm")

#=========================================

# 2. juvenile mortality rate

rm(list=ls())
graphics.off()

# raw data
df <- as_tibble(read.csv('../data/TraitData.csv')) 
dv <- df %>% rename(temp = interactor1temp, species = interactor1, zj = standardisedtraitvalue) %>%
  select(species, temp, standardisedtraitname, zj) %>% 
  filter(standardisedtraitname == 'zj', zj != 'NA') %>%
  filter(species != 'Amblyseius womersleyi' & species != 'Clavigralla tomentosicollis') %>%
  mutate(temp = as.numeric(temp))

# load in predictions
ZetaJPredictions <- as_tibble(read_csv('../data/ZetaJPredictions.csv'))

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

save_plot(ZetaJPlot, file="../results/SI/ZetaJFits.pdf", 
          base_height=18,base_width = 14, base_asp = 0.75,units="cm")


#===============================================

# 3. adult mortality rate

rm(list=ls())
graphics.off()

# raw data
df <- as_tibble(read.csv('../data/TraitData.csv')) 
dv <- df %>% rename(temp = interactor1temp, species = interactor1, z = standardisedtraitvalue) %>%
  select(species, temp, standardisedtraitname, z) %>% 
  filter(standardisedtraitname == 'z', z != 'NA') %>%
  mutate(temp = as.numeric(temp))

#load in predictions
ZetaPredictions <- as_tibble(read_csv('../data/ZetaPredictions.csv'))

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

save_plot(ZetaPlot, file="../results/SI/ZetaFits.pdf", 
          base_height=18,base_width = 14, base_asp = 0.75,units="cm")

#======================================================

# 4. fecundity 

rm(list=ls())
graphics.off()

# load in raw data
df <- as_tibble(read.csv('../data/TraitData.csv')) 
dv <- df %>% select(interactor1, interactor1temp, standardisedtraitname, standardisedtraitvalue) %>%
  rename(species = interactor1, temp = interactor1temp, bmax = standardisedtraitvalue) %>%
  filter(standardisedtraitname == 'bmax', bmax != 'NA') %>%
  mutate(temp = as.numeric(temp))


# load in predictions
BetaPredictions <- as_tibble(read_csv('../data/BetaPredictions.csv'))

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


save_plot(BetaPlot, file="../results/SI/BetaFits.pdf", 
          base_height=18,base_width = 14, base_asp = 0.75,units="cm")

#==================================================

# 5.  Population growth rate (r_m) curves 

rm(list=ls())
graphics.off()

df   <- as_tibble(read.csv('../data/r_mCalcs.csv', header = TRUE))

# truncate data for plotting 

rmData <- df %>%
  mutate_at(vars(c(rm_opt)), 
            ~ifelse(rm_opt < -0.001, -0.001, .)) %>%
  mutate_at(vars(c(rm_optLwr)), 
            ~ifelse(rm_optLwr < -0.001, -0.001, .)) %>%
  mutate_at(vars(c(rm_optUpr)), 
            ~ifelse(rm_optUpr < -0.001, -0.001, .)) %>%
  filter(rm_opt > -0.001)

# plot r_m TPCs

rmPlot <- ggplot()+
  geom_line(aes(temp, rm_opt), rmData)+
  facet_wrap(~species, ncol = 4)+
  scale_y_continuous(expression(plain(paste(" Maximal population growth rate ("~italic(r[m])~")"))),
                     limits=c(-0.001,0.6),
                     expand = c(0.01, 0),
                     breaks=seq(0,0.5, by=0.1))+
  theme_bw()+
  geom_ribbon(aes(temp, ymin=rm_optLwr, ymax=rm_optUpr), rmData, fill="#004225",alpha=0.3)+
  theme(text = element_text(size=8))+theme(strip.text = element_text(face = "italic"))+
  labs(x=expression(plain(paste(" Temperature, ",degree,"C"))))+
  theme(legend.position = 'none')

save_plot(rmPlot, file="../results/SI/rmTPCs.pdf", 
          base_height=15, base_width = 17.5, base_asp = 1, units="cm")



#===================================================
#     4. Rearing Temperature 

# alpha  

RearDF <- as_tibble(read_csv('../data/TraitData.csv')) %>% 
          select(standardisedtraitname,interactor1,interactor1growthtemp) %>%
          filter(interactor1growthtemp != 'not stated' & 
                 interactor1growthtemp != 'various' & interactor1growthtemp != 'NA' ) %>%
          rename(trait = standardisedtraitname, species = interactor1, RearTemp = interactor1growthtemp) %>%
          mutate(RearTemp = as.numeric(RearTemp))

RearDF %>% distinct(species) %>% print(n=60)

RearAlpha <- RearDF %>% filter(trait == '1/alpha') %>% 
             group_by(species) %>% 
             summarise(avg = mean(RearTemp)) %>%
             arrange(avg) %>%
             mutate(Rcurve_ID = case_when(species == 'Aedes albopictus' ~ '1',
                            species == 'Anthonomus grandis' ~ '2',
                            species == 'Paracoccus marginatus' ~ '3',
                            species == 'Acyrthosiphon pisum' ~ '4',
                            species == 'Harmonia axyridis' ~ '10',
                            species == 'Tribolium castaneum' ~ '6',
                            species == 'Aedes krombeini' ~ '7',
                            species == 'Tetraneura nigriabdominalis' ~ '12',
                            species == 'Stethorus punctillum' ~ '14',
                            species == 'Tetranychus mcdanieli' ~ '15',
                            species == 'Tetranychus urticae' ~ '16',
                            species == 'Harmonia axyridis' ~ '10',
                            species == 'Planococcus citri' ~ '11',
                            species == 'Muscidifurax zaraptor' ~ '9',
                            species == 'Aphis nasturtii' ~ '8',
                            species == 'Rhopalosiphum maidis' ~ '5',
                            species == 'Anopheles gambiae' ~ '20',
                            species == 'Amblyseius womersleyi' ~ '17',
                            species == 'Macrocentrus iridescens' ~ '18',
                            species == 'Otiorhynchus sulcatus' ~ '19',
                            species == 'Drosophila suzukii' ~ '20',
                            species == 'Gastrolina depressa' ~ '21',
                            species == 'Aubeonymus mariaefranciscae' ~ '22',
                            species == 'Iphiseius degenerans' ~ '23',
                            species == 'Amblyseius swirskii' ~ '24',
                            species == 'Macrosiphum euphorbia' ~ '25',
                            species == 'Myzus persicae' ~ '26',
                            #species == 'Tetranychus evansi' ~ '27',
                            species == 'Helicoverpa armigera' ~ '28',
                            species == 'Antestiopsis thunbergii' ~ '29',
                            species == 'Monochamus leuconotus' ~ '30',
                            species == 'Kampimodromus aberrans' ~ '31',
                            species == 'Phenacoccus solenopsis' ~ '32',
                            species == 'Leptinotarsa decemlineata' ~ '33',
                            species == 'Thrips hawaiiensis' ~ '34')) %>%
                            arrange(Rcurve_ID) %>% rename(Rspecies = species, RTemp = avg)
  
alpha   <- as_tibble(read.csv('../data/alpha_Tpks_AllParams.csv', header = TRUE)) %>%
           filter(param == 'topt') %>% 
           select(species, estimate, conf_lower, conf_upper, trait) %>%
           mutate(curve_ID = case_when(species == 'Aedes albopictus' ~ '1',
                               species == 'Anthonomus grandis' ~ '2',
                               species == 'Paracoccus marginatus' ~ '3',
                               species == 'Acyrthosiphon pisum' ~ '4',
                               species == 'Harmonia axyridis' ~ '10',
                               species == 'Tribolium castaneum' ~ '6',
                               species == 'Aedes krombeini' ~ '7',
                               species == 'Tetraneura nigriabdominalis' ~ '12',
                               species == 'Stethorus punctillum' ~ '14',
                               species == 'Tetranychus mcdanieli' ~ '15',
                               species == 'Tetranychus urticae' ~ '16',
                               species == 'Harmonia axyridis' ~ '10',
                               species == 'Planococcus citri' ~ '11',
                               species == 'Muscidifurax zaraptor' ~ '9',
                               species == 'Aphis nasturtii' ~ '8',
                               species == 'Rhopalosiphum maidis' ~ '5',
                               species == 'Anopheles gambiae' ~ '20',
                               species == 'Amblyseius womersleyi' ~ '17',
                               species == 'Macrocentrus iridescens' ~ '18',
                               species == 'Otiorhynchus sulcatus' ~ '19',
                               species == 'Drosophila suzukii' ~ '20',
                               species == 'Gastrolina depressa' ~ '21',
                               species == 'Aubeonymus mariaefranciscae' ~ '22',
                               species == 'Iphiseius degenerans' ~ '23',
                               species == 'Amblyseius swirskii' ~ '24',
                               species == 'Macrosiphum euphorbia' ~ '25',
                               species == 'Myzus persicae' ~ '26',
                               #species == 'Tetranychus evansi' ~ '27',
                               species == 'Helicoverpa armigera' ~ '28',
                               species == 'Antestiopsis thunbergii' ~ '29',
                               species == 'Monochamus leuconotus' ~ '30',
                               species == 'Kampimodromus aberrans' ~ '31',
                               species == 'Phenacoccus solenopsis' ~ '32',
                               species == 'Leptinotarsa decemlineata' ~ '33',
                               species == 'Thrips hawaiiensis' ~ '34')) %>%
              arrange(curve_ID) %>% filter(curve_ID != 'NA')

  
RearAlpha <- bind_cols(alpha,RearAlpha)


RAlp_model <- lm(estimate ~ RTemp, RearAlpha); summary(RAlp_model)

Rear_alp <- 
  
  ggplot(RearAlpha, aes(x=RTemp, y=estimate)) +
  geom_smooth(method = 'lm', colour = '#636363', size=0.15, fill = NA)+
  geom_linerange(aes(ymin =conf_lower , ymax = conf_upper),size= 0.1,col="#000000") +
  geom_point(size = 2.5,stroke=0.2, col = '#000000', shape=24, fill ='#e66101')+
  theme_bw()+
  scale_x_continuous(expression(plain(paste("Rearing Temperature (",degree,"C)"))),
                    limits =c(20,30.5),
                    expand = c(0, 0),
                    breaks=seq(21,29, by=1))+
  scale_y_continuous(expression(paste(italic(T[pk]), " of Juvenile Development Time (", italic(alpha),")")),
                     limits =c(21,38),
                     expand = c(0, 0),
                     breaks=seq(22,36, by=2))+
  theme(legend.position = 'none',legend.text = element_text(size = 10))+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = 20.5, y = 37,label = "A"), 
            parse = TRUE, size = 4, colour = "black")


#===================================================
# Zeta J

zetaJ   <- as_tibble(read.csv('../data/zj_Tpks_AllParams.csv', header = TRUE)) %>%
  filter(param == 'topt') %>% 
  select(species, estimate, conf_lower, conf_upper, trait) %>%
  mutate(curve_ID = case_when(species == 'Aedes albopictus' ~ '1',
                              species == 'Anthonomus grandis' ~ '2',
                              species == 'Paracoccus marginatus' ~ '3',
                              species == 'Acyrthosiphon pisum' ~ '4',
                              species == 'Harmonia axyridis' ~ '10',
                              species == 'Tribolium castaneum' ~ '6',
                              #species == 'Aedes krombeini' ~ '7',
                              species == 'Tetraneura nigriabdominalis' ~ '12',
                              species == 'Stethorus punctillum' ~ '14',
                              species == 'Tetranychus mcdanieli' ~ '15',
                              species == 'Tetranychus urticae' ~ '16',
                              species == 'Planococcus citri' ~ '11',
                              species == 'Muscidifurax zaraptor' ~ '9',
                              species == 'Aphis nasturtii' ~ '8',
                              species == 'Rhopalosiphum maidis' ~ '5',
                              species == 'Anopheles gambiae' ~ '20',
                              #species == 'Amblyseius womersleyi' ~ '17',
                              #species == 'Macrocentrus iridescens' ~ '18',
                              #species == 'Otiorhynchus sulcatus' ~ '19',
                              #species == 'Drosophila suzukii' ~ '20',
                              #species == 'Gastrolina depressa' ~ '21',
                              species == 'Aubeonymus mariaefranciscae' ~ '22',
                              #species == 'Iphiseius degenerans' ~ '23',
                              #species == 'Amblyseius swirskii' ~ '24',
                              #species == 'Macrosiphum euphorbia' ~ '25',
                              #species == 'Myzus persicae' ~ '26',
                              #species == 'Tetranychus evansi' ~ '27',
                              species == 'Helicoverpa armigera' ~ '28',
                              species == 'Antestiopsis thunbergii' ~ '29',
                              #species == 'Monochamus leuconotus' ~ '30',
                              #species == 'Kampimodromus aberrans' ~ '31',
                              #species == 'Phenacoccus solenopsis' ~ '32',
                              #species == 'Leptinotarsa decemlineata' ~ '33',
                              species == 'Thrips hawaiiensis' ~ '34')) %>%
  arrange(curve_ID) %>% filter(curve_ID != 'NA')

zetaJ %>% distinct(species,curve_ID) %>% print(n=30)

RearDF <- as_tibble(read_csv('../data/TraitData.csv')) %>% 
  select(standardisedtraitname,interactor1,interactor1growthtemp) %>%
  filter(interactor1growthtemp != 'not stated' & 
           interactor1growthtemp != 'various' & interactor1growthtemp != 'NA' ) %>%
  rename(trait = standardisedtraitname, species = interactor1, RearTemp = interactor1growthtemp) %>%
  mutate(RearTemp = as.numeric(RearTemp))

RearZetaJ <- RearDF %>% filter(trait == 'zj') %>% 
  group_by(species) %>% 
  summarise(avg = mean(RearTemp)) %>%
  arrange(avg) %>%
  mutate(Rcurve_ID = case_when(species == 'Aedes albopictus' ~ '1',
                               species == 'Anthonomus grandis' ~ '2',
                               species == 'Paracoccus marginatus' ~ '3',
                               species == 'Acyrthosiphon pisum' ~ '4',
                               species == 'Harmonia axyridis' ~ '10',
                               species == 'Tribolium castaneum' ~ '6',
                               #species == 'Aedes krombeini' ~ '7',
                               species == 'Tetraneura nigriabdominalis' ~ '12',
                               species == 'Stethorus punctillum' ~ '14',
                               species == 'Tetranychus mcdanieli' ~ '15',
                               species == 'Tetranychus urticae' ~ '16',
                               species == 'Planococcus citri' ~ '11',
                               species == 'Muscidifurax zaraptor' ~ '9',
                               species == 'Aphis nasturtii' ~ '8',
                               species == 'Rhopalosiphum maidis' ~ '5',
                               species == 'Anopheles gambiae' ~ '20',
                               #species == 'Amblyseius womersleyi' ~ '17',
                               #species == 'Macrocentrus iridescens' ~ '18',
                               #species == 'Otiorhynchus sulcatus' ~ '19',
                               #species == 'Drosophila suzukii' ~ '20',
                               #species == 'Gastrolina depressa' ~ '21',
                               species == 'Aubeonymus mariaefranciscae' ~ '22',
                               #species == 'Iphiseius degenerans' ~ '23',
                               #species == 'Amblyseius swirskii' ~ '24',
                               #species == 'Macrosiphum euphorbia' ~ '25',
                               #species == 'Myzus persicae' ~ '26',
                               #species == 'Tetranychus evansi' ~ '27',
                               species == 'Helicoverpa armigera' ~ '28',
                               species == 'Antestiopsis thunbergii' ~ '29',
                               #species == 'Monochamus leuconotus' ~ '30',
                               #species == 'Kampimodromus aberrans' ~ '31',
                               #species == 'Phenacoccus solenopsis' ~ '32',
                               #species == 'Leptinotarsa decemlineata' ~ '33',
                               species == 'Thrips hawaiiensis' ~ '34')) %>%
  arrange(Rcurve_ID) %>% rename(Rspecies = species, RTemp = avg) %>%
  filter(Rcurve_ID != 'NA')

RearZetaJ %>% distinct(Rspecies,Rcurve_ID) %>% print(n=60)

RearZetaJ <- bind_cols(zetaJ,RearZetaJ)

RzetaJ_model <- lm(estimate ~ RTemp, RearZetaJ); summary(RzetaJ_model)

Rear_zj <- 
  
  ggplot(RearZetaJ, aes(x=RTemp, y=estimate)) +
  geom_smooth(method = 'lm', colour = '#636363', size=0.15, fill = NA)+
  geom_linerange(aes(ymin =conf_lower , ymax = conf_upper),size= 0.1,col="#000000") +
  geom_point(size = 2.5,stroke=0.2, col = '#000000', shape = 21, fill="#1f78b4")+
  theme_bw()+
  scale_x_continuous(expression(plain(paste("Rearing Temperature (",degree,"C)"))),
                     limits =c(21.5,30.75),
                     expand = c(0, 0),
                     breaks=seq(21,32, by=1))+
  scale_y_continuous(expression(paste(italic(T[pk]), " of Juvenile Mortality Rate (", italic(z[J]),")")),
                     limits =c(7.5,35),
                     expand = c(0, 0),
                     breaks=seq(10,30, by=5))+
  theme(legend.position = 'none',legend.text = element_text(size = 10))+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = 22, y = 33,label = "D"), 
            parse = TRUE, size = 4, colour = "black")


#===================================================
# Zeta 

RearDF <- as_tibble(read_csv('../data/TraitData.csv')) %>% 
  select(standardisedtraitname,interactor1,interactor1growthtemp) %>%
  filter(interactor1growthtemp != 'not stated' & 
           interactor1growthtemp != 'various' & interactor1growthtemp != 'NA' ) %>%
  rename(trait = standardisedtraitname, species = interactor1, RearTemp = interactor1growthtemp) %>%
  mutate(RearTemp = as.numeric(RearTemp))

RearZeta <- RearDF %>% filter(trait == 'z') %>% 
  group_by(species) %>% 
  summarise(avg = mean(RearTemp)) %>%
  arrange(avg) %>%
  mutate(Rcurve_ID = case_when(species == 'Plutella xylostella' ~ '1',
                               species == 'Aedes albopictus' ~ '32',
                               species == 'Paracoccus marginatus' ~ '3',
                               species == 'Acyrthosiphon pisum' ~ '4',
                               species == 'Tribolium castaneum' ~ '6',
                               species == 'Aedes krombeini' ~ '7',
                               species == 'Tetraneura nigriabdominalis' ~ '12',
                               species == 'Stethorus punctillum' ~ '14',
                               species == 'Tetranychus mcdanieli' ~ '15',
                               species == 'Muscidifurax zaraptor' ~ '9',
                               species == 'Aphis nasturtii' ~ '2',
                               species == 'Rhopalosiphum maidis' ~ '5',
                               species == 'Anopheles gambiae' ~ '20',
                               species == 'Anthonomus grandis' ~ '27',
                               species == 'Helicoverpa armigera' ~ '28',
                               species == 'Antestiopsis thunbergii' ~ '29',
                               species == 'Monochamus leuconotus' ~ '30',
                               species == 'Culex pipiens' ~ '31',
                               species == 'Phenacoccus solenopsis' ~ '32',
                               species == 'Culex quinquefasciatus' ~ '33',
                               species == 'Thrips hawaiiensis' ~ '34')) %>%
  arrange(Rcurve_ID) %>% rename(Rspecies = species, RTemp = avg) %>%
  filter(Rcurve_ID != 'NA')

RearZeta %>% distinct(Rspecies) %>% print(n=60)


zeta  <- as_tibble(read.csv('../data/z_Tpks_AllParams.csv', header = TRUE)) %>%
  filter(param == 'topt') %>% 
  select(species, estimate, conf_lower, conf_upper, trait) %>%
  mutate(curve_ID = case_when(species == 'Plutella xylostella' ~ '1',
                               species == 'Aedes albopictus' ~ '32',
                               species == 'Paracoccus marginatus' ~ '3',
                               species == 'Acyrthosiphon pisum' ~ '4',
                               species == 'Tribolium castaneum' ~ '6',
                               species == 'Aedes krombeini' ~ '7',
                               species == 'Tetraneura nigriabdominalis' ~ '12',
                               species == 'Stethorus punctillum' ~ '14',
                               species == 'Tetranychus mcdanieli' ~ '15',
                               species == 'Muscidifurax zaraptor' ~ '9',
                               species == 'Aphis nasturtii' ~ '2',
                               species == 'Rhopalosiphum maidis' ~ '5',
                               species == 'Anopheles gambiae' ~ '20',
                               species == 'Anthonomus grandis' ~ '27',
                               species == 'Helicoverpa armigera' ~ '28',
                               species == 'Antestiopsis thunbergii' ~ '29',
                               species == 'Monochamus leuconotus' ~ '30',
                               species == 'Culex pipiens' ~ '31',
                               species == 'Phenacoccus solenopsis' ~ '32',
                               species == 'Culex quinquefasciatus' ~ '33',
                               species == 'Thrips hawaiiensis' ~ '34')) %>%
  arrange(curve_ID) %>% 
  filter(curve_ID != 'NA')

  
zeta %>% distinct(species) %>% print(n=30)

RearZeta <- bind_cols(zeta,RearZeta)

Rzeta_model <- lm(estimate ~ RTemp, RearZeta); summary(Rzeta_model)

Rear_z <- 
  ggplot(RearZeta, aes(x=RTemp, y=estimate)) +
  geom_smooth(method = 'lm', colour = '#636363', size=0.15, fill = NA)+
  geom_linerange(aes(ymin =conf_lower , ymax = conf_upper),size= 0.1,col="#000000") +
  geom_point(size = 2.5,stroke=0.2, col = '#000000', shape=22, fill ="#a6cee3")+
  theme_bw()+
  scale_x_continuous(expression(plain(paste("Rearing Temperature (",degree,"C)"))),
                     limits =c(22.5,30.75),
                     expand = c(0, 0),
                     breaks=seq(21,32, by=1))+
  scale_y_continuous(expression(paste(italic(T[pk]), " of Adult Mortality Rate (", italic(z),")")),
                     limits =c(5,28),
                     expand = c(0, 0),
                     breaks=seq(10,25, by=5))+
  theme(legend.position = 'none',legend.text = element_text(size = 10))+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = 23, y = 26.5,label = "C"), 
            parse = TRUE, size = 4, colour = "black")

#======================================================
# bmax 

RearDF <- as_tibble(read_csv('../data/TraitData.csv')) %>% 
  select(standardisedtraitname,interactor1,interactor1growthtemp) %>%
  filter(interactor1growthtemp != 'not stated' & 
           interactor1growthtemp != 'various' & interactor1growthtemp != 'NA' ) %>%
  rename(trait = standardisedtraitname, species = interactor1, RearTemp = interactor1growthtemp) %>%
  mutate(RearTemp = as.numeric(RearTemp))

RearBeta <- RearDF %>% filter(trait == 'bmax') %>% 
  group_by(species) %>% 
  summarise(avg = mean(RearTemp)) %>%
  arrange(avg) %>%
  mutate(Rcurve_ID = case_when(species == 'Plutella xylostella' ~ '1',
                               species == 'Drosophila suzukii' ~ '2',
                               species == 'Paracoccus marginatus' ~ '3',
                               species == 'Acyrthosiphon pisum' ~ '4',
                               species == 'Tetraneura nigriabdominalis' ~ '5',
                               species == 'Stethorus punctillum' ~ '6',
                               species == 'Tetranychus mcdanieli' ~ '7',
                               species == 'Muscidifurax zaraptor' ~ '8',
                               species == 'Aphis nasturtii' ~ '9',
                               species == 'Rhopalosiphum maidis' ~ '10',
                               species == 'Hylobius transversovittatus' ~ '11',
                               species == 'Helicoverpa armigera' ~ '12',
                               species == 'Antestiopsis thunbergii' ~ '13',
                               species == 'Monochamus leuconotus' ~ '14',
                               species == 'Culex pipiens' ~ '15',
                               species == 'Phenacoccus solenopsis' ~ '16',
                               species == 'Culex quinquefasciatus' ~ '17',
                               species == 'Thrips hawaiiensis' ~ '18',
                               species == 'Iphiseius degenerans' ~ '19',
                               species == 'Kampimodromus aberrans' ~ '20',
                               species == 'Planococcus citri' ~ '21',
                               species == 'Amblyseius womersleyi' ~ '22',    
                               species == 'Anthonomus grandis'   ~ '23',      
                               species == 'Callosobruchus analis'  ~ '24',    
                               species == 'Callosobruchus chinensis' ~ '25',  
                               species == 'Callosobruchus maculatus'  ~ '26', 
                               species == 'Callosobruchus rhodesianus' ~ '27',
                               species == 'Aedes krombeini' ~ '29')) %>%
  arrange(Rcurve_ID) %>% rename(Rspecies = species, RTemp = avg) %>%
  filter(Rcurve_ID != 'NA')

RearBeta %>% distinct(Rspecies, Rcurve_ID) %>% print(n=60)


beta  <- as_tibble(read.csv('../data/bmax_Tpks_AllParams.csv', header = TRUE)) %>%
  filter(param == 'topt') %>% 
  select(species, estimate, conf_lower, conf_upper, trait) %>%
  mutate(curve_ID = case_when(species == 'Plutella xylostella' ~ '1',
                               species == 'Drosophila suzukii' ~ '2',
                               species == 'Paracoccus marginatus' ~ '3',
                               species == 'Acyrthosiphon pisum' ~ '4',
                               species == 'Tetraneura nigriabdominalis' ~ '5',
                               species == 'Stethorus punctillum' ~ '6',
                               species == 'Tetranychus mcdanieli' ~ '7',
                               species == 'Muscidifurax zaraptor' ~ '8',
                               species == 'Aphis nasturtii' ~ '9',
                               species == 'Rhopalosiphum maidis' ~ '10',
                               species == 'Hylobius transversovittatus' ~ '11',
                               species == 'Helicoverpa armigera' ~ '12',
                               species == 'Antestiopsis thunbergii' ~ '13',
                               species == 'Monochamus leuconotus' ~ '14',
                               species == 'Culex pipiens' ~ '15',
                               species == 'Phenacoccus solenopsis' ~ '16',
                               species == 'Culex quinquefasciatus' ~ '17',
                               species == 'Thrips hawaiiensis' ~ '18',
                               species == 'Iphiseius degenerans' ~ '19',
                               species == 'Kampimodromus aberrans' ~ '20',
                               species == 'Planococcus citri' ~ '21',
                               species == 'Amblyseius womersleyi' ~ '22',    
                               species == 'Anthonomus grandis'   ~ '23',      
                               species == 'Callosobruchus analis'  ~ '24',    
                               species == 'Callosobruchus chinensis' ~ '25',  
                               species == 'Callosobruchus maculatus'  ~ '26', 
                               species == 'Callosobruchus rhodesianus' ~ '27',
                               species == 'Aedes krombeini' ~ '29')) %>%
  arrange(curve_ID) %>% 
  filter(curve_ID != 'NA')


beta %>% distinct(species) %>% print(n=30)

RearBeta <- bind_cols(beta,RearBeta)

Rbeta_model <- lm(estimate ~ RTemp, RearBeta); summary(Rbeta_model)

Rear_beta <- 
  
  ggplot(RearBeta, aes(x=RTemp, y=estimate)) +
  geom_smooth(method = 'lm', colour = '#000000', size=0.15, fill = NA)+
  geom_linerange(aes(ymin =conf_lower , ymax = conf_upper),size= 0.1,col="#000000") +
  geom_point(size = 2.5,stroke=0.2, col = '#000000', fill="#fdb863", shape=23)+
  theme_bw()+
  scale_x_continuous(expression(plain(paste("Rearing Temperature (",degree,"C)"))),
                     limits =c(21.5,27.25),
                     expand = c(0, 0),
                     breaks=seq(22,27, by=1))+
  scale_y_continuous(expression(paste(italic(T[pk]), " of Fecundity Rate (", italic(b[max]),")")),
                     limits =c(19,36),
                     expand = c(0, 0),
                     breaks=seq(20,35, by=5))+
  theme(legend.position = 'none',legend.text = element_text(size = 10))+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = 22, y = 35, label = "B"), 
            parse = TRUE, size = 4, colour = "black")


FigSI <- ggarrange(Rear_alp,Rear_beta,Rear_z,Rear_zj,nrow=2,ncol=2)


ggsave("../results/SI/RearTemp.pdf",FigSI, width = 15, height =15, 
       units = "cm",device = cairo_pdf)

