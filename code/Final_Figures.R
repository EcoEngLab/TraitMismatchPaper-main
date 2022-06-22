#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#                   All Figures                    #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

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

# Figure 3: Tpks + activation energy plot 

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

#order by development rate alpha
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
                               expression(plain(paste("Juvenile Development Time (",italic(alpha),")")))),
                    values = c("#1f78b4","#a6cee3","#fdb863","#e66101"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow = 1,ncol =4,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(labels = c(expression(plain(paste("Juvenile Mortality Rate (",italic(z[J]),")"))),
                                 expression(plain(paste("Adult Mortality Rate (",italic(z),")"))),
                                 expression(plain(paste("Fecundity (",italic(b[max]),")"))),
                                 expression(plain(paste("Juvenile Development Time (",italic(alpha),")")))),
                      values = c("#1f78b4","#a6cee3","#fdb863","#e66101"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=1,ncol=4,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  scale_shape_manual(labels = c(expression(plain(paste("Juvenile Mortality Rate (",italic(z[J]),")"))),
                                expression(plain(paste("Adult Mortality Rate (",italic(z),")"))),
                                expression(plain(paste("Fecundity (",italic(b[max]),")"))),
                                expression(plain(paste("Juvenile Development Time (",italic(alpha),")")))),
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

activation_e <-  filter(topt, trait == 'juvenile development rate' | trait == 'fecundity rate', param =="e")  

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
                       trait = case_when(trait == "fecundity rate" ~ "Peak Fecundity",
                                         trait == "juvenile development rate" ~ "Development Time"))
                                         

TraitLevel <- c("Peak Fecundity", "Development Time" )
activation_e$trait  <- as.character(activation_e$trait)
activation_e$trait  <- factor(activation_e$trait, levels = TraitLevel)
activation_e$species <- as.character(activation_e$species)
activation_e$species <- factor(activation_e$species, levels=SPorder)
activation_e <- na.omit(activation_e)


fig3b <- ggplot(activation_e, aes(estimate, species, shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(xmin = conf_lower, xmax = conf_upper),width=0.35,size=0.35) +
  geom_point(size = 2.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("Activation Energy (",italic(E),")"))),
                     limits =c(0,6),
                     expand = c(0, 0),
                     breaks=seq(0,6, by=1))+
  scale_fill_manual(values = c("#fdb863","#e66101"),
                    name=expression(bold("")))+
  scale_colour_manual(values = c("#fdb863","#e66101"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=1,ncol=2,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  scale_shape_manual(values = c(23,24),
                     name=expression(bold("")),
                     guide = guide_legend(nrow = 1,ncol =2,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  theme(legend.position = 'none',legend.text = element_text(size = 8.5),
        axis.text.y = element_markdown(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  theme(legend.margin=margin(t = -0.4, unit='cm'))+
  geom_vline(xintercept=0.65, lty ='dashed')+
  theme(axis.text.y=element_blank())+
  theme(plot.margin=margin(l=-4,unit="cm"))+
  theme(aspect.ratio=2)+
  ggtitle("B")+
  theme(plot.title=element_text(face="bold", size = 15, vjust = -1))

fig3b

plotMain <- fig3a+theme(legend.position="none")+fig3b+theme(legend.position="none")
# plotMain
legend <- get_legend(fig3a+theme(legend.position = "bottom"))


# plotMain <- plot_grid(plotMain, NULL, rel_widths = c(1,0))

fig3 <- plot_grid(plotMain, legend,ncol = 1, rel_heights = c(1, 0.05))

# 
save_plot(fig3, file="../results/Fig3.pdf",
          base_height=15,base_asp=1.5, units="cm")

#=======================================================================
# Figure 5: rm_opt plots 

rm(list=ls())
graphics.off()

# plot mass corrected value of r_m at Topt against the peak temperature for r_m
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

#rmTpk_model <- lm(rm_massCor ~ rmTpk, rm_data); summary(rmTpk_model)

#prepare data for plotting
rm_data <- rm_data %>% 
           mutate(rm_massCor = log(rm_opt/mass^coef(rm_model)[2]))
           
#plot rm_opt in 1/kT, correcting for mass
MassCorrectedrm_opt <- 
  rm_data %>%
  ggplot(aes(x = rmTpk, y = rm_massCor)) +
  geom_smooth(method = 'lm', size=0.3, col='#636363', fill="#004225")+
  scale_y_continuous(expression(plain(paste("ln(",italic(r[paste(m,",", opt)])," / ",
                                            italic(M^-0.16),"))"))),
                     limits=c(-4.5,-0.5),
                     expand = c(0.01, 0),
                     breaks=seq(-4,-1, by=1))+
  scale_x_continuous(expression(plain(paste(italic("T"[pk])))),
                     limits=c(19,34.5),
                     expand = c(0.01, 0),
                     breaks=seq(20,34, by=2))+
  geom_point(aes(shape = species, fill = species), size=2,stroke=0.25)+
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
  theme(text = element_text(size=9),
        legend.position = 'none',
        legend.text = element_text(size = 4.25, face = 'italic'),
        legend.background = element_rect(colour = "black", size = 0.125), 
        legend.margin=margin(t = 0.05, b = 0.1, r=0.1,l=0.1, unit='cm'),
        legend.key.size = unit(0.2, 'cm')) +
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.5,vjust=1.4,
                label = "A"),size = 5, colour = "black")+
  theme(axis.title.y = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

MassCorrectedrm_opt

#save_plot(MassCorrectedrm_opt, file="../results/MassCorrectedrm_opt_Tpk.pdf", 
#           base_height=9,base_width = 9, base_asp = 0.75,units="cm")


#=======================================================
# sum of trait Tpks versus mass-corrected r_m, opt

Species <- rm_data$species

alphaTpks <- as_tibble(read.csv('../data/alpha_Tpks_AllParams.csv')) %>% 
  filter(param=='topt') %>%
  filter(species %in% Species)%>%
  select(param,species,estimate, conf_lower, conf_upper,trait) %>%
  mutate(estimate = as.numeric(estimate), 
         conf_lower = as.numeric(conf_lower), conf_upper = as.numeric(conf_upper))

zjTpks <- as_tibble(read.csv('../data/zj_Tpks_AllParams.csv')) %>% 
  filter(param=='topt') %>%
  filter(species %in% Species)%>%
  select(param,species,estimate, conf_lower, conf_upper,trait) %>%
  mutate(estimate = as.numeric(estimate), 
         conf_lower = as.numeric(conf_lower), 
         conf_upper = as.numeric(conf_upper))


zTpks <- as_tibble(read.csv('../data/z_Tpks_AllParams.csv')) %>% 
  filter(param=='topt') %>%
  filter(species %in% Species)%>%
  select(param,species,estimate, conf_lower, conf_upper,trait) %>%
  mutate(estimate = as.numeric(estimate), 
         conf_lower = as.numeric(conf_lower), 
         conf_upper = as.numeric(conf_upper))

bmaxTpks <- as_tibble(read.csv('../data/bmax_Tpks_AllParams.csv')) %>% 
  filter(param=='topt') %>%
  filter(species %in% Species)%>%
  select(param,species,estimate, conf_lower, conf_upper,trait) %>%
  mutate(estimate = as.numeric(estimate), 
         conf_lower = as.numeric(conf_lower), 
         conf_upper = as.numeric(conf_upper))


AllTpks <- bind_rows(alphaTpks,zjTpks,zTpks,bmaxTpks)

#write_csv(AllTpks, 'AllTpkParams.csv')

# Calculate of Tpks

OptSum <- AllTpks %>% group_by(species) %>%
  summarise(sum = sum(estimate))

SumPlot <- left_join(rm_data, OptSum) 
           
SumTpks_plot <-
  SumPlot %>%
  ggplot(aes(x = sum, y = rm_massCor))+
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="#004225")+
  scale_y_continuous(expression(plain(paste("Log mass-corrected ",italic(r[m])," at its optimal temperature (ln(",italic(r[paste(m,",", opt)])," / ",
                                            italic(M^-0.16),"))"))),
                     limits=c(-3.8,-0.5),
                     expand = c(0.01, 0),
                     breaks=seq(-3,-1, by=1))+
  scale_x_continuous(expression(plain(paste("Sum of ",italic("T"[pk]),"'s"))))+
  geom_point(aes(shape=species, 
                 fill=species),
             size=2,
             stroke=0.25)+
  scale_shape_manual(values = c(21,21,21,21,
                                22,22,22,22,
                                23,23,23,23,
                                24,24,24,24),
                     name=expression(bold("")),
                     guide = guide_legend(nrow=16,ncol=1,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  scale_fill_manual(values = c("#EEEEEE","#CCCCCC","#999999","#333333",
                               "#EEEEEE","#CCCCCC","#999999","#333333",
                               "#EEEEEE","#CCCCCC","#999999","#333333",
                               "#EEEEEE","#CCCCCC","#999999","#333333"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow=16,ncol=1,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(values = c("#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=16,ncol=1,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  theme_bw()+
  theme(text=element_text(family="Times"))+
  theme(text = element_text(size=9),
        legend.position = 'none',
        legend.text = element_text(size = 4.5, face = 'italic'),
        legend.background = element_rect(colour = "black", size = 0.125), 
        legend.margin=margin(t = 0.05, b = 0.1, r=0.1,l=0.1, unit='cm'),
        legend.key.size = unit(0.3, 'cm'))+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.5,vjust=1.4,
                label = "C"),size = 5, colour = "black")+
  theme(axis.title.y = element_text(hjust=0.1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

SumTpks_plot

#save_plot(SumTpks_plot, file="../results/MassCorrectedrm_optSumTpks.pdf", 
#          base_height=9,base_width = 9, base_asp = 0.75,units="cm")


#=======================================================================

# Relationship between r_m opt and 1/alpha Tpk

alphaMass <- as_tibble(read_csv('../data/a_pksT_pksMass.csv')) %>% 
  select(species, a_pk, a_pkLwr, a_pkUpr, mass) %>%
  mutate(curve_ID = case_when(species == 'Anoplophora glabripennis' ~ '1',
                              species == 'Halyomorpha halys' ~ '2',       
                              species == 'Aedes aegypti' ~ '3',
                              species == 'Anthonomus grandis' ~ '4',
                              species == 'Paracoccus marginatus' ~ '5',
                              species == 'Acyrthosiphon pisum' ~ '6',
                              species == 'Aphis gossypii' ~ '7',
                              species == 'Bemisia tabaci' ~ '8',
                              species == 'Tetraneura nigriabdominalis' ~ '9',
                              species == 'Stethorus punctillum' ~ '10',
                              species == 'Tetranychus mcdanieli' ~ '11',
                              species == 'Muscidifurax zaraptor' ~ '12',
                              species == 'Aphis nasturtii' ~ '13',
                              species == 'Rhopalosiphum maidis' ~ '14',
                              species == 'Thrips hawaiiensis' ~ '15',
                              species == 'Helicoverpa armigera' ~ '16')) %>%
  filter(curve_ID != 'NA') %>% 
  rename(alphaspecies = species, alphamass = mass, alphacurve = curve_ID ) %>%
  arrange(alphacurve)

alpharm_data <- bind_cols(rm_data, alphaMass) %>% 
  select(-alphaspecies,-alphamass, -alphacurve)

# linear model (note the allometry is linear in log-log scale)
alpharm_model <- lm(log(rm_opt/mass^-0.16) ~ log(a_pk/mass^-0.265), data = alpharm_data)
summary(alpharm_model)
anova(alpharm_model)

#plot a_pk vs r_opt

MassCorrectedrm_optAlphaTpk <-
  alpharm_data %>%
  ggplot(aes(x=log(a_pk/mass^-0.265), y = log(rm_opt/mass^-0.16)))+
  geom_point()+
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="#004225")+
  scale_y_continuous(expression(plain(paste("ln(",italic(r[paste(m,",", opt)])," / ",
                                            italic(M^-0.16),"))"))),
                     limits=c(-4.5,-0.5),
                     expand = c(0.01, 0),
                     breaks=seq(-4,-1, by=1))+
  scale_x_continuous(expression(plain(paste("ln(", italic(1/alpha[pk])~")/",
                                            italic(M^-0.265),")"))))+
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
                                           title.hjust=0.4))+
  theme_bw()+
  theme(text=element_text(family="Times"))+
  theme(text = element_text(size=9),
        legend.position = 'none',
        legend.text = element_text(size = 5, face = 'italic'),
        legend.background = element_rect(colour = "white", size = 0.125), 
        legend.margin=margin(t = 0.01, b = 0.1, r=0.1,l=0.1, unit='cm'),
        legend.key.size = unit(0.3, 'cm'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.5,vjust=1.4,
                label = "B"),size = 5, colour = "black")+
  theme(axis.title.y = element_blank())
  
        
#save_plot(MassCorrectedrm_optAlphaTpk, file="../results/MassCorrectedrm_optAlphaTpk.pdf", 
#          base_height=5,base_width = 5, base_asp = 0.75,units="cm")

MassCorrectedrm_optAlphaTpk


# plotMain
legend <- get_legend(MassCorrectedrm_optAlphaTpk+theme(legend.position = c(0.5,0.7)))

Fig5 <- plot_grid(MassCorrectedrm_opt+MassCorrectedrm_optAlphaTpk+
                    SumTpks_plot+legend, nrow = 1, rel_widths = c(1,2,1,1))
                    

save_plot(Fig5, file="../results/Fig5.pdf", 
          base_height=10,base_width = 11, base_asp = 0.75,units="cm")


#==========================================================
#   Figure 6: Analysis of 'hotter is better' theory

# Fig. 6a -- 1/alpha

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
  geom_smooth(method = 'lm', colour = '#636363', size=0.15, fill="#004225")+
  scale_y_continuous(expression(plain(paste("ln((", italic(1/alpha[pk])~")/",italic(M^-0.265),")"))))+
  scale_x_continuous(expression(plain(paste(italic("T"[pk])))),
                     limits =c(16.5,36),
                     expand = c(0, 0),
                     breaks=seq(20,35, by=5))+
  theme_bw()+
  theme(text=element_text(family="Times", size=8), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_point(size = 1, col="#000000",stroke=0.1, shape = 24, fill="#e66101")

MassCorrectedApkTpk

#save_plot(MassCorrectedApkTpk, file="../results/MassCorrectedApkTpk.pdf", 
#          base_height=5,base_width = 6, base_asp = 0.75,units="cm")

#============================

# Fig. 6b -- zj

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
  geom_smooth(method = 'lm', colour = '#636363', size=0.15, fill = '#004225')+
  scale_y_continuous(expression(plain(paste("ln(", italic(z[J][pk])~"/",italic(M^-0.193),")"))))+
  scale_x_continuous(expression(plain(paste(italic("T"[pk])))))+
  theme_bw()+
  theme(text=element_text(family="Times", size=8))+
  geom_point(size = 1, col="#000000",stroke=0.1, shape = 21, fill="#1f78b4")+
  theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())

MassCorrectedzjpkTpk

#save_plot(MassCorrectedzjpkTpk, file="../results/MassCorrectedzjpkTpk.pdf", 
#          base_height=5,base_width = 6, base_asp = 0.75,units="cm")


#============================

# Fig. 6c -- z

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
  geom_smooth(method = 'lm', colour = '#636363', size=0.15, fill='#004225')+
  scale_y_continuous(expression(plain(paste("ln(", italic(z[pk])~"/",italic(M^-0.124),")"))))+
  scale_x_continuous(expression(plain(paste(italic("T"[pk])))),
                     limits=c(8,27),
                     expand = c(0, 0),
                     breaks=seq(10,25, by=5))+
  theme_bw()+
  theme(text=element_text(family="Times", size=8))+
  geom_point(size = 1, col="#000000",stroke=0.1, shape = 22, fill="#a6cee3")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


MassCorrectedzpkTpk

#save_plot(MassCorrectedzpkTpk, file="../results/MassCorrectedzpkTpk.pdf", 
#          base_height=5,base_width = 6, base_asp = 0.75,units="cm")

# Fig. 6c -- bmax

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
  geom_smooth(method = 'lm', colour = '#636363', size=0.15, fill='#004225')+
  scale_y_continuous(expression(plain(paste("ln(", italic(b[max])~"/",italic(M^0.08),")"))))+
  scale_x_continuous(expression(plain(paste(italic("T"[pk])))))+
  theme_bw()+
  theme(text=element_text(family="Times", size=8))+
  geom_point(size = 1, col="#000000",stroke=0.1, shape = 23, fill="#fdb863")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

MassCorrectedBmaxTpk

#save_plot(MassCorrectedBmaxTpk, file="../results/MassCorrectedBmaxTpk.pdf", 
#          base_height=5,base_width = 6, base_asp = 0.75,units="cm")


#========= plot hotter-is-better panel 

Fig6 <- (MassCorrectedApkTpk+MassCorrectedBmaxTpk)/(MassCorrectedzpkTpk+MassCorrectedzjpkTpk); Fig6

save_plot(Fig6, file="../results/Fig6.pdf", 
          base_height=10,base_width = 12, base_asp = 0.75,units="cm")

