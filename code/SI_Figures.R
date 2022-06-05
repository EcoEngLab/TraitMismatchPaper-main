
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

setwd("~/Dropbox/TraitTesting/data") # n.b. working directory is data directory

#======================================

# 1. Variance  #######

AllTpks <- as_tibble(read_csv('AllTpkParams.csv'))
rm_data <- as_tibble(read_csv('rm_optSizeScaling.csv'))

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
  mutate(rm_massCor_lwr = log(rm_optLwr/mass^coef(rm_model)[2]),
         rm_massCor     = log(rm_opt/mass^coef(rm_model)[2]),
         rm_massCor_upr = log(rm_optUpr/mass^coef(rm_model)[2])) %>%
  mutate_at(vars(c(rm_massCor_lwr)), 
            ~ifelse(rm_massCor_lwr == 'NaN', -2.8765094, .))

#Calculate variance and Sum

OptVar <- AllTpks %>% group_by(species) %>%
  summarise(variance = var(estimate))
OptSum <- AllTpks %>% group_by(species) %>%
  summarise(sum = sum(estimate))


str(VarPlot)

###### 2. Variance Plot #######
VarPlot <- left_join(rm_data, OptVar, by="species")
VarLM <-   lm(VarPlot$rm_massCor~VarPlot$variance)

summary(VarLM)

VarLMData <- predict(VarLM, interval = "confidence")
VarPlot <- cbind(VarPlot,VarLMData)

Var_rm <- 
  ggplot(VarPlot, aes(x=variance, y=rm_massCor, colour = species, fill = species)) +
  geom_linerange(aes(ymin = rm_massCor_lwr, ymax = rm_massCor_upr),col="#000000", size = 0.2) +
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
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill=NA)+
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
  theme(text = element_text(size=8),
        legend.position = 'right',
        legend.text = element_text(size = 6, face = 'italic'),
        legend.background = element_rect(colour = "black", size = 0.125), 
        legend.margin=margin(t = 0.05, b = 0.1, r=0.1,l=0.1, unit='cm'),
        legend.key.size = unit(0.3, 'cm'))
  
  
save_plot(Var_rm, file="../results/SI1.pdf", 
          base_height=9,base_width = 12, base_asp = 0.75,units="cm")



####### 2. Activation Energy historgram ########

rm(list=ls())
graphics.off()

alpha   <- as_tibble(read.csv('alpha_Tpks_AllParams.csv', header = TRUE))
zj      <- as_tibble(read.csv('zj_Tpks_AllParams.csv', header = TRUE))
z       <- as_tibble(read.csv('z_Tpks_AllParams.csv', header = TRUE))
bmax    <- as_tibble(read.csv('bmax_Tpks_AllParams.csv', header = TRUE))

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
  labs(x="Activation Energy (E)", y="Density")+
  geom_histogram(position="identity", binwidth=1)+
  geom_density()+
  facet_grid(.~trait)


ggsave("../results/eHist.pdf",eplot, width = 9, height = 9, units = "cm")








































#Read data
rmpredictions <- read.csv("rmprediction.csv")

rmPeak <- rmpredictions %>% group_by(species) %>%
  summarise(which(rmpredictions$estimate == max(estimate)))

rmPeakdata <- rmpredictions[rmPeak$`which(rmpredictions$estimate == max(estimate))`, ]
rmPeakdata <- rmPeakdata %>% rename(rmax = estimate)


Species <- rmPeakdata$species

alpha  <- as_tibble(read.csv('alpha_Tpks.csv', header = TRUE)) %>%
  filter(param=="topt") %>%
  filter(species %in% Species)%>%
  select(species,estimate, conf_lower, conf_upper,trait)

bpk  <- as_tibble(read.csv('bpk_Tpks.csv', header = TRUE)) %>%
  filter(param=="topt") %>%
  filter(species %in% Species)%>%
  select(species,estimate, conf_lower, conf_upper,trait)

z  <- as_tibble(read.csv('z_Tpks.csv', header = TRUE)) %>%
  filter(param=="topt") %>%
  filter(species %in% Species)%>%
  select(species,estimate, conf_lower, conf_upper,trait)

zj  <- as_tibble(read.csv('zj_Tpks.csv', header = TRUE)) %>%
  filter(param=="topt") %>%
  filter(species %in% Species)%>%
  select(species,estimate, conf_lower, conf_upper,trait)

AllOpt<- rbind(alpha,bpk,z,zj)
#Calculate variance and Sum

OptVar <- AllOpt %>% group_by(species) %>%
  summarise(variance = var(estimate))
OptSum <- AllOpt %>% group_by(species) %>%
  summarise(sum = sum(estimate))




##changing order of df for plotting
colors <- c("#d9d9d9","#bdbdbd","#969696","#737373","#525252","#252525")
OrderCol <- colors[order(rmPeakdata$rmax)]
SPorder <- rmPeakdata$species[order(rmPeakdata$rmax)]













## 1. Fitness with alpha
AlpPlot <- left_join(rmPeakdata, alpha, by="species")
AlpLM <- lm(AlpPlot$rmax~AlpPlot$estimate)

AlpLMData <- predict(AlpLM, interval = "confidence")
AlpPlot <- cbind(AlpPlot,AlpLMData)

AlpPlot$species <- as.character(AlpPlot$species)
AlpPlot$species <- factor(AlpPlot$species, levels=SPorder)

Alp_rmax <- ggplot(AlpPlot, aes(x=estimate, y=rmax, col=factor(species, levels=SPorder))) +
  geom_errorbar(aes(ymin = rmax_lwr, ymax = rmax_upr),col="#000000") +
  geom_point(size = 2.5,stroke=0.2)+
  theme_bw(base_size = 12.5) +
  # theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste(italic(T)[pk], " of Juvenile Development Rate (", italic(alpha),")"))),
                     limits =c(26,35),
                     expand = c(0, 0),
                     breaks=seq(26,35, by=1))+
  scale_y_continuous(expression(plain(paste(italic(r[m])," at Optimal Temperature (",italic(r[paste(m,",", opt)]), ")"))),
                     limits =c(0,0.35),
                     expand = c(0, 0),
                     breaks=seq(0,1, by=0.1))+
  scale_colour_manual(values = colors,
                      name=expression(bold("species")),
                      labels =SPorder)+
  theme(legend.position = 'none',legend.text = element_text(size = 10),
        axis.text.y = element_text(face = 'italic'),aspect.ratio = 1)+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.5,vjust=1.4,
                label = "A"),size = 5, colour = "black")+
  theme(legend.title = element_blank())+
  geom_line(aes(estimate,fit),AlpPlot,size=0.35,col="#636363")+
  theme(axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10))

Alp_rmax

###### 2. Variance Plot #######
VarPlot <- left_join(rmPeakdata, OptVar, by="species")
VarLM <- lm(VarPlot$rmax~VarPlot$variance)

VarLMData <- predict(VarLM, interval = "confidence")
VarPlot <- cbind(VarPlot,VarLMData)

Var_rmax<- ggplot(VarPlot, aes(x=variance, y=rmax, col=factor(species, levels=SPorder))) +
  geom_errorbar(aes(ymin = rmax_lwr, ymax = rmax_upr),col="#000000") +
  geom_point(size = 2.5, stroke=0.2)+
  theme_bw(base_size = 12.5) +
  # theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("Variance of ", italic(T)[pk]))),
                     limits =c(20,90),
                     expand = c(0, 0),
                     breaks=seq(20,90, by=10))+
  scale_y_continuous(expression(plain(paste(italic(r[m])," at Optimal Temperature (",italic(r[paste(m,",", opt)]), ")"))),
                     limits =c(0,0.31),
                     expand = c(0, 0),
                     breaks=seq(0,1, by=0.1))+
  scale_colour_manual(values = colors,
                      name=expression(bold("species")),
                      labels =SPorder)+
  theme(legend.position = 'none',legend.text = element_text(size = 10),
        axis.text.y = element_text(face = 'italic'),aspect.ratio = 1)+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.5,vjust=1.4,
                label = "B"),size = 5, colour = "black")+
  theme(legend.title = element_blank())+
  geom_line(aes(variance,fit),VarPlot,size=0.35,col="#636363")+
  theme(axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10))
Var_rmax

si1 <- ggarrange(Alp_rmax,Var_rmax,nrow=1,ncol=2,
                 common.legend = T, legend="right")

ggsave("../results/SI1.pdf",si1, width = 20, height =10, 
       units = "cm",device = cairo_pdf)



####### 2. Activation Energy historgram ########

rm(list=ls())
graphics.off()

##Read in data ##
alpha  <- as_tibble(read.csv('alpha_Tpks_AllParams.csv', header = TRUE))
zj     <- as_tibble(read.csv('zj_Tpks_AllParams.csv', header = TRUE))
z      <- as_tibble(read.csv('z_Tpks_AllParams.csv', header = TRUE))
bpk    <- as_tibble(read.csv('bpk_Tpks_AllParams.csv', header = TRUE))

toptz  <- rbind(alpha,zj,z,bpk)

### activation energy###
topte <- filter(toptz, param =="e" )
write.csv(topte, file="../data/e_Estimate.csv")

topte <- mutate(topte, trait = case_when(trait =="adult mortality rate" ~ "Adult Mortality Rate",
                                         trait == "fecundity" ~ "Peak Fecundity",
                                         trait == "juvenile development rate" ~ "Development Rate",
                                         trait == "juvenile mortality rate" ~ "Juvenile Mortality Rate"))


eplot <- ggplot(topte, aes(x=estimate))+
  theme_bw()+
  labs(x="Activation Energy (E)", y="Density")+
  geom_histogram(position="identity", binwidth=1)+
  geom_density()+
  facet_grid(.~trait)
# eplot

ggsave("../results/eHist.pdf",eplot, width = 20, height = 10, units = "cm")




#########3. Rearing Temperature ##########
# 

RearDF <- as_tibble(read_csv('TraitData.csv')) %>% 
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
  
alpha   <- as_tibble(read.csv('alpha_Tpks_AllParams.csv', header = TRUE)) %>%
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
  geom_point(size = 2.5,stroke=0.2, col = '#000000', shape=24, fill ='#FFFFFF')+
  theme_bw(base_size = 12.5)+
  scale_x_continuous(expression(plain(paste("Rearing Temperature (",degree,"C)"))),
                    limits =c(20,30.5),
                    expand = c(0, 0),
                    breaks=seq(21,29, by=1))+
  scale_y_continuous(expression(paste(italic(T)[pk], " of Juvenile Development Time (", italic(alpha),")")),
                     limits =c(21,38),
                     expand = c(0, 0),
                     breaks=seq(22,36, by=2))+
  theme(legend.position = 'none',legend.text = element_text(size = 10))
  


+
  theme(text=element_text(family="Times"))



+
  geom_text(aes(x = 23.5, y = 37,label = "A"), 
            parse = TRUE, size = 4, colour = "black")















# Tpks for all species

alpha   <- as_tibble(read.csv('alpha_Tpks_AllParams.csv', header = TRUE))
zj      <- as_tibble(read.csv('zj_Tpks_AllParams.csv', header = TRUE))
z       <- as_tibble(read.csv('z_Tpks_AllParams.csv', header = TRUE))
bmax    <- as_tibble(read.csv('bmax_Tpks_AllParams.csv', header = TRUE))

topt  <- rbind(alpha,zj,z,bmax)
toptz <- subset(topt, topt$param=="topt")

toptz$species <- as.factor(toptz$species)
toptz$trait   <- as.factor(toptz$trait)


toptz  <- rbind(alpha,zj,z,bpk)


toptz <-toptz %>% filter(param=="topt") %>%
  left_join(RearDF) %>%
  mutate(RearTemp = as.numeric(RearTemp))

#Traits Plots
##Alpha
alp <- subset(toptz, toptz$trait=="juvenile development rate")
alp <- drop_na(alp, "RearTemp")
alpLM <- lm(alp$estimate~ alp$RearTemp)

alpLMData <- predict(alpLM, interval = "confidence")
alp <- cbind(alp,alpLMData)


Rear_alp <- ggplot(alp, aes(x=RearTemp, y=estimate)) +
  geom_errorbar(aes(ymin =conf_lower , ymax = conf_upper),width = 0.2,col="#000000") +
  geom_point(size = 2.5,stroke=0.2)+
  theme_bw(base_size = 12.5) +
  # theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("Rearing Temperature (C",degree,")"))),
                     limits =c(23,28),
                     expand = c(0, 0),
                     breaks=seq(0,50, by=1)
  )+
  scale_y_continuous(expression(plain(paste(italic(T)[pk], " of Juvenile Development Rate (", italic(alpha),")"))),
                     limits =c(23,38),
                     expand = c(0, 0),
                     breaks=seq(0,50, by=2)
  )+
  theme(legend.position = 'none',legend.text = element_text(size = 10),
        axis.text.y = element_text(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = 23.5, y = 37,label = "A"), 
            parse = TRUE, size = 6, colour = "black") +
  geom_line(aes(RearTemp,fit),alp,size=0.5,col="#636363")
# geom_abline(intercept = 0, slope=1, linetype="dashed")

# Rear_alp


##Fecundity###
Bopt <- subset(toptz, toptz$trait=="fecundity")
Bopt <- drop_na(Bopt, "RearTemp")
BoptLM <- lm(Bopt$estimate~ Bopt$RearTemp)

BoptLMData <- predict(BoptLM, interval = "confidence")
Bopt <- cbind(Bopt,BoptLMData)


Rear_Bopt <- ggplot(Bopt, aes(x=RearTemp, y=estimate)) +
  geom_errorbar(aes(ymin =conf_lower , ymax = conf_upper),width = 0.2,col="#000000") +
  geom_point(size = 2.5,stroke=0.2)+
  theme_bw(base_size = 12.5) +
  # theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("Rearing Temperature (C",degree,")"))),
                     limits =c(23,28),
                     expand = c(0, 0),
                     breaks=seq(0,50, by=1)
  )+
  scale_y_continuous(expression(plain(paste(italic(T)[pk], " of Fecundity (", italic(B[max]),")"))),
                     limits =c(20,32),
                     expand = c(0, 0),
                     breaks=seq(0,50, by=2)
  )+
  theme(legend.position = 'none',legend.text = element_text(size = 10),
        axis.text.y = element_text(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = 23.5, y = 31.25,label = "B"),
            parse = TRUE, size = 6, colour = "black") +
  geom_line(aes(RearTemp,fit),Bopt,size=0.5,col="#636363")

# Rear_Bopt

##Adult Mortality Rate ###
Z <- subset(toptz, toptz$trait=="adult mortality rate")
Z <- drop_na(Z, "RearTemp")
ZLM <- lm(Z$estimate~ Z$RearTemp)

ZLMData <- predict(ZLM, interval = "confidence")
Z <- cbind(Z,ZLMData)


Rear_Z <- ggplot(Z, aes(x=RearTemp, y=estimate)) +
  geom_errorbar(aes(ymin =conf_lower , ymax = conf_upper),width = 0.2,col="#000000") +
  geom_point(size = 2.5,stroke=0.2)+
  theme_bw(base_size = 12.5) +
  # theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("Rearing Temperature (C",degree,")"))),
                     limits =c(23,28),
                     expand = c(0, 0),
                     breaks=seq(0,50, by=1)
  )+
  scale_y_continuous(expression(plain(paste(italic(T)[pk], " of Adult Mortality Rate (",italic(z),")"))),
                     limits =c(13,32),
                     expand = c(0, 0),
                     breaks=seq(0,50, by=2)
  )+
  theme(legend.position = 'none',legend.text = element_text(size = 10),
        axis.text.y = element_text(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = 23.5, y = 30.75,label = "C"),
            parse = TRUE, size = 6, colour = "black") +
  geom_line(aes(RearTemp,fit),Z,size=0.5,col="#636363")

# Rear_Z


##Juvi Mortality##
Zj <- subset(toptz, toptz$trait=="juvenile mortality rate")
Zj <- subset(toptz, toptz$trait=="fecundity")
Zj <- drop_na(Zj, "RearTemp")
ZjLM <- lm(Zj$estimate~ Zj$RearTemp)

ZjLMData <- predict(ZjLM, interval = "confidence")
Zj <- cbind(Zj,ZjLMData)


Rear_Zj <- ggplot(Zj, aes(x=RearTemp, y=estimate)) +
  geom_errorbar(aes(ymin =conf_lower , ymax = conf_upper),width = 0.2,col="#000000") +
  geom_point(size = 2.5,stroke=0.2)+
  theme_bw(base_size = 12.5) +
  # theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("Rearing Temperature (C",degree,")"))),
                     limits =c(23,28),
                     expand = c(0, 0),
                     breaks=seq(0,50, by=1)
  )+
  scale_y_continuous(expression(plain(paste(italic(T)[pk], " of Juvenile Mortality Rate (", italic(z[j]),")"))),
                     limits =c(20,32),
                     expand = c(0, 0),
                     breaks=seq(0,50, by=2)
  )+
  theme(legend.position = 'none',legend.text = element_text(size = 10),
        axis.text.y = element_text(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = 23.5, y = 31.25,label = "D"),
            parse = TRUE, size = 6, colour = "black") +
  geom_line(aes(RearTemp,fit),Zj,size=0.5,col="#636363")

# Rear_Zj


fig6 <- ggarrange(Rear_alp,Rear_Bopt,Rear_Z, Rear_Zj,nrow=2,ncol=2)
# fig6

ggsave("../results/RearTemp.pdf",fig6, width = 20, height =20, 
       units = "cm",device = cairo_pdf)

