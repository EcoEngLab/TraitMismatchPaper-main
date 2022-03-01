#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#       Mismatch manuscript figures 1c and 4       #
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

rm(list=ls())
graphics.off()

# setwd("~/Dropbox/ph_thesis/Topt_paper/data")


#read in the trait TPC data for Fig 4
alpha  <- as_tibble(read.csv('alpha_Tpks.csv', header = TRUE))
zj     <- as_tibble(read.csv('zj_Tpks.csv', header = TRUE))
z      <- as_tibble(read.csv('z_Tpks.csv', header = TRUE))
bpk    <- as_tibble(read.csv('bpk_Tpks.csv', header = TRUE))

toptz  <- rbind(alpha,zj,z,bpk)
  
toptz$species <- as.factor(toptz$species)
toptz$trait   <- as.factor(toptz$trait)

levels(toptz$species)

toptz$trait <- fct_relevel(toptz$trait, "juvenile mortality rate")

toptz <- toptz %>% filter(param== "topt")

# remove species with TPCs for only one trait
toptz <- toptz %>% filter(species!= "Trichogramma sp. nr. Lutea"
                    & species!= "Trichogramma bruni "
                    & species!= "Culex annulirostris"
                    & species!= "Macrocentrus iridescens")

# relevel so that zj goes from low to high Tpk
toptz$species <- fct_relevel(toptz$species, "Tetraneura nigriabdominalis") 
toptz$species <- fct_relevel(toptz$species, "Rhopalosiphum maidis", after = 7)
toptz$species <- fct_relevel(toptz$species, "Aphis nasturtii", after = Inf)
toptz$species <- fct_relevel(toptz$species, "Bemisia tabaci", after = 13)
toptz$species <- fct_relevel(toptz$species, "Aphis gossypii", after = 9)
toptz$species <- fct_relevel(toptz$species, "Aedes albopictus", after = Inf)
toptz$species <- fct_relevel(toptz$species, "Aedes krombeini", after = 9)
toptz$species <- fct_relevel(toptz$species, "Stethorus punctillum", after = Inf)
toptz$species <- fct_relevel(toptz$species, "Aedes aegypti", after = Inf)
toptz$species <- fct_relevel(toptz$species, "Anopheles gambiae s.s.", after = 10)

# All traits

fig4a <- ggplot(toptz, aes(estimate, species, shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(xmin = conf_lower, xmax = conf_upper),width=0.35,size=0.35) +
  geom_point(size = 2.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("", italic(T)[pk]))),
                     limits =c(10,37),
                     expand = c(0, 0),
                     breaks=seq(12,36, by=4))+
  scale_fill_manual(values = c("#1f78b4","#a6cee3","#fdb863","#e66101"),
                      name=expression(bold("")),
                    labels=expression(italic(z[J]),italic(z),italic(b[pk]),italic("\u03B1")),
                      guide = guide_legend(nrow = 1,ncol =4,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  scale_colour_manual(values = c("#1f78b4","#a6cee3","#fdb863","#e66101"),
                    name=expression(bold("")),
                    labels=expression(italic(z[J]),italic(z),italic(b[pk]),italic("\u03B1")),
                    guide = guide_legend(nrow=1,ncol=4,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_shape_manual(values = c(21,22,23,24),
                      name=expression(bold("")),
                     labels=expression(italic(z[J]),italic(z),italic(b[pk]),italic("\u03B1")),
                      guide = guide_legend(nrow = 1,ncol =4,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  theme(legend.position = c(0.5,-0.1),legend.text = element_text(size = 8.5),
        axis.text.y = element_text(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = 12, y = 15.5,label = "A"), 
            parse = TRUE, size = 6, colour = "black")+
  theme(legend.margin=margin(t = -0.5,b=-0.3, unit='cm'))
  

#ggsave("~/Dropbox/ph_thesis/Topt_paper/results/fig3a.pdf",fig4a, width = 20, height = 20, units = "cm")

#%%%%%%%%%%%%%%%%%%%%%%%% dataset for zj minus z (fig 4c(?))

zjz <- zj  %>% filter(trait == "juvenile mortality rate")
zjz <- zjz %>% filter(param == "topt")
zjz <- zjz %>% filter(species!= "Stethorus punctillum" # remove species that don't have zj or z
                      & species!= "Aphis gossypii" 
                      & species!= "Aedes krombeini"
                      & species!= "Telenomus isis")

zjz <- zjz %>% rename(zj=estimate,zj_lwr=conf_lower,zj_upr=conf_upper)

#%%

zz <- z  %>% filter(trait == "adult mortality rate")
zz <- zz %>% filter(param == "topt")
zz <- zz %>% filter(species!= "Stethorus punctillum"
                    & species!= "Aphis gossypii"
                    & species!= "Aedes krombeini"
                    & species!= "Telenomus isis")

zz <- zz %>% rename(z=estimate,
                    z_lwr=conf_lower,
                    z_upr=conf_upper,specz=species,trtz=trait)

#%%%%%%%%%%%%%%%%%%%%%%%%

# organise data for plotting 
mortz <- cbind(zjz,zz)
mortz <- select(mortz,trait, zj,zj_lwr,zj_upr,species,trtz,z,z_lwr,z_upr)

# subtract z from zj
mortz$zjminusz     <- mortz$zj-mortz$z
mortz$zjminusz_lwr <- mortz$zj_lwr-mortz$z_lwr
mortz$zjminusz_upr <- mortz$zj_upr-mortz$z_upr

mortz$species <- as.factor(mortz$species)

# order data 
mortz$species <- fct_relevel(mortz$species, "Rhopalosiphum maidis", after = 10)
mortz$species <- fct_relevel(mortz$species, "Muscidifurax zaraptor", after = 9)
mortz$species <- fct_relevel(mortz$species, "Anthonomus grandis", after = 8)
mortz$species <- fct_relevel(mortz$species, "Paracoccus marginatu", after = 7)
mortz$species <- fct_relevel(mortz$species, "Anopheles gambiae s.s.", after = 6)
mortz$species <- fct_relevel(mortz$species, "Bemisia tabaci", after = 5)
mortz$species <- fct_relevel(mortz$species, "Aphis nasturtii", after = 4)
mortz$species <- fct_relevel(mortz$species, "Tetranychus mcdanieli", after = 3)
mortz$species <- fct_relevel(mortz$species, "Tetraneura nigriabdominalis", after = Inf) 


fig4c <- ggplot(mortz, aes(zjminusz, species, shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(xmin = zjminusz_lwr, xmax = zjminusz_upr),width=0.25,size=0.35,col="#000000") +
  geom_point(size = 2.5, col="#000000",stroke=0.2)+
  theme_bw(base_size = 12.5) +
  theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste(" Temperature Mismatch, ",degree,"C (",italic("\u0394")~
                                              italic(T)[pk]^italic(z[J]~"\u2212"~italic(z))~")"))),
                                              limits =c(-16,18),
                                              expand = c(0, 0),
                                              breaks=seq(-12,16, by=4))+
  scale_fill_manual(values = c("#bababa"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow = 4,ncol =1 ,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5,
                                         label.position = "bottom"))+
  scale_colour_manual(values = c("#bababa"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow = 4,ncol =1 ,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5,
                                           label.position = "bottom"))+
  scale_shape_manual(values = c(21),
                     name=expression(bold("")),
                     guide = guide_legend(nrow = 4,ncol =1 ,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5,
                                          label.position = "bottom"))+
  theme(legend.position = 'none',legend.text = element_text(size = 10),
        axis.text.y = element_text(face = 'italic'))+
  geom_vline(aes(xintercept = 0), linetype = 2)+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = -12, y = 10.5,label = "C"), 
            parse = TRUE, size = 6, colour = "black")
  

#%%%%%%%%%%%%%%%%%%%%%%%% alpha minus zj (fig 4b)
azz <- alpha
azz <- azz  %>% filter(trait == "juvenile development rate")
azz <- azz  %>% filter(param == "topt")
azz <- azz  %>% filter(species!= "Stethorus punctillum"
                       & species!="Aphis gossypii" 
                       & species!= "Aedes krombeini" 
                       & species!= "Telenomus isis"
                       & species!= "Planococcus citri"
                       & species!= "Culex annulirostris"
                       & species!= "Trichogramma bruni "
                       & species!= "Macrocentrus iridescens"
                       & species!= "Trichogramma sp. nr. Lutea")

azz <- azz %>% rename(alpha=estimate,
                    alpha_lwr=conf_lower,
                    alpha_upr=conf_upper,specz=species,trtz=trait)

juvz <- cbind(zjz,azz)
juvz <- select(juvz,trait,zj,zj_lwr,zj_upr,species,trtz,alpha,alpha_lwr,alpha_upr)

juvz$alphaminuszj         <- juvz$alpha-juvz$zj
juvz$alphaminuszj_lwr     <- juvz$alpha_lwr-juvz$zj_lwr
juvz$alphaminuszj_upr     <- juvz$alpha_upr-juvz$zj_upr

juvz$species <- as.factor(juvz$species)

juvz <- juvz %>% arrange(juvz, alphaminuszj)

juvz$species <- fct_relevel(juvz$species, "Anthonomus grandis", after = 10)
juvz$species <- fct_relevel(juvz$species, "Tetranychus mcdanieli", after = 9)
juvz$species <- fct_relevel(juvz$species, "Aphis nasturtii", after = 8)
juvz$species <- fct_relevel(juvz$species, "Aedes albopictus", after = 7)
juvz$species <- fct_relevel(juvz$species, "Bemisia tabaci", after = 5)
juvz$species <- fct_relevel(juvz$species, "Paracoccus marginatu", after = 4)
juvz$species <- fct_relevel(juvz$species, "Anopheles gambiae s.s.", after = 3)
juvz$species <- fct_relevel(juvz$species, "Muscidifurax zaraptor", after = 2)
juvz$species <- fct_relevel(juvz$species, "Tetraneura nigriabdominalis", after = 1)
juvz$species <- fct_relevel(juvz$species, "Rhopalosiphum maidis") 
juvz$species <- fct_relevel(juvz$species, "Aedes aegypti", after = Inf) 


fig4b <- ggplot(juvz, aes(alphaminuszj, species, shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(xmin = alphaminuszj_lwr, xmax = alphaminuszj_upr),width=0.25,size=0.35,col="#000000") +
  geom_point(size = 2.5, col="#000000",stroke=0.2)+
  theme_bw(base_size = 12.5) +
  theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste(" Temperature Mismatch, ",degree,"C (",italic("\u0394")~
                                            italic(T)[pk]^italic("\u03B1 \u2212"~italic(z[J]))~")"))),
                     limits =c(0,25),
                     expand = c(0, 0),
                     breaks=seq(0,25, by=3))+
  scale_fill_manual(values = c("#bababa"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow = 4,ncol =1 ,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5,
                                         label.position = "bottom"))+
  scale_colour_manual(values = c("#bababa"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow = 4,ncol =1 ,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5,
                                           label.position = "bottom"))+
  scale_shape_manual(values = c(21),
                     name=expression(bold("")),
                     guide = guide_legend(nrow = 4,ncol =1 ,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5,
                                          label.position = "bottom"))+
  theme(legend.position = 'none',legend.text = element_text(size = 10),
        axis.text.y = element_text(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = 3, y = 10.5,label = "B"), 
            parse = TRUE, size = 6, colour = "black")

p1 <-  fig4a | fig4b / fig4c; p1
 
ggsave("../results/Fig4.pdf",
       p1, width = 24, height = 17, units = "cm",device = cairo_pdf)


#%%%%%%%%%%%%%%%%%%%%%%%% 

# regression data for fig 1c, d and e(?)

alpha  <- as_tibble(read.csv('alpha_Tpks.csv', header = TRUE))
zj     <- as_tibble(read.csv('zj_Tpks.csv', header = TRUE))
z      <- as_tibble(read.csv('z_Tpks.csv', header = TRUE))
bpk    <- as_tibble(read.csv('bpk_Tpks.csv', header = TRUE))


alphaTopts         <- alpha %>% filter(param == "topt")
alphaBpks          <- alpha %>% filter(param == "rmax")
alphaBpks$temp     <- alphaTopts$estimate # for predictor in linear regression model 
alphaBpks$temp_lwr <- alphaTopts$conf_lower
alphaBpks$temp_upr <- alphaTopts$conf_upper

bpkTopts         <- bpk %>% filter(param == "topt")
bpkBpks          <- bpk %>% filter(param == "rmax")
bpkBpks$temp     <- bpkTopts$estimate 
bpkBpks$temp_lwr <- bpkTopts$conf_lower
bpkBpks$temp_upr <- bpkTopts$conf_upper


zjz <- rbind(zj,z)
zjzTopts <- filter(zjz, param == "topt")
zjzBpks <-  filter(zjz, param == "rmax")

zjzBpks$conf_lower <- 1/zjzBpks$conf_lower # invert for rates
zjzBpks$estimate   <- 1/zjzBpks$estimate
zjzBpks$conf_upper <- 1/zjzBpks$conf_upper
zjzBpks$temp       <- zjzTopts$estimate    # for predictor in linear regression model 
zjzBpks$temp_lwr   <- zjzTopts$conf_lower
zjzBpks$temp_upr   <- zjzTopts$conf_upper

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# hotter is better linear regression models

moddv <- lm(log(estimate)~temp, data = alphaBpks)
anova(moddv)
summary(moddv)

preddv <- predict(moddv, interval = "confidence")

dv <- cbind(alphaBpks,preddv)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zj <- zjzBpks %>% filter(trait == "juvenile mortality rate")

modzj <- lm(log(estimate)~temp, data = zj)
anova(modzj)
summary(modzj)

predzj <- predict(modzj, interval = "confidence")

zj <- cbind(zj,predzj)

zj <- zj %>% select(species,conf_lower,estimate,conf_upper,trait,temp_lwr,temp,temp_upr,lwr,fit,upr)
zj <- zj %>% rename(Bpk_lwr=conf_lower,Bpk=estimate,Bpk_upr=conf_upper,linear_lwr=lwr,linearfit=fit,linear_upr=upr)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z <- zjzBpks %>% filter(trait == "adult mortality rate")

modz <- lm(log(estimate)~temp, data = z)
anova(modz)
summary(modz)

predz <- predict(modz, interval = "confidence")

z <- cbind(z,predz)

z <- z %>% select(species,conf_lower,estimate,conf_upper,trait,temp_lwr,temp,temp_upr,lwr,fit,upr)
z <- z %>% rename(Bpk_lwr=conf_lower,Bpk=estimate,Bpk_upr=conf_upper,linear_lwr=lwr,linearfit=fit,linear_upr=upr)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modbpk <- lm(log(estimate)~temp, data = bpkBpks)
anova(modbpk)
summary(modbpk)

predbpk <- predict(modbpk, interval = "confidence")

bpkBpks <- cbind(bpkBpks,predbpk)

# The bpk model might be improved with additional data. I have data that was 
# not included because it was not suitable for the rmax model i.e. it was not possible  
# to standardise these data to a rate because adult lifespan was not included.
# however it is not necessary to standardise these data if they're only being used 
# for the regression analysis

# The bpk model has not been added to regression fig below


##Plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig1c <- ggplot(dv, aes(temp, log(estimate), shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(ymin = log(conf_lower), ymax = log(conf_upper)),width=0.15,size=0.15) +
  geom_errorbar(aes(xmin = temp_lwr, xmax = temp_upr),width=0.05,size=0.15) +
  geom_point(size = 1.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  scale_y_continuous(expression(plain(paste("Juvenile Development Rate (", alpha,")"))),
                     limits =c(-4,-1),
                     expand = c(0, 0),
                     breaks=seq(-4,0, by=1))+
  scale_x_continuous(expression(plain(paste("", italic(T)[pk]))),
                     limits =c(23,37),
                     expand = c(0, 0),
                     breaks=seq(24,36, by=2))+
  scale_fill_manual(values = c("#e66101"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow=1,ncol=1,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(values = c("#e66101"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=1,ncol=1,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  scale_shape_manual(values = c(24),
                     name=expression(bold("")),
                     guide = guide_legend(nrow=1,ncol =1,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  # theme(legend.position = c(0.3,0.86), 
  #       legend.background = element_rect(fill=alpha("#FFFFFF",1), colour = "#636363", size = 0.1),
  #       legend.text = element_text(size = 6))+
  theme(text=element_text(family="Times",size=16),
        axis.title = element_text(face="bold"))+
  geom_line(aes(temp, fit),dv,size=0.35,col="#636363")+
  geom_ribbon(aes(ymin=lwr, ymax=upr),dv,alpha = 0.25,show.legend = NA, col="#e66101",fill="#e66101",lwd=0.1)+
  theme(legend.margin=margin(t = -0.2, b = 0.1,r=0.1,l = 0.1, unit='cm'))+
  # annotate("text", x = 28, y = -3.85,label = "slope = 0.12 ± 0.10 (95% CI)",
  #          alpha = 1, family="Times", size = 2)+
  # annotate("text", x = 33, y = -3.85,label = "paste(italic(R) ^ 2, \" = 0.25\")", parse = TRUE,
  #          alpha = 1, family="Times", size = 2)+
  # theme(legend.key.size = unit(0.1, 'cm'))+
  theme(legend.position = "none")+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.5,vjust=1.4,label = "C"),
            size = 5, colour = "black")+
  theme(plot.margin=margin(t=0.5,unit="cm"))


fig1c


##Fecundity
fig1d <- ggplot(bpkBpks, aes(temp, log(estimate), shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(ymin = log(conf_lower), ymax = log(conf_upper)),width=0.15,size=0.15) +
  geom_errorbar(aes(xmin = temp_lwr, xmax = temp_upr),width=0.05,size=0.15) +
  geom_point(size = 1.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  scale_y_continuous(expression(plain(paste("Fecundity (", beta,")"))),
                     limits =c(0,4),
                     expand = c(0, 0),
                     breaks=seq(-4,4, by=1))+
  scale_x_continuous(expression(plain(paste("", italic(T)[pk]))),
                     limits =c(22,31),
                     expand = c(0, 0),
                     breaks=seq(22,36, by=2))+
  scale_fill_manual(values = c("#fdb863"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow=1,ncol=1,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(values = c("#fdb863"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=1,ncol=1,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  scale_shape_manual(values = c(23),
                     name=expression(bold("")),
                     guide = guide_legend(nrow=1,ncol =1,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  # theme(legend.position = c(0.2,0.9), 
  #       legend.background = element_rect(fill=alpha("#FFFFFF",1), colour = "#636363", size = 0.1),
  #       legend.text = element_text(size = 10))+
  theme(text=element_text(family="Times",size=16),
        axis.title = element_text(face="bold"))+
  geom_line(aes(temp, fit),bpkBpks,size=0.35,col="#636363")+
  geom_ribbon(aes(ymin=lwr, ymax=upr),bpkBpks,alpha = 0.25,show.legend = NA, col="#e69201",fill="#e69201",lwd=0.1)+
  theme(legend.margin=margin(t = -0.2, b = 0.1,r=0.1,l = 0.1, unit='cm'))+
  #fecundity:
  # annotate("text", x = 24, y = 0.5,label = "slope = -0.08 ± 0.23 (95% CI)",
  #          alpha = 1, family="Times", size = 2)+
  # annotate("text", x = 28, y = 0.5,label = "paste(italic(R) ^ 2, \" = 0.08\")", parse = TRUE,
  #          alpha = 1, family="Times", size = 2)+
  theme(legend.position = "none")+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.5,vjust=1.4,
                label = "D"),size = 5, colour = "black")+
  theme(plot.margin=margin(t=0.5,unit="cm"))

fig1d


##Adult Mortality
z <- rename(z, estimate=Bpk, conf_lower=Bpk_lwr, conf_upper = Bpk_upr)


fig1e <- ggplot(z, aes(temp, log(estimate), shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(ymin = log(conf_lower), ymax = log(conf_upper)),width=0.15,size=0.15) +
  geom_errorbar(aes(xmin = temp_lwr, xmax = temp_upr),width=0.05,size=0.15) +
  geom_point(size = 1.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  scale_y_continuous(expression(plain(paste("Adult Mortality Rate (", italic(z),")"))),
                     limits =c(-5,-2),
                     expand = c(0, 0),
                     breaks=seq(-10,0, by=1))+
  scale_x_continuous(expression(plain(paste("", italic(T)[pk]))),
                     limits =c(10,31),
                     expand = c(0, 0),
                     breaks=seq(0,36, by=2))+
  scale_fill_manual(values = c("#a6cee3"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow=1,ncol=1,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(values = c("#a6cee3"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=1,ncol=1,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  scale_shape_manual(values = c(22),
                     name=expression(bold("")),
                     guide = guide_legend(nrow=1,ncol =1,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  # theme(legend.position = c(0.3,0.86), 
  #       legend.background = element_rect(fill=alpha("#FFFFFF",1), colour = "#636363", size = 0.1),
  #       legend.text = element_text(size = 6))+
  theme(text=element_text(family="Times",size=16),
        axis.title = element_text(face="bold"))+
  geom_line(aes(temp, linearfit),z,size=0.35,col="#636363")+
  geom_ribbon(aes(ymin=linear_lwr, ymax=linear_upr),z,alpha = 0.5,show.legend = NA, col="#a6cee3",fill="#a6cee3",lwd=0.1)+
  theme(legend.margin=margin(t = -0.2, b = 0.1,r=0.1,l = 0.1, unit='cm'))+
  # annotate("text", x = 28, y = -3.85,label = "slope = 0.12 ± 0.10 (95% CI)",
  #          alpha = 1, family="Times", size = 2)+
  # annotate("text", x = 33, y = -3.85,label = "paste(italic(R) ^ 2, \" = 0.25\")", parse = TRUE,
  #          alpha = 1, family="Times", size = 2)+
  # theme(legend.key.size = unit(0.1, 'cm'))+
  theme(legend.position = "none")+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.5,vjust=1.4,
                label = "E"),size = 5, colour = "black")
  # theme(plot.margin=margin(l=3,t=0.5,unit="cm"))
fig1e

##Juvenile Mortality
zj <- rename(zj, estimate=Bpk, conf_lower=Bpk_lwr, conf_upper = Bpk_upr)


fig1f <- ggplot(zj, aes(temp, log(estimate), shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(ymin = log(conf_lower), ymax = log(conf_upper)),width=0.15,size=0.15) +
  geom_errorbar(aes(xmin = temp_lwr, xmax = temp_upr),width=0.05,size=0.15) +
  geom_point(size = 1.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  scale_y_continuous(expression(plain(paste("Juvenile Mortality Rate (", italic(z)[j],")"))),
                     limits =c(-9,-3),
                     expand = c(0, 0),
                     breaks=seq(-10,0, by=1))+
  scale_x_continuous(expression(plain(paste("", italic(T)[pk]))),
                     limits =c(10,35),
                     expand = c(0, 0),
                     breaks=seq(0,36, by=2))+
  scale_fill_manual(values = c("#1f78b4"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow=1,ncol=1,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(values = c("#1f78b4"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=1,ncol=1,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  scale_shape_manual(values = c(21),
                     name=expression(bold("")),
                     guide = guide_legend(nrow=1,ncol =1,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  # theme(legend.position = c(0.3,0.86), 
  #       legend.background = element_rect(fill=alpha("#FFFFFF",1), colour = "#636363", size = 0.1),
  #       legend.text = element_text(size = 6))+
  theme(text=element_text(family="Times",size=16),
        axis.title = element_text(face="bold"))+
  geom_line(aes(temp, linearfit),zj,size=0.35,col="#636363")+
  geom_ribbon(aes(ymin=linear_lwr, ymax=linear_upr),zj,alpha = 0.25,show.legend = NA, col="#1f78b4",fill="#1f78b4",lwd=0.1)+  theme(legend.margin=margin(t = -0.2, b = 0.1,r=0.1,l = 0.1, unit='cm'))+
  # annotate("text", x = 28, y = -3.85,label = "slope = 0.12 ± 0.10 (95% CI)",
  #          alpha = 1, family="Times", size = 2)+
  # annotate("text", x = 33, y = -3.85,label = "paste(italic(R) ^ 2, \" = 0.25\")", parse = TRUE,
  #          alpha = 1, family="Times", size = 2)+
  # theme(legend.key.size = unit(0.1, 'cm'))+
  theme(legend.position = "none")+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.5,vjust=1.4,
                label = "F"),size = 5, colour = "black")
  # theme(plot.margin=margin(l=-4,t=0.5,unit="cm"))
  




# fig1f

##plots
fig1 <- plot_grid(fig1c, fig1d,fig1e,fig1f,align="hv",nrow=2,ncol=2)
#load legend from mismatch plot
# load("../results/legend.rda")

# fig1 <- plot_grid(plotMain, legend, ncol = 1, rel_heights = c(1, 0.1))

save_plot(fig1, file="../results/Fig1.pdf",
          base_height=20,base_asp=1.25, units="cm")


