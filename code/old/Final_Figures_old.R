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

setwd("")

###Note: working directory is data directory


##### Hotter is Better plots #########

alpha  <- as_tibble(read.csv('alpha_Tpks_AllParams.csv', header = TRUE))
zj     <- as_tibble(read.csv('zj_Tpks_AllParams.csv', header = TRUE))
z      <- as_tibble(read.csv('z_Tpks_AllParams.csv', header = TRUE))
bpk    <- as_tibble(read.csv('bpk_Tpks_AllParams.csv', header = TRUE))


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


##Plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig1c <- ggplot(dv, aes(temp, log(estimate), shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(ymin = log(conf_lower), ymax = log(conf_upper)),width=0.15,size=0.15) +
  geom_errorbar(aes(xmin = temp_lwr, xmax = temp_upr),width=0.05,size=0.15) +
  geom_point(size = 1.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  scale_y_continuous(expression(plain(paste("Log Juvenile Development Rate (",italic(alpha),")"))),
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
  theme(text=element_text(family="Times", size=17),
        axis.title = element_text(face="bold"))+
  theme(axis.text = element_text(size=18))+
  geom_line(aes(temp, fit),dv,size=0.35,col="#636363")+
  geom_ribbon(aes(ymin=lwr, ymax=upr),dv,alpha = 0.25,show.legend = NA, col="#e66101",fill="#e66101",lwd=0.1)+
  theme(legend.margin=margin(t = -0.2, b = 0.1,r=0.1,l = 0.1, unit='cm'))+
  theme(legend.position = "none")+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.5,vjust=1.4,label = "C"),
            size = 5, colour = "black")+
  theme(plot.margin=margin(t=0.5,unit="cm"))


# fig1c


##Fecundity
fig1d <- ggplot(bpkBpks, aes(temp, log(estimate), shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(ymin = log(conf_lower), ymax = log(conf_upper)),width=0.15,size=0.15) +
  geom_errorbar(aes(xmin = temp_lwr, xmax = temp_upr),width=0.05,size=0.15) +
  geom_point(size = 1.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  scale_y_continuous(expression(plain(paste("Log Fecundity (", italic(b[max]),")"))),
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
  theme(text=element_text(family="Times",size=17),
        axis.title = element_text(face="bold"))+
  theme(axis.text = element_text(size=18))+
  geom_line(aes(temp, fit),bpkBpks,size=0.35,col="#636363")+
  geom_ribbon(aes(ymin=lwr, ymax=upr),bpkBpks,alpha = 0.25,show.legend = NA, col="#e69201",fill="#e69201",lwd=0.1)+
  theme(legend.margin=margin(t = -0.2, b = 0.1,r=0.1,l = 0.1, unit='cm'))+
  theme(legend.position = "none")+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.5,vjust=1.4,
                label = "D"),size = 5, colour = "black")+
  theme(plot.margin=margin(t=0.5,unit="cm"))

# fig1d


##Adult Mortality
z <- rename(z, estimate=Bpk, conf_lower=Bpk_lwr, conf_upper = Bpk_upr)


fig1e <- ggplot(z, aes(temp, log(estimate), shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(ymin = log(conf_lower), ymax = log(conf_upper)),width=0.15,size=0.15) +
  geom_errorbar(aes(xmin = temp_lwr, xmax = temp_upr),width=0.05,size=0.15) +
  geom_point(size = 1.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  scale_y_continuous(expression(plain(paste("Log Adult Mortality Rate (", italic(z),")"))),
                     limits =c(-5,-2),
                     expand = c(0, 0),
                     breaks=seq(-10,0, by=1))+
  scale_x_continuous(expression(plain(paste("", italic(T)[pk]))),
                     limits =c(8,31),
                     expand = c(0, 0),
                     breaks=seq(0,36, by=4))+
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
  theme(text=element_text(family="Times",size=17),
        axis.title = element_text(face="bold"))+
  theme(axis.text = element_text(size=18))+
  geom_line(aes(temp, linearfit),z,size=0.35,col="#636363")+
  geom_ribbon(aes(ymin=linear_lwr, ymax=linear_upr),z,alpha = 0.5,show.legend = NA, col="#a6cee3",fill="#a6cee3",lwd=0.1)+
  theme(legend.margin=margin(t = -0.2, b = 0.1,r=0.1,l = 0.1, unit='cm'))+
  theme(legend.position = "none")+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.5,vjust=1.4,
                label = "E"),size = 5, colour = "black")
# fig1e

##Juvenile Mortality
zj <- rename(zj, estimate=Bpk, conf_lower=Bpk_lwr, conf_upper = Bpk_upr)


fig1f <- ggplot(zj, aes(temp, log(estimate), shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(ymin = log(conf_lower), ymax = log(conf_upper)),width=0.15,size=0.15) +
  geom_errorbar(aes(xmin = temp_lwr, xmax = temp_upr),width=0.05,size=0.15) +
  geom_point(size = 1.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  scale_y_continuous(expression(plain(paste("Log Juvenile Mortality Rate (", italic(z[J]),")"))),
                     limits =c(-9,-3),
                     expand = c(0, 0),
                     breaks=seq(-10,0, by=1))+
  scale_x_continuous(expression(plain(paste("", italic(T)[pk]))),
                     limits =c(10,35),
                     expand = c(0, 0),
                     breaks=seq(0,36, by=4))+
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
  theme(text=element_text(family="Times",size=19))+
  theme(axis.title = element_text(face="bold"))+
  theme(axis.text = element_text(size=18))+
  geom_line(aes(temp, linearfit),zj,size=0.35,col="#636363")+
  geom_ribbon(aes(ymin=linear_lwr, ymax=linear_upr),zj,alpha = 0.25,show.legend = NA, col="#1f78b4",fill="#1f78b4",lwd=0.1)+  theme(legend.margin=margin(t = -0.2, b = 0.1,r=0.1,l = 0.1, unit='cm'))+
  theme(legend.position = "none")+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.5,vjust=1.4,
                label = "F"),size = 5, colour = "black")


# fig1f

##plots
fig1 <- plot_grid(fig1c,NULL, fig1d,fig1e,NULL,fig1f,
                  align="hv",nrow=2,ncol=3,rel_widths = c(1,0.15,1),scale = 1.1)+
  theme(plot.margin=unit(c(1,1,1,1), "cm"))

fig1 <- ggarrange(fig1c,fig1d,fig1e,fig1f, nrow=2,ncol = 2, align="hv")+
  theme(plot.margin=unit(c(1,1,1,1), "cm"))

save_plot(fig1, file="../results/Fig1.pdf",
          base_height=20,base_asp=1.25, units="cm")


##### 2. Tpk + activation energy mistmatch plot #########

rm(list=ls())
graphics.off()

#Topt mismatch

alpha   <- as_tibble(read.csv('alpha_Tpks_AllParams.csv', header = TRUE))
zj      <- as_tibble(read.csv('zj_Tpks_AllParams.csv', header = TRUE))
z       <- as_tibble(read.csv('z_Tpks_AllParams.csv', header = TRUE))
bmax    <- as_tibble(read.csv('bmax_Tpks_AllParams.csv', header = TRUE))

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

#order by developement alpha
alp <- subset(toptz, toptz$trait=="juvenile development rate")
SPorder <- alp$species[order(alp$estimate)]
toptz$species <- factor(toptz$species, levels=SPorder)



# All traits

fig5a <- ggplot(toptz, aes(estimate, species, shape=trait, colour=trait,fill=trait)) +
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
  #theme(text=element_text(family="Times"))+
  theme(legend.margin=margin(t = -0.4, unit='cm'))+
  # coord_fixed(ratio = 1.5)
  theme(aspect.ratio=2)+
  ggtitle("A")+
  theme(plot.title=element_text(face="bold", size = 15, vjust = -1))

fig5a

###e mismatch ###

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


fig5b <- ggplot(activation_e, aes(estimate, species, shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(xmin = conf_lower, xmax = conf_upper),width=0.35,size=0.35) +
  geom_point(size = 2.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("Activation Energy (",italic(E),")"))),
                     limits =c(-10,38),
                     expand = c(0, 0),
                     breaks=seq(-20,36, by=4))+
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
  #theme(text=element_text(family="Times"))+
  theme(legend.margin=margin(t = -0.4, unit='cm'))+
  theme(axis.text.y=element_blank())+
  theme(plot.margin=margin(l=-4,unit="cm"))+
  theme(aspect.ratio=2)+
  ggtitle("B")+
  theme(plot.title=element_text(face="bold", size = 15, vjust = -1))

fig5b


plotMain <- fig5a+theme(legend.position="none")+fig5b+theme(legend.position="none")
# plotMain
legend <- get_legend(fig5a+theme(legend.position = "bottom"))


# plotMain <- plot_grid(plotMain, NULL, rel_widths = c(1,0))

fig5 <- plot_grid(plotMain, legend,ncol = 1, rel_heights = c(1, 0.05))

# 
save_plot(fig5, file="../results/Fig5.pdf",
          base_height=15,base_asp=1.5, units="cm")


##### 3. Pop growth rate, max pop growth rate with sum of traits #########

##reading data
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

##Fig 5a: Rm Curves##
fig6a <-  ggplot(rmpredictions, aes(temp,estimate,colour=factor(species, levels=SPorder)))+
  scale_x_continuous(expression(plain(paste(" Temperature (",degree,"C)"))))+
  scale_y_continuous(expression(plain(paste("Maximal population growth rate (",italic(r)[m], ")"))),
                     limits=c(-0.001,0.255),
                     expand = c(0.01, 0),
                     breaks=seq(0,0.25, by=0.05))+
  theme_bw(base_size = 12)+
  geom_line(size=0.85)+
  # geom_hline(aes(yintercept = 0), linetype = 2, show.legend = FALSE)+
  scale_colour_manual(values = colors,
                      name=expression(bold("species")),
                      labels =SPorder,
                      guide = guide_legend(nrow = 6,ncol =1 ,
                                           direction = "horizontal",
                                           title.position = "top",
                                           title.hjust=0.5))+
  theme(text=element_text(family="Times"))+
  theme(legend.text=element_text(family="Times",face = 'italic', size = 7), 
        legend.position = 'none')+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.5,vjust=1.4,
                label = "A"),size = 5, colour = "black")+
  theme(legend.title = element_blank())+
  theme(legend.margin=margin(t = -0.2, b = 0.1,r=0.1,l = 0.1, unit='cm'),
        aspect.ratio = 1)+
  geom_point(data= rmPeakdata, aes(x=temp, y=rmax), size=3)

fig6a

#####Sum of Topts #### 
SumPlot <- left_join(rmPeakdata, OptSum)
SumLM   <- lm(SumPlot$rmax~ SumPlot$sum)

SumLMData <- predict(SumLM, interval = "confidence")
SumPlot   <- cbind(SumPlot,SumLMData)


Sum_rmax<- ggplot(SumPlot, aes(x=sum, y=rmax, col=factor(species, levels=SPorder))) +
  geom_errorbar(aes(ymin = rmax_lwr, ymax = rmax_upr),col="#000000") +
  geom_point(size = 2.5, stroke=0.2)+
  theme_bw(base_size = 12.5) +
  # theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("Sum of ", italic(T[pk]),"'s"))),
                     limits =c(90,105),
                     expand = c(0, 0),
                     breaks=seq(20,120, by=5))+
  scale_y_continuous(expression(plain(paste(italic(r[m])," at Optimal Temperature (",italic(r[paste(m,",", opt)]), ")"))),
                     limits =c(0,0.31),
                     expand = c(0, 0),
                     breaks=seq(0,1, by=0.1))+
  scale_colour_manual(values = colors,
                      name=expression(bold("species")),
                      labels =SPorder)+
  theme(legend.position = 'none',legend.text = element_text(size = 10),
        aspect.ratio = 1)+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.5,vjust=1.4,
                label = "B"),size = 5, colour = "black")+
  geom_line(aes(sum,fit),SumPlot,size=0.35,col="#636363")

Sum_rmax

##save

fig6 <- ggarrange(fig6a,Sum_rmax,nrow=1,ncol=2,
                  common.legend = T, legend="right")
fig6
ggsave("../results/Fig6.pdf",fig6, width = 20, height =10, 
       units = "cm",device = cairo_pdf)

