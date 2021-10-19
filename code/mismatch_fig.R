#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#              Mismatch figures                    #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

require('tidyverse')
require('data.table')
require('rTPC')
require('car')
require('grDevices')
require('Cairo')
require('ggplot2')
require('patchwork')

rm(list=ls())
graphics.off()

# setwd("~/Dropbox/ph_thesis/Topt_paper/data")

#read in the Topt data

alpha  <- as_tibble(read.csv('alpha_Tpks.csv', header = TRUE))
alpha  <- alpha %>% filter(X <= 38)
zj     <- as_tibble(read.csv('zj_Tpks.csv', header = TRUE))
z      <- as_tibble(read.csv('z_Tpks.csv', header = TRUE))
bpk    <- as_tibble(read.csv('bpk_Tpks.csv', header = TRUE))

toptz  <- rbind(alpha,zj,z,bpk)
toptz <- subset(toptz, toptz$param=="topt")
  
toptz$species <- as.factor(toptz$species)
toptz$trait   <- as.factor(toptz$trait)

levels(toptz$species)

toptz$trait <- fct_relevel(toptz$trait, "juvenile mortality rate")

toptz <- toptz %>% filter(param!= "rmax")

# 
# toptz$species <- fct_relevel(toptz$species, "Aedes aegypti", after = Inf)
# toptz$species <- fct_relevel(toptz$species, "Helicoverpa armigera", after = Inf)
# toptz$species <- fct_relevel(toptz$species, "Anopheles gambiae s.s.", after = Inf)
# toptz$species <- fct_relevel(toptz$species, "Tetraneura nigriabdominalis") 
# toptz$species <- fct_relevel(toptz$species, "Rhopalosiphum maidis", after = 7)
# toptz$species <- fct_relevel(toptz$species, "Liposcelis bostrychophila", after = 9)
# toptz$species <- fct_relevel(toptz$species, "Aphis nasturtii", after = Inf)
# toptz$species <- fct_relevel(toptz$species, "Bemisia tabaci", after = 12)
# toptz$species <- fct_relevel(toptz$species, "Aphis gossypii", after = 9)
# toptz$species <- fct_relevel(toptz$species, "Aedes albopictus", after = Inf)
# toptz$species <- fct_relevel(toptz$species, "Aedes krombeini", after = 9)

#order by developement alpha
alp <- subset(toptz, toptz$trait == "juvenile development rate")
SPorder <- alp$species[order(alp$estimate)]
toptz$species <- as.character(toptz$species)
toptz$species <- factor(toptz$species, levels=SPorder)

#Oder by funnel shape --> max to min max diff
Diff <- toptz %>% group_by(species) %>%
        summarise(diff = max(estimate)-min(estimate))
SPorder <- Diff$species[order(Diff$diff)]
toptz$species <- as.character(toptz$species)
toptz$species <- factor(toptz$species, levels=SPorder)

# All traits

fig3a <- ggplot(toptz, aes(estimate, species, shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(xmin = conf_lower, xmax = conf_upper),width=0.35,size=0.35) +
  geom_point(size = 2.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("", italic(T)[opt]))),
                     limits =c(10,37),
                     expand = c(0, 0),
                     breaks=seq(12,36, by=4))+
  scale_fill_manual(values = c("#1f78b4","#a6cee3","#fdb863","#e66101"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow = 1,ncol =4,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
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
        axis.text.y = element_text(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  theme(legend.margin=margin(t = -0.4, unit='cm'))
  
fig3a
ggsave("../results/Fig4.pdf",fig3a, width = 20, height = 15, units = "cm")

#%%%%%%%%%%%%%%%%%%%%%%%% zj versus z 
zjz <- zj  %>% filter(trait == "juvenile mortality rate")
zjz <- zjz %>% filter(param == "topt")
zjz <- zjz %>% filter(species!= "Stethorus punctillum")
zjz <- zjz %>% filter(species!= "Aphis gossypii") 
zjz <- zjz %>% filter(species!= "Aedes krombeini")
zjz <- zjz %>% filter(species!= "Telenomus isis")
zjz <- zjz %>% filter(species!= "Liposcelis bostrychophila")

zjz <- zjz %>% rename(zj=estimate,zj_lwr=conf_lower,zj_upr=conf_upper)

zz <- z  %>% filter(trait == "adult mortality rate")
zz <- zz %>% filter(param == "topt")
zz <- zz %>% filter(species!= "Stethorus punctillum")
zz <- zz %>% filter(species!= "Aphis gossypii") 
zz <- zz %>% filter(species!= "Aedes krombeini")
zz <- zz %>% filter(species!= "Telenomus isis")

zz <- zz %>% rename(z=estimate,
                    z_lwr=conf_lower,
                    z_upr=conf_upper,specz=species,trtz=trait)

#%%%%%%%%%%%%%%%%%%%%%%%%

mortz <- cbind(zjz,zz)
mortz <- select(mortz,trait, zj,zj_lwr,zj_upr,species,trtz,z,z_lwr,z_upr)

mortz$zjminusz     <- mortz$zj-mortz$z
mortz$zjminusz_lwr <- mortz$zj_lwr-mortz$z_lwr
mortz$zjminusz_upr <- mortz$zj_upr-mortz$z_upr

mortz$species <- as.factor(mortz$species)

mortz$species <- fct_relevel(mortz$species, "Tetraneura nigriabdominalis", after = Inf) 
mortz$species <- fct_relevel(mortz$species, "Tetranychus mcdanieli", after = 3)
mortz$species <- fct_relevel(mortz$species, "Helicoverpa armigera", after = 2)
mortz$species <- fct_relevel(mortz$species, "Anopheles gambiae s.s.", after = 7)
mortz$species <- fct_relevel(mortz$species, "Paracoccus marginatu", after = 8)
mortz$species <- fct_relevel(mortz$species, "Anthonomus grandis", after = 8)
mortz$species <- fct_relevel(mortz$species, "Bemisia tabaci", after = 5)
mortz$species <- fct_relevel(mortz$species, "Aphis nasturtii", after = 4)


fig3b <- ggplot(mortz, aes(zjminusz, species, shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(xmin = zjminusz_lwr, xmax = zjminusz_upr),width=0.25,size=0.35,col="#000000") +
  geom_point(size = 2.5, col="#000000",stroke=0.2)+
  theme_bw(base_size = 12.5) +
  theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste(" Temperature Mismatch, ",degree,"C (",
                                italic(z[J]),", ",italic(T)[opt],
                                            " \u2212 ",italic(z),", ",italic(T)[opt],")"))),
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
  geom_text(aes(x = -12, y = 13,label = "B"), 
            parse = TRUE, size = 6, colour = "black")
  

p1 <- fig3a+fig3b

ggsave("~/Dropbox/ph_thesis/Topt_paper/results/Traitmismatch.pdf",p1, width = 30, height = 16.5, units = "cm",device = cairo_pdf)


#%%%%%%%%%%%%%%%%%%%%%%%% rm across temps

alpha <- as_tibble(read.csv('dvrate_topts.csv', header = TRUE))

trtT  <- rbind(alpha,zj,z)
tempz <- trtT %>% filter(param == "topt")
trtT  <- trtT %>% filter(param == "rmax")

trtT_zjz            <- trtT %>% filter(trait!="juvenile development rate")
trtT_zjz$conf_lower <- 1/trtT_zjz$conf_lower
trtT_zjz$estimate   <- 1/trtT_zjz$estimate
trtT_zjz$conf_upper <- 1/trtT_zjz$conf_upper

trtT <- trtT %>% filter(trait == 'juvenile development rate')
trtT <- rbind(trtT,trtT_zjz)
trtT$temp <- tempz$estimate
trtT$temp_lwr <- tempz$conf_lower
trtT$temp_upr <- tempz$conf_upper

trtT$trait <- fct_relevel(trtT$trait, "juvenile mortality rate")

traitBpk <- trtT %>% filter(trait!= "juvenile development rate") # subset for fig 4


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dv   <- trtT %>% filter(trait == "juvenile development rate")

moddv <- lm(log(estimate)~temp, data = dv)
anova(moddv)
summary(moddv)

preddv <- predict(moddv, interval = "confidence")

dv <- cbind(dv,preddv)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zj <- trtT %>% filter(trait == "juvenile mortality rate")

modzj <- lm(log(estimate)~temp, data = zj)
anova(modzj)
summary(modzj)

predzj <- predict(modzj, interval = "confidence")

zj <- cbind(zj,predzj)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z <- trtT %>% filter(trait == "adult mortality rate")

modz <- lm(log(estimate)~temp, data = z)
anova(modz)
summary(modz)

predz <- predict(modz, interval = "confidence")

z <- cbind(z,predz)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%

si_1 <- ggplot(traitBpk, aes(temp, log(estimate), shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(ymin = log(conf_lower), ymax = log(conf_upper)),width=0.15,size=0.15) +
  geom_errorbar(aes(xmin = temp_lwr, xmax = temp_upr),width=0.1,size=0.15) +
  geom_point(size = 1.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  scale_y_continuous(expression(plain(paste("ln (", italic(B)[pk],")"))),
                     limits =c(-8.4,-1),
                     expand = c(0, 0),
                     breaks=seq(-8,0, by=2))+
  scale_x_continuous(expression(plain(paste("", italic(T)[opt]))),
                     limits =c(10.5,33.5),
                     expand = c(0, 0),
                     breaks=seq(12,36, by=4))+
  scale_fill_manual(values = c("#1f78b4","#a6cee3"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow=1,ncol=2,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(values = c("#1f78b4","#a6cee3"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=1,ncol=2,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  scale_shape_manual(values = c(21,22),
                     name=expression(bold("")),
                     guide = guide_legend(nrow=1,ncol=2,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  theme(legend.position = c(0.5,0.89), 
        legend.background = element_rect(fill=alpha("#FFFFFF",1), colour = "#636363", size = 0.1),
        legend.text = element_text(size = 6))+
        theme(text=element_text(family="Times"))+
  geom_line(aes(temp, fit),z,size=0.35,col="#636363")+
  geom_ribbon(aes(ymin=lwr, ymax=upr),z,alpha = 0.25,show.legend = NA, col="#a6cee3",fill="#a6cee3",lwd=0.1)+
  geom_line(aes(temp, fit),zj,size=0.35,col="#636363")+
  geom_ribbon(aes(ymin=lwr, ymax=upr),zj,alpha = 0.25,show.legend = NA, col="#1f78b4",fill="#1f78b4",lwd=0.1)+
  theme(legend.margin=margin(t = -0.2, b = 0.1,r=0.1,l = 0.1, unit='cm'))+
  annotate("text", x = 24, y = -3,label = "slope = −0.04 ± 0.08 (95% CI)",
           alpha = 1, family="Times", size = 2)+
  annotate("text", x = 31, y = -3,label = "paste(italic(r) ^ 2, \" = 0.09\")", parse = TRUE,
           alpha = 1, family="Times", size = 2)+
  annotate("text", x = 16, y = -6,label = "slope = −0.09 ± 0.10 (95% CI)",
           alpha = 1, family="Times", size = 2)+
  annotate("text", x = 12.5, y = -6.55,label = "paste(italic(r) ^ 2, \" = 0.19\")", parse = TRUE,
           alpha = 1, family="Times", size = 2)+
  theme(legend.key.size = unit(0.1, 'cm'))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%


si_2 <- ggplot(dv, aes(temp, log(estimate), shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(ymin = log(conf_lower), ymax = log(conf_upper)),width=0.15,size=0.15) +
  geom_errorbar(aes(xmin = temp_lwr, xmax = temp_upr),width=0.05,size=0.15) +
  geom_point(size = 1.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  scale_y_continuous(expression(plain(paste("ln (", italic(B)[pk],")"))),
                     limits =c(-4,-1),
                     expand = c(0, 0),
                     breaks=seq(-4,0, by=1))+
  scale_x_continuous(expression(plain(paste("", italic(T)[opt]))),
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
  theme(legend.position = c(0.3,0.89), 
        legend.background = element_rect(fill=alpha("#FFFFFF",1), colour = "#636363", size = 0.1),
        legend.text = element_text(size = 6))+
  theme(text=element_text(family="Times"))+
  geom_line(aes(temp, fit),dv,size=0.35,col="#636363")+
  geom_ribbon(aes(ymin=lwr, ymax=upr),dv,alpha = 0.25,show.legend = NA, col="#e66101",fill="#e66101",lwd=0.1)+
  theme(legend.margin=margin(t = -0.2, b = 0.1,r=0.1,l = 0.1, unit='cm'))+
  annotate("text", x = 28, y = -3.85,label = "slope = 0.03 ± 0.11 (95% CI)",
           alpha = 1, family="Times", size = 2)+
  annotate("text", x = 33, y = -3.85,label = "paste(italic(r) ^ 2, \" = 0.01\")", parse = TRUE,
           alpha = 1, family="Times", size = 2)+
  theme(legend.key.size = unit(0.1, 'cm'))
  
si <- si_1 + si_2; si

ggsave("~/Dropbox/ph_thesis/Topt_paper/results/SI_weakhotbetter.pdf",si, width = 15, height =8, 
       units = "cm",device = cairo_pdf)



