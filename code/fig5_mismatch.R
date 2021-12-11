#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#              Mismatch figures                    #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

require('tidyverse')
require('data.table')
require('car')
require('grDevices')
require('Cairo')
require('ggplot2')
require('ggpubr')
require('patchwork')
require('cowplot')


rm(list=ls())
graphics.off()

#Topt mismatch

alpha  <- as_tibble(read.csv('alpha_Tpks_AllParams.csv', header = TRUE))
# alpha  <- alpha %>% filter(X <= 38)
zj     <- as_tibble(read.csv('zj_Tpks_AllParams.csv', header = TRUE))
z      <- as_tibble(read.csv('z_Tpks_AllParams.csv', header = TRUE))
bpk    <- as_tibble(read.csv('bpk_Tpks_AllParams.csv', header = TRUE))

topt  <- rbind(alpha,zj,z,bpk)
toptz <- subset(topt, topt$param=="topt")
# toptz <- rename(toptz, trait=Trait)  

toptz$species <- as.factor(toptz$species)
toptz$trait   <- as.factor(toptz$trait)

levels(toptz$species)

toptz$trait <- fct_relevel(toptz$trait, "juvenile mortality rate")

toptz <- toptz %>% filter(param!= "rmax")

SpCount <- table(toptz$species)
RmSp <- names(SpCount)[which(SpCount==1)]

toptz <- filter(toptz, !(species %in% RmSp))


#order by developement alpha
alp <- subset(toptz, toptz$trait=="juvenile development rate")
SPorder <- alp$species[order(alp$estimate)]
toptz$species <- as.character(toptz$species)
toptz$species <- factor(toptz$species, levels=SPorder)


# All traits

fig5a <- ggplot(toptz, aes(estimate, species, shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(xmin = conf_lower, xmax = conf_upper),width=0.35,size=0.35) +
  geom_point(size = 2.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("", italic(T)[pk]))),
                     limits =c(9,37),
                     expand = c(0, 0),
                     breaks=seq(12,36, by=4))+
  scale_fill_manual(labels = c(expression(plain(paste("Juvenile Mortality Rate (",z[j],")"))),
                               expression(plain(paste("Adult Mortality Rate (",z,")"))),
                               expression(plain(paste("Fecundity (",beta,")"))),
                               expression(plain(paste("Juvenile Development Rate (",alpha,")")))),
                    values = c("#1f78b4","#a6cee3","#fdb863","#e66101"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow = 1,ncol =4,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(labels = c(expression(plain(paste("Juvenile Mortality Rate (",z[j],")"))),
                                 expression(plain(paste("Adult Mortality Rate (",z,")"))),
                                 expression(plain(paste("Fecundity (",beta,")"))),
                                 expression(plain(paste("Juvenile Development Rate (",alpha,")")))),
                      values = c("#1f78b4","#a6cee3","#fdb863","#e66101"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=1,ncol=4,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  scale_shape_manual(labels = c(expression(plain(paste("Juvenile Mortality Rate (",z[j],")"))),
                                expression(plain(paste("Adult Mortality Rate (",z,")"))),
                                expression(plain(paste("Fecundity (",beta,")"))),
                                expression(plain(paste("Juvenile Development Rate (",alpha,")")))),
                     values = c(21,22,23,24),
                     name=expression(bold("")),
                     guide = guide_legend(nrow = 1,ncol =4,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  theme(legend.position = c(0.35,-0.1),legend.text = element_text(size = 10),
        axis.text.y = element_text(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  theme(legend.margin=margin(t = -0.4, unit='cm'))+
  coord_fixed(ratio = 1.5)

# fig5a


###e mismatch ###
topte <- filter(topt, param =="e" )
topte <- mutate(topte, trait = case_when(trait =="adult mortality rate" ~ "Adult Mortality Rate",
                                         trait == "fecundity" ~ "Peak Fecundity",
                                         trait == "juvenile development rate" ~ "Development Rate",
                                         trait == "juvenile mortality rate" ~ "Juvenile Mortality Rate"))

TraitLevel <- c("Juvenile Mortality Rate", "Adult Mortality Rate","Peak Fecundity","Development Rate" )
topte$trait <- as.character(topte$trait)
topte$trait <- factor(topte$trait, levels = TraitLevel)
topte$species <- as.character(topte$species)
topte$species <- factor(topte$species, levels=SPorder)
topte <- na.omit(topte)


##remove species with 1 trait
SpCount <- table(topte$species)
# RmSp <- names(SpCount)[which(SpCount==1)]
# topte <- filter(topte, !(species %in% RmSp))


fig5b <- ggplot(topte, aes(estimate, species, shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(xmin = conf_lower, xmax = conf_upper),width=0.35,size=0.35) +
  geom_point(size = 2.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("Activation Energy (",e,")"))),
                     limits =c(0,11),
                     expand = c(0, 0),
                     breaks=seq(-20,36, by=2))+
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
        axis.text.y = element_text(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  theme(legend.margin=margin(t = -0.4, unit='cm'))+
  theme(axis.text.y=element_blank())+
  theme(plot.margin=margin(l=-4,unit="cm"))+
  coord_fixed(ratio = 1.5)

# fig5b

plotMain <- plot_grid(fig5a+theme(legend.position="none")
                  ,fig5b+theme(legend.position="none")
                  , ncol=2, align="h")
legend <- get_legend(fig5a+theme(legend.position = "bottom"))

#custom legend for fig1
# legend <- get_legend(fig5a+theme(legend.position = "bottom")
#                      +guides(fill=guide_legend(nrow=2,byrow=TRUE),
#                              colour=guide_legend(nrow=2,byrow=TRUE),
#                              shape=guide_legend(nrow=2,byrow=TRUE)))
# save(legend, file="../results/legend.rda")

fig5 <- plot_grid(plotMain, legend, ncol = 1, rel_heights = c(1, 0.1))

# 
save_plot(fig5, file="../results/Fig5.pdf",
          base_height=15,base_asp=1.5, units="cm")


