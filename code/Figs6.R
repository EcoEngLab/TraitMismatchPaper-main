## Figure 6 a b: R max with mistmatch
rm(list=ls())

library(ggplot2)
library(ggpubr)
library(tidyverse)

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
  scale_y_continuous(expression(plain(paste(" Population growth rate ("~italic(r)[m]~")"))),
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
        legend.position = c(-0.4,0.68))+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.5,vjust=1.4,
                label = "A"),size = 5, colour = "black")+
  theme(legend.title = element_blank())+
  theme(legend.margin=margin(t = -0.2, b = 0.1,r=0.1,l = 0.1, unit='cm'),
        aspect.ratio = 1)+
  geom_point(data= rmPeakdata, aes(x=temp, y=rmax), size=3)

fig6a

#####Sum of Topts #### 
SumPlot <- left_join(rmPeakdata, OptSum)
SumLM <- lm(SumPlot$rmax~ SumPlot$sum)

SumLMData <- predict(SumLM, interval = "confidence")
SumPlot <- cbind(SumPlot,SumLMData)


Sum_rmax<- ggplot(SumPlot, aes(x=sum, y=rmax, col=factor(species, levels=SPorder))) +
  geom_errorbar(aes(ymin = rmax_lwr, ymax = rmax_upr),col="#000000") +
  geom_point(size = 2.5, stroke=0.2)+
  theme_bw(base_size = 12.5) +
  # theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("Sum of ", italic(T)[pk]))),
                     limits =c(90,105),
                     expand = c(0, 0),
                     breaks=seq(20,120, by=5))+
  scale_y_continuous(expression(plain(paste("Peak Population growth rate ("~italic(r)[m]~")"))),
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
  geom_line(aes(sum,fit),SumPlot,size=0.35,col="#636363")
Sum_rmax

##save

fig6 <- ggarrange(fig6a,Sum_rmax,nrow=1,ncol=2,
                  common.legend = T, legend="right")
fig6
ggsave("../results/Fig6.pdf",fig6, width = 20, height =10, 
       units = "cm",device = cairo_pdf)


######Plots with fitness#######

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
  scale_x_continuous(expression(plain(paste("Peak Juvenile Development Rate (", italic(T)[pk]^alpha,")"))),
                     limits =c(26,35),
                     expand = c(0, 0),
                     breaks=seq(26,35, by=1))+
  scale_y_continuous(expression(plain(paste("Peak Population growth rate ("~italic(r)[m]~")"))),
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
  theme(axis.title.x = element_text(size=10))

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
  scale_y_continuous(expression(plain(paste("Peak Population growth rate ("~italic(r)[m]~")"))),
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
  geom_line(aes(variance,fit),VarPlot,size=0.35,col="#636363")
Var_rmax

si1 <- ggarrange(Alp_rmax,Var_rmax,nrow=1,ncol=2,
                  common.legend = T, legend="right")

ggsave("../results/SI1.pdf",si1, width = 20, height =10, 
       units = "cm",device = cairo_pdf)
