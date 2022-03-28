########Figures for Supplementary Information ##########


##Packages##
library(ggplot2)
library(tidyverse)


##########1. Varaince + Development rate against fitness #######
rm(list=ls())
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
##Write Rear DF
RearDF <- data.frame(matrix(ncol=4,nrow=0))

#Entering rearing temperatures manually
RearDF<- rbind(RearDF, c("Anopheles gambiae s.s.", 26,25,27))
RearDF<- rbind(RearDF, c("Telenomus isis", 25,24,26))
RearDF<- rbind(RearDF, c("Triticum aestivum L.", 30,NA,NA))
RearDF<- rbind(RearDF, c("Aedes aegypti", 25,NA,NA))
RearDF<- rbind(RearDF, c("Aedes albopictus", 27.5,26.5,28.5))
RearDF<- rbind(RearDF, c("Anthonomus grandis", 27,26,28))
RearDF<- rbind(RearDF, c("Aphis gossypii", 26,24,28))
RearDF<- rbind(RearDF, c("Bactrocera  correcta", 27.5,25,30))
RearDF<- rbind(RearDF, c("Corythucha ciliata", 26,25.5,26.5))
RearDF<- rbind(RearDF, c("Culex annulirostris", 27,26,28))
RearDF<- rbind(RearDF, c("Pseudaletia sequax", 21,20,22))
RearDF<- rbind(RearDF, c("Lepinotus reticulatus", 30,NA,NA))
RearDF<- rbind(RearDF, c("Macrocentrus iridescens", 25,NA,NA))
RearDF<- rbind(RearDF, c("Moina macrocopa", 25,NA,NA))
RearDF<- rbind(RearDF, c("Planococcus citri", 24,22,26))
RearDF<- rbind(RearDF, c("Sitona discoideus", 24,22,26))
RearDF<- rbind(RearDF, c("Stethorus punctillum", 24,NA,NA))
RearDF<- rbind(RearDF, c("Telenomus chrysopae", 24,23,25))
RearDF<- rbind(RearDF, c("Telenomus lobatus", 24,23,25))
RearDF<- rbind(RearDF, c("Tetraneura nigriabdominalis ", 25,NA,NA))
RearDF<- rbind(RearDF, c("Theocolax elegans",25,NA,NA))
RearDF<- rbind(RearDF, c("Helicoverpa armigera",25,NA,NA))
RearDF<- rbind(RearDF, c("Anopheles gambiae",26,25,27))
RearDF<- rbind(RearDF, c("Aphis nasturtii", 24, 23, 25))
RearDF<- rbind(RearDF, c("Muscidifurax zaraptor", 25, NA, NA))
RearDF<- rbind(RearDF, c("Paracoccus marginatus", 25, 24, 26))
RearDF<- rbind(RearDF, c("Rhopalosiphum maidis ", 25, NA, NA))
RearDF<- rbind(RearDF, c("Tetranychus mcdanieli", 24, NA, NA))
RearDF<- rbind(RearDF, c("Trichogramma bruni", 16.6,10.1,23))

names(RearDF) <- c("species","RearTemp","RearTemp_upper", "RearTemp_lower")

write.csv(RearDF, "RearTemp.csv")


##### plots? ####
alpha  <- as_tibble(read.csv('alpha_Tpks.csv', header = TRUE))
zj     <- as_tibble(read.csv('zj_Tpks.csv', header = TRUE))
z      <- as_tibble(read.csv('z_Tpks.csv', header = TRUE))
bpk    <- as_tibble(read.csv('bpk_Tpks.csv', header = TRUE))

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
  scale_y_continuous(expression(plain(paste("Peak Juvenile Development Rate (", italic(T)[pk]^alpha,")"))),
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
  scale_y_continuous(expression(plain(paste("Peak Fecundity (", italic(T)[pk]^B[pk],")"))),
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
  scale_y_continuous(expression(plain(paste("Peak Adult Mortality Rate (", italic(T)[pk]^Z,")"))),
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
  scale_y_continuous(expression(plain(paste("Peak Juvenile Mortality Rate (", italic(T)[pk]^Z[j],")"))),
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










