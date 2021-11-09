### Extract rearing temperatures for each species, plot it
rm(list=ls())
library(tidyverse)

# setwd("~/Documents/TraitMismatchPaper-main/data/")

##Get publication list
VecTrait <- read.csv("old/Final_Traitofinterest_ph.csv")
# VecTrait <- read.csv("mismatch_dat.csv")
Refs <- unique(VecTrait$citation)
print(Refs)

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
  scale_y_continuous(expression(plain(paste("Peak Juvenile Development Rate (", italic(T)[opt]^alpha,")"))),
                     limits =c(23,38),
                     expand = c(0, 0),
                     breaks=seq(0,50, by=1)
                     )+
  theme(legend.position = 'none',legend.text = element_text(size = 10),
        axis.text.y = element_text(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = 27.5, y = 37,label = "A"), 
            parse = TRUE, size = 6, colour = "black") +
  geom_line(aes(RearTemp,fit),alp,size=0.5,col="#636363")
  # geom_abline(intercept = 0, slope=1, linetype="dashed")

Rear_alp


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
  scale_y_continuous(expression(plain(paste("Peak Fecundity (", italic(T)[opt]^B[pk],")"))),
                     limits =c(20,32),
                     expand = c(0, 0),
                     breaks=seq(0,50, by=1)
  )+
  theme(legend.position = 'none',legend.text = element_text(size = 10),
        axis.text.y = element_text(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = 27.5, y = 31.5,label = "B"),
  parse = TRUE, size = 6, colour = "black") +
  geom_line(aes(RearTemp,fit),Bopt,size=0.5,col="#636363")

Rear_Bopt

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
  scale_y_continuous(expression(plain(paste("Peak Adult Mortality Rate (", italic(T)[opt]^Z,")"))),
                     limits =c(13,32),
                     expand = c(0, 0),
                     breaks=seq(0,50, by=1)
  )+
  theme(legend.position = 'none',legend.text = element_text(size = 10),
        axis.text.y = element_text(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = 27.5, y = 31,label = "C"),
            parse = TRUE, size = 6, colour = "black") +
  geom_line(aes(RearTemp,fit),Z,size=0.5,col="#636363")

Rear_Z


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
  scale_y_continuous(expression(plain(paste("Peak Juvenile Mortality Rate (", italic(T)[opt]^Z[j],")"))),
                     limits =c(20,32),
                     expand = c(0, 0),
                     breaks=seq(0,50, by=1)
  )+
  theme(legend.position = 'none',legend.text = element_text(size = 10),
        axis.text.y = element_text(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  geom_text(aes(x = 27.5, y = 31.5,label = "D"),
            parse = TRUE, size = 6, colour = "black") +
  geom_line(aes(RearTemp,fit),Zj,size=0.5,col="#636363")

Rear_Zj


fig6 <- ggarrange(Rear_alp,Rear_Bopt,Rear_Z, Rear_Zj,nrow=2,ncol=2)
fig6

ggsave("../results/Fig6.pdf",fig6, width = 20, height =20, 
       units = "cm",device = cairo_pdf)



