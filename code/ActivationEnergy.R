###Get histograms of all paramters

rm(list=ls())
graphics.off()

##Packages##
library(ggplot2)
library(tidyverse)

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
  labs(x="Activation Energy")+
  geom_histogram(position="identity", binwidth=1)+
  geom_density()+
  facet_grid(.~trait)
eplot

ggsave("../results/eHist.pdf",eplot, width = 20, height = 10, units = "cm")

###mistmatch fig

TraitLevel <- c("Juvenile Mortality Rate", "Adult Mortality Rate","Peak Fecundity","Development Rate" )
topte$trait <- as.character(topte$trait)
topte$trait <- factor(topte$trait, levels = TraitLevel)

##remove species with 1 trait
SpCount <- table(topte$species)
RmSp <- names(SpCount)[which(SpCount==1)]
toptz <- filter(toptz, !(species %in% RmSp))


fig3a <- ggplot(topte, aes(estimate, species, shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(xmin = conf_lower, xmax = conf_upper),width=0.35,size=0.35) +
  geom_point(size = 2.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain("Activation Energy")),
                     limits =c(0,8),
                     expand = c(0, 0),
                     breaks=seq(-20,36, by=4))+
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
ggsave("../results/eMismatch.pdf",fig3a, width = 20, height = 14, units = "cm")


###Get activation energy old fashioned way using LM##
df <- read.csv("mismatch_dat.csv") %>%
  rename(temp = ambienttemp,species = interactor1, trait=standardisedtraitname) %>%
  select(curve_ID, species, temp, trait, standardisedtraitvalue, stage)%>%
  rename(rate = standardisedtraitvalue)%>%
  mutate(rate= ifelse(rate==0, 10^-6, rate))%>%
  mutate(trait=case_when(trait=="Female Mortality Rate"~ "adult mortality rate",
                         trait=="Juvenile Development Rate"~"juvenile development rate",
                         trait=="Juvenile Mortality Rate" ~"juvenile mortality rate",
                         trait=="Fecundity Rate"~"fecundity"))

##get index
df <- df %>%
  group_by(species,trait) %>%
  mutate(curve_ID = cur_group_id()) %>%
  ungroup() %>%
  arrange(curve_ID)

#all unique curves
# Curves <- unique(df$curve_ID)

GetActivation <- function(df){
  ##get topt
  # browser()
  toptDF <- filter(toptz, toptz$trait==df$trait[1]&toptz$species==df$species[1])
  topt <- toptDF$estimate[which(toptDF$param=="topt")]
  e <- toptDF$estimate[which(toptDF$param=="e")]
  
  if(any(df$trait[1]%in% c("juvenile development rate","fecundity" ))){
    #normal sharpe schollfield
    CalcDF <- filter(df, df$temp<topt)
  }else{
    #reverse sharpe
    # df$rate <- 1/df$rate
    CalcDF <- filter(df, df$temp>topt)

  }

  Boltz <- 0.00008617
  CalcDF$InvT <- 1/((CalcDF$temp+273.15)*Boltz)
  out <- lm(log(CalcDF$rate)~CalcDF$InvT)
  NewE <- -coef(out)[2]
  #browser()
  plot(df$temp,df$rate)
  plot(CalcDF$temp,CalcDF$rate)
  plot(CalcDF$InvT, log(CalcDF$rate))
  
  if(is.null(e)){e <- NA}
    
  return(c(df$species[1],df$trait[1],e,NewE))
}

##Calc E
DFList <- split(df, f=df$curve_ID)
# GetActivation(DFList[[2]])
out <- lapply(DFList, GetActivation)
eDF <- data.frame(do.call(rbind,out))
# plot(eDF$V3,eDF$CalcDF.InvT)
# lines(1:10,1:10)

plot(eDF$V3, eDF$CalcDF.InvT, xlab="Old e", ylab="New e")

abline(a=1:20,b=1:20)

