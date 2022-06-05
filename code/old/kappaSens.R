# load packages
require('tidyverse')
require('patchwork')

rm(list=ls())
graphics.off()

setwd("~/Desktop")

#read in the trait data


rmpredictions <- as_tibble(read.csv('kappaSensAnalysis.csv')) %>%
                 mutate(species = factor(species, levels = c("Rhopalosiphum maidis", "Muscidifurax zaraptor",
                                                             "Tetraneura nigriabdominalis", "Aphis nasturtii",
                                                             "Paracoccus marginatu","Anthonomus grandis")))

levels(rmpredictions$species)

rmpredictions$species <- fct_relevel(rmpredictions$species, "Anthonomus grandis")
rmpredictions$species <- fct_relevel(rmpredictions$species, "Rhopalosiphum maidis", after = Inf)
rmpredictions$species <- fct_relevel(rmpredictions$species, "Aphis nasturtii", after = 1)
rmpredictions$species <- fct_relevel(rmpredictions$species, "Muscidifurax zaraptor", after = 3)



rmestimates <-  
  
  ggplot()+
  geom_line(aes(temp,estimate, colour=species), rmpredictions, size=0.6, lty =4)+
  geom_line(aes(temp,ksmallestimate, colour=species), rmpredictions, size=0.85, lty =1)+
  
  scale_x_continuous(expression(plain(paste(" Temperature, ",degree,"C"))))+
  scale_y_continuous(expression(plain(paste(" Population growth rate ("~italic(r)[m]~")"))),
                     limits=c(-0.001,0.255),
                     expand = c(0.01, 0),
                     breaks=seq(0,0.25, by=0.05))+
  theme_bw(base_size = 12)+
  geom_hline(aes(yintercept = 0), linetype = 2, show.legend = FALSE)+
  scale_colour_manual(values = c("#252525","#525252","#737373","#969696","#d9d9d9", "#bdbdbd"),
                      name=expression(bold("species")),
                      labels = c("R. maidis", "M. zaraptor",
                                 "T. nigriabdominalis", "A. nasturtii",
                                 "P. marginatu","A. grandis"),
                      guide = guide_legend(nrow = 6,ncol =1 ,
                                           direction = "horizontal",
                                           title.position = "top",
                                           title.hjust=0.5))+
  theme(text=element_text(family="Times"))+
  theme(legend.text=element_text(family="Times",face = 'italic', size = 8), 
        legend.position = c(0.225,0.65))+
  geom_text(aes(x = 4, y = 0.245,label = "A"), 
            parse = TRUE, size = 6, colour = "black")+
  theme(legend.title = element_blank())+
  theme(legend.margin=margin(t = -0.2, b = 0.1,r=0.1,l = 0.1, unit='cm'))

ggsave("../results/rmestimates.pdf",
       rmestimates, width = 10, height = 10, units = "cm",device = cairo_pdf)