
############### Analyses of Trait TPCs ################

rm(list=ls())
graphics.off()
unlink("results/*") # remove all existing results files
unlink("results/Figs/*") # remove all existing results plots
unlink("results/SI/*")

require('tidyverse')
require('ggtext')
require('cowplot')
require('patchwork')
require('Cairo')


# load in raw data
df <- as_tibble(read.csv('data/TraitData.csv')) 
dv <- df %>% select(interactor1, interactor1temp, standardisedtraitname, standardisedtraitvalue) %>%
  rename(species = interactor1, temp = interactor1temp, alpha = standardisedtraitvalue) %>%
  filter(standardisedtraitname == '1/alpha', alpha != 'NA') %>%
  mutate(temp = as.numeric(temp))

head(dv)


AlphaPredictions <- as_tibble(read_csv('data/AlphaPredictions.csv')) %>% 
  mutate(alpha = 1/alpha, alphaLwr = 1/alphaLwr, alphaUpr = 1/alphaUpr)

head(AlphaPredictions)

AlphaPlot <- ggplot(AlphaPredictions) +
  geom_line(aes(temp, alpha)) +
  geom_point(aes(temp, alpha), dv, size = 0.75, alpha =0.3) +
  facet_wrap(~species, scales = 'free_y', ncol = 4) +
  theme_bw()+
  geom_ribbon(aes(temp, ymin=alphaLwr, ymax=alphaUpr), AlphaPredictions, fill="#e66101",alpha=0.3,
              inherit.aes = T)+
  theme(text = element_text(size=6, family='Times'), strip.text = element_text(face = "italic"))+
  labs(y=expression(italic(1/alpha)), x=expression(plain(paste(" Temperature, ",degree,"C"))))+
  theme(legend.position = 'none')

AlphaPlot

save_plot(AlphaPlot, file="results/SI/AlphaFits.pdf", base_height=22,base_width = 14, base_asp = 0.75,units="cm")

# fecundity 

dv <- df %>% select(interactor1, interactor1temp, standardisedtraitname, standardisedtraitvalue) %>%
  rename(species = interactor1, temp = interactor1temp, bmax = standardisedtraitvalue) %>%
  filter(standardisedtraitname == 'bmax', bmax != 'NA') %>%
  mutate(temp = as.numeric(temp))

head(dv)


# load in predictions
bmaxPredictions <- as_tibble(read_csv('data/BetaPredictions.csv'))

head(bmaxPredictions)

bmaxPlot <- ggplot(bmaxPredictions) +
  geom_line(aes(temp, bmax)) +
  geom_point(aes(temp, bmax), dv, size = 0.75, alpha =0.3) +
  facet_wrap(~species, scales = 'free_y', ncol = 4) +
  theme_bw() +
  geom_ribbon(aes(temp, ymin=bmaxLwr, ymax=bmaxUpr), bmaxPredictions, fill="#fdb863",alpha=0.3,
              inherit.aes = T)+
  theme(text = element_text(size=6, family='Times'))+
  theme(strip.text = element_text(face = "italic"))+
  labs(y=expression(italic(b)[max]), x=expression(plain(paste(" Temperature, ",degree,"C"))))+
  theme(legend.position = 'none'); bmaxPlot


save_plot(bmaxPlot, file="results/SI/bmaxFits.pdf",  
          base_height=18,base_width = 14, base_asp = 0.75,units="cm")



# zj

dv <- df %>% rename(temp = interactor1temp, species = interactor1, zj = standardisedtraitvalue) %>%
  select(species, temp, standardisedtraitname, zj) %>% 
  filter(standardisedtraitname == 'zj', zj != 'NA') %>%
  mutate(temp = as.numeric(temp))

head(dv)

# load in predictions
zJPredictions <- as_tibble(read_csv('data/zJPredictions.csv'))

# truncate fits for plotting 
zJPredictions <- zJPredictions %>%
  mutate_at(vars(c(zj)), 
            ~ifelse(zj > 0.2, 0.2, .)) %>%
  mutate_at(vars(c(zjLwr)), 
            ~ifelse(zjLwr > 0.2, 0.2, .)) %>%
  mutate_at(vars(c(zjUpr)), 
            ~ifelse(zjUpr > 0.2, 0.2, .)) 

zJPredictions <- zJPredictions %>% filter(zj < 0.2)

head(zJPredictions)

zJPlot <- ggplot(zJPredictions) +
  geom_line(aes(temp, zj)) +
  geom_point(aes(temp, zj), dv, size = 0.75, alpha =0.3) +
  scale_y_continuous(limits=c(-0.001,0.2))+
  facet_wrap(~species, ncol = 4) +
  theme_bw() +
  geom_ribbon(aes(temp, ymin=zjLwr, ymax=zjUpr), zJPredictions, fill="#1f78b4",alpha=0.3,
              inherit.aes = F)+
  theme(text = element_text(size=6, family='Times'))+theme(strip.text = element_text(face = "italic"))+
  labs(y=expression(italic(z[J])), x=expression(plain(paste(" Temperature, ",degree,"C"))))+
  theme(legend.position = 'none')

zJPlot

save_plot(zJPlot, file="results/SI/zJFits.pdf", 
          base_height=20,base_width = 14, base_asp = 0.75,units="cm")



# z


dv <- df %>% rename(temp = interactor1temp, species = interactor1, z = standardisedtraitvalue) %>%
  select(species, temp, standardisedtraitname, z) %>% 
  filter(standardisedtraitname == 'z', z != 'NA') %>%
  mutate(temp = as.numeric(temp))

head(dv)

#load in predictions
zPredictions <- as_tibble(read_csv('data/zPredictions.csv'))

# truncate for plotting 
zPredictions <- zPredictions %>%
  mutate_at(vars(c(z)), 
            ~ifelse(z > 0.2, 0.2, .)) %>%
  mutate_at(vars(c(zLwr)), 
            ~ifelse(zLwr > 0.2, 0.2, .)) %>%
  mutate_at(vars(c(zUpr)), 
            ~ifelse(zUpr > 0.2, 0.2, .)) 

zPredictions <- zPredictions %>% filter(z < 0.2)

head(zPredictions)

zPlot <- ggplot(zPredictions) +
  geom_line(aes(temp, z)) +
  geom_point(aes(temp, z), dv, size = 0.75, alpha =0.3) +
  facet_wrap(~species, ncol = 4)+
  scale_y_continuous(limits=c(-0.001,0.2))+
  theme_bw() +
  geom_ribbon(aes(temp, ymin=zLwr, ymax=zUpr), zPredictions, fill="#a6cee3",alpha=0.5,
              inherit.aes = T)+
  theme(text = element_text(size=6, family='Times'))+theme(strip.text = element_text(face = "italic"))+
  labs(y=expression(italic(z)), x=expression(plain(paste(" Temperature, ",degree,"C"))))+
  theme(legend.position = 'none'); zPlot

save_plot(zPlot, file="results/SI/zFits.pdf", 
          base_height=20,base_width = 14, base_asp = 0.75,units="cm")



####### Distributions of Trait $T_{pk}$s and activation energies ####### 

# Tpks for all species

alpha   <- as_tibble(read.csv('data/alpha_Tpks_AllParams.csv', header = TRUE))
zj      <- as_tibble(read.csv('data/zj_Tpks_AllParams.csv', header = TRUE))
z       <- as_tibble(read.csv('data/z_Tpks_AllParams.csv', header = TRUE))
bmax    <- as_tibble(read.csv('data/bmax_Tpks_AllParams.csv', header = TRUE))

topt  <- bind_rows(alpha,zj,z,bmax) %>%
         mutate(species = as.factor(species), trait = as.factor(trait))

toptz <- topt %>% filter(param == 'topt') 

levels(toptz$species)
levels(toptz$trait)

toptz$trait <- fct_relevel(toptz$trait, "juvenile mortality rate")

toptz <- toptz %>% filter(param!= "rmax")

SpCount <- table(toptz$species)
RmSp    <- names(SpCount)[which(SpCount==1)]

toptz <- filter(toptz, !(species %in% RmSp)) %>%
         filter(species != 'NA')

head(toptz)


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
toptz$species[which(toptz$species=="Scapsipedus icipe")] <- "**Scapsipedus icipe***"
toptz$species[which(toptz$species=="Hyperaspis notata")] <- "**Hyperaspis notata***"
toptz$species[which(toptz$species=="Toxoptera citricida")] <- "**Toxoptera citricida***"
toptz$species[which(toptz$species=="Ophraella communa")] <- "**Ophraella communa***"
toptz$species[which(toptz$species=="Bemisia argentifolii")] <- "**Bemisia argentifolii***"
toptz$species[which(toptz$species=="Dactylopius austrinus")] <- "**Dactylopius austrinus***"
toptz$species[which(toptz$species=="Aulacorthum solani")] <- "**Aulacorthum solani***"

toptz <- toptz %>% filter(species != 'Plutella xylostella')

#order by development rate alpha
alp <-     subset(toptz, toptz$trait=="juvenile development rate")
SPorder <- alp$species[order(alp$estimate)]
toptz$species <- factor(toptz$species, levels=SPorder)

T_pksPlot <- ggplot(toptz, aes(estimate, species, shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(xmin = conf_lower, xmax = conf_upper),width=0.35,size=0.35) +
  geom_point(size = 2.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("", italic(T)[pk]))),
                     limits =c(5,40),
                     expand = c(0, 0),
                     breaks=seq(10,35, by=5))+
  scale_fill_manual(labels = c(expression(plain(paste("juvenile mortality rate (",italic(z[J]),")"))),
                               expression(plain(paste("adult mortality rate (",italic(z),")"))),
                               expression(plain(paste("fecundity (",italic(b[max]),")"))),
                               expression(plain(paste("juvenile development time (",italic(alpha),")")))),
                    values = c("#1f78b4","#a6cee3","#fdb863","#e66101"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow = 1,ncol =4,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(labels = c(expression(plain(paste("juvenile mortality rate (",italic(z[J]),")"))),
                                 expression(plain(paste("adult mortality rate (",italic(z),")"))),
                                 expression(plain(paste("fecundity (",italic(b[max]),")"))),
                                 expression(plain(paste("juvenile development time (",italic(alpha),")")))),
                      values = c("#1f78b4","#a6cee3","#fdb863","#e66101"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=1,ncol=4,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  scale_shape_manual(labels = c(expression(plain(paste("juvenile mortality rate (",italic(z[J]),")"))),
                                expression(plain(paste("adult mortality rate (",italic(z),")"))),
                                expression(plain(paste("fecundity (",italic(b[max]),")"))),
                                expression(plain(paste("juvenile development time (",italic(alpha),")")))),
                     values = c(21,22,23,24),
                     name=expression(bold("")),
                     guide = guide_legend(nrow = 1,ncol =4,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  theme(legend.position = c(0.35,-0.1),legend.text = element_text(size = 10),
        axis.text.y = element_markdown(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  theme(legend.margin=margin(t = -0.4, unit='cm'))+
  # coord_fixed(ratio = 1.5)
  theme(aspect.ratio=2)+
  ggtitle("A")+
  theme(plot.title=element_text(face="bold", size = 15, vjust = -1))

T_pksPlot


E <-  filter(topt, trait == 'juvenile development rate' | trait == 'fecundity rate', param =="e")  
Eh <- filter(topt, trait == 'juvenile mortality rate' | trait == 'adult mortality rate', param =="eh") %>%
  mutate(param=replace(param, param=='eh', 'e')) %>%
  filter(estimate > -8 & estimate < 50)



activation_e <- bind_rows(E, Eh) %>%
                mutate(trait = factor(trait, levels = c('juvenile development rate',
                                                        'juvenile mortality rate', 
                                                        'adult mortality rate',
                                                        'fecundity rate')))

eplot <- ggplot(activation_e, aes(x=estimate))+
  theme_bw()+
  scale_x_continuous(expression(plain(paste("Activation Energy (",italic(E),")"))))+
  labs(y="Density")+
  geom_density(fill = "lightblue")+
  facet_wrap(.~trait, scales = "free")+
  theme(text = element_text(family = 'Times', size = 8))

eplot

save_plot(eplot, file="results/SI/eDists.pdf", base_height=7,base_width=8,base_asp=0.5, units="cm")




# Mean activation energies 

activation_e <-  filter(topt, trait == 'juvenile development rate' |
                          trait == 'fecundity rate', param =="e")  

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
activation_e$species[which(activation_e$species=="Scapsipedus icipe")] <- "**Scapsipedus icipe***"
activation_e$species[which(activation_e$species=="Hyperaspis notata")] <- "**Hyperaspis notata***"
activation_e$species[which(activation_e$species=="Toxoptera citricida")] <- "**Toxoptera citricida***"
activation_e$species[which(activation_e$species=="Ophraella communa")] <- "**Ophraella communa***"
activation_e$species[which(activation_e$species=="Bemisia argentifolii")] <- "**Bemisia argentifolii***"
activation_e$species[which(activation_e$species=="Dactylopius austrinus")] <- "**Dactylopius austrinus***"
activation_e$species[which(activation_e$species=="Aulacorthum solani")] <- "**Aulacorthum solani***"


activation_e <- activation_e %>% 
                mutate(trait = factor(trait, levels = c('fecundity rate','juvenile development rate'))) %>%
                mutate(species = factor(species, levels = SPorder))
                                        
                                    
activation_e <- na.omit(activation_e)

EsPlot <- ggplot(activation_e, aes(estimate, species, shape=trait, colour=trait,fill=trait)) +
  geom_errorbar(aes(xmin = conf_lower, xmax = conf_upper),width=0.35,size=0.35) +
  geom_point(size = 2.5, col="#000000",stroke=0.1)+
  theme_bw(base_size = 12.5) +
  theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("Activation Energy (",italic(E),")"))),
                     limits =c(0,6),
                     expand = c(0, 0),
                     breaks=seq(0,6, by=1))+
  scale_fill_manual(values = c("#fdb863","#e66101"),
                    name=expression(bold("")))+
  scale_colour_manual(values = c("#fdb863","#e66101"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=1,ncol=2,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  scale_shape_manual(values = c(23,24),
                     name=expression(bold("")),
                     guide = guide_legend(nrow = 1,ncol =2,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  theme(legend.position = 'none',legend.text = element_text(size = 8.5),
        axis.text.y = element_markdown(face = 'italic'))+
  theme(text=element_text(family="Times"))+
  theme(legend.margin=margin(t = -0.4, unit='cm'))+
  geom_vline(xintercept=0.65, lty ='dashed')+
  theme(axis.text.y=element_blank())+
  theme(plot.margin=margin(l=-4,unit="cm"))+
  theme(aspect.ratio=2)+
  ggtitle("B")+
  theme(plot.title=element_text(face="bold", size = 15, vjust = -1))

EsPlot

plotMain <- T_pksPlot + theme(legend.position="none") +
  EsPlot + theme(legend.position="none")

legend <- get_legend(T_pksPlot + theme(legend.position = "bottom"))

Tpks_Es_Plot <- plot_grid(plotMain, legend,ncol = 1, rel_heights = c(1, 0.05))

Tpks_Es_Plot

save_plot(Tpks_Es_Plot, file="results/Figs/Tpks_Es_Plot.pdf", base_height=15,base_asp=1.5, units="cm")


#_________________________________


speciesmass <- as_tibble(read_csv('data/TraitData.csv')) %>%
  select(originaltraitname, originaltraitvalue, interactor1) %>%
  filter(originaltraitname == 'body size' & originaltraitvalue != 'NA')


sizeMeans <- 
  speciesmass %>% 
  group_by(interactor1) %>% 
  summarise(avg = mean(originaltraitvalue), sd = sd(originaltraitvalue)) %>%
  arrange(avg)

head(sizeMeans)

write_csv(sizeMeans, 'data/sizeMeans.csv')

distinct(sizeMeans, interactor1) %>% print(n=70)

bodyMass <- sizeMeans %>% 
  rename(species = interactor1) %>%
  mutate(curve_ID = case_when(species == 'Aedes albopictus' ~ '1',
                              species == 'Aedes aegypti' ~ '2',
                              species == 'Anthonomus grandis' ~ '3',
                              species == 'Paracoccus marginatus' ~ '4',
                              species == 'Aphis gossypii' ~ '6',
                              species == 'Harmonia axyridis' ~ '7',
                              species == 'Tribolium castaneum' ~ '8',
                              species == 'Aedes krombeini' ~ '9',
                              species == 'Bemisia tabaci' ~ '10',
                              species == 'Tetraneura nigriabdominalis' ~ '11',
                              species == 'Stethorus punctillum' ~ '12',
                              species == 'Tetranychus mcdanieli' ~ '13',
                              species == 'Tetranychus urticae' ~ '14',
                              species == 'Clavigralla tomentosicollis' ~ '15',
                              species == 'Planococcus citri' ~ '16',
                              species == 'Muscidifurax zaraptor' ~ '17',
                              species == 'Aphis nasturtii' ~ '18',
                              species == 'Rhopalosiphum maidis' ~ '19',
                              species == 'Anopheles gambiae' ~ '20',
                              species == 'Anoplophora glabripennis' ~ '21',
                              species == 'Amblyseius womersleyi' ~ '22',
                              species == 'Trichogramma sp. nr. Lutea' ~ '23',
                              species == 'Trichogramma bruni' ~ '24',
                              species == 'Culex annulirostris' ~ '25',
                              species == 'Macrocentrus iridescens' ~ '26',
                              species == 'Otiorhynchus sulcatus' ~ '27',
                              species == 'Drosophila suzukii' ~ '28',
                              species == 'Gastrolina depressa' ~ '29',
                              species == 'Laricobius nigrinus' ~ '30',
                              species == 'Aubeonymus mariaefranciscae' ~ '31',
                              species == 'Iphiseius degenerans' ~ '32',
                              species == 'Amblyseius swirskii' ~ '33',
                              species == 'Macrosiphum euphorbia' ~ '34',
                              species == 'Myzus persicae' ~ '35',
                              species == 'Tetranychus evansi' ~ '36',
                              species == 'Helicoverpa armigera' ~ '37',
                              species == 'Antestiopsis thunbergii' ~ '38',
                              species == 'Monochamus leuconotus' ~ '39',
                              species == 'Kampimodromus aberrans' ~ '40',
                              species == 'Phenacoccus solenopsis' ~ '41',
                              species == 'Leptinotarsa decemlineata' ~ '42',
                              species == 'Halyomorpha halys' ~ '43',
                              species == 'Muscidifurax raptorellus' ~ '44',
                              species == 'Thrips hawaiiensis' ~ '45',
                              species == 'Scapsipedus icipe' ~ '46',
                              species == 'Hyperaspis notata' ~ '47',
                              species == 'Toxoptera citricida' ~ '48',
                              species == 'Ophraella communa' ~ '49',
                              species == 'Bemisia argentifolii' ~ '50',
                              species == 'Dactylopius austrinus' ~ '51',
                              species == 'Aulacorthum solani' ~ '52')) %>%
  arrange(curve_ID) %>% 
  rename(massspecies = species, masscurve_ID = curve_ID) %>% 
  filter(masscurve_ID != 'NA')


Tc <- as_tibble(read.csv('data/alpha_Tpks_AllParams.csv')) %>%
  mutate(curve_ID = case_when(species == 'Aedes albopictus' ~ '1',
                              species == 'Aedes aegypti' ~ '2',
                              species == 'Anthonomus grandis' ~ '3',
                              species == 'Paracoccus marginatus' ~ '4',
                              species == 'Aphis gossypii' ~ '6',
                              species == 'Harmonia axyridis' ~ '7',
                              species == 'Tribolium castaneum' ~ '8',
                              species == 'Aedes krombeini' ~ '9',
                              species == 'Bemisia tabaci' ~ '10',
                              species == 'Tetraneura nigriabdominalis' ~ '11',
                              species == 'Stethorus punctillum' ~ '12',
                              species == 'Tetranychus mcdanieli' ~ '13',
                              species == 'Tetranychus urticae' ~ '14',
                              species == 'Clavigralla tomentosicollis' ~ '15',
                              species == 'Planococcus citri' ~ '16',
                              species == 'Muscidifurax zaraptor' ~ '17',
                              species == 'Aphis nasturtii' ~ '18',
                              species == 'Rhopalosiphum maidis' ~ '19',
                              species == 'Anopheles gambiae' ~ '20',
                              species == 'Anoplophora glabripennis' ~ '21',
                              species == 'Amblyseius womersleyi' ~ '22',
                              species == 'Trichogramma sp. nr. Lutea' ~ '23',
                              species == 'Trichogramma bruni' ~ '24',
                              species == 'Culex annulirostris' ~ '25',
                              species == 'Macrocentrus iridescens' ~ '26',
                              species == 'Otiorhynchus sulcatus' ~ '27',
                              species == 'Drosophila suzukii' ~ '28',
                              species == 'Gastrolina depressa' ~ '29',
                              species == 'Laricobius nigrinus' ~ '30',
                              species == 'Aubeonymus mariaefranciscae' ~ '31',
                              species == 'Iphiseius degenerans' ~ '32',
                              species == 'Amblyseius swirskii' ~ '33',
                              species == 'Macrosiphum euphorbia' ~ '34',
                              species == 'Myzus persicae' ~ '35',
                              species == 'Tetranychus evansi' ~ '36',
                              species == 'Helicoverpa armigera' ~ '37',
                              species == 'Antestiopsis thunbergii' ~ '38',
                              species == 'Monochamus leuconotus' ~ '39',
                              species == 'Kampimodromus aberrans' ~ '40',
                              species == 'Phenacoccus solenopsis' ~ '41',
                              species == 'Leptinotarsa decemlineata' ~ '42',
                              species == 'Halyomorpha halys' ~ '43',
                              species == 'Muscidifurax raptorellus' ~ '44',
                              species == 'Thrips hawaiiensis' ~ '45',
                              species == 'Scapsipedus icipe' ~ '46',
                              species == 'Hyperaspis notata' ~ '47',
                              species == 'Toxoptera citricida' ~ '48',
                              species == 'Ophraella communa' ~ '49',
                              species == 'Bemisia argentifolii' ~ '50',
                              species == 'Dactylopius austrinus' ~ '51',
                              species == 'Aulacorthum solani' ~ '52')) %>%
  arrange(curve_ID) 

Tc %>% distinct(species, curve_ID) %>% print(n=55)

# Extract relevant data 
a_pk <- Tc %>% filter(param == 'rmax') %>% 
  rename(a_pk = estimate, a_pkLwr = conf_lower,a_pkUpr = conf_upper) %>% 
  select(a_pk, a_pkLwr, a_pkUpr, species, curve_ID)

T_pk <- Tc %>% filter(param == 'topt') %>% 
  rename(T_pk = estimate, T_pkLwr = conf_lower,T_pkUpr = conf_upper) %>% 
  select(T_pk, T_pkLwr, T_pkUpr)

bodyMass <- bodyMass %>% select(avg, masscurve_ID)

Tc <- bind_cols(a_pk, T_pk, bodyMass)

Tc <- Tc %>% rename(mass = avg) %>% select(-masscurve_ID)

write_csv(Tc, 'data/a_pksT_pksMass.csv')


# fit linear model (note that alpha_pk is linear in log-log scale)
TcLm <- lm(log(a_pk) ~ T_pk, data = Tc); summary(TcLm)

Tc %>%
  ggplot(aes(x=T_pk,y = log(a_pk)))+
  geom_point()+
  geom_smooth(method = 'lm')



a_data <- as_tibble(read.csv('data/a_pksT_pksMass.csv')) # load in data 
a_data$kT <- 1/(8.617333262145 * 10^-5 * (a_data$T_pk+273.15))
head(a_data)


#plot uncorrected data in log-log scale
a_data %>%
  ggplot(aes(x=log(mass), y = log(a_pk)))+
  geom_point()+
  geom_smooth(method = 'lm')

#plot a_pk vs Tpk
a_data %>%
  ggplot(aes(x = T_pk, y = a_pk)) +
  geom_point()+
  geom_smooth(method = 'lm')


# linear model (note the allometry is linear in log-log scale)
a_model <- lm(log(a_pk) ~ log(mass) + kT, data = a_data)
summary(a_model)
coef(a_model)
save(a_model,file="results/a_MTE_model.Rdata")


a_cf <-  confint(a_model,level = .95)
anova(a_model)
coef(a_model)[2]

#plot a_pk in 1/kT, correcting for mass
MassCorrectedApkTpk <- a_data %>%
  ggplot(aes(x = T_pk, y = log(a_pk/mass^coef(a_model)[2]))) +
  geom_smooth(method = 'lm', colour = '#636363', size=0.15) +
  scale_y_continuous(expression(plain(paste("ln((", italic(1/alpha[','][' '][pk])~")/",
                                            italic(M^-0.27),")"))))+
  scale_x_continuous(expression(plain(paste(italic("T"[pk])))),
                     limits =c(17,36.5),
                     expand = c(0, 0),
                     breaks=seq(20,35, by=5))+
  theme_bw()+
  theme(text=element_text(family="Times"), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank())+
  geom_point(size = 1, col="#000000",stroke=0.1, shape = 24, fill="#e66101")+
  annotate("text", x = 18.5, y = -1.1, label = 'A', size=5, family='Times')

MassCorrectedApkTpk


#===============================================================

# 4. Peak fecundity (bmax)

bodyMass <- as_tibble(read.csv('data/sizeMeans.csv')) %>% 
  rename(species = interactor1) %>% 
  mutate(curve_ID = case_when(species == 'Aedes aegypti' ~ '2',
                              species == 'Anthonomus grandis' ~ '3',
                              species == 'Paracoccus marginatus' ~ '4',
                              species == 'Acyrthosiphon pisum' ~ '5',
                              species == 'Aphis gossypii' ~ '6',
                              species == 'Aedes krombeini' ~ '9',
                              species == 'Bemisia tabaci' ~ '10',
                              species == 'Tetraneura nigriabdominalis' ~ '11',
                              species == 'Stethorus punctillum' ~ '12',
                              species == 'Tetranychus mcdanieli' ~ '13',
                              species == 'Clavigralla tomentosicollis' ~ '15',
                              species == 'Planococcus citri' ~ '16',
                              species == 'Muscidifurax zaraptor' ~ '17',
                              species == 'Aphis nasturtii' ~ '18',
                              species == 'Rhopalosiphum maidis' ~ '19',
                              species == 'Anoplophora glabripennis' ~ '21',
                              species == 'Amblyseius womersleyi' ~ '22',
                              species == 'Trichogramma sp. nr. Lutea' ~ '23',
                              species == 'Trichogramma bruni' ~ '24',
                              species == 'Drosophila suzukii' ~ '28',
                              species == 'Iphiseius degenerans' ~ '32',
                              species == 'Tetranychus evansi' ~ '36',
                              species == 'Helicoverpa armigera' ~ '37',
                              species == 'Antestiopsis thunbergii' ~ '38',
                              species == 'Monochamus leuconotus' ~ '39',
                              species == 'Kampimodromus aberrans' ~ '40',
                              species == 'Phenacoccus solenopsis' ~ '41',
                              species == 'Halyomorpha halys' ~ '43',
                              species == 'Muscidifurax raptorellus' ~ '44',
                              species == 'Thrips hawaiiensis' ~ '45',
                              species == 'Hylobius transversovittatus' ~ '46',
                              species == 'Callosobruchus maculatus' ~ '47',   
                              species == 'Callosobruchus chinensis' ~ '48',   
                              species == 'Callosobruchus analis' ~ '49',      
                              species == 'Callosobruchus rhodesianus' ~ '50',
                              species == 'Sepedon spinipes' ~ '51',
                              species == 'Plutella xylostella' ~ '52',
                              species == 'Scapsipedus icipe' ~ '53',
                              species == 'Hyperaspis notata' ~ '54',
                              species == 'Toxoptera citricida' ~ '55',
                              species == 'Ophraella communa' ~ '56',
                              species == 'Bemisia argentifolii' ~ '57',
                              species == 'Dactylopius austrinus' ~ '58',
                              species == 'Aulacorthum solani' ~ '59')) %>%
  arrange(curve_ID) %>% filter(curve_ID != 'NA') %>%
  rename(massspecies = species, masscurve_ID = curve_ID) 

distinct(bodyMass, massspecies) %>% print(n=55)


bmax <- as_tibble(read.csv('data/bmax_Tpks_AllParams.csv')) %>%
  mutate(curve_ID = case_when(species == 'Aedes aegypti' ~ '2',
                              species == 'Anthonomus grandis' ~ '3',
                              species == 'Paracoccus marginatus' ~ '4',
                              species == 'Acyrthosiphon pisum' ~ '5',
                              species == 'Aphis gossypii' ~ '6',
                              species == 'Aedes krombeini' ~ '9',
                              species == 'Bemisia tabaci' ~ '10',
                              species == 'Tetraneura nigriabdominalis' ~ '11',
                              species == 'Stethorus punctillum' ~ '12',
                              species == 'Tetranychus mcdanieli' ~ '13',
                              species == 'Clavigralla tomentosicollis' ~ '15',
                              species == 'Planococcus citri' ~ '16',
                              species == 'Muscidifurax zaraptor' ~ '17',
                              species == 'Aphis nasturtii' ~ '18',
                              species == 'Rhopalosiphum maidis' ~ '19',
                              species == 'Anoplophora glabripennis' ~ '21',
                              species == 'Amblyseius womersleyi' ~ '22',
                              species == 'Trichogramma sp. nr. Lutea' ~ '23',
                              species == 'Trichogramma bruni' ~ '24',
                              species == 'Drosophila suzukii' ~ '28',
                              species == 'Iphiseius degenerans' ~ '32',
                              species == 'Tetranychus evansi' ~ '36',
                              species == 'Helicoverpa armigera' ~ '37',
                              species == 'Antestiopsis thunbergii' ~ '38',
                              species == 'Monochamus leuconotus' ~ '39',
                              species == 'Kampimodromus aberrans' ~ '40',
                              species == 'Phenacoccus solenopsis' ~ '41',
                              species == 'Halyomorpha halys' ~ '43',
                              species == 'Muscidifurax raptorellus' ~ '44',
                              species == 'Thrips hawaiiensis' ~ '45',
                              species == 'Hylobius transversovittatus' ~ '46',
                              species == 'Callosobruchus maculatus' ~ '47',   
                              species == 'Callosobruchus chinensis' ~ '48',   
                              species == 'Callosobruchus analis' ~ '49',      
                              species == 'Callosobruchus rhodesianus' ~ '50',
                              species == 'Sepedon spinipes' ~ '51',
                              species == 'Plutella xylostella' ~ '52',
                              species == 'Scapsipedus icipe' ~ '53',
                              species == 'Hyperaspis notata' ~ '54',
                              species == 'Toxoptera citricida' ~ '55',
                              species == 'Ophraella communa' ~ '56',
                              species == 'Bemisia argentifolii' ~ '57',
                              species == 'Dactylopius austrinus' ~ '58',
                              species == 'Aulacorthum solani' ~ '59')) %>%
  arrange(curve_ID) %>% filter(curve_ID != 'NA')

distinct(bmax, species) %>% print(n=50)


b_pk <- bmax %>% filter(param == 'rmax') %>% 
  rename(bmax = estimate, bmaxLwr = conf_lower, bmaxUpr = conf_upper) %>% 
  select(bmax, bmaxLwr, bmaxUpr, species, curve_ID)

T_pk <- bmax %>% filter(param == 'topt') %>% 
  rename(T_pk = estimate, T_pkLwr = conf_lower,T_pkUpr = conf_upper) %>% 
  select(T_pk, T_pkLwr, T_pkUpr)

bodyMass <- bodyMass %>% select(avg, masscurve_ID)

bmax_data <- bind_cols(b_pk, T_pk, bodyMass)

bmax_data <- bmax_data %>% rename(mass = avg)  %>% select(-masscurve_ID)

write_csv(bmax_data, 'data/bmaxT_pksMass.csv')

bmax_data %>%
  ggplot(aes(x=T_pk,y = log(bmax)))+
  geom_point()+
  geom_smooth(method = 'lm')


bmax_data <- as_tibble(read.csv('data/bmaxT_pksMass.csv')) # load in data 
bmax_data$kT <- 1/(8.617333262145 * 10^-5 * (bmax_data$T_pk+273.15))
head(bmax_data)


#plot uncorrected data in log-log scale
bmax_data %>%
  ggplot(aes(x=log(mass), y = log(bmax)))+
  geom_point()+
  geom_smooth(method = 'lm')


#plot b_pk vs T
bmax_data %>%
  ggplot(aes(x = T_pk, y = bmax)) +
  geom_point()+
  geom_smooth(method = 'lm')


# linear model (note the allometry is linear in log-log scale)
bmax_model <- lm(log(bmax) ~ log(mass) + kT, data = bmax_data)
summary(bmax_model)
coef(bmax_model)

save(bmax_model,file="results/bmax_MTE_model.Rdata")

bmax_cf <-  confint(bmax_model,level = .95)
anova(bmax_model)

#plot bmax in 1/kT, correcting for mass

MassCorrectedBmaxTpk <- 
  bmax_data %>%
  ggplot(aes(x = T_pk, y = log(bmax/mass^coef(bmax_model)[2]))) +
  geom_smooth(method = 'lm', colour = '#636363', size=0.15)+
  scale_y_continuous(expression(plain(paste("ln(", italic(b[max][','][' '][pk])~"/",
                                            italic(M^0.12),")"))))+
  scale_x_continuous(expression(plain(paste(italic("T"[pk])))),
                     limits =c(19,35),
                     expand = c(0, 0),
                     breaks=seq(20,35, by=5))+  
  theme_bw()+
  theme(text=element_text(family="Times"), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  geom_point(size = 1, col="#000000",stroke=0.1, shape = 23, fill="#fdb863")+
  annotate("text", x = 20.25, y = 3.45, label = 'D', size=5, family='Times')

MassCorrectedBmaxTpk



#±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±

# 2. Juvenile mortality rate (z_j)

bodyMass <- as_tibble(read.csv('data/sizeMeans.csv')) %>% 
  rename(species = interactor1) %>% 
  mutate(curve_ID = case_when(species == 'Thrips hawaiiensis' ~ '1',
                              species == 'Aedes aegypti' ~ '2',
                              species == 'Anthonomus grandis' ~ '3',
                              species == 'Paracoccus marginatus' ~ '4',
                              species == 'Acyrthosiphon pisum' ~ '5',
                              species == 'Aphis gossypii' ~ '6',
                              species == 'Harmonia axyridis' ~ '7',
                              species == 'Tribolium castaneum' ~ '8',
                              species == 'Tetranychus mcdanieli' ~ '9',
                              species == 'Bemisia tabaci' ~ '10',
                              species == 'Tetraneura nigriabdominalis' ~ '11',
                              species == 'Amblyseius womersleyi' ~ '13',
                              species == 'Stethorus punctillum' ~ '12',
                              species == 'Tetranychus urticae' ~ '14',
                              species == 'Planococcus citri' ~ '16',
                              species == 'Muscidifurax zaraptor' ~ '17',
                              species == 'Aphis nasturtii' ~ '18',
                              species == 'Rhopalosiphum maidis' ~ '19',
                              species == 'Anopheles gambiae' ~ '20',
                              species == 'Anoplophora glabripennis' ~ '21',
                              species == 'Culex annulirostris' ~ '25',
                              species == 'Laricobius nigrinus' ~ '30',
                              species == 'Aubeonymus mariaefranciscae' ~ '31',
                              species == 'Helicoverpa armigera' ~ '37',
                              species == 'Halyomorpha halys' ~ '43',
                              species == 'Muscidifurax raptorellus' ~ '44',
                              species == 'Aedes albopictus' ~ '45',
                              species == 'Scapsipedus icipe' ~ '46',
                              species == 'Hyperaspis notata' ~ '47',
                              species == 'Toxoptera citricida' ~ '48',
                              species == 'Ophraella communa' ~ '49',
                              species == 'Bemisia argentifolii' ~ '50',
                              species == 'Dactylopius austrinus' ~ '51',
                              species == 'Aulacorthum solani' ~ '52')) %>%
  arrange(curve_ID) %>% 
  rename(massspecies = species, masscurve_ID = curve_ID) %>% 
  filter(masscurve_ID != 'NA')

bodyMass %>% distinct(massspecies,masscurve_ID) %>% print(n=50)

zjPk <- as_tibble(read.csv('data/zj_Tpks_AllParams.csv')) %>%
  mutate(curve_ID = case_when(species == 'Thrips hawaiiensis' ~ '1',
                              species == 'Aedes aegypti' ~ '2',
                              species == 'Anthonomus grandis' ~ '3',
                              species == 'Paracoccus marginatus' ~ '4',
                              species == 'Acyrthosiphon pisum' ~ '5',
                              species == 'Aphis gossypii' ~ '6',
                              species == 'Harmonia axyridis' ~ '7',
                              species == 'Tribolium castaneum' ~ '8',
                              species == 'Tetranychus mcdanieli' ~ '9',
                              species == 'Bemisia tabaci' ~ '10',
                              species == 'Tetraneura nigriabdominalis' ~ '11',
                              species == 'Amblyseius womersleyi' ~ '13',
                              species == 'Stethorus punctillum' ~ '12',
                              species == 'Tetranychus urticae' ~ '14',
                              species == 'Planococcus citri' ~ '16',
                              species == 'Muscidifurax zaraptor' ~ '17',
                              species == 'Aphis nasturtii' ~ '18',
                              species == 'Rhopalosiphum maidis' ~ '19',
                              species == 'Anopheles gambiae' ~ '20',
                              species == 'Anoplophora glabripennis' ~ '21',
                              species == 'Culex annulirostris' ~ '25',
                              species == 'Laricobius nigrinus' ~ '30',
                              species == 'Aubeonymus mariaefranciscae' ~ '31',
                              species == 'Helicoverpa armigera' ~ '37',
                              species == 'Halyomorpha halys' ~ '43',
                              species == 'Muscidifurax raptorellus' ~ '44',
                              species == 'Aedes albopictus' ~ '45',
                              species == 'Scapsipedus icipe' ~ '46',
                              species == 'Hyperaspis notata' ~ '47',
                              species == 'Toxoptera citricida' ~ '48',
                              species == 'Ophraella communa' ~ '49',
                              species == 'Bemisia argentifolii' ~ '50',
                              species == 'Dactylopius austrinus' ~ '51',
                              species == 'Aulacorthum solani' ~ '52')) %>%
  arrange(curve_ID) 

zjPk %>% distinct(species,curve_ID) %>% print(n=50)

zj_pk <- zjPk %>% filter(param == 'rmax') %>% 
  rename(zjpk = estimate, zjpkLwr = conf_lower, zjpkUpr = conf_upper) %>% 
  select(zjpk, zjpkLwr, zjpkUpr, species, curve_ID) %>% mutate(zjpk = as.numeric(zjpk)) %>%
  mutate(zjpk = 1/zjpk, zjpkLwr = 1/zjpkLwr, zjpkUpr = 1/zjpkUpr)

T_pk <- zjPk %>% filter(param == 'topt') %>% 
  rename(T_pk = estimate, T_pkLwr = conf_lower,T_pkUpr = conf_upper) %>% 
  select(T_pk, T_pkLwr, T_pkUpr)

bodyMass <- bodyMass %>% select(avg, masscurve_ID)

zj_data <- bind_cols(zj_pk, T_pk, bodyMass)

zj_data <- zj_data %>% rename(mass = avg) %>% select(-masscurve_ID)

write_csv(zj_data, 'data/zj_pksT_pksMass.csv')

zj_model <- lm(log(zjpk) ~ T_pk, data = zj_data); summary(zj_model)

zj_data %>%
  ggplot(aes(x=T_pk,y = log(zjpk)))+
  geom_point()+
  geom_smooth(method = 'lm')

zj_data <- as_tibble(read.csv('data/zj_pksT_pksMass.csv')) # load in data 
zj_data$kT <- 1/(8.617333262145 * 10^-5 * (zj_data$T_pk+273.15))
head(zj_data)


#plot uncorrected data in log-log scale
zj_data %>%
  ggplot(aes(x=log(mass), y = log(zjpk)))+
  geom_point()+
  geom_smooth(method = 'lm')


#plot zj_pk vs Tpk
zj_data %>%
  ggplot(aes(x = T_pk, y = zjpk)) +
  geom_point()+
  geom_smooth(method = 'lm')


# linear model (note the allometry is linear in log-log scale)
zj_model <- lm(log(zjpk) ~ log(mass) + kT, data = zj_data)
summary(zj_model)
coef(zj_model)

save(zj_model,file="results/zj_MTE_model.Rdata")

b_cf <-  confint(zj_model,level = .95)
anova(zj_model)


#plot a_pk in 1/kT, correcting for mass
MassCorrectedzjpkTpk <- 
  zj_data %>%
  ggplot(aes(x = T_pk, y = log(zjpk/mass^coef(zj_model)[2])))+
  geom_smooth(method = 'lm', colour = '#636363', size=0.15)+
  scale_y_continuous(expression(plain(paste("ln(", italic(z[J][','][' '][pk])~"/",
                                                            italic(M^-0.22),")"))))+
  scale_x_continuous(expression(plain(paste(italic("T"[pk])))),
                     limits =c(5,34),
                     expand = c(0, 0),
                     breaks=seq(10,30, by=5))+
  theme_bw()+
  theme(text=element_text(family="Times"), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  geom_point(size = 1, col="#000000",stroke=0.1, shape = 21, fill="#1f78b4")+
  annotate("text", x = 7.5, y = -3, label = 'B', size=5, family='Times')

MassCorrectedzjpkTpk



#==================================================================

# 3. Adult mortality rate (z)

bodyMass <- as_tibble(read.csv('data/sizeMeans.csv')) %>% 
  rename(species = interactor1) %>% 
  mutate(curve_ID = case_when(species == 'Culex pipiens' ~ '1',
                              species == 'Plutella xylostella' ~ '2',
                              species == 'Thrips hawaiiensis' ~ '3',
                              species == 'Phenacoccus solenopsis' ~ '4',
                              species == 'Antestiopsis thunbergii' ~ '5',
                              species == 'Culex quinquefasciatus' ~ '6',
                              species == 'Halyomorpha halys' ~ '7',
                              species == 'Monochamus leuconotus' ~ '8',
                              species == 'Anthonomus grandis' ~ '9',
                              species == 'Paracoccus marginatus' ~ '10',
                              species == 'Aphis gossypii' ~ '12',
                              species == 'Tribolium castaneum' ~ '14',
                              species == 'Tetranychus mcdanieli' ~ '15',
                              species == 'Bemisia tabaci' ~ '16',
                              species == 'Tetraneura nigriabdominalis' ~ '17',
                              species == 'Stethorus punctillum' ~ '18',
                              species == 'Aedes krombeini' ~ '19',
                              species == 'Muscidifurax zaraptor' ~ '21',
                              species == 'Aphis nasturtii' ~ '22',
                              species == 'Rhopalosiphum maidis' ~ '23',
                              species == 'Anopheles gambiae' ~ '24',
                              species == 'Anoplophora glabripennis' ~ '25',
                              species == 'Helicoverpa armigera' ~ '26',
                              species == 'Aedes albopictus' ~ '27',
                              species == 'Trichogramma bruni' ~ '28',
                              species == 'Trichogramma sp. nr. Lutea' ~ '29',
                              species == 'Aedes aegypti' ~ '30',
                              species == 'Scapsipedus icipe' ~ '46',
                              species == 'Hyperaspis notata' ~ '47',
                              species == 'Toxoptera citricida' ~ '48',
                              species == 'Ophraella communa' ~ '49',
                              species == 'Bemisia argentifolii' ~ '50',
                              species == 'Dactylopius austrinus' ~ '51',
                              species == 'Aulacorthum solani' ~ '52')) %>%
  arrange(curve_ID) %>% 
  rename(massspecies = species, masscurve_ID = curve_ID) %>% 
  filter(masscurve_ID != 'NA')

distinct(bodyMass, massspecies,masscurve_ID) %>% print(n=50)

zPk <- as_tibble(read.csv('data/z_Tpks_AllParams.csv')) %>%
  mutate(curve_ID = case_when(species == 'Culex pipiens' ~ '1',
                              species == 'Plutella xylostella' ~ '2',
                              species == 'Thrips hawaiiensis' ~ '3',
                              species == 'Phenacoccus solenopsis' ~ '4',
                              species == 'Antestiopsis thunbergii' ~ '5',
                              species == 'Culex quinquefasciatus' ~ '6',
                              species == 'Halyomorpha halys' ~ '7',
                              species == 'Monochamus leuconotus' ~ '8',
                              species == 'Anthonomus grandis' ~ '9',
                              species == 'Paracoccus marginatus' ~ '10',
                              species == 'Aphis gossypii' ~ '12',
                              species == 'Tribolium castaneum' ~ '14',
                              species == 'Tetranychus mcdanieli' ~ '15',
                              species == 'Bemisia tabaci' ~ '16',
                              species == 'Tetraneura nigriabdominalis' ~ '17',
                              species == 'Stethorus punctillum' ~ '18',
                              species == 'Aedes krombeini' ~ '19',
                              species == 'Muscidifurax zaraptor' ~ '21',
                              species == 'Aphis nasturtii' ~ '22',
                              species == 'Rhopalosiphum maidis' ~ '23',
                              species == 'Anopheles gambiae' ~ '24',
                              species == 'Anoplophora glabripennis' ~ '25',
                              species == 'Helicoverpa armigera' ~ '26',
                              species == 'Aedes albopictus' ~ '27',
                              species == 'Trichogramma bruni' ~ '28',
                              species == 'Trichogramma sp. nr. Lutea' ~ '29',
                              species == 'Aedes aegypti' ~ '30',
                              species == 'Scapsipedus icipe' ~ '46',
                              species == 'Hyperaspis notata' ~ '47',
                              species == 'Toxoptera citricida' ~ '48',
                              species == 'Ophraella communa' ~ '49',
                              species == 'Bemisia argentifolii' ~ '50',
                              species == 'Dactylopius austrinus' ~ '51',
                              species == 'Aulacorthum solani' ~ '52')) %>%
  arrange(curve_ID) 

distinct(zPk, species) %>% print(n=50)


z_pk <- zPk %>% filter(param == 'rmax') %>% 
  rename(zpk = estimate, zpkLwr = conf_lower, zpkUpr = conf_upper) %>% 
  select(zpk, zpkLwr, zpkUpr, species, curve_ID) %>% 
  mutate(zpk = 1/zpk, zpkLwr = 1/zpkLwr, zpkUpr = 1/zpkUpr)

T_pk <- zPk %>% filter(param == 'topt') %>% 
  rename(T_pk = estimate, T_pkLwr = conf_lower,T_pkUpr = conf_upper) %>% 
  select(T_pk, T_pkLwr, T_pkUpr)

bodyMass <- bodyMass %>% select(avg, masscurve_ID)

z_data <- bind_cols(z_pk, T_pk, bodyMass)

z_data <- z_data %>% rename(mass = avg) %>% select(-masscurve_ID)

write_csv(z_data, 'data/z_pksT_pksMass.csv')


z_model <- lm(log(zpk) ~ T_pk, data = z_data); summary(z_model)

z_data %>%
  ggplot(aes(x=T_pk,y = log(zpk)))+
  geom_point()+
  geom_smooth(method = 'lm')


z_data <- as_tibble(read.csv('data/z_pksT_pksMass.csv')) # load in data 
z_data$kT <- 1/(8.617333262145 * 10^-5 * (z_data$T_pk+273.15))
head(z_data)


#plot uncorrected data in log-log scale
z_data %>%
  ggplot(aes(x=log(mass), y = log(zpk)))+
  geom_point()+
  geom_smooth(method = 'lm')

#plot z_pk vs Tpk
z_data %>%
  ggplot(aes(x = T_pk, y = zpk)) +
  geom_point()+
  geom_smooth(method = 'lm')

# linear model (note the allometry is linear in log-log scale)
z_model <- lm(log(zpk) ~ log(mass) + kT, data = z_data)
summary(z_model)
coef(z_model)

save(z_model,file="results/z_MTE_model.Rdata")

z_cf <-  confint(z_model,level = .95)
anova(z_model)


#plot a_pk in 1/kT, correcting for mass
MassCorrectedzpkTpk <- 
  z_data %>%
  ggplot(aes(x = T_pk, y = log(zpk/mass^coef(z_model)[2]))) +
  geom_smooth(method = 'lm', colour = '#636363', size=0.15)+
  scale_y_continuous(expression(plain(paste("ln(", italic(z[','][' '][pk])~"/",italic(M^-0.09),")"))))+
  scale_x_continuous(expression(plain(paste(italic("T"[pk])))),
                     limits=c(6.25,27),
                     expand = c(0, 0),
                     breaks=seq(10,25, by=5))+  
  theme_bw()+
  theme(text=element_text(family="Times"), panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
  geom_point(size = 1, col="#000000",stroke=0.1, shape = 22, fill="#a6cee3")+
  annotate("text", x = 8, y = -2.45, label = 'C', size=5, family='Times')

MassCorrectedzpkTpk


#========= plot hotter-is-better panel 

p1 <- (MassCorrectedApkTpk+MassCorrectedzjpkTpk)/(MassCorrectedzpkTpk+MassCorrectedBmaxTpk); p1

save_plot(p1, file="results/SI/Traits_HotterIsBetter.pdf", 
          base_height=10,base_width = 12, base_asp = 0.75,units="cm")



# Maximum population growth rate ($r_m$) calculations


# Create argument containing target species
spps <- "case_when(species == 'Anoplophora glabripennis' ~ '1', 
        species == 'Halyomorpha halys' ~ '2', 
        species == 'Aedes aegypti' ~ '3',
        species == 'Anthonomus grandis' ~ '4',
        species == 'Paracoccus marginatus' ~ '5',
        species == 'Aphis gossypii' ~ '6',
        species == 'Bemisia tabaci' ~ '7',
        species == 'Tetraneura nigriabdominalis' ~ '8',
        species == 'Stethorus punctillum' ~ '9',
        species == 'Tetranychus mcdanieli' ~ '10',
        species == 'Muscidifurax zaraptor' ~ '11',
        species == 'Aphis nasturtii' ~ '12',
        species == 'Rhopalosiphum maidis' ~ '13',
        species == 'Thrips hawaiiensis' ~ '14',
        species == 'Helicoverpa armigera' ~ '15',
        species == 'Scapsipedus icipe' ~ '16',
        species == 'Hyperaspis notata' ~ '17',
        species == 'Toxoptera citricida' ~ '18',
        species == 'Ophraella communa' ~ '19',
        species == 'Bemisia argentifolii' ~ '20',
        species == 'Dactylopius austrinus' ~ '21',
        species == 'Aulacorthum solani' ~ '22')"

# Read in the trait data

alpha <- as_tibble(read.csv('data/AlphaPredictions.csv')) %>% 
  select(species, temp, alpha, alphaLwr, alphaUpr) %>%
  mutate(temp = as.numeric(temp)) %>%
  filter(species == 'Anoplophora glabripennis' |
           species == 'Halyomorpha halys' |
           species == 'Aedes aegypti'|
           species == 'Anthonomus grandis' |
           species == 'Paracoccus marginatus' |
           species == 'Aphis gossypii' |
           species == 'Bemisia tabaci' |
           species == 'Tetraneura nigriabdominalis' |
           species == 'Stethorus punctillum' |
           species == 'Tetranychus mcdanieli' |
           species == 'Muscidifurax zaraptor' |
           species == 'Aphis nasturtii' |
           species == 'Rhopalosiphum maidis' |
           species == 'Thrips hawaiiensis' |
           species == 'Helicoverpa armigera' |
           species == 'Scapsipedus icipe' |
           species == 'Hyperaspis notata' |
           species == 'Toxoptera citricida'|
           species == 'Ophraella communa' |
           species == 'Bemisia argentifolii'|
           species == 'Dactylopius austrinus'|
           species == 'Aulacorthum solani') %>%
  mutate(curve_ID = eval(parse(text=spps))) %>%
  arrange(curve_ID)

alpha %>% distinct(species, curve_ID) %>% print(n=50)


zj  <- as_tibble(read.csv('data/zjPredictions.csv'))  %>% 
  select(species, temp, zj, zjLwr, zjUpr) %>%
  mutate(temp = as.numeric(temp)) %>%
  filter(species == 'Anoplophora glabripennis' |
           species == 'Halyomorpha halys' |
           species == 'Aedes aegypti'|
           species == 'Anthonomus grandis' |
           species == 'Paracoccus marginatus' |
           species == 'Aphis gossypii' |
           species == 'Bemisia tabaci' |
           species == 'Tetraneura nigriabdominalis' |
           species == 'Stethorus punctillum' |
           species == 'Tetranychus mcdanieli' |
           species == 'Muscidifurax zaraptor' |
           species == 'Aphis nasturtii' |
           species == 'Rhopalosiphum maidis' |
           species == 'Thrips hawaiiensis' |
           species == 'Helicoverpa armigera' |
           species == 'Scapsipedus icipe' |
           species == 'Hyperaspis notata' |
           species == 'Toxoptera citricida'|
           species == 'Ophraella communa' |
           species == 'Bemisia argentifolii'|
           species == 'Dactylopius austrinus'|
           species == 'Aulacorthum solani') %>%
  mutate(curve_ID = eval(parse(text=spps))) %>%
  arrange(curve_ID)


#________________

z  <- as_tibble(read.csv('data/zPredictions.csv'))  %>% 
  select(species, temp, z, zLwr, zUpr) %>%
  mutate(temp = as.numeric(temp)) %>%
  filter(species == 'Anoplophora glabripennis' |
           species == 'Halyomorpha halys' |
           species == 'Aedes aegypti'|
           species == 'Anthonomus grandis' |
           species == 'Paracoccus marginatus' |
           species == 'Aphis gossypii' |
           species == 'Bemisia tabaci' |
           species == 'Tetraneura nigriabdominalis' |
           species == 'Stethorus punctillum' |
           species == 'Tetranychus mcdanieli' |
           species == 'Muscidifurax zaraptor' |
           species == 'Aphis nasturtii' |
           species == 'Rhopalosiphum maidis' |
           species == 'Thrips hawaiiensis' |
           species == 'Helicoverpa armigera' |
           species == 'Scapsipedus icipe' |
           species == 'Hyperaspis notata' |
           species == 'Toxoptera citricida'|
           species == 'Ophraella communa' |
           species == 'Bemisia argentifolii'|
           species == 'Dactylopius austrinus'|
           species == 'Aulacorthum solani') %>%
  mutate(curve_ID = eval(parse(text=spps))) %>%
  arrange(curve_ID)

#_____________________________

bmax <- as_tibble(read.csv('data/BetaPredictions.csv'))  %>% 
  select(species, temp, bmax, bmaxLwr, bmaxUpr) %>%
  mutate(temp = as.numeric(temp)) %>%
  filter(species == 'Anoplophora glabripennis' |
           species == 'Halyomorpha halys' |
           species == 'Aedes aegypti'|
           species == 'Anthonomus grandis' |
           species == 'Paracoccus marginatus' |
           species == 'Aphis gossypii' |
           species == 'Bemisia tabaci' |
           species == 'Tetraneura nigriabdominalis' |
           species == 'Stethorus punctillum' |
           species == 'Tetranychus mcdanieli' |
           species == 'Muscidifurax zaraptor' |
           species == 'Aphis nasturtii' |
           species == 'Rhopalosiphum maidis' |
           species == 'Thrips hawaiiensis' |
           species == 'Helicoverpa armigera' |
           species == 'Scapsipedus icipe' |
           species == 'Hyperaspis notata' |
           species == 'Toxoptera citricida'|
           species == 'Ophraella communa' |
           species == 'Bemisia argentifolii'|
           species == 'Dactylopius austrinus'|
           species == 'Aulacorthum solani') %>%
  mutate(curve_ID = eval(parse(text=spps))) %>%
  arrange(curve_ID)

df1 <- inner_join(alpha,bmax)

df2 <- inner_join(z, zj)

df <-  inner_join(df1, df2)

df <- df %>% mutate(kappa = 0.01)

#_______________

# Define parameters

zj    <- df$zj
alpha <- df$alpha
z     <- df$z
bmax  <- df$bmax
k     <- df$kappa

# Calculate rmax
df <- df %>% mutate(rm_opt = (((k+z)*((log(bmax/(k+z)))-(alpha*zj)))/(alpha*(k+z)+1)))

# lower
zj_lwr    <- df$zjLwr
alpha_lwr <- df$alphaLwr
z_lwr     <- df$zLwr
bmax_lwr  <- df$bmaxLwr
k         <- df$kappa

df <- df %>% mutate(rm_optLwr = (((k+z_lwr)*((log(bmax_lwr/(k+z_lwr)))-(alpha_lwr*zj_lwr)))/(alpha_lwr*(k+z_lwr)+1)))

# upper 
zj_upr    <- df$zjUpr
alpha_upr <- df$alphaUpr
z_upr     <- df$zUpr
bmax_upr  <- df$bmaxUpr
k         <- df$kappa

df <- df %>% mutate(rm_optUpr = (((k+z_upr)*((log(bmax_upr/(k+z_upr)))-(alpha_upr*zj_upr)))/(alpha_upr*(k+z_upr)+1)))

write_csv(df, 'results/r_mCalcs.csv')


bodyMass <- as_tibble(read.csv("data/sizeMeans.csv")) %>%
  rename(species = interactor1, mass = avg, mass_sd = sd) %>%
  mutate(curve_ID = eval(parse(text=spps))) %>%
  filter(curve_ID != 'NA') %>% arrange(curve_ID)

head(bodyMass)


bodyMass <- bodyMass %>% select(species, mass, curve_ID) 

df <- inner_join(df, bodyMass)

write_csv(df, "results/rm_optSizeScaling.csv")


df <- df %>%
  mutate_at(vars(c(rm_opt)), 
            ~ifelse(rm_opt < -0.001, -0.001, .)) %>%
  mutate_at(vars(c(rm_optLwr)), 
            ~ifelse(rm_optLwr < -0.001, -0.001, .)) %>%
  mutate_at(vars(c(rm_optUpr)), 
            ~ifelse(rm_optUpr < -0.001, -0.001, .)) %>%
  filter(rm_opt > -0.001)


rmPlot <- ggplot()+
  geom_line(aes(temp, rm_opt), df)+
  facet_wrap(~species, ncol = 4, scales = "free_y")+
  scale_y_continuous(expression(plain(paste(" Maximal population growth rate ("~italic(r[m])~")"))))+
  theme_bw()+
  theme(strip.text = element_text(face = "italic"),
        text=element_text(family="Times"))+
  labs(x=expression(plain(paste(" Temperature, ",degree,"C"))))+
  theme(legend.position = 'none')

rmPlot

save_plot(rmPlot, file="results/SI/rmTPCs.pdf", 
          base_height=16, base_width = 20, base_asp = 1, units="cm")

#______________________________

rm_data <- read_csv('results/rm_optSizeScaling.csv') %>% 
  group_by(species) %>% 
  slice(which.max(rm_opt)) %>%
  arrange(curve_ID)

rm_data$kT <- 1/(8.617333262145 * 10^-5 * (rm_data$temp+273.15))

rm_dataA <- rm_data %>% filter(species != 'Anoplophora glabripennis') # remove outlier from regression

#plot uncorrected data in log-log scale
rm_dataA %>%
  ggplot(aes(x=log(mass), y = log(rm_opt)))+
  geom_point()+
  geom_smooth(method = 'lm')

#plot rm_opt vs T
rm_dataA %>%
  ggplot(aes(x = temp, y = rm_opt)) +
  geom_point()+
  geom_smooth(method = 'lm')

# linear model (note the allometry is linear in log-log scale)
rm_model <- lm(log(rm_opt) ~ log(mass) + kT, data = rm_dataA)
summary(rm_model)
coef(rm_model)

save(rm_model,file= "results/rm_MTE_model.Rdata")

cf <-  confint(rm_model,level = .95)
anova(rm_model)

r <- lm(log(rm_opt/mass^coef(rm_model)[2]) ~ temp, data = rm_dataA)
summary(r)

#plot rm_opt in 1/kT, correcting for mass
MassCorr_rm_opt <- 
  rm_dataA %>%
  ggplot(aes(x = temp, y = log(rm_opt/mass^coef(rm_model)[2]))) +
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="gray")+
  scale_y_continuous(expression(plain(paste("Log mass-corrected ",italic(r[m]),
                                            " at optimal temperature (ln(",
                                            italic(r[paste(m,",", opt)])," / ",
                                            italic(M^-0.10),"))"))),
                     limits=c(-2.8,-1.25),
                     expand = c(0.01, 0),
                     breaks=seq(-2.5,-1.5, by=0.5))+
  scale_x_continuous(expression(paste(italic(T[pk]))),
                     limits=c(20.25,33.6),
                     expand = c(0, 0),
                     breaks=seq(21,33, by=3))+
  geom_point(aes(shape = species, fill = species), size=2,stroke=0.25)+
  scale_shape_manual(values = c(21,21,21,21,21,21,
                                22,22,22,22,22,22,
                                23,23,23,23,23,23,
                                24,24,24),
                     name=expression(bold("")),
                     guide = guide_legend(nrow=7,ncol=3,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  scale_fill_manual(values = c("#EEEEEE","#CCCCCC","#999999",'#666666',"#333333","#000000",
                               "#EEEEEE","#CCCCCC","#999999",'#666666',"#333333","#000000",
                               "#EEEEEE","#CCCCCC","#999999",'#666666',"#333333","#000000",
                               "#EEEEEE","#CCCCCC","#999999"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow=7,ncol=3,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(values = c("#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333","#333333","#333333","#333333"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=7,ncol=3,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  theme_bw()+
  theme(text = element_text(size=14, family = 'Times'),
        legend.position = 'none',
        legend.text = element_text(size = 4.5, face = 'italic'),
        legend.background = element_rect(colour = "black", size = 0.125), 
        legend.margin=margin(t = 0.1, b = 0.1, r=0.1,l=0.1, unit='cm'),
        legend.key.size = unit(0.2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text(aes(x = -Inf, y = Inf, hjust = -0.7, vjust=2, family='Times',
                label = "B"), size = 5, colour = "black")+
  theme(axis.title.y = element_text(size = 10))

MassCorr_rm_opt

#________________________________________________________________
# Sum of trait Tpks versus mass corrected r_m, opt

Species <- rm_data$species

alphaTpks <- as_tibble(read.csv('data/alpha_Tpks_AllParams.csv')) %>% 
  filter(param=='topt') %>%
  filter(species %in% Species)%>%
  select(param,species,estimate, conf_lower, conf_upper,trait) %>%
  mutate(estimate = as.numeric(estimate), 
         conf_lower = as.numeric(conf_lower), conf_upper = as.numeric(conf_upper))

zjTpks <- as_tibble(read.csv('data/zj_Tpks_AllParams.csv')) %>% 
  filter(param=='topt') %>%
  filter(species %in% Species)%>%
  select(param,species,estimate, conf_lower, conf_upper,trait) %>%
  mutate(estimate = as.numeric(estimate), 
         conf_lower = as.numeric(conf_lower), 
         conf_upper = as.numeric(conf_upper))


zTpks <- as_tibble(read.csv('data/z_Tpks_AllParams.csv')) %>% 
  filter(param=='topt') %>%
  filter(species %in% Species)%>%
  select(param,species,estimate, conf_lower, conf_upper,trait) %>%
  mutate(estimate = as.numeric(estimate), 
         conf_lower = as.numeric(conf_lower), 
         conf_upper = as.numeric(conf_upper))

bmaxTpks <- as_tibble(read.csv('data/bmax_Tpks_AllParams.csv')) %>% 
  filter(param=='topt') %>%
  filter(species %in% Species)%>%
  select(param,species,estimate, conf_lower, conf_upper,trait) %>%
  mutate(estimate = as.numeric(estimate), 
         conf_lower = as.numeric(conf_lower), 
         conf_upper = as.numeric(conf_upper))

AllTpks <- bind_rows(alphaTpks,zjTpks,zTpks,bmaxTpks)

head(AllTpks)

# write_csv(AllTpks, '../results/AllTpkParams.csv')

load("results/rm_MTE_model.Rdata") # Load linear model fitted in previous step (04_rmCalcs)


#prepare data for plotting
rm_data <- rm_data %>% 
  mutate(rm_massCor = rm_opt/mass^coef(rm_model)[2])

#Calculate variance and Sum

OptVar <- AllTpks %>% group_by(species) %>%
  summarise(variance = var(estimate))
OptSum <- AllTpks %>% group_by(species) %>%
  summarise(sum = sum(estimate))

PlotData <- left_join(rm_data, OptVar, by="species")
PlotData <- left_join(PlotData, OptSum, by="species")
head(PlotData)

write_csv(PlotData,'results/rm_optSizeScaling_w_tpksums.csv')

# load simulation results, filtering by only r_opt values of rows with predicted optimal order of T_pk's 

sim_results_opt <- read.csv("results/sim_results.csv") %>% filter(Opt_T_pk_order == 1) 

sim_results_notopt <- read.csv("results/sim_results.csv") %>% filter(Opt_T_pk_order == 0) 

head(sim_results_opt)

SumTpks_plot <- ggplot() +
  geom_point(data = sim_results_notopt, aes(x = T_pk_sums, y = r_m_opt), color="lightyellow2", size = 3) +
  geom_point(data = sim_results_opt, aes(x = T_pk_sums, y = r_m_opt), color="gray", size = 3)+
  geom_point(data = PlotData, aes(x = sum, y = rm_massCor, shape=species, fill=species), 
             size=2, stroke=0.25)+ 
  scale_y_continuous(expression(plain(paste("Mass-corrected ",italic(r[m]),
                                            " at its optimal temperature"))))+
  scale_x_continuous(expression(plain(paste("Sum of ",italic("T"[pk]),"'s")))) +
  theme_bw()+
  theme(text=element_text(family="Times")) + 
  scale_shape_manual(values = c(21,21,21,21,21,21,
                                22,22,22,22,22,22,
                                23,23,23,23,23,23,
                                24,24,24,24),
                     name=expression(bold("")),
                     guide = guide_legend(nrow=7,ncol=3,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  scale_fill_manual(values = c("#EEEEEE","#CCCCCC","#999999",'#666666',"#333333","#000000",
                               "#EEEEEE","#CCCCCC","#999999",'#666666',"#333333","#000000",
                               "#EEEEEE","#CCCCCC","#999999",'#666666',"#333333","#000000",
                               "#EEEEEE","#CCCCCC","#999999",'#666666'),
                    name=expression(bold("")),
                    guide = guide_legend(nrow=7,ncol=3,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(values = c("#333333","#333333","#333333",
                                 "#333333","#333333","#333333",
                                 "#333333","#333333","#333333",
                                 "#333333","#333333","#333333",
                                 "#333333","#333333","#333333",
                                 "#333333", "#333333","#333333",
                                 "#333333","#333333","#333333","#333333"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=7,ncol=3,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+ 
  theme(text = element_text(size=14),
        legend.position = 'none',
        legend.text = element_text(size = 4.5, face = 'italic'),
        legend.background = element_rect(colour = "black", size = 0.125), 
        legend.margin=margin(t = 0.05, b = 0.1, r=0.1,l=0.1, unit='cm'),
        legend.key.size = unit(0.3, 'cm'))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  theme(axis.title.y = element_text(size=12))+
  geom_text(aes(x = -Inf, y = Inf, hjust = -0.7, vjust=2, family='Times',
                label = "A"), size = 5, colour = "black")

SumTpks_plot



# Create argument containing target species


# Read in the trait data
alphaMass <- as_tibble(read.csv('data/a_pksT_pksMass.csv')) %>% 
  select(species, a_pk, a_pkLwr, a_pkUpr) %>%
  filter(species == 'Anoplophora glabripennis' |
           species == 'Halyomorpha halys' |
           species == 'Aedes aegypti'|
           species == 'Anthonomus grandis' |
           species == 'Paracoccus marginatus' |
           species == 'Aphis gossypii' |
           species == 'Bemisia tabaci' |
           species == 'Tetraneura nigriabdominalis' |
           species == 'Stethorus punctillum' |
           species == 'Tetranychus mcdanieli' |
           species == 'Muscidifurax zaraptor' |
           species == 'Aphis nasturtii' |
           species == 'Rhopalosiphum maidis' |
           species == 'Thrips hawaiiensis' |
           species == 'Helicoverpa armigera' |
           species == 'Scapsipedus icipe' |
           species == 'Hyperaspis notata' |
           species == 'Toxoptera citricida'|
           species == 'Ophraella communa' |
           species == 'Bemisia argentifolii'|
           species == 'Dactylopius austrinus'|
           species == 'Aulacorthum solani') %>%
  mutate(curve_ID = eval(parse(text=spps))) %>%
  mutate(curve_ID = as.numeric(curve_ID)) %>%
  arrange(curve_ID)

alphaMass

alpharm_data <- inner_join(rm_data, alphaMass)
alpharm_data

# linear model (note the allometry is linear in log-log scale)
alpharm_model <- lm(log(rm_opt/mass^-0.1) ~ log(a_pk/mass^-0.27), data = alpharm_data)
summary(alpharm_model)
anova(alpharm_model)


## plot $a_pk$ vs $r_opt$

MassCorr_rm_opt_vs_AlphaTpk <-
  alpharm_data %>%
  ggplot(aes(x=log(a_pk/mass^-0.27), y = log(rm_opt/mass^-0.10)))+
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="gray")+
  scale_y_continuous(expression(plain(paste("ln(",italic(r[m][','][' '][opt])," / ",
                                            italic(M^-0.1),"))"))),
                     limits=c(-4,-0.75),
                     expand = c(0.01, 0),
                     breaks=seq(-4,-1, by=1))+
  scale_x_continuous(expression(plain(paste("ln(", italic(1/alpha[','][' '][pk])~")/",
                                            italic(M^-0.27),")"))),
                     limits=c(-4.35,-1.35),
                     expand = c(0, 0),
                     breaks=seq(-4,-2, by=1))+
  geom_point(aes(shape=species, fill=species), size=2, stroke=0.25)+
  scale_shape_manual(values = c(21,21,21,21,21,21,
                                22,22,22,22,22,22,
                                23,23,23,23,23,23,
                                24,24,24,24),
                     name=expression(bold("")),
                     guide = guide_legend(nrow=11,ncol=2,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  scale_fill_manual(values = c("#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF", "#EEEEEE","#CCCCCC","#999999"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow=11,ncol=2,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(values = c("#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=11,ncol=2,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  theme_bw()+
  theme(text=element_text(family="Times"))+
  theme(text = element_text(size=14),
        legend.text = element_text(size = 9, face = 'italic'),
        legend.background = element_rect(colour = "white", size = 0.125), 
        legend.margin=margin(t = 0.1, b = 0.1, r=0.1,l=0.1, unit='cm'),
        legend.key.size = unit(0.4, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.7, vjust=2, family='Times',
                label = "C"),size = 5, colour = "black")+
  theme(axis.title.y = element_blank(), legend.position = 'none')
  
MassCorr_rm_opt_vs_AlphaTpk

legend <- get_legend(MassCorr_rm_opt_vs_AlphaTpk + theme(legend.position = c(0.515,0.525)))

rm_calcs_Plot <- (SumTpks_plot+legend)/(MassCorr_rm_opt+MassCorr_rm_opt_vs_AlphaTpk)

rm_calcs_Plot

save_plot(rm_calcs_Plot, file="results/Figs/SumTpks.pdf", 
          base_height=18, base_width = 19.5, base_asp = 1, units="cm")


#___________________________________________________


# bmax (fecundity)

betaMass <- as_tibble(read_csv('data/bmaxT_pksMass.csv')) %>%
  select(species, bmax, bmaxLwr, bmaxUpr) %>%
  filter(species == 'Anoplophora glabripennis' |
           species == 'Halyomorpha halys' |
           species == 'Aedes aegypti'|
           species == 'Anthonomus grandis' |
           species == 'Paracoccus marginatus' |
           species == 'Aphis gossypii' |
           species == 'Bemisia tabaci' |
           species == 'Tetraneura nigriabdominalis' |
           species == 'Stethorus punctillum' |
           species == 'Tetranychus mcdanieli' |
           species == 'Muscidifurax zaraptor' |
           species == 'Aphis nasturtii' |
           species == 'Rhopalosiphum maidis' |
           species == 'Thrips hawaiiensis' |
           species == 'Helicoverpa armigera' |
           species == 'Scapsipedus icipe' |
           species == 'Hyperaspis notata' |
           species == 'Toxoptera citricida'|
           species == 'Ophraella communa' |
           species == 'Bemisia argentifolii'|
           species == 'Dactylopius austrinus'|
           species == 'Aulacorthum solani') %>%
  mutate(curve_ID = eval(parse(text=spps))) %>%
  mutate(curve_ID = as.numeric(curve_ID)) %>%
  rename(b_maxpk = bmax, b_maxpkLwr = bmaxLwr, b_maxpkUpr = bmaxUpr) %>%
  arrange(curve_ID)

betaMass

betarm_data <- inner_join(rm_data, betaMass)
betarm_data


# linear model (note the allometry is linear in log-log scale)
betarm_model <- lm(log(rm_opt/mass^-0.1) ~ log(bmax/mass^0.12), data = betarm_data)
summary(betarm_model)
anova(betarm_model)


MassCorr_rm_opt_vs_BetaTpk <-
  betarm_data %>%
  ggplot(aes(x=log(bmax/mass^0.12), y = log(rm_opt/mass^-0.1)))+
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="gray")+
  scale_y_continuous(expression(plain(paste("ln(",italic(r[m][','][' '][opt])," / ",
                                            italic(M^-0.1),"))"))),
                     limits=c(-4,-0.75),
                     expand = c(0.01, 0),
                     breaks=seq(-4,-1, by=1))+
  scale_x_continuous(expression(plain(paste("ln(", italic(b[max][','][' '][pk])~")/",
                                            italic(M^0.12),")"))))+
  geom_point(aes(shape=species, fill=species), size=2, stroke=0.25)+
  scale_shape_manual(values = c(21,21,21,21,21,21,
                                22,22,22,22,22,22,
                                23,23,23,23,23,23,
                                24,24,24,24),
                     name=expression(bold("")),
                     guide = guide_legend(nrow=11,ncol=2,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  scale_fill_manual(values = c("#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF", "#EEEEEE","#CCCCCC","#999999"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow=11,ncol=2,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(values = c("#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=11,ncol=2,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  theme_bw()+
  theme(text=element_text(family="Times"))+
  theme(text = element_text(size=14),
        legend.position = 'none',
        legend.text = element_text(size = 9, face = 'italic'),
        legend.background = element_rect(colour = "white", size = 0.125), 
        legend.margin=margin(t = 0.1, b = 0.1, r=0.1,l=0.1, unit='cm'),
        legend.key.size = unit(0.4, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.7, vjust=2,
                label = "C"),size = 5, colour = "black", family='Times')

MassCorr_rm_opt_vs_BetaTpk

#_________________________________________

# z (adult mortality rate)

zMass <- as_tibble(read_csv('data/z_pksT_pksMass.csv')) %>% 
  select(species, zpk, zpkLwr, zpkUpr) %>%
  filter(species == 'Anoplophora glabripennis' |
           species == 'Halyomorpha halys' |
           species == 'Aedes aegypti'|
           species == 'Anthonomus grandis' |
           species == 'Paracoccus marginatus' |
           species == 'Aphis gossypii' |
           species == 'Bemisia tabaci' |
           species == 'Tetraneura nigriabdominalis' |
           species == 'Stethorus punctillum' |
           species == 'Tetranychus mcdanieli' |
           species == 'Muscidifurax zaraptor' |
           species == 'Aphis nasturtii' |
           species == 'Rhopalosiphum maidis' |
           species == 'Thrips hawaiiensis' |
           species == 'Helicoverpa armigera' |
           species == 'Scapsipedus icipe' |
           species == 'Hyperaspis notata' |
           species == 'Toxoptera citricida'|
           species == 'Ophraella communa' |
           species == 'Bemisia argentifolii'|
           species == 'Dactylopius austrinus'|
           species == 'Aulacorthum solani') %>%
  mutate(curve_ID = eval(parse(text=spps))) %>%
  mutate(curve_ID = as.numeric(curve_ID)) %>%
  arrange(curve_ID)

zMass

zrm_data <- inner_join(rm_data, zMass)

head(zrm_data)


# linear model (note the allometry is linear in log-log scale)
zrm_model <- lm(log(rm_opt/mass^-0.1) ~ log(zpk/mass^-0.09), data = zrm_data)
summary(zrm_model)
anova(zrm_model)



## plot $z$ vs $r_opt$

MassCorr_rm_opt_vs_zTpk <-
  zrm_data %>%
  ggplot(aes(x=log(zpk/mass^-0.09), y = log(rm_opt/mass^-0.1)))+
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="gray")+
  scale_y_continuous(expression(plain(paste("ln(",italic(r[m][','][' '][opt])," / ",
                                            italic(M^-0.1),"))"))),
                     limits=c(-4,-0.75),
                     expand = c(0.01, 0),
                     breaks=seq(-4,-1, by=1))+
  scale_x_continuous(expression(plain(paste("ln(", italic(z[','][' '][pk])~")/",
                                            italic(M^-0.09),")"))))+
  geom_point(aes(shape=species, fill=species), size=2, stroke=0.25)+
  scale_shape_manual(values = c(21,21,21,21,21,21,
                                22,22,22,22,22,22,
                                23,23,23,23,23,23,
                                24,24,24,24),
                     name=expression(bold("")),
                     guide = guide_legend(nrow=11,ncol=2,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  scale_fill_manual(values = c("#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF", "#EEEEEE","#CCCCCC","#999999"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow=11,ncol=2,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(values = c("#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=11,ncol=2,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  theme_bw()+
  theme(text=element_text(family="Times"))+
  theme(text = element_text(size=14),
        legend.position = 'none',
        legend.text = element_text(size = 9, face = 'italic'),
        legend.background = element_rect(colour = "white", size = 0.125), 
        legend.margin=margin(t = 0.1, b = 0.1, r=0.1,l=0.1, unit='cm'),
        legend.key.size = unit(0.4, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text(aes(x = -Inf, y = Inf,hjust = -0.7, vjust=2,
                label = "B"),size = 5, colour = "black", family = 'Times') 


MassCorr_rm_opt_vs_zTpk


# $z_J$ (juvenile mortality rate)

zJMass <- as_tibble(read_csv('data/zj_pksT_pksMass.csv')) %>% 
  select(species, zjpk, zjpkLwr, zjpkUpr) %>%
  filter(species == 'Anoplophora glabripennis' |
           species == 'Halyomorpha halys' |
           species == 'Aedes aegypti'|
           species == 'Anthonomus grandis' |
           species == 'Paracoccus marginatus' |
           species == 'Aphis gossypii' |
           species == 'Bemisia tabaci' |
           species == 'Tetraneura nigriabdominalis' |
           species == 'Stethorus punctillum' |
           species == 'Tetranychus mcdanieli' |
           species == 'Muscidifurax zaraptor' |
           species == 'Aphis nasturtii' |
           species == 'Rhopalosiphum maidis' |
           species == 'Thrips hawaiiensis' |
           species == 'Helicoverpa armigera' |
           species == 'Scapsipedus icipe' |
           species == 'Hyperaspis notata' |
           species == 'Toxoptera citricida'|
           species == 'Ophraella communa' |
           species == 'Bemisia argentifolii'|
           species == 'Dactylopius austrinus'|
           species == 'Aulacorthum solani') %>%
  mutate(curve_ID = eval(parse(text=spps))) %>%
  mutate(curve_ID = as.numeric(curve_ID)) %>%
  arrange(curve_ID)

zJMass

zJrm_data <- inner_join(rm_data, zJMass)

head(zJrm_data)

# linear model (note the allometry is linear in log-log scale)
zJrm_model <- lm(log(rm_opt/mass^-0.1) ~ log(zjpk/mass^-0.22), data = zJrm_data)
summary(zJrm_model)
anova(zJrm_model)

MassCorr_rm_opt_vs_zJTpk <-
  zJrm_data %>%
  ggplot(aes(x=log(zjpk/mass^-0.22), y = log(rm_opt/mass^-0.1)))+
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="gray")+
  scale_y_continuous(expression(plain(paste("ln(",italic(r[m][','][' '][opt])," / ",
                                            italic(M^-0.1),"))"))),
                     limits=c(-4,-0.75),
                     expand = c(0.01, 0),
                     breaks=seq(-4,-1, by=1))+
  scale_x_continuous(expression(plain(paste("ln(", italic(z[J][','][' '][pk])~")/",
                                            italic(M^-0.22),")"))))+
  geom_point(aes(shape=species, fill=species), size=2, stroke=0.25)+
  scale_shape_manual(values = c(21,21,21,21,21,21,
                                22,22,22,22,22,22,
                                23,23,23,23,23,23,
                                24,24,24,24),
                     name=expression(bold("")),
                     guide = guide_legend(nrow=11,ncol=2,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  scale_fill_manual(values = c("#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF", "#EEEEEE","#CCCCCC","#999999"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow=11,ncol=2,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(values = c("#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=11,ncol=2,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  theme_bw()+
  theme(text=element_text(family="Times"))+
  theme(text = element_text(size=14),
        legend.position = 'none',
        legend.text = element_text(size = 9, face = 'italic'),
        legend.background = element_rect(colour = "white", size = 0.125), 
        legend.margin=margin(t = 0.1, b = 0.1, r=0.1,l=0.1, unit='cm'),
        legend.key.size = unit(0.4, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text(aes(x = -Inf, y = Inf, hjust = -0.7, vjust=2,
                label = "A"), size = 5, colour = "black", family='Times')


MassCorr_rm_opt_vs_zJTpk


legend <- get_legend(MassCorr_rm_opt_vs_zJTpk + theme(legend.position = c(0.515,0.525)))

rm_vs_OtherTraitPks <- MassCorr_rm_opt_vs_zJTpk+legend+
  MassCorr_rm_opt_vs_zTpk+
  MassCorr_rm_opt_vs_BetaTpk+plot_layout(nrow = 2, byrow = TRUE)

rm_vs_OtherTraitPks


save_plot(rm_vs_OtherTraitPks, file= "results/SI/rm_vs_OtherTraitPks.pdf",
          base_height=16, base_width = 17.5, base_asp = 1, units="cm")



VarLM <-   lm(PlotData$rm_massCor ~ PlotData$variance)

summary(VarLM)


Var_rm <- ggplot(PlotData, aes(x=variance, y=rm_massCor, colour = species, fill = species)) +
  scale_x_continuous(expression(plain(paste("Variance of ", italic(T)[pk],'s'))),
                     limits =c(19,126),
                     expand = c(0, 0),
                     breaks=seq(20,120, by=30))+
  scale_y_continuous(expression(plain(paste("Mass-corrected ",italic(r[m]),
                                            " at optimal temperature (ln(",italic(r[m][','][' '][opt]),
                                            " / ", italic(M^-0.1),"))"))))+
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="gray")+
  geom_point(aes(shape=species, 
                 fill=species), size=2,stroke=0.25)+
  scale_shape_manual(values = c(21,21,21,21,21,21,
                                22,22,22,22,22,22,
                                23,23,23,23,23,23,
                                24,24,24,24),
                     name=expression(bold("")),
                     guide = guide_legend(nrow=11,ncol=2,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  scale_fill_manual(values = c("#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF", "#EEEEEE","#CCCCCC","#999999"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow=11,ncol=2,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(values = c("#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=11,ncol=2,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  theme_bw()+
  theme(text=element_text(family="Times"))+
  theme(text = element_text(size=11),
        legend.position = 'none',
        legend.text = element_text(size = 9, face = 'italic'),
        legend.background = element_rect(colour = "black", size = 0.125), 
        legend.margin=margin(t = 0.05, b = 0.1, r=0.1,l=0.1, unit='cm'),
        legend.key.size = unit(0.3, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=9, hjust=1))+
  geom_text(aes(x = -Inf, y = Inf, hjust = -0.7, vjust=2,
                label = "A"),size = 5, colour = "black", family='Times')


Var_rm



Sum_vs_Var_Tpks_plot <- PlotData %>%
  ggplot(aes(x = sum, y = variance))+
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="gray")+
  scale_y_continuous(expression(plain(paste("Variance of ",italic("T"[pk]),"s"))))+
  scale_x_continuous(expression(plain(paste("Sum of ",italic("T"[pk]),"s"))))+
  geom_point(aes(shape=species, fill=species), size=2, stroke=0.25)+
  scale_shape_manual(values = c(21,21,21,21,21,21,
                                22,22,22,22,22,22,
                                23,23,23,23,23,23,
                                24,24,24,24),
                     name=expression(bold("")),
                     guide = guide_legend(nrow=11,ncol=2,
                                          direction = "vertical",
                                          title.position = "top",
                                          title.hjust=0.5))+
  scale_fill_manual(values = c("#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF","#EEEEEE","#CCCCCC","#999999",'#666666',"#000000",
                               "#FFFFFF", "#EEEEEE","#CCCCCC","#999999"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow=11,ncol=2,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_colour_manual(values = c("#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333","#333333","#333333","#333333","#333333","#333333","#333333",
                                 "#333333"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow=11,ncol=2,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  theme_bw()+
  theme(text=element_text(family="Times"))+
  theme(text = element_text(size=11),
        legend.position = 'none',
        legend.text = element_text(size = 4.5, face = 'italic'),
        legend.background = element_rect(colour = "white", size = 0.125), 
        legend.margin=margin(t = 0.05, b = 0.1, r=0.1,l=0.1, unit='cm'),
        legend.key.size = unit(0.3, 'cm'))+
  geom_text(aes(x = -Inf, y = Inf, hjust = -0.7, vjust=2,
                label = "B"), size = 5, colour = "black", family='Times')+
  theme(axis.title.y = element_text(hjust=0.5), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

Sum_vs_Var_Tpks_plot


legendz <- get_legend(Sum_vs_Var_Tpks_plot + theme(legend.position = c(0.515,0.525), 
                                                   legend.text = element_text(size = 8, face = 'italic')))

rm_vs_Var_and_Sum <- (Var_rm+Sum_vs_Var_Tpks_plot+legendz)
rm_vs_Var_and_Sum



save_plot(rm_vs_Var_and_Sum, file= "results/SI/rm_vs_Var_and_Sum.pdf",
          base_height=8.5, base_width = 25, base_asp = 1, units="cm")




##### Relationship of $T_{pk}$ of $\alpha$ vs latitude #####

alphaLat <- as_tibble(read_csv('data/TraitData.csv')) %>%
  rename(species = interactor1) %>%
  filter(standardisedtraitname == '1/alpha' & latitude != 'NA') %>%
  distinct(species, latitude)

alphaMass <- as_tibble(read_csv('data/a_pksT_pksMass.csv'))  %>% 
  select(-curve_ID) 


head(alphaLat)

head(alphaMass)

alphaLat_data <- alphaMass %>% inner_join(alphaLat) %>%
  mutate(a_pkmassCor = a_pk/mass^-0.27)

head(alphaLat_data)


a_T_pklat_plot <-
  alphaLat_data %>%
  ggplot(aes(x= abs(latitude), y = T_pk))+
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="gray")+
  scale_y_continuous(expression(paste(italic(T[pk]), " of juvenile development time (", italic(alpha),")")))+
  scale_x_continuous(expression(plain(paste("Absolute latitude"))))+
  geom_linerange(aes(ymin =T_pkLwr , ymax = T_pkUpr),size= 0.1,col="#000000")+
  theme_bw()+
  geom_point(size = 2.5,stroke=0.2, col = '#000000', shape=24, fill ='#e66101')+
  theme(text=element_text(family="Times"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text(aes(x = 0, y = 37.5,label = "A"), 
            parse = TRUE, size = 4, colour = "black", family='Times')

a_T_pklat_plot

a_T_pklatmodel <- lm( T_pk ~ latitude, data = alphaLat_data)
summary(a_T_pklatmodel)
anova(a_T_pklatmodel)

##### ## Relationship of $T_{pk}$ of $b_{max}$ vs latitude #####

bmaxLat <- as_tibble(read_csv('data/TraitData.csv')) %>%
  rename(species = interactor1) %>%
  filter(standardisedtraitname == 'bmax' & latitude != 'NA') %>%
  distinct(species, latitude)

bmaxMass <- as_tibble(read_csv('data/bmaxT_pksMass.csv'))  %>% 
  select(-curve_ID) 


head(bmaxLat)

head(bmaxMass)

bmaxLat_data <- bmaxMass %>% inner_join(bmaxLat) %>%
  mutate(bmaxmassCor = bmax/mass^0.12)

head(bmaxLat_data)


bmax_T_pklat_plot <-
  bmaxLat_data %>%
  ggplot(aes(x= abs(latitude), y = T_pk))+
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="gray")+
  geom_linerange(aes(ymin =T_pkLwr , ymax = T_pkUpr),size= 0.1,col="#000000")+
  scale_y_continuous(expression(paste(italic(T[pk]), " of fecundity rate (", italic(b[max]),")")))+
  scale_x_continuous(expression(plain(paste("Absolute latitude"))))+
  geom_point(size = 2.5,stroke=0.2, col = '#000000', fill="#fdb863", shape=23)+
  theme_bw()+
  theme(text=element_text(family="Times"))+
  theme(legend.position = 'none',legend.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text(aes(x = 0, y = 35, label = "D"), 
            parse = TRUE, size = 4, colour = "black",family='Times')


bmax_T_pklat_plot

bmax_T_pklatmodel <- lm( T_pk ~ latitude, data = bmaxLat_data)
summary(bmax_T_pklatmodel)
anova(bmax_T_pklatmodel)


##### Relationship of $T_{pk}$ of $z_J$ vs latitude #####

zJLat <- as_tibble(read_csv('data/TraitData.csv')) %>%
  rename(species = interactor1) %>%
  filter(standardisedtraitname == 'zj' & latitude != 'NA') %>%
  distinct(species, latitude)

zJMass <- as_tibble(read_csv('data/zj_pksT_pksMass.csv'))  %>% 
  select(-curve_ID)

head(zJLat)

head(zJMass)

zJLat_data <- zJMass %>% inner_join(zJLat) %>%
  mutate(z_JmassCor = zjpk/mass^-0.22)

head(zJLat_data)

zJ_T_pklatmodel <- lm( T_pk ~ latitude, data = zJLat_data)
summary(zJ_T_pklatmodel)
anova(zJ_T_pklatmodel)


zJ_T_pklat_plot <-
  zJLat_data %>%
  ggplot(aes(x= abs(latitude), y = T_pk)) +
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="gray")+
  scale_y_continuous(expression(paste(italic(T[pk]), " of juvenile mortality rate (", italic(z[J]),")")))+
  scale_x_continuous(expression(plain(paste("Absolute latitude"))))+
  geom_linerange(aes(ymin =T_pkLwr , ymax = T_pkUpr),size= 0.1,col="#000000")+
  geom_point(size = 2.5,stroke=0.2, col = '#000000', shape = 21, fill="#1f78b4")+
  theme_bw()+
  theme(legend.position = 'none',legend.text = element_text(size = 10))+
  theme(text=element_text(family="Times"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text(aes(x = 0, y = 33,label = "B"), 
            parse = TRUE, size = 4, colour = "black",family='Times')

zJ_T_pklat_plot



##### Relationship of $T_{pk}$ of $z$ vs latitude #####

zLat <- as_tibble(read_csv('data/TraitData.csv')) %>%
  rename(species = interactor1) %>%
  filter(standardisedtraitname == 'z' & latitude != 'NA') %>%
  distinct(species, latitude)

zMass <- as_tibble(read_csv('data/z_pksT_pksMass.csv'))  %>% 
  select(-curve_ID)


head(zLat)

head(zMass)

zLat_data <- zMass %>% inner_join(zLat) %>%
  mutate(zmassCor = zpk/mass^-0.17)

head(zLat_data)


z_T_pklatmodel <- lm( T_pk ~ latitude, data = zLat_data)
summary(z_T_pklatmodel)
anova(z_T_pklatmodel)


z_T_pklat_plot <-
  zLat_data %>%
  ggplot(aes(x= abs(latitude), y = T_pk))+
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="gray")+
  scale_y_continuous(expression(paste(italic(T[pk]), " of adult mortality rate (", italic(z),")")))+
  scale_x_continuous(expression(plain(paste("Absolute latitude"))))+
  geom_linerange(aes(ymin =T_pkLwr , ymax = T_pkUpr),size= 0.1,col="#000000")+
  geom_point(size = 2.5,stroke=0.2, col = '#000000', shape=22, fill ="#a6cee3")+
  theme_bw()+
  theme(legend.position = 'none',legend.text = element_text(size = 10))+
  theme(text=element_text(family="Times"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text(aes(x = 0, y = 28,label = "C"), 
            parse = TRUE, size = 4, colour = "black", family='Times')

z_T_pklat_plot


p1 <- a_T_pklat_plot+zJ_T_pklat_plot+z_T_pklat_plot+bmax_T_pklat_plot

p1

save_plot(p1, file= "results/SI/Traits_Lat_Adapt.pdf",
          base_height=15, base_width = 17, base_asp = 1, units="cm")



#### Relationship of $T_{pk}$ of $\alpha$ vs Rearing Temperature #### 

RearDF <- as_tibble(read_csv('data/TraitData.csv')) %>% 
  select(standardisedtraitname,interactor1,interactor1growthtemp) %>%
  filter(interactor1growthtemp != 'not stated' & 
         interactor1growthtemp != 'various' & 
         interactor1growthtemp != 'NA' ) %>%
  rename(trait = standardisedtraitname, species = interactor1, RearTemp = interactor1growthtemp) %>%
  mutate(RearTemp = as.numeric(RearTemp))

RearAlpha <- RearDF %>% filter(trait == '1/alpha') %>% 
    group_by(species) %>% 
    summarise(avg = mean(RearTemp)) %>%
    arrange(avg) %>% 
  mutate(curve_ID = case_when(species ==  'Aulacorthum solani' ~ '1',        
  species ==  'Macrosiphum euphorbia' ~ '2',       
  species == 'Myzus persicae' ~ '3',             
  species == 'Drosophila suzukii' ~ '4',         
  species == 'Harmonia axyridis' ~ '5',          
  species == 'Otiorhynchus sulcatus'  ~ '6',     
  species == 'Aphis nasturtii'  ~ '7',      
  species == 'Gastrolina depressa' ~ '8',        
  species == 'Stethorus punctillum'~ '9',       
  species == 'Tetranychus mcdanieli'~ '10',      
  species == 'Thrips hawaiiensis'  ~ '11',       
  species == 'Toxoptera citricida' ~ '12',       
  species == 'Aedes krombeini' ~ '13',          
  species == 'Amblyseius swirskii' ~ '14',    
  species == 'Antestiopsis thunbergii' ~ '15',
  species == 'Helicoverpa armigera' ~ '16',
  species == 'Iphiseius degenerans' ~ '17',      
  species == 'Kampimodromus aberrans' ~ '18',    
  species == 'Leptinotarsa decemlineata' ~ '19',  
  species == 'Macrocentrus iridescens' ~ '20',    
  species == 'Monochamus leuconotus' ~ '21',      
  species == 'Muscidifurax zaraptor' ~ '22',      
  species == 'Paracoccus marginatus' ~ '23',      
  species == 'Phenacoccus solenopsis' ~ '24',    
  species == 'Rhopalosiphum maidis' ~ '25',       
  species == 'Tetraneura nigriabdominalis'~ '26',
  species == 'Aedes albopictus' ~ '27',        
  species == 'Anopheles gambiae' ~ '28',          
  species == 'Aubeonymus mariaefranciscae' ~ '29',
  species == 'Planococcus citri' ~ '30',          
  species == 'Amblyseius womersleyi' ~ '31',      
  species == 'Ophraella communa'~ '32',          
  species == 'Tetranychus urticae' ~ '33',        
  species == 'Anthonomus grandis'~ '34',         
  species == 'Hyperaspis notata' ~ '35',          
  species == 'Tribolium castaneum'~ '36',
  species == 'Scapsipedus icipe'~ '37')) %>%
    arrange(curve_ID) %>% 
    rename(RTemp = avg)
  
head(RearAlpha)


alpha <- as_tibble(read.csv('data/alpha_Tpks_AllParams.csv', header = TRUE)) %>%
  filter(param == 'topt') %>% 
  select(species, estimate, conf_lower, conf_upper, trait) %>%
  mutate(curve_ID = case_when(species ==  'Aulacorthum solani' ~ '1',        
                               species ==  'Macrosiphum euphorbia' ~ '2',       
                               species == 'Myzus persicae' ~ '3',             
                               species == 'Drosophila suzukii' ~ '4',         
                               species == 'Harmonia axyridis' ~ '5',          
                               species == 'Otiorhynchus sulcatus'  ~ '6',     
                               species == 'Aphis nasturtii'  ~ '7',      
                               species == 'Gastrolina depressa' ~ '8',        
                               species ==  'Stethorus punctillum'~ '9',       
                               species == 'Tetranychus mcdanieli'~ '10',      
                               species == 'Thrips hawaiiensis'  ~ '11',       
                               species == 'Toxoptera citricida' ~ '12',       
                               species == 'Aedes krombeini' ~ '13',          
                               species == 'Amblyseius swirskii' ~ '14',    
                               species == 'Antestiopsis thunbergii' ~ '15',
                               species == 'Helicoverpa armigera' ~ '16',
                               species == 'Iphiseius degenerans' ~ '17',      
                               species == 'Kampimodromus aberrans' ~ '18',    
                               species == 'Leptinotarsa decemlineata' ~ '19',  
                               species == 'Macrocentrus iridescens' ~ '20',    
                               species == 'Monochamus leuconotus' ~ '21',      
                               species == 'Muscidifurax zaraptor' ~ '22',      
                               species == 'Paracoccus marginatus' ~ '23',      
                               species == 'Phenacoccus solenopsis' ~ '24',    
                               species == 'Rhopalosiphum maidis' ~ '25',       
                               species == 'Tetraneura nigriabdominalis'~ '26',
                               species == 'Aedes albopictus' ~ '27',        
                               species == 'Anopheles gambiae' ~ '28',          
                               species == 'Aubeonymus mariaefranciscae' ~ '29',
                               species == 'Planococcus citri' ~ '30',          
                               species == 'Amblyseius womersleyi' ~ '31',      
                               species == 'Ophraella communa'~ '32',          
                               species == 'Tetranychus urticae' ~ '33',        
                               species == 'Anthonomus grandis'~ '34',         
                               species == 'Hyperaspis notata' ~ '35',          
                               species == 'Tribolium castaneum'~ '36',
                               species == 'Scapsipedus icipe'~ '37')) %>%
                arrange(curve_ID) %>% filter(curve_ID != 'NA')


RearAlpha <- inner_join(alpha,RearAlpha)

head(RearAlpha)

RAlp_model <- lm(estimate ~ RTemp, RearAlpha); summary(RAlp_model)


Rear_alp <- ggplot(RearAlpha, aes(x=RTemp, y=estimate)) +
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="gray")+
  geom_linerange(aes(ymin =conf_lower , ymax = conf_upper),size= 0.1,col="#000000") +
  geom_point(size = 2.5,stroke=0.2, col = '#000000', shape=24, fill ='#e66101')+
  theme_bw()+
  scale_x_continuous(expression(plain(paste("Rearing Temperature (",degree,"C)"))),
                     limits =c(19.8,30.2),
                     expand = c(0, 0),
                     breaks=seq(20,30, by=2))+
  scale_y_continuous(expression(paste(italic(T[pk]), " of juvenile development time (", italic(alpha),")")),
                     limits =c(21,39),
                     expand = c(0, 0),
                     breaks=seq(21,39, by=3))+
  theme(legend.position = 'none',legend.text = element_text(size = 10))+
  theme(text=element_text(family="Times"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text(aes(x = 20.5, y = 38,label = "A"), 
            parse = TRUE, size = 4, colour = "black", family='Times')

Rear_alp


Rearbmax <- RearDF %>% filter(trait == 'bmax') %>% 
  group_by(species) %>% 
  summarise(avg = mean(RearTemp)) %>%
  arrange(avg) %>%
  mutate(curve_ID = case_when(species == 'Aulacorthum solani' ~ '1',
                               species == 'Plutella xylostella' ~ '3',
                               species == 'Aphis nasturtii' ~ '4',
                               species == 'Stethorus punctillum' ~ '5',
                               species == 'Thrips hawaiiensis' ~ '6',
                               species == 'Toxoptera citricida' ~ '7',
                               species == 'Acyrthosiphon pisum' ~ '8',
                               species == 'Aedes krombeini' ~ '9',
                               species == 'Antestiopsis thunbergii' ~ '10',
                               species == 'Helicoverpa armigera' ~ '11',
                               species == 'Hylobius transversovittatus' ~ '12',
                               species == 'Iphiseius degenerans' ~ '13',
                               species == 'Kampimodromus aberrans' ~ '14',
                               species == 'Monochamus leuconotus' ~ '15',
                               species == 'Muscidifurax zaraptor' ~ '16',
                               species == 'Paracoccus marginatus' ~ '17',
                               species == 'Phenacoccus solenopsis' ~ '18',
                               species == 'Rhopalosiphum maidis' ~ '19',
                               species == 'Tetraneura nigriabdominalis' ~ '20',
                               species == 'Tetranychus mcdanieli' ~ '21',
                               species == 'Planococcus citri' ~ '22',    
                               species == 'Amblyseius womersleyi' ~ '23',      
                               species == 'Ophraella communa'  ~ '24',    
                               species == 'Anthonomus grandis' ~ '25',  
                               species == 'Callosobruchus analis'  ~ '26', 
                               species == 'Callosobruchus chinensis' ~ '27',
                               species == 'Callosobruchus maculatus' ~ '28',
                               species == 'Callosobruchus rhodesianus' ~ '29',
                               species == 'Hyperaspis notata' ~ '30',
                               species == 'Scapsipedus icipe' ~ '31')) %>%
  arrange(curve_ID) %>% rename(species = species, RTemp = avg) %>%
  filter(curve_ID != 'NA')


bmax  <- as_tibble(read.csv('data/bmax_Tpks_AllParams.csv', header = TRUE)) %>%
  filter(param == 'topt') %>% 
  select(species, estimate, conf_lower, conf_upper, trait) %>%
  mutate(curve_ID = case_when(species == 'Aulacorthum solani' ~ '1',
                              species == 'Plutella xylostella' ~ '3',
                              species == 'Aphis nasturtii' ~ '4',
                              species == 'Stethorus punctillum' ~ '5',
                              species == 'Thrips hawaiiensis' ~ '6',
                              species == 'Toxoptera citricida' ~ '7',
                              species == 'Acyrthosiphon pisum' ~ '8',
                              species == 'Aedes krombeini' ~ '9',
                              species == 'Antestiopsis thunbergii' ~ '10',
                              species == 'Helicoverpa armigera' ~ '11',
                              species == 'Hylobius transversovittatus' ~ '12',
                              species == 'Iphiseius degenerans' ~ '13',
                              species == 'Kampimodromus aberrans' ~ '14',
                              species == 'Monochamus leuconotus' ~ '15',
                              species == 'Muscidifurax zaraptor' ~ '16',
                              species == 'Paracoccus marginatus' ~ '17',
                              species == 'Phenacoccus solenopsis' ~ '18',
                              species == 'Rhopalosiphum maidis' ~ '19',
                              species == 'Tetraneura nigriabdominalis' ~ '20',
                              species == 'Tetranychus mcdanieli' ~ '21',
                              species == 'Planococcus citri' ~ '22',    
                              species == 'Amblyseius womersleyi' ~ '23',      
                              species == 'Ophraella communa'  ~ '24',    
                              species == 'Anthonomus grandis' ~ '25',  
                              species == 'Callosobruchus analis'  ~ '26', 
                              species == 'Callosobruchus chinensis' ~ '27',
                              species == 'Callosobruchus maculatus' ~ '28',
                              species == 'Callosobruchus rhodesianus' ~ '29',
                              species == 'Hyperaspis notata' ~ '30',
                              species == 'Scapsipedus icipe' ~ '31')) %>%
  arrange(curve_ID) %>% 
  filter(curve_ID != 'NA')

Rearbmax <- inner_join(bmax, Rearbmax)

Rbmax_model <- lm(estimate ~ RTemp, Rearbmax); summary(Rbmax_model)


Rear_bmax <- ggplot(Rearbmax, aes(x=RTemp, y=estimate)) +
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="gray") +
  geom_linerange(aes(ymin =conf_lower , ymax = conf_upper),size= 0.1,col="#000000") +
  geom_point(size = 2.5,stroke=0.2, col = '#000000', fill="#fdb863", shape=23)+
  theme_bw()+
  scale_x_continuous(expression(plain(paste("Rearing Temperature (",degree,"C)"))),
                     limits =c(19.8,30.2),
                     expand = c(0, 0),
                     breaks=seq(20,30, by=2))+
  scale_y_continuous(expression(paste(italic(T[pk]), " of fecundity rate (", italic(b[max]),")")),
                     limits =c(17.5,38),
                     expand = c(0, 0),
                     breaks=seq(20,35, by=5))+
  theme(legend.position = 'none',legend.text = element_text(size = 10))+
  theme(text=element_text(family="Times"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text(aes(x = 20.23, y = 37, label = "B"), 
            parse = TRUE, size = 4, colour = "black", family='Times')

Rear_bmax


RearZetaJ <- RearDF %>% filter(trait == 'zj') %>% 
  group_by(species) %>% 
  summarise(avg = mean(RearTemp)) %>%
  arrange(avg) %>% 
  mutate(curve_ID = case_when(species == 'Aulacorthum solani' ~ '1',
                               species == 'Harmonia axyridis' ~ '2',
                               species == 'Aphis nasturtii' ~ '3',
                               species == 'Stethorus punctillum' ~ '4',
                               species == 'Tetranychus mcdanieli' ~ '5',
                               species == 'Thrips hawaiiensis' ~ '6',
                               species == 'Toxoptera citricida' ~ '7',
                               species == 'Acyrthosiphon pisum' ~ '8',
                               species == 'Aedes albopictus' ~ '9',
                               species == 'Helicoverpa armigera' ~ '10',
                               species == 'Muscidifurax zaraptor' ~ '11',
                               species == 'Paracoccus marginatus' ~ '12',
                               species == 'Rhopalosiphum maidis' ~ '13',
                               species == 'Tetraneura nigriabdominalis' ~ '14',
                               species == 'Planococcus citri' ~ '15',
                               species == 'Anopheles gambiae' ~ '16',
                               species == 'Aubeonymus mariaefranciscae' ~ '17',
                               species == 'Amblyseius womersleyi' ~ '18',
                               species == 'Ophraella communa' ~ '19',
                               species == 'Anthonomus grandis' ~ '20',
                               species == 'Hyperaspis notata' ~ '21',
                               species == 'Tetranychus urticae' ~ '22',
                               species == 'Scapsipedus icipe' ~ '23',
                               species == 'Tribolium castaneum' ~ '24')) %>%
  arrange(curve_ID) %>% rename(species = species, RTemp = avg) %>%
  filter(curve_ID != 'NA')

zetaJ <- as_tibble(read.csv('data/zj_Tpks_AllParams.csv', header = TRUE)) %>%
  filter(param == 'topt') %>% 
  select(species, estimate, conf_lower, conf_upper, trait) %>%
  mutate(curve_ID = case_when(species == 'Aulacorthum solani' ~ '1',
                              species == 'Harmonia axyridis' ~ '2',
                              species == 'Aphis nasturtii' ~ '3',
                              species == 'Stethorus punctillum' ~ '4',
                              species == 'Tetranychus mcdanieli' ~ '5',
                              species == 'Thrips hawaiiensis' ~ '6',
                              species == 'Toxoptera citricida' ~ '7',
                              species == 'Acyrthosiphon pisum' ~ '8',
                              species == 'Aedes albopictus' ~ '9',
                              species == 'Helicoverpa armigera' ~ '10',
                              species == 'Muscidifurax zaraptor' ~ '11',
                              species == 'Paracoccus marginatus' ~ '12',
                              species == 'Rhopalosiphum maidis' ~ '13',
                              species == 'Tetraneura nigriabdominalis' ~ '14',
                              species == 'Planococcus citri' ~ '15',
                              species == 'Anopheles gambiae' ~ '16',
                              species == 'Aubeonymus mariaefranciscae' ~ '17',
                              species == 'Amblyseius womersleyi' ~ '18',
                              species == 'Ophraella communa' ~ '19',
                              species == 'Anthonomus grandis' ~ '20',
                              species == 'Hyperaspis notata' ~ '21',
                              species == 'Tetranychus urticae' ~ '22',
                              species == 'Scapsipedus icipe' ~ '23',
                              species == 'Tribolium castaneum' ~ '24')) %>%
                              arrange(curve_ID) %>% filter(curve_ID != 'NA')

 zetaJ %>% distinct(species,curve_ID) %>% print(n=30)
 RearZetaJ %>% distinct(species) %>% print(n=30)

 RearZetaJ <- inner_join(zetaJ,RearZetaJ)
 
 head(RearZetaJ)
 
 RzetaJ_model <- lm(estimate ~ RTemp, RearZetaJ); summary(RzetaJ_model)

 
 Rear_zj <- ggplot(RearZetaJ, aes(x=RTemp, y=estimate)) +
   geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="gray")+
   geom_linerange(aes(ymin =conf_lower , ymax = conf_upper),size= 0.1,col="#000000") +
   geom_point(size = 2.5,stroke=0.2, col = '#000000', shape = 21, fill="#1f78b4")+
   theme_bw()+
   scale_x_continuous(expression(plain(paste("Rearing Temperature (",degree,"C)"))),
                      limits =c(19.8,30.2),
                      expand = c(0, 0),
                      breaks=seq(20,30, by=2))+
   scale_y_continuous(expression(paste(italic(T[pk]), " of juvenile mortality rate (", italic(z[J]),")")),
                      limits =c(5,35),
                      expand = c(0, 0),
                      breaks=seq(10,30, by=5))+
   theme(legend.position = 'none',legend.text = element_text(size = 10))+
   theme(text=element_text(family="Times"), panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())+
   geom_text(aes(x = 20.5, y = 33.5,label = "D"), 
             parse = TRUE, size = 4, colour = "black", family= 'Times')
 
 Rear_zj

 
 RearZeta <- RearDF %>% filter(trait == 'z') %>% 
   group_by(species) %>% 
   summarise(avg = mean(RearTemp)) %>%
   arrange(avg) %>%
   mutate(curve_ID = case_when(species == 'Aulacorthum solani' ~ '1',
                               species == 'Plutella xylostella' ~ '2',
                               species == 'Aphis nasturtii' ~ '3',
                               species == 'Stethorus punctillum' ~ '4',
                               species == 'Thrips hawaiiensis' ~ '5',
                               species == 'Toxoptera citricida' ~ '6',
                               species == 'Aedes krombeini' ~ '7',
                               species == 'Antestiopsis thunbergii' ~ '8',
                               species == 'Helicoverpa armigera' ~ '9',
                               species == 'Monochamus leuconotus' ~ '10',
                               species == 'Muscidifurax zaraptor' ~ '11',
                               species == 'Paracoccus marginatus' ~ '12',
                               species == 'Phenacoccus solenopsis' ~ '13',
                               species == 'Rhopalosiphum maidis' ~ '14',
                               species == 'Tetraneura nigriabdominalis' ~ '15',
                               species == 'Tetranychus mcdanieli' ~ '16',
                               species == 'Aedes albopictus' ~ '17',
                               species == 'Anopheles gambiae' ~ '18',
                               species == 'Ophraella communa' ~ '19',
                               species == 'Anthonomus grandis' ~ '20',
                               species == 'Culex pipiens' ~ '21',
                               species == 'Culex quinquefasciatus' ~ '22',
                               species == 'Hyperaspis notata' ~ '23',
                               species == 'Scapsipedus icipe' ~ '24',
                               species == 'Tribolium castaneum' ~ '25')) %>%
   arrange(curve_ID) %>% rename(RTemp = avg) %>%
   filter(curve_ID != 'NA')
 
RearZeta %>% distinct(species) %>% print(n=60)

 zeta  <- as_tibble(read.csv('data/z_Tpks_AllParams.csv', header = TRUE)) %>%
   filter(param == 'topt') %>% 
   select(species, estimate, conf_lower, conf_upper, trait) %>%
   mutate(curve_ID = case_when(species == 'Aulacorthum solani' ~ '1',
                               species == 'Plutella xylostella' ~ '2',
                               species == 'Aphis nasturtii' ~ '3',
                               species == 'Stethorus punctillum' ~ '4',
                               species == 'Thrips hawaiiensis' ~ '5',
                               species == 'Toxoptera citricida' ~ '6',
                               species == 'Aedes krombeini' ~ '7',
                               species == 'Antestiopsis thunbergii' ~ '8',
                               species == 'Helicoverpa armigera' ~ '9',
                               species == 'Monochamus leuconotus' ~ '10',
                               species == 'Muscidifurax zaraptor' ~ '11',
                               species == 'Paracoccus marginatus' ~ '12',
                               species == 'Phenacoccus solenopsis' ~ '13',
                               species == 'Rhopalosiphum maidis' ~ '14',
                               species == 'Tetraneura nigriabdominalis' ~ '15',
                               species == 'Tetranychus mcdanieli' ~ '16',
                               species == 'Aedes albopictus' ~ '17',
                               species == 'Anopheles gambiae' ~ '18',
                               species == 'Ophraella communa' ~ '19',
                               species == 'Anthonomus grandis' ~ '20',
                               species == 'Culex pipiens' ~ '21',
                               species == 'Culex quinquefasciatus' ~ '22',
                               species == 'Hyperaspis notata' ~ '23',
                               species == 'Scapsipedus icipe' ~ '24',
                               species == 'Tribolium castaneum' ~ '25')) %>%
   arrange(curve_ID) %>% 
   filter(curve_ID != 'NA')
 
 
zeta %>% distinct(species,curve_ID) %>% print(n=30)

RearZeta <- inner_join(zeta,RearZeta) 

Rzeta_model <- lm(estimate ~ RTemp, RearZeta); summary(Rzeta_model)

Rear_z <- ggplot(RearZeta, aes(x=RTemp, y=estimate)) +
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="gray") +
  geom_linerange(aes(ymin =conf_lower , ymax = conf_upper),size= 0.1,col="#000000") +
  geom_point(size = 2.5,stroke=0.2, col = '#000000', shape=22, fill ="#a6cee3")+
  theme_bw()+
  scale_x_continuous(expression(plain(paste("Rearing Temperature (",degree,"C)"))),
                     limits =c(19.8,30.2),
                     expand = c(0, 0),
                     breaks=seq(20,32, by=2))+
  scale_y_continuous(expression(paste(italic(T[pk]), " of adult mortality rate (", italic(z),")")),
                     limits =c(5,30.5),
                     expand = c(0, 0),
                     breaks=seq(10,30, by=5))+
  theme(legend.position = 'none',legend.text = element_text(size = 10))+
  theme(text=element_text(family="Times"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text(aes(x = 20.5, y = 29,label = "C"), 
            parse = TRUE, size = 4, colour = "black", family='Times')

Rear_z  


p1 <- Rear_alp+Rear_bmax+Rear_z+Rear_zj

save_plot(p1, file="results/SI/Traits_RearTemp_Adapt.pdf", 
          base_height=15,base_width = 17, base_asp = 0.75,units="cm")


##### Relationship between Latitude and Rearing Temperature #####

RearTemp_vs_Lat_Data <- as_tibble(read_csv('data/TraitData.csv')) %>% 
  select(standardisedtraitname, interactor1, interactor1temp, interactor1growthtemp, latitude) %>%
  filter(interactor1growthtemp != 'not stated' & interactor1growthtemp != 'various' & 
           interactor1growthtemp != 'NA' & latitude != 'NA') %>%
  rename(trait = standardisedtraitname, species = interactor1, RearTemp = interactor1growthtemp) %>%
  mutate(RearTemp = as.numeric(RearTemp), RearTemp = as.numeric(RearTemp))

head(RearTemp_vs_Lat_Data)


RearTemp_vs_Lat <- ggplot(RearTemp_vs_Lat_Data, aes(x=abs(latitude), y=RearTemp)) +
  geom_point(size = 1.2,stroke=0.1, col = '#000000')+
  geom_smooth(method = 'lm', colour = '#636363', size=0.3, fill="gray") +
  theme_bw()+
  scale_y_continuous(expression(plain(paste("Rearing Temperature (",degree,"C)")))) +
  scale_x_continuous(expression(plain(paste("Absolute latitude"))))+
  theme(legend.position = 'none', legend.text = element_text(size = 10))+
  theme(text=element_text(family="Times", size = 8), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

RearTemp_vs_Lat

save_plot(RearTemp_vs_Lat, file="results/SI/RearTemp_vs_Latitude.pdf", 
          base_height=8,base_width =9, base_asp = 1,units="cm")

RearTemp_vs_Lat_model <- lm(RearTemp ~ abs(latitude), RearTemp_vs_Lat_Data); summary(RearTemp_vs_Lat_model)



 