# Load libraries
require(tidyverse)
require(patchwork) # easy way of creating panels
require(car)
require(forcats)
require(cowplot)

setwd("~/Dropbox/TraitTesting/data")

speciesmass <- as_tibble(read_csv('TraitData.csv')) %>%
  select(originaltraitname, originaltraitvalue, interactor1) %>%
  filter(originaltraitname == 'body size' & originaltraitvalue != 'NA')
      
sizeMeans <- 
  speciesmass %>% 
  group_by(interactor1) %>% 
  summarise(avg = mean(originaltraitvalue), sd = sd(originaltraitvalue)) %>%
  arrange(avg)
  
write_csv(sizeMeans, 'sizeMeans.csv')

# plot means to get idea of size variation across species 
sizeMeans %>%
  mutate(interactor1 = fct_reorder(interactor1,avg)) %>%
  arrange(avg) %>%
  ggplot(aes(x=avg, y=interactor1, colour = interactor1))+
  geom_point()

  
  
  





