#__________________________________________________#
#              Mismatch plot                       #
#__________________________________________________#

require('dplyr')
require('ggplot2')
require('nls.multstart')
require('broom')
require('tidyverse')
require('rTPC')
require('data.table')
require('car')
require('boot')
require('patchwork')
require('minpack.lm')
require("tidyr")
require('purrr')

# update.packages(ask = FALSE)

rm(list=ls())
graphics.off()

#setwd("/home/primuser/Documents/VByte/VecMismatchPaper1/code/")
setwd("~/Dropbox/ph_thesis/Topt_paper")

#read in the Topt data

mnuniquefinalmoz <- as_tibble(read.csv('eulerparammnfin.csv', header = TRUE))

mnuniquefinalmoz$species <- as.factor(mnuniquefinalmoz$species)
mnuniquefinalmoz$param   <- as.factor(mnuniquefinalmoz$param)
mnuniquefinalmoz$stage   <- as.factor(mnuniquefinalmoz$stage)

mnuniquefinalmoz$stage <- relevel(mnuniquefinalmoz$stage, "juvenile")

ggplot(mnuniquefinalmoz, aes(estimate, species, shape=param, colour=param,fill=param)) +
  geom_errorbar(aes(xmin = conf_lower, xmax = conf_upper),width=0.5) +
  geom_point(size = 3.5, col="#000000") + 
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank())+
  scale_x_continuous(expression(plain(paste("", italic(T)[opt]))),
                     limits =c(10,38),
                     expand = c(0, 0),
                     breaks=seq(12,36, by=2))+
  scale_fill_manual(values = c("#39558CFF","#FDE725FF","#FDE725FF","#39558CFF"),
                      name=expression(bold("")),
                      guide = guide_legend(nrow = 1,ncol =4 ,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  scale_colour_manual(values = c("#39558CFF","#FDE725FF","#FDE725FF","#39558CFF"),
                    name=expression(bold("")),
                    guide = guide_legend(nrow = 1,ncol =4 ,
                                         direction = "vertical",
                                         title.position = "top",
                                         title.hjust=0.5))+
  scale_shape_manual(values = c(21,22,23,24),
                      name=expression(bold("")),
                      guide = guide_legend(nrow = 1,ncol =4 ,
                                           direction = "vertical",
                                           title.position = "top",
                                           title.hjust=0.5))+
  theme(legend.position = 'top')




