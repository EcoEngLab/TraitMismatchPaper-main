library(dplyr)
library(ggplot2)
library(nls.multstart)
library(broom)
library(tidyverse)
library(rTPC)
library(data.table)
library(car)
library(boot)
library(patchwork)
library(minpack.lm)
library(tidyr)
library(purrr)
library(rTPC)


#install.packages('rTPC')
#update.packages(ask = FALSE)

rm(list=ls())
graphics.off()

setwd("/home/primuser/Documents/VByte/VecMismatchPaper1/code/")

#take a look at the different models available
get_model_names()

#read in the trait data
final_trait_data <- read.csv('../data/Final_Traitofinterest.csv')

#remove completely irrelevant columns 
df <- final_trait_data[,colSums(is.na(final_trait_data))<nrow(final_trait_data)]

df <- df %>%
  select('originalid', 'originaltraitname', 'originaltraitunit', 'originaltraitvalue', 'interactor1','interactor1stage', 'ambienttemp', 'citation')


d<- df %>%
  rename(temp = ambienttemp,
         rate = originaltraitvalue)

d$rate <- d$rate + 10^-6



d <-mutate(d, inversetraitvalue = 1/rate)


#change heading to correct variable
d <-  mutate(d, adjustedtraitname = ifelse(originaltraitname == 'Mortality Rate' , 'zj',
                                           ifelse(originaltraitname == 'Egg development time' , 'a',
                                                  ifelse(originaltraitname == 'Generation Time' , 'a',
                                                         ifelse(originaltraitname == 'Survival Rate' , 'z',
                                                                ifelse(originaltraitname == 'Development time' , 'a',
                                                                       ifelse(originaltraitname == 'Development Time' , 'a',
                                                                              ifelse(originaltraitname == 'Development Rate' , 'a',
                                                                                     ifelse(originaltraitname == 'Survivorship' , 'z',
                                                                                            ifelse(originaltraitname == 'Longevity' , 'z',
                                                                                                   ifelse(originaltraitname == 'Survival Time' , 'z',
                                                                                                          ifelse(originaltraitname == 'Percentage Survival' , 'z',
                                                                                                                 ifelse(originaltraitname == 'Oviposition Rate' , 'bpk',
                                                                                                                        ifelse(originaltraitname == 'Juvenile survival ' , 'zj',
                                                                                                                               ifelse(originaltraitname == 'Fecundity Rate' , 'bpk',
                                                                                                                                      ifelse(originaltraitname == 'Fecundity' , 'bpk',
                                                                                                                                             originaltraitname))))))))))))))))

#concatenate to get curves
d$cit <- gsub("([A-Za-z]+).*", "\\1", d$citation)
#d$conc <- paste(d$interactor1, "_", d$originaltraitname, "_", d$cit)
d$conc <- paste(d$interactor1, '_', d$interactor1stage)
d$conc <- as.character(d$conc)
d$conc <- gsub(" ", "", d$conc, fixed = TRUE)
numberofcurves <- unique(d$conc)
#write.csv(d, '../data/finaltraitofinterestmn.csv')
d$conc <- paste(d$interactor1, '_', d$interactor1stage)
#____________________________________________________#
#               Development rate                     #
#____________________________________________________#



#real rate
d_rate <- filter(d, adjustedtraitname=="a")
d_rate <- filter(d_rate, originaltraitunit=="1/day")
d_rate <- select(d_rate,temp,rate,conc)
d_rate <- as_tibble(d_rate)
d_rate$temp <- as.numeric(d_rate$temp)
d_rate$rate <- as.numeric(d_rate$rate)

#real rate
d_rate32 <- filter(d, adjustedtraitname=="a")
d_rate32 <- filter(d_rate32, originaltraitunit=="1/days")
d_rate32 <- select(d_rate32,temp,rate,conc)
d_rate32 <- as_tibble(d_rate32)
d_rate32$temp <- as.numeric(d_rate32$temp)
d_rate32$rate <- as.numeric(d_rate32$rate)

#just day
d_time <- filter(d, adjustedtraitname=="a")
d_time <- filter(d_time, originaltraitunit=="days")
d_time$rate <- 1/d_time$rate
d_time <- select(d_time,temp,rate,conc)
d_time <- as_tibble(d_time)
d_time$temp <- as.numeric(d_time$temp)
d_time$rate <- as.numeric(d_time$rate)

#day but caps
d_time2 <- filter(d, adjustedtraitname=="a")
d_time2 <- filter(d_time2, originaltraitunit=="Days")
d_time2$rate <- 1/d_time2$rate
d_time2 <- select(d_time2,temp,rate,conc)
d_time2 <- as_tibble(d_time2)
d_time2$temp <- as.numeric(d_time2$temp)
d_time2$rate <- as.numeric(d_time2$rate)

#100 per day
d_rate_hunned <- filter(d, adjustedtraitname=="a")
d_rate_hunned <- filter(d_rate_hunned, originaltraitunit=="100/day")
d_rate_hunned$rate <- d_rate_hunned$rate/100
d_rate_hunned <- select(d_rate_hunned,temp,rate,conc)
d_rate_hunned <- as_tibble(d_rate_hunned)
d_rate_hunned$temp <- as.numeric(d_rate_hunned$temp)
d_rate_hunned$rate <- as.numeric(d_rate_hunned$rate)

#day per individual
d_time3 <- filter(d, adjustedtraitname=="a")
d_time3 <- filter(d_time3, originaltraitunit=="day / 1 individual")
d_time3$rate <- 1/d_time3$rate
d_time3 <- select(d_time3,temp,rate,conc)
d_time3 <- as_tibble(d_time3)
d_time3$temp <- as.numeric(d_time3$temp)
d_time3$rate <- as.numeric(d_time3$rate)

#real rate but with a bullshit label
d_rate2 <- filter(d, adjustedtraitname=="a")
d_rate2 <- filter(d_rate2, originaltraitunit=="event / (1 individual * 24 hour)")
d_rate2 <- select(d_rate2,temp,rate,conc)
d_rate2 <- as_tibble(d_rate2)
d_rate2$temp <- as.numeric(d_rate2$temp)
d_rate2$rate <- as.numeric(d_rate2$rate)



d_rate <- rbind(d_rate,d_time, d_rate_hunned, d_rate2, d_time2, d_time3, d_rate32)


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(d_rate$temp, d_rate$rate, model_name = 'sharpeschoolhigh_1981')

d_fits <- nest(d_rate, data = c(temp, rate)) %>%
  mutate(sharpeschoolhigh = map(data, ~nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                                     data = .x,
                                                     iter = c(3,3,3,3),
                                                     start_lower = start_vals - 10,
                                                     start_upper = start_vals + 10,
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)))

summary(d_fits$sharpeschoolhigh[[1]])

d_preds <- mutate(d_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(sharpeschoolhigh)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(conc, model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


glimpse(d_preds)

# plot panel

ggplot(d_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, d_rate) +
  facet_wrap(~conc, scales = 'free_y', ncol = 6) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")



# stack models and calculate extra params
d_params <- pivot_longer(d_fits, names_to = 'model_name', values_to = 'fit', c(sharpeschoolhigh)) %>%
  mutate(params = map(fit, calc_params)) %>%
  select(conc,model_name, params) %>%
  unnest(params)

d_params



