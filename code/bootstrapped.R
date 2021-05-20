require('ggplot2')
require('nls.multstart')
require('broom')
require('tidyverse')
require('rTPC')
require('dplyr')
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
setwd("/home/primuser/Documents/VByte/VecMismatchPaper1/code/")

#take a look at the different models available
get_model_names()

#read in the trait data
final_trait_data <- read.csv('../data/Final_Traitofinterest.csv')

#filter out sets that less than the required parameters for schoolfield-high ()
final_trait_data <- dplyr::filter(final_trait_data,
                                  originaltraitname != 'Adult survival' &
                                    originalid != 'MSS0059' &
                                    originaltraitname != 'Adult longevity (female, bloodfed)' &
                                    originaltraitname != 'Adult longevity (male)' &
                                    originaltraitname != 'Adult survival (female, bloodfed)' &
                                    originaltraitname != 'Adult survival (male)' )
#remove completely irrelevant columns 
df <- final_trait_data[,colSums(is.na(final_trait_data))<nrow(final_trait_data)]

#filter to single species and trait
df2 <- dplyr::filter(df, originalid == 'csm7I')


df1 <- df %>%
  select('originalid', 'originaltraitname', 'originaltraitunit', 'originaltraitvalue', 'interactor1', 'ambienttemp')



#visualize
ggplot(df2, aes(ambienttemp, originaltraitvalue))+
  geom_point()+
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Development Rate',
       title = 'Development Rate across temperatures for Aedes albopictus')


# choose model
mod = 'sharpschoolhigh_1981'

#mutate the titles because I am lazy and I don't want to change the formula
df2 <- df2 %>%
  mutate(temp = ambienttemp,
         rate = originaltraitvalue)

# fit Sharpe-Schoolfield model
d_fit <- nest(df2, data = c(temp, rate)) %>%
  mutate(sharpeschoolhigh = map(data, ~nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                                     data = .x,
                                                     iter = c(3,3,3,3),
                                                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') - 10,
                                                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') + 10,
                                                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))),
         # predict over that data,
         preds =  map2(sharpeschoolhigh, new_data, ~augment(.x, newdata = .y)))


# 
# maptest <- nest(df2, data = c(temp, rate)) %>% map(data, ~nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
#                          data = .,
#                          iter = c(3,3,3,3),
#                          start_lower = get_start_vals(.$temp, .$rate, model_name = 'sharpeschoolhigh_1981') - 10,
#                          start_upper = get_start_vals(.$temp, .$rate, model_name = 'sharpeschoolhigh_1981') + 10,
#                          lower = get_lower_lims(.$temp, .$rate, model_name = 'sharpeschoolhigh_1981'),
#                          upper = get_upper_lims(.$temp, .$rate, model_name = 'sharpeschoolhigh_1981'),
#                          supp_errors = 'Y',
#                          convergence_count = FALSE))


# unnest predictions
d_preds <- select(d_fit, preds) %>%
  unnest(preds)

# plot data and predictions
ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'blue') +
  geom_point(aes(temp, rate), d, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Growth rate',
       title = 'Growth rate across temperatures')

# fit model
fit <- nls_multstart(originaltraitvalue~sharpeschoolhigh_1981(temp = ambienttemp, r_tref,e,eh,th, tref = 15),
                     data = df2,
                     iter = 500,
                     start_lower = start_vals - 10,
                     start_upper = start_vals + 10,
                     lower = low_lims,
                     upper = upper_lims,
                     supp_errors = 'Y')
fit


# calculate additional traits
calc_params(fit) %>%
  # round for easy viewing
  mutate_all(round, 2)


# predict new data
new_data <- data.frame(ambienttemp = seq(min(df2$ambienttemp), max(df2$ambienttemp), 0.5))
preds <- augment(fit, newdata = new_data)

#visualize with new line
ggplot(df2, aes(ambienttemp, originaltraitvalue))+
  geom_point()+
  geom_line(aes(ambienttemp, .fitted), preds, col = 'blue') +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Development Rate',
       title = 'Development Rate across Temperatures for Aedes albopictus')

unique(df$originalid)


#visualize all together
ggplot(df, aes(ambienttemp, originaltraitvalue))+
  geom_point()+
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Trait Data')+
  facet_wrap(df$originaltraitname, scales = 'free' )


#mutate some things
df10<-  mutate(df1, modifiedtraitvalue = ifelse(originaltraitname == 'Mortality Rate' , 1/originaltraitvalue,
                                                ifelse(originaltraitname == 'Egg development time' , 1/originaltraitvalue,
                                                       ifelse(originaltraitname == 'Generation time' , 1/originaltraitvalue,
                                                              ifelse(originaltraitname == 'Survival Rate' , 1/originaltraitvalue,
                                                                     ifelse(originaltraitname == 'Development time' , 1/originaltraitvalue,
                                                                            ifelse(originaltraitname == 'Development Time' , 1/originaltraitvalue,
                                                                                   ifelse(originaltraitname == 'Survivorship' , 1/originaltraitvalue,
                                                                                          ifelse(originaltraitname == 'Longevity' , 1/originaltraitvalue,
                                                                                                 originaltraitvalue)))))))))


#visualize all together with modified data
ggplot(df10, aes(ambienttemp, modifiedtraitvalue))+
  geom_point()+
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Trait Data')+
  facet_wrap(df10$originaltraitname, scales = 'free' )

#test a loop mean of temperature
shoobop <- list()
for (x in unique(df$originalid)){
  # Subset
  dfsub <- subset(df, df$originalid == x)
  # Apply function
  mean(dfsub$ambienttemp)
  shoobop[[x]] <- print(paste0(x, ": ", mean(dfsub$ambienttemp)))
}

#base loop for what needs to be flipped
modelfit <- list()
parameters <- list()
for (x in unique(df$originalid)){
  # Subset
  dfsub <- subset(df, df$originalid == x)
  # get start vals
  start_vals <- get_start_vals(dfsub$ambienttemp, dfsub$originaltraitvalue, model_name = 'sharpeschoolhigh_1981')
  
  # get limits
  low_lims <- get_lower_lims(dfsub$ambienttemp, dfsub$originaltraitvalue, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(dfsub$ambienttemp, dfsub$originaltraitvalue, model_name = 'sharpeschoolhigh_1981')
  
  modelfit[[x]] <-  nls_multstart(originaltraitvalue~sharpeschoolhigh_1981(temp = ambienttemp, r_tref,e,eh,th, tref = 15),
                                  data = dfsub,
                                  iter = 500,
                                  start_lower = start_vals - 10,
                                  start_upper = start_vals + 10,
                                  lower = low_lims,
                                  upper = upper_lims,
                                  supp_errors = 'Y')
  parameters[[x]] <- calc_params(modelfit[[x]])
  
}



#loop with flips incorporated
modelfit2 <- list()
parameters2 <- list()
preds2 <- list()
for (x in unique(df10$originalid)){
  # Subset
  dfsub <- subset(df10, df10$originalid == x)
  # get start vals
  start_vals <- get_start_vals(dfsub$ambienttemp, dfsub$modifiedtraitvalue, model_name = 'sharpeschoolhigh_1981')
  
  # get limits
  low_lims <- get_lower_lims(dfsub$ambienttemp, dfsub$modifiedtraitvalue, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(dfsub$ambienttemp, dfsub$modifiedtraitvalue, model_name = 'sharpeschoolhigh_1981')
  
  modelfit2[[x]] <-  nls_multstart(modifiedtraitvalue~sharpeschoolhigh_1981(temp = ambienttemp, r_tref,e,eh,th, tref = 15),
                                   data = dfsub,
                                   iter = 500,
                                   start_lower = start_vals - 10,
                                   start_upper = start_vals + 10,
                                   lower = low_lims,
                                   upper = upper_lims,
                                   supp_errors = 'Y')
  parameters2[[x]] <- calc_params(modelfit2[[x]])
  
  #stuff a prediction loop into this because why not?
  new_data <- data.frame(ambienttemp = seq(min(dfsub$ambienttemp), max(dfsub$ambienttemp), 0.5))
  preds2[[x]] <- augment(modelfit2[[x]], newdata = new_data)
}


#visualize
ggplot(df10, aes(ambienttemp, originaltraitvalue))+
  geom_point()+
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Trait Data')+
  facet_wrap(df10$originalid, scales = 'free' )

#visualize with new line across all
ggplot(df10, aes(ambienttemp, originaltraitvalue))+
  geom_point()+
  # geom_line(aes(ambienttemp, preds2$), preds, col = 'blue') +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Trait Data')+
  facet_wrap(df10$originalid, scales = 'free' )

#build and clean final dataframe (104 work well)
paramdf <- df10 %>%
  select('originalid', 'interactor1', 'originaltraitname')
paramdf <- paramdf %>% distinct((originalid), .keep_all = TRUE)
boundparams <- bind_rows(parameters2, .id = "originalid")
boundparams <- merge(paramdf, boundparams)
drops <- c("(originalid)")
boundparams <- boundparams[ , !(names(boundparams) %in% drops)]


#mutate some things
boundparams<-  mutate(boundparams, adjustedtraitname = ifelse(originaltraitname == 'Mortality Rate' , 'zj',
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
                                                                                                                                                  ifelse(originaltraitname == 'Fecundity Rate' , 'k',
                                                                                                                                                         ifelse(originaltraitname == 'Fecundity' , 'k',
                                                                                                                                                                originaltraitname))))))))))))))))



#plot topt
ggplot(boundparams, aes(topt, interactor1, color = adjustedtraitname))+
  geom_pointrange(aes(xmin=topt-breadth, xmax=topt+breadth))+
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC) (Thermal Optimum)',
       y = 'Insect Species')

#plot activation energy
ggplot(boundparams, aes(e, interactor1, color = originaltraitname))+
  geom_point()+
  theme_bw(base_size = 12) +
  labs(x = 'activation energy',
       y = 'Insect Species')