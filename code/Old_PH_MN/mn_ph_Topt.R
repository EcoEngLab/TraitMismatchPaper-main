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
  select('originalid', 'originaltraitname', 'originaltraitunit', 'originaltraitvalue', 'interactor1', 'ambienttemp', 'citation')


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
                                                                                                                               ifelse(originaltraitname == 'Fecundity Rate' , 'k',
                                                                                                                                      ifelse(originaltraitname == 'Fecundity' , 'k',
                                                                                                                                             originaltraitname))))))))))))))))

#concatenate to get curves
d$cit <- gsub("([A-Za-z]+).*", "\\1", d$citation)
#d$conc <- paste(d$interactor1, "_", d$originaltraitname, "_", d$cit)
d$conc <- paste(d$interactor1)
d$conc <- as.character(d$conc)
d$conc <- gsub(" ", "", d$conc, fixed = TRUE)
numberofcurves <- unique(d$conc)


#____________________________________________________#
#               Development rate                     #
#____________________________________________________#

d_rate <- filter(d, adjustedtraitname=="a")
d_rate <- filter(d_rate, originaltraitname=="Development Rate")
d_rate <- filter(d_rate, rate <= 1)
d_rate <- select(d_rate,temp,rate,conc)
d_rate <- as_tibble(d_rate)
d_rate$temp <- as.numeric(d_rate$temp)
d_rate$rate <- as.numeric(d_rate$rate)

d_time <- filter(d, adjustedtraitname=="a")
d_time <- filter(d_time, rate >= 1)
d_time$rate <- 1/d_time$rate
d_time <- select(d_time,temp,rate,conc)
d_time <- as_tibble(d_time)
d_time$temp <- as.numeric(d_time$temp)
d_time$rate <- as.numeric(d_time$rate)

d_rate <- rbind(d_rate,d_time)


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


# sample(1:2, 20, replace=TRUE)


#_____________________________________________________________________________#
#                         Bootstrapping by species                            #
#_____________________________________________________________________________#

d1 <- filter(d_rate, conc=="Aedesalbopictus")

fit_nlsLM1 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                               data = d1,
                               start = coef(d_fits$sharpeschoolhigh[[1]]),
                               lower = get_lower_lims(d1$temp, d1$rate, model_name = 'sharpeschoolhigh_1981'),
                               upper = get_upper_lims(d1$temp, d1$rate, model_name = 'sharpeschoolhigh_1981'),
                               weights = rep(1, times = nrow(d1)))

# bootstrap using case resampling
boot1 <- Boot(fit_nlsLM1, method = 'residual')

# look at the data
head(boot1$t)


boot1_preds <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d1$temp), max(d1$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot1_conf_preds <- group_by(boot1_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

dp1 <- filter(d_preds, conc=="Aedesalbopictus")

# plot bootstrapped CIs
  
p1 <- ggplot() +
  geom_line(aes(temp, .fitted), dp1, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d1, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')

# calculate params with CIs

extra_params1 <- calc_params(fit_nlsLM1) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params1 <- Boot(fit_nlsLM1, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM1)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params1 <- left_join(ci_extra_params1, extra_params1)
#> Joining, by = "param"

#ggplot(ci_extra_params1, aes(param, estimate)) +
#geom_point(size = 1) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = ' dev rate TPC; using residual resampling')

# topt

topt1 <- as_tibble(ci_extra_params1[2,])
topt1$species <- as.character("Aedesalbopictus")

ggplot(topt1, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())


#____________________________________
  
  d2 <- filter(d_rate, conc=="Anthonomusgrandis")
  
  fit_nlsLM2 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                  data = d2,
                                  start = coef(d_fits$sharpeschoolhigh[[2]]),
                                  lower = get_lower_lims(d2$temp, d2$rate, model_name = 'sharpeschoolhigh_1981'),
                                  upper = get_upper_lims(d2$temp, d2$rate, model_name = 'sharpeschoolhigh_1981'),
                                  weights = rep(1, times = nrow(d2)))
  
  # bootstrap using case resampling
  boot2 <- Boot(fit_nlsLM2, method = 'residual')
  
  # look at the data
  head(boot2$t)
  
  
  boot2_preds <- boot2$t %>%
    as.data.frame() %>%
    drop_na() %>%
    mutate(iter = 1:n()) %>%
    group_by_all() %>%
    do(data.frame(temp = seq(min(d2$temp), max(d2$temp), length.out = 100))) %>%
    ungroup() %>%
    mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))
  
  # calculate bootstrapped confidence intervals
  boot2_conf_preds <- group_by(boot2_preds, temp) %>%
    summarise(conf_lower = quantile(pred, 0.025),
              conf_upper = quantile(pred, 0.975)) %>%
    ungroup()
  
  
  dp2 <- filter(d_preds, conc=="Anthonomusgrandis")
  
# plot bootstrapped CIs
 p2 <- ggplot() +
    geom_line(aes(temp, .fitted), dp2, col = 'blue') +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds, fill = 'blue', alpha = 0.3) +
    geom_point(aes(temp, rate), d2, size = 2, alpha = 0.5) +
    theme_bw(base_size = 12) +
    labs(x = 'Temperature (?C)',
         y = 'dev rate')
 
# calculate params with CIs
 
 extra_params2 <- calc_params(fit_nlsLM2) %>%
   pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
 
 ci_extra_params2 <- Boot(fit_nlsLM2, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM2)), R = 200, method = 'residual') %>%
   confint(., method = 'bca') %>%
   as.data.frame() %>%
   rename(conf_lower = 1, conf_upper = 2) %>%
   rownames_to_column(., var = 'param') %>%
   mutate(method = 'residual bootstrap')
 
 ci_extra_params2 <- left_join(ci_extra_params2, extra_params2)
 #> Joining, by = "param"
 
 #ggplot(ci_extra_params2, aes(param, estimate)) +
 #geom_point(size = 2) +
 #geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
 #theme_bw() +
 #facet_wrap(~param, scales = 'free') +
 #scale_x_discrete('') +
 #labs(title = 'Calculation of confidence intervals for extra parameters',
 #subtitle = 'dev rate TPC; using residual resampling')
 
 # topt
 
 topt2 <- as_tibble(ci_extra_params2[2,])
 topt2$species <- as.character("Anthonomusgrandis")
 
 ggplot(topt2, aes(estimate, species)) +
   geom_point(size = 4) + 
   geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
   theme_bw(base_size = 12) +
   facet_wrap(~species, scales = 'free') +
   scale_x_continuous('') +
   labs(title = 'calculation of Topt with CIs',
        subtitle = 'dev rate TPC; using residual resampling')+
   theme(axis.title.y = element_blank())
 
  
  
#_______________________________________________
  
  d3 <- filter(d_rate, conc=="Aphisgossypii")
  
  fit_nlsLM3 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                  data = d3,
                                  start = coef(d_fits$sharpeschoolhigh[[3]]),
                                  lower = get_lower_lims(d3$temp, d3$rate, model_name = 'sharpeschoolhigh_1981'),
                                  upper = get_upper_lims(d3$temp, d3$rate, model_name = 'sharpeschoolhigh_1981'),
                                  weights = rep(1, times = nrow(d3)))
  
  # bootstrap using case resampling
  boot3 <- Boot(fit_nlsLM3, method = 'residual')
  
  boot3_preds <- boot3$t %>%
    as.data.frame() %>%
    drop_na() %>%
    mutate(iter = 1:n()) %>%
    group_by_all() %>%
    do(data.frame(temp = seq(min(d3$temp), max(d3$temp), length.out = 100))) %>%
    ungroup() %>%
    mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))
  
  # calculate bootstrapped confidence intervals
  boot3_conf_preds <- group_by(boot3_preds, temp) %>%
    summarise(conf_lower = quantile(pred, 0.025),
              conf_upper = quantile(pred, 0.975)) %>%
    ungroup()
  
  
  dp3 <- filter(d_preds, conc=="Aphisgossypii")
  
# plot bootstrapped CIs
p3 <-  ggplot() +
    geom_line(aes(temp, .fitted), dp3, col = 'blue') +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot3_conf_preds, fill = 'blue', alpha = 0.3) +
    geom_point(aes(temp, rate), d3, size = 2, alpha = 0.5) +
    theme_bw(base_size = 12) +
    labs(x = 'Temperature (?C)',
         y = 'dev rate')

# calculate params with CIs

extra_params3 <- calc_params(fit_nlsLM3) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params3 <- Boot(fit_nlsLM3, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM3)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params3 <- left_join(ci_extra_params3, extra_params3)
#> Joining, by = "param"

#ggplot(ci_extra_params3, aes(param, estimate)) +
#geom_point(size = 2) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt3 <- as_tibble(ci_extra_params3[2,])
topt3$species <- as.character("Aphisgossypii")

ggplot(topt3, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())



#_______________________________________________

d4 <- filter(d_rate, conc=="Macrocentrusiridescens")

fit_nlsLM4 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                data = d4,
                                start = coef(d_fits$sharpeschoolhigh[[4]]),
                                lower = get_lower_lims(d4$temp, d4$rate, model_name = 'sharpeschoolhigh_1981'),
                                upper = get_upper_lims(d4$temp, d4$rate, model_name = 'sharpeschoolhigh_1981'),
                                weights = rep(1, times = nrow(d4)))

# bootstrap using case resampling
boot4 <- Boot(fit_nlsLM4, method = 'residual')

boot4_preds <- boot4$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d4$temp), max(d4$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot4_conf_preds <- group_by(boot4_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp4 <- filter(d_preds, conc=="Macrocentrusiridescens")

# plot bootstrapped CIs
p4 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp4, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot4_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d4, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')

# calculate params with CIs

extra_params4 <- calc_params(fit_nlsLM4) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params4 <- Boot(fit_nlsLM4, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM4)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params4 <- left_join(ci_extra_params4, extra_params4)
#> Joining, by = "param"

#ggplot(ci_extra_params4, aes(param, estimate)) +
  #geom_point(size = 4) +
  #geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  #theme_bw() +
  #facet_wrap(~param, scales = 'free') +
  #scale_x_discrete('') +
  #labs(title = 'Calculation of confidence intervals for extra parameters',
       #subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt4 <- as_tibble(ci_extra_params4[2,])
topt4$species <- as.character("Macrocentrusiridescens")

ggplot(topt4, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())
  



#_______________________________________________

d5 <- filter(d_rate, conc=="Stethoruspunctillum")

fit_nlsLM5 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                data = d5,
                                start = coef(d_fits$sharpeschoolhigh[[5]]),
                                lower = get_lower_lims(d5$temp, d5$rate, model_name = 'sharpeschoolhigh_1981'),
                                upper = get_upper_lims(d5$temp, d5$rate, model_name = 'sharpeschoolhigh_1981'),
                                weights = rep(1, times = nrow(d5)))

# bootstrap using case resampling
boot5 <- Boot(fit_nlsLM5, method = 'residual')

boot5_preds <- boot5$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d5$temp), max(d5$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot5_conf_preds <- group_by(boot5_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp5 <- filter(d_preds, conc=="Stethoruspunctillum")

# plot bootstrapped CIs
p5 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp5, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot5_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d5, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')


# calculate params with CIs

extra_params5 <- calc_params(fit_nlsLM5) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params5 <- Boot(fit_nlsLM5, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM5)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params5 <- left_join(ci_extra_params5, extra_params5)
#> Joining, by = "param"

#ggplot(ci_extra_params5, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = ' dev rate TPC; using residual resampling')

# topt

topt5 <- as_tibble(ci_extra_params5[2,])
topt5$species <- as.character("Stethoruspunctillum")

ggplot(topt5, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())


#_______________________________________________

d6 <- filter(d_rate, conc=="Telenomusisis")

fit_nlsLM6 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                data = d6,
                                start = coef(d_fits$sharpeschoolhigh[[6]]),
                                lower = get_lower_lims(d6$temp, d6$rate, model_name = 'sharpeschoolhigh_1981'),
                                upper = get_upper_lims(d6$temp, d6$rate, model_name = 'sharpeschoolhigh_1981'),
                                weights = rep(1, times = nrow(d6)))

# bootstrap using case resampling
boot6 <- Boot(fit_nlsLM6, method = 'residual')

boot6_preds <- boot6$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d6$temp), max(d6$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot6_conf_preds <- group_by(boot6_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp6 <- filter(d_preds, conc=="Telenomusisis")

# plot bootstrapped CIs
p6 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp6, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot6_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d6, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')


# calculate params with CIs

extra_params6 <- calc_params(fit_nlsLM6) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params6 <- Boot(fit_nlsLM6, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM6)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params6 <- left_join(ci_extra_params6, extra_params6)
#> Joining, by = "param"

#ggplot(ci_extra_params6, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt6 <- as_tibble(ci_extra_params6[2,])
topt6$species <- as.character("Telenomusisis")

ggplot(topt6, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())


#_______________________________________________

d7 <- filter(d_rate, conc=="Tetranychusmcdanieli")

fit_nlsLM7 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                data = d7,
                                start = coef(d_fits$sharpeschoolhigh[[7]]),
                                lower = get_lower_lims(d7$temp, d7$rate, model_name = 'sharpeschoolhigh_1981'),
                                upper = get_upper_lims(d7$temp, d7$rate, model_name = 'sharpeschoolhigh_1981'),
                                weights = rep(1, times = nrow(d7)))

# bootstrap using case resampling
boot7 <- Boot(fit_nlsLM7, method = 'residual')

boot7_preds <- boot7$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d7$temp), max(d7$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot7_conf_preds <- group_by(boot7_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp7 <- filter(d_preds, conc=="Tetranychusmcdanieli")

# plot bootstrapped CIs
p7 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp7, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot7_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d7, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')


# calculate params with CIs

extra_params7 <- calc_params(fit_nlsLM7) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params7 <- Boot(fit_nlsLM7, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM7)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params7 <- left_join(ci_extra_params7, extra_params7)
#> Joining, by = "param"

#ggplot(ci_extra_params7, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt7 <- as_tibble(ci_extra_params7[2,])
topt7$species <- as.character("Tetranychusmcdanieli")

ggplot(topt7, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())


#_______________________________________________


d8 <- filter(d_rate, conc=="Aedescamptorhynchus")

fit_nlsLM8 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                data = d8,
                                start = coef(d_fits$sharpeschoolhigh[[8]]),
                                lower = get_lower_lims(d8$temp, d8$rate, model_name = 'sharpeschoolhigh_1981'),
                                upper = get_upper_lims(d8$temp, d8$rate, model_name = 'sharpeschoolhigh_1981'),
                                weights = rep(1, times = nrow(d8)))

# bootstrap using case resampling
boot8 <- Boot(fit_nlsLM8, method = 'residual')

boot8_preds <- boot8$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d8$temp), max(d8$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot8_conf_preds <- group_by(boot8_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.985)) %>%
  ungroup()


dp8 <- filter(d_preds, conc=="Aedescamptorhynchus")

# plot bootstrapped CIs
p8 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp8, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot8_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d8, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')


# calculate params with CIs

extra_params8 <- calc_params(fit_nlsLM8) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params8 <- Boot(fit_nlsLM8, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM8)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params8 <- left_join(ci_extra_params8, extra_params8)
#> Joining, by = "param"

#ggplot(ci_extra_params8, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt8 <- as_tibble(ci_extra_params8[2,])
topt8$species <- as.character("Aedescamptorhynchus")

ggplot(topt8, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())

#________________________________________________________

# ******** doesn't fit *********** #

d9 <- filter(d_rate, conc=="Aedesnotoscriptus")

fit_nlsLM9 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                data = d9,
                                start = coef(d_fits$sharpeschoolhigh[[9]]),
                                lower = get_lower_lims(d9$temp, d9$rate, model_name = 'sharpeschoolhigh_1981'),
                                upper = get_upper_lims(d9$temp, d9$rate, model_name = 'sharpeschoolhigh_1981'),
                                weights = rep(1, times = nrow(d9)))

# bootstrap using residual resampling
boot9 <- Boot(fit_nlsLM9, method = 'residual')

boot9_preds <- boot9$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d9$temp), max(d9$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1991(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot9_conf_preds <- group_by(boot9_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.995)) %>%
  ungroup()


dp9 <- filter(d_preds, conc=="Aedesnotoscriptus")

# plot bootstrapped CIs
p9 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp9, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot9_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d9, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')


# calculate params with CIs

extra_params9 <- calc_params(fit_nlsLM9) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params9 <- Boot(fit_nlsLM9, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM9)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params9 <- left_join(ci_extra_params9, extra_params9)
#> Joining, by = "param"

#ggplot(ci_extra_params9, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt9 <- as_tibble(ci_extra_params9[2,])
topt9$species <- as.character("Aedesnotoscriptus")

ggplot(topt9, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())

#________________________________________________________________

d10 <- filter(d_rate, conc=="Bactroceracorrecta")

fit_nlsLM10 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                data = d10,
                                start = coef(d_fits$sharpeschoolhigh[[10]]),
                                lower = get_lower_lims(d10$temp, d10$rate, model_name = 'sharpeschoolhigh_1981'),
                                upper = get_upper_lims(d10$temp, d10$rate, model_name = 'sharpeschoolhigh_1981'),
                                weights = rep(1, times = nrow(d10)))

# bootstrap using case resampling
boot10 <- Boot(fit_nlsLM10, method = 'residual')

boot10_preds <- boot10$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d10$temp), max(d10$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot10_conf_preds <- group_by(boot10_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp10 <- filter(d_preds, conc=="Bactroceracorrecta")

# plot bootstrapped CIs
p10 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp10, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot10_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d10, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')

# calculate params with CIs

extra_params10 <- calc_params(fit_nlsLM10) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params10 <- Boot(fit_nlsLM10, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM10)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params10 <- left_join(ci_extra_params10, extra_params10)
#> Joining, by = "param"

#ggplot(ci_extra_params10, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt10 <- as_tibble(ci_extra_params10[2,])
topt10$species <- as.character("Bactroceracorrecta")

ggplot(topt10, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())

#________________________________________________________________

d11 <- filter(d_rate, conc=="Corythuchaciliata")

fit_nlsLM11 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                 data = d11,
                                 start = coef(d_fits$sharpeschoolhigh[[11]]),
                                 lower = get_lower_lims(d11$temp, d11$rate, model_name = 'sharpeschoolhigh_1981'),
                                 upper = get_upper_lims(d11$temp, d11$rate, model_name = 'sharpeschoolhigh_1981'),
                                 weights = rep(1, times = nrow(d11)))

# bootstrap using case resampling
boot11 <- Boot(fit_nlsLM11, method = 'residual')

boot11_preds <- boot11$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d11$temp), max(d11$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot11_conf_preds <- group_by(boot11_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp11 <- filter(d_preds, conc=="Corythuchaciliata")

# plot bootstrapped CIs
p11 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp11, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot11_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d11, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')

# calculate params with CIs

extra_params11 <- calc_params(fit_nlsLM11) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params11 <- Boot(fit_nlsLM11, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM11)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params11 <- left_join(ci_extra_params11, extra_params11)
#> Joining, by = "param"

#ggplot(ci_extra_params11, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt11 <- as_tibble(ci_extra_params11[2,])
topt11$species <- as.character("Corythuchaciliata")

ggplot(topt11, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())
#________________________________________________________________

d12 <- filter(d_rate, conc=="Culexannulirostris")

fit_nlsLM12 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                 data = d12,
                                 start = coef(d_fits$sharpeschoolhigh[[12]]),
                                 lower = get_lower_lims(d12$temp, d12$rate, model_name = 'sharpeschoolhigh_1981'),
                                 upper = get_upper_lims(d12$temp, d12$rate, model_name = 'sharpeschoolhigh_1981'),
                                 weights = rep(1, times = nrow(d12)))

# bootstrap using case resampling
boot12 <- Boot(fit_nlsLM12, method = 'residual')

boot12_preds <- boot12$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d12$temp), max(d12$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot12_conf_preds <- group_by(boot12_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp12 <- filter(d_preds, conc=="Culexannulirostris")

# plot bootstrapped CIs
p12 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp12, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot12_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d12, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')

# calculate params with CIs

extra_params12 <- calc_params(fit_nlsLM12) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params12 <- Boot(fit_nlsLM12, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM12)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params12 <- left_join(ci_extra_params12, extra_params12)
#> Joining, by = "param"

#ggplot(ci_extra_params12, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt12 <- as_tibble(ci_extra_params12[2,])
topt12$species <- as.character("Culexannulirostris")

ggplot(topt12, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())
#________________________________________________________________

d13 <- filter(d_rate, conc=="Cydiapomonella")

fit_nlsLM13 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                 data = d13,
                                 start = coef(d_fits$sharpeschoolhigh[[13]]),
                                 lower = get_lower_lims(d13$temp, d13$rate, model_name = 'sharpeschoolhigh_1981'),
                                 upper = get_upper_lims(d13$temp, d13$rate, model_name = 'sharpeschoolhigh_1981'),
                                 weights = rep(1, times = nrow(d13)))

# bootstrap using case resampling
boot13 <- Boot(fit_nlsLM13, method = 'residual')

boot13_preds <- boot13$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d13$temp), max(d13$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot13_conf_preds <- group_by(boot13_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp13 <- filter(d_preds, conc=="Cydiapomonella")

# plot bootstrapped CIs
p13 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp13, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot13_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d13, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')

# calculate params with CIs

extra_params13 <- calc_params(fit_nlsLM13) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params13 <- Boot(fit_nlsLM13, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM13)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params13 <- left_join(ci_extra_params13, extra_params13)
#> Joining, by = "param"

#ggplot(ci_extra_params13, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt13 <- as_tibble(ci_extra_params13[2,])
topt13$species <- as.character("Cydiapomonella")

ggplot(topt13, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())

d_rate$conc

##############################################################################
###############             d14 does not work ################################
# #________________________________________________________________
# 
# d14 <- filter(d_rate, conc=="Euplectrusronnai")
# 
# fit_nlsLM14 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
#                                  data = d14,
#                                  start = coef(d_fits$sharpeschoolhigh[[14]]),
#                                  lower = get_lower_lims(d14$temp, d14$rate, model_name = 'sharpeschoolhigh_1981'),
#                                  upper = get_upper_lims(d14$temp, d14$rate, model_name = 'sharpeschoolhigh_1981'),
#                                  weights = rep(1, times = nrow(d14)))
# 
# # bootstrap using case resampling
# boot14 <- Boot(fit_nlsLM14, method = 'residual')
# 
# boot14_preds <- boot14$t %>%
#   as.data.frame() %>%
#   drop_na() %>%
#   mutate(iter = 1:n()) %>%
#   group_by_all() %>%
#   do(data.frame(temp = seq(min(d14$temp), max(d14$temp), length.out = 140))) %>%
#   ungroup() %>%
#   mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))
# 
# # calculate bootstrapped confidence intervals
# boot14_conf_preds <- group_by(boot14_preds, temp) %>%
#   summarise(conf_lower = quantile(pred, 0.025),
#             conf_upper = quantile(pred, 0.975)) %>%
#   ungroup()
# 
# 
# dp14 <- filter(d_preds, conc=="Euplectrusronnai")
# 
# # plot bootstrapped CIs
# p14 <-  ggplot() +
#   geom_line(aes(temp, .fitted), dp14, col = 'blue') +
#   geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot14_conf_preds, fill = 'blue', alpha = 0.3) +
#   geom_point(aes(temp, rate), d14, size = 2, alpha = 0.5) +
#   theme_bw(base_size = 12) +
#   labs(x = 'Temperature (?C)',
#        y = 'dev rate')
# 
# # calculate params with CIs
# 
# extra_params14 <- calc_params(fit_nlsLM14) %>%
#   pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
# 
# ci_extra_params14 <- Boot(fit_nlsLM14, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM14)), R = 200, method = 'residual') %>%
#   confint(., method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'residual bootstrap')
# 
# ci_extra_params14 <- left_join(ci_extra_params14, extra_params14)
# #> Joining, by = "param"
# 
# #ggplot(ci_extra_params14, aes(param, estimate)) +
# #geom_point(size = 4) +
# #geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
# #theme_bw() +
# #facet_wrap(~param, scales = 'free') +
# #scale_x_discrete('') +
# #labs(title = 'Calculation of confidence intervals for extra parameters',
# #subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')
# 
# # topt
# 
# topt14 <- as_tibble(ci_extra_params14[2,])
# topt14$species <- as.character("Euplectrusronnai")
# 
# ggplot(topt14, aes(estimate, species)) +
#   geom_point(size = 4) + 
#   geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
#   theme_bw(base_size = 12) +
#   facet_wrap(~species, scales = 'free') +
#   scale_x_continuous('') +
#   labs(title = 'calculation of Topt with CIs',
#        subtitle = 'dev rate TPC; using residual resampling')+
#   theme(axis.title.y = element_blank())
# 
# d_rate$conc

#________________________________________________________________

d15 <- filter(d_rate, conc=="Glyptapantelesmuesebecki")

fit_nlsLM15 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                 data = d15,
                                 start = coef(d_fits$sharpeschoolhigh[[15]]),
                                 lower = get_lower_lims(d15$temp, d15$rate, model_name = 'sharpeschoolhigh_1981'),
                                 upper = get_upper_lims(d15$temp, d15$rate, model_name = 'sharpeschoolhigh_1981'),
                                 weights = rep(1, times = nrow(d15)))

# bootstrap using case resampling
boot15 <- Boot(fit_nlsLM15, method = 'residual')

boot15_preds <- boot15$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d15$temp), max(d15$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot15_conf_preds <- group_by(boot15_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp15 <- filter(d_preds, conc=="Glyptapantelesmuesebecki")

# plot bootstrapped CIs
p15 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp15, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot15_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d15, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')

# calculate params with CIs

extra_params15 <- calc_params(fit_nlsLM15) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params15 <- Boot(fit_nlsLM15, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM15)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params15 <- left_join(ci_extra_params15, extra_params15)
#> Joining, by = "param"

#ggplot(ci_extra_params15, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt15 <- as_tibble(ci_extra_params15[2,])
topt15$species <- as.character("Glyptapantelesmuesebecki")

ggplot(topt15, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())

d_rate$conc
#________________________________________________________________

d16 <- filter(d_rate, conc=="Lepinotusreticulatus")

fit_nlsLM16 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                 data = d16,
                                 start = coef(d_fits$sharpeschoolhigh[[16]]),
                                 lower = get_lower_lims(d16$temp, d16$rate, model_name = 'sharpeschoolhigh_1981'),
                                 upper = get_upper_lims(d16$temp, d16$rate, model_name = 'sharpeschoolhigh_1981'),
                                 weights = rep(1, times = nrow(d16)))

# bootstrap using case resampling
boot16 <- Boot(fit_nlsLM16, method = 'residual')

boot16_preds <- boot16$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d16$temp), max(d16$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot16_conf_preds <- group_by(boot16_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp16 <- filter(d_preds, conc=="Lepinotusreticulatus")

# plot bootstrapped CIs
p16 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp16, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot16_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d16, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')

# calculate params with CIs

extra_params16 <- calc_params(fit_nlsLM16) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params16 <- Boot(fit_nlsLM16, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM16)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params16 <- left_join(ci_extra_params16, extra_params16)
#> Joining, by = "param"

#ggplot(ci_extra_params16, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt16 <- as_tibble(ci_extra_params16[2,])
topt16$species <- as.character("Lepinotusreticulatus")

ggplot(topt16, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())

d_rate$conc
#________________________________________________________________

d17 <- filter(d_rate, conc=="Planococcuscitri")

fit_nlsLM17 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                 data = d17,
                                 start = coef(d_fits$sharpeschoolhigh[[17]]),
                                 lower = get_lower_lims(d17$temp, d17$rate, model_name = 'sharpeschoolhigh_1981'),
                                 upper = get_upper_lims(d17$temp, d17$rate, model_name = 'sharpeschoolhigh_1981'),
                                 weights = rep(1, times = nrow(d17)))

# bootstrap using case resampling
boot17 <- Boot(fit_nlsLM17, method = 'residual')

boot17_preds <- boot17$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d17$temp), max(d17$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot17_conf_preds <- group_by(boot17_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp17 <- filter(d_preds, conc=="Planococcuscitri")

# plot bootstrapped CIs
p17 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp17, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot17_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d17, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')

# calculate params with CIs

extra_params17 <- calc_params(fit_nlsLM17) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params17 <- Boot(fit_nlsLM17, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM17)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params17 <- left_join(ci_extra_params17, extra_params17)
#> Joining, by = "param"

#ggplot(ci_extra_params17, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt17 <- as_tibble(ci_extra_params17[2,])
topt17$species <- as.character("Planococcuscitri")

ggplot(topt17, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())

d_rate$conc
#________________________________________________________________

d18 <- filter(d_rate, conc=="Procambarusclarkii")

fit_nlsLM18 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                 data = d18,
                                 start = coef(d_fits$sharpeschoolhigh[[18]]),
                                 lower = get_lower_lims(d18$temp, d18$rate, model_name = 'sharpeschoolhigh_1981'),
                                 upper = get_upper_lims(d18$temp, d18$rate, model_name = 'sharpeschoolhigh_1981'),
                                 weights = rep(1, times = nrow(d18)))

# bootstrap using case resampling
boot18 <- Boot(fit_nlsLM18, method = 'residual')

boot18_preds <- boot18$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d18$temp), max(d18$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot18_conf_preds <- group_by(boot18_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp18 <- filter(d_preds, conc=="Procambarusclarkii")

# plot bootstrapped CIs
p18 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp18, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot18_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d18, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')

# calculate params with CIs

extra_params18 <- calc_params(fit_nlsLM18) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params18 <- Boot(fit_nlsLM18, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM18)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params18 <- left_join(ci_extra_params18, extra_params18)
#> Joining, by = "param"

#ggplot(ci_extra_params18, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt18 <- as_tibble(ci_extra_params18[2,])
topt18$species <- as.character("Procambarusclarkii")

ggplot(topt18, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())

d_rate$conc
#________________________________________________________________

d19 <- filter(d_rate, conc=="Sitonadiscoideus")

fit_nlsLM19 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                 data = d19,
                                 start = coef(d_fits$sharpeschoolhigh[[19]]),
                                 lower = get_lower_lims(d19$temp, d19$rate, model_name = 'sharpeschoolhigh_1981'),
                                 upper = get_upper_lims(d19$temp, d19$rate, model_name = 'sharpeschoolhigh_1981'),
                                 weights = rep(1, times = nrow(d19)))

# bootstrap using case resampling
boot19 <- Boot(fit_nlsLM19, method = 'residual')

boot19_preds <- boot19$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d19$temp), max(d19$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot19_conf_preds <- group_by(boot19_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp19 <- filter(d_preds, conc=="Sitonadiscoideus")

# plot bootstrapped CIs
p19 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp19, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot19_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d19, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')

# calculate params with CIs

extra_params19 <- calc_params(fit_nlsLM19) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params19 <- Boot(fit_nlsLM19, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM19)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params19 <- left_join(ci_extra_params19, extra_params19)
#> Joining, by = "param"

#ggplot(ci_extra_params19, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt19 <- as_tibble(ci_extra_params19[2,])
topt19$species <- as.character("Sitonadiscoideus")

ggplot(topt19, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())

d_rate$conc
#________________________________________________________________

d20 <- filter(d_rate, conc=="Telenomuschrysopae")

fit_nlsLM20 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                 data = d20,
                                 start = coef(d_fits$sharpeschoolhigh[[20]]),
                                 lower = get_lower_lims(d20$temp, d20$rate, model_name = 'sharpeschoolhigh_1981'),
                                 upper = get_upper_lims(d20$temp, d20$rate, model_name = 'sharpeschoolhigh_1981'),
                                 weights = rep(1, times = nrow(d20)))

# bootstrap using case resampling
boot20 <- Boot(fit_nlsLM20, method = 'residual')

boot20_preds <- boot20$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d20$temp), max(d20$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot20_conf_preds <- group_by(boot20_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp20 <- filter(d_preds, conc=="Telenomuschrysopae")

# plot bootstrapped CIs
p20 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp20, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot20_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d20, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')

# calculate params with CIs

extra_params20 <- calc_params(fit_nlsLM20) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params20 <- Boot(fit_nlsLM20, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM20)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params20 <- left_join(ci_extra_params20, extra_params20)
#> Joining, by = "param"

#ggplot(ci_extra_params20, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt20 <- as_tibble(ci_extra_params20[2,])
topt20$species <- as.character("Telenomuschrysopae")

ggplot(topt20, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())

d_rate$conc
#________________________________________________________________

d21 <- filter(d_rate, conc=="Tetraneuranigriabdominalis")

fit_nlsLM21 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                 data = d21,
                                 start = coef(d_fits$sharpeschoolhigh[[21]]),
                                 lower = get_lower_lims(d21$temp, d21$rate, model_name = 'sharpeschoolhigh_1981'),
                                 upper = get_upper_lims(d21$temp, d21$rate, model_name = 'sharpeschoolhigh_1981'),
                                 weights = rep(1, times = nrow(d21)))

# bootstrap using case resampling
boot21 <- Boot(fit_nlsLM21, method = 'residual')

boot21_preds <- boot21$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d21$temp), max(d21$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot21_conf_preds <- group_by(boot21_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp21 <- filter(d_preds, conc=="Tetraneuranigriabdominalis")

# plot bootstrapped CIs
p21 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp21, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot21_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d21, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')

# calculate params with CIs

extra_params21 <- calc_params(fit_nlsLM21) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params21 <- Boot(fit_nlsLM21, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM21)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params21 <- left_join(ci_extra_params21, extra_params21)
#> Joining, by = "param"

#ggplot(ci_extra_params21, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt21 <- as_tibble(ci_extra_params21[2,])
topt21$species <- as.character("Tetraneuranigriabdominalis")

ggplot(topt21, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())

d_rate$conc
#________________________________________________________________

d22 <- filter(d_rate, conc=="Theocolaxelegans")

fit_nlsLM22 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                 data = d22,
                                 start = coef(d_fits$sharpeschoolhigh[[22]]),
                                 lower = get_lower_lims(d22$temp, d22$rate, model_name = 'sharpeschoolhigh_1981'),
                                 upper = get_upper_lims(d22$temp, d22$rate, model_name = 'sharpeschoolhigh_1981'),
                                 weights = rep(1, times = nrow(d22)))

# bootstrap using case resampling
boot22 <- Boot(fit_nlsLM22, method = 'residual')

boot22_preds <- boot22$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d22$temp), max(d22$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot22_conf_preds <- group_by(boot22_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp22 <- filter(d_preds, conc=="Theocolaxelegans")

# plot bootstrapped CIs
p22 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp22, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot22_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d22, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')

# calculate params with CIs

extra_params22 <- calc_params(fit_nlsLM22) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params22 <- Boot(fit_nlsLM22, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM22)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params22 <- left_join(ci_extra_params22, extra_params22)
#> Joining, by = "param"

#ggplot(ci_extra_params22, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt22 <- as_tibble(ci_extra_params22[2,])
topt22$species <- as.character("Theocolaxelegans")

ggplot(topt22, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())

d_rate$conc
#_______________________________________________

d23 <- filter(d_rate, conc=="Trichogrammabruni")

fit_nlsLM23 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                data = d23,
                                start = coef(d_fits$sharpeschoolhigh[[23]]),
                                lower = get_lower_lims(d23$temp, d23$rate, model_name = 'sharpeschoolhigh_1981'),
                                upper = get_upper_lims(d23$temp, d23$rate, model_name = 'sharpeschoolhigh_1981'),
                                weights = rep(1, times = nrow(d23)))

# bootstrap using case resampling
boot23 <- Boot(fit_nlsLM23, method = 'residual')

boot23_preds <- boot23$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d23$temp), max(d23$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot23_conf_preds <- group_by(boot23_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp23 <- filter(d_preds, conc=="Trichogrammabruni")

# plot bootstrapped CIs
p23 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp23, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot23_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d23, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')


# calculate params with CIs

extra_params23 <- calc_params(fit_nlsLM23) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params23 <- Boot(fit_nlsLM23, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM23)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params23 <- left_join(ci_extra_params23, extra_params23)
#> Joining, by = "param"

#ggplot(ci_extra_params23, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt23 <- as_tibble(ci_extra_params23[2,])
topt23$species <- as.character("Trichogrammabruni")

ggplot(topt23, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())

d_rate$conc


###############################################################################
######################## d24 does not work #####################################
###############################################################################
# #_______________________________________________
# 
# d24 <- filter(d_rate, conc=="Trichogrammasp.nr.Lutea")
# 
# fit_nlsLM24 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
#                                  data = d24,
#                                  start = coef(d_fits$sharpeschoolhigh[[24]]),
#                                  lower = get_lower_lims(d24$temp, d24$rate, model_name = 'sharpeschoolhigh_1981'),
#                                  upper = get_upper_lims(d24$temp, d24$rate, model_name = 'sharpeschoolhigh_1981'),
#                                  weights = rep(1, times = nrow(d24)))
# 
# # bootstrap using case resampling
# boot24 <- Boot(fit_nlsLM24, method = 'residual')
# 
# boot24_preds <- boot24$t %>%
#   as.data.frame() %>%
#   drop_na() %>%
#   mutate(iter = 1:n()) %>%
#   group_by_all() %>%
#   do(data.frame(temp = seq(min(d24$temp), max(d24$temp), length.out = 100))) %>%
#   ungroup() %>%
#   mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))
# 
# # calculate bootstrapped confidence intervals
# boot24_conf_preds <- group_by(boot24_preds, temp) %>%
#   summarise(conf_lower = quantile(pred, 0.025),
#             conf_upper = quantile(pred, 0.975)) %>%
#   ungroup()
# 
# 
# dp24 <- filter(d_preds, conc=="Trichogrammasp.nr.Lutea")
# 
# # plot bootstrapped CIs
# p24 <-  ggplot() +
#   geom_line(aes(temp, .fitted), dp24, col = 'blue') +
#   geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot24_conf_preds, fill = 'blue', alpha = 0.3) +
#   geom_point(aes(temp, rate), d24, size = 2, alpha = 0.5) +
#   theme_bw(base_size = 12) +
#   labs(x = 'Temperature (?C)',
#        y = 'dev rate')
# 
# 
# # calculate params with CIs
# 
# extra_params24 <- calc_params(fit_nlsLM24) %>%
#   pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
# 
# ci_extra_params24 <- Boot(fit_nlsLM24, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM24)), R = 200, method = 'residual') %>%
#   confint(., method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'residual bootstrap')
# 
# ci_extra_params24 <- left_join(ci_extra_params24, extra_params24)
# #> Joining, by = "param"
# 
# #ggplot(ci_extra_params24, aes(param, estimate)) +
# #geom_point(size = 4) +
# #geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
# #theme_bw() +
# #facet_wrap(~param, scales = 'free') +
# #scale_x_discrete('') +
# #labs(title = 'Calculation of confidence intervals for extra parameters',
# #subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')
# 
# # topt
# 
# topt24 <- as_tibble(ci_extra_params24[2,])
# topt24$species <- as.character("Trichogrammasp.nr.Lutea")
# 
# ggplot(topt24, aes(estimate, species)) +
#   geom_point(size = 4) + 
#   geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
#   theme_bw(base_size = 12) +
#   facet_wrap(~species, scales = 'free') +
#   scale_x_continuous('') +
#   labs(title = 'calculation of Topt with CIs',
#        subtitle = 'dev rate TPC; using residual resampling')+
#   theme(axis.title.y = element_blank())
# 
# d_rate$conc

#_______________________________________________

d25 <- filter(d_rate, conc=="Trichogrammasp.nr.Mwanzai")

fit_nlsLM25 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                 data = d25,
                                 start = coef(d_fits$sharpeschoolhigh[[25]]),
                                 lower = get_lower_lims(d25$temp, d25$rate, model_name = 'sharpeschoolhigh_1981'),
                                 upper = get_upper_lims(d25$temp, d25$rate, model_name = 'sharpeschoolhigh_1981'),
                                 weights = rep(1, times = nrow(d25)))

# bootstrap using case resampling
boot25 <- Boot(fit_nlsLM25, method = 'residual')

boot25_preds <- boot25$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d25$temp), max(d25$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot25_conf_preds <- group_by(boot25_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp25 <- filter(d_preds, conc=="Trichogrammasp.nr.Mwanzai")

# plot bootstrapped CIs
p25 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp25, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot25_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d25, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')


# calculate params with CIs

extra_params25 <- calc_params(fit_nlsLM25) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params25 <- Boot(fit_nlsLM25, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM25)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params25 <- left_join(ci_extra_params25, extra_params25)
#> Joining, by = "param"

#ggplot(ci_extra_params25, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt25 <- as_tibble(ci_extra_params25[2,])
topt25$species <- as.character("Trichogrammasp.nr.Mwanzai")

ggplot(topt25, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())

d_rate$conc

#_______________________________________________

d26 <- filter(d_rate, conc=="Trichogrammabruni")

fit_nlsLM26 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                 data = d26,
                                 start = coef(d_fits$sharpeschoolhigh[[26]]),
                                 lower = get_lower_lims(d26$temp, d26$rate, model_name = 'sharpeschoolhigh_1981'),
                                 upper = get_upper_lims(d26$temp, d26$rate, model_name = 'sharpeschoolhigh_1981'),
                                 weights = rep(1, times = nrow(d26)))

# bootstrap using case resampling
boot26 <- Boot(fit_nlsLM26, method = 'residual')

boot26_preds <- boot26$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d26$temp), max(d26$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot26_conf_preds <- group_by(boot26_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp26 <- filter(d_preds, conc=="Trichogrammabruni")

# plot bootstrapped CIs
p26 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp26, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot26_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d26, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')


# calculate params with CIs

extra_params26 <- calc_params(fit_nlsLM26) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params26 <- Boot(fit_nlsLM26, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM26)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params26 <- left_join(ci_extra_params26, extra_params26)
#> Joining, by = "param"

#ggplot(ci_extra_params26, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt26 <- as_tibble(ci_extra_params26[2,])
topt26$species <- as.character("Trichogrammabruni")

ggplot(topt26, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())

d_rate$conc
#####Telenomuslobatus
#_______________________________________________

d27 <- filter(d_rate, conc=="Telenomuslobatus")

fit_nlsLM27 <- minpack.lm::nlsLM(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                 data = d27,
                                 start = coef(d_fits$sharpeschoolhigh[[27]]),
                                 lower = get_lower_lims(d27$temp, d27$rate, model_name = 'sharpeschoolhigh_1981'),
                                 upper = get_upper_lims(d27$temp, d27$rate, model_name = 'sharpeschoolhigh_1981'),
                                 weights = rep(1, times = nrow(d27)))

# bootstrap using case resampling
boot27 <- Boot(fit_nlsLM27, method = 'residual')

boot27_preds <- boot27$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(d27$temp), max(d27$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = sharpeschoolhigh_1981(temp, r_tref, e, eh, th, tref = 15))

# calculate bootstrapped confidence intervals
boot27_conf_preds <- group_by(boot27_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


dp27 <- filter(d_preds, conc=="Telenomuslobatus")

# plot bootstrapped CIs
p27 <-  ggplot() +
  geom_line(aes(temp, .fitted), dp27, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot27_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), d27, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')


# calculate params with CIs

extra_params27 <- calc_params(fit_nlsLM27) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params27 <- Boot(fit_nlsLM27, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM27)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params27 <- left_join(ci_extra_params27, extra_params27)
#> Joining, by = "param"

#ggplot(ci_extra_params27, aes(param, estimate)) +
#geom_point(size = 4) +
#geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
#theme_bw() +
#facet_wrap(~param, scales = 'free') +
#scale_x_discrete('') +
#labs(title = 'Calculation of confidence intervals for extra parameters',
#subtitle = 'Macrocentrusiridescens dev rate TPC; using residual resampling')

# topt

topt27 <- as_tibble(ci_extra_params27[2,])
topt27$species <- as.character("Telenomuslobatus")

ggplot(topt27, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  facet_wrap(~species, scales = 'free') +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())

d_rate$conc

#________________________________________________________________

topts <- rbind(topt1,topt2,topt3,topt4,topt5,topt6,topt7,topt8,topt10, topt11, topt12, topt13, topt15, topt16,
               topt17, topt18, topt19, topt20, topt21, topt22, topt23)
#________________________________________________________________

tpklab <- expression(italic(T)[opt])

ggplot(topts, aes(estimate, species, fill=species))+
  geom_errorbar(aes(xmin = conf_lower, xmax = conf_upper),width=0.25,size=0.1,col="#000000") +
  geom_point(size = 3, shape=21,stroke=1,col="#000000") + 
  theme_bw(base_size = 12) +
  scale_x_continuous(tpklab) +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme_bw(base_size = 18)+
  theme(legend.position = "none")+
  theme(axis.title.y = element_blank())


