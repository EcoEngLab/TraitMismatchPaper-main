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

#setwd("C:/Users/pjh215/Dropbox/ph_thesis/Topt_paper")
setwd("/home/primuser/Documents/VByte/VecMismatchPaper1/code/")

#read in the trait data

df <- as_tibble(read.csv('../data/Final_Traitofinterest_ph.csv'))

df <- df %>% filter(mismatch_keep =="yes")

df <- df %>%
  select('originaltraitname','originaltraitvalue','originaltraitunit', 
         'standardisedtraitname','standardisedtraitvalue','standardisedtraitunit',
         'interactor1', 'ambienttemp', 'citation')


df<- df %>%
  rename(temp = ambienttemp,species = interactor1)


#_____________________________________________________#
#               The mosquitoes                        #
#_____________________________________________________#

##                Aedes aegypti

ae.ae <- df %>% filter(species == "Aedes aegypti", standardisedtraitname != "NA")

#- Development Rate

ae.dv <- ae.ae %>% filter(standardisedtraitname == "Juvenile Development Rate")

ae.dv <- ae.dv %>% select(standardisedtraitvalue,temp)
ae.dv <- ae.dv %>% rename(rate = standardisedtraitvalue)

ae.dv$rate <- ae.dv$rate + 10^-6


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(ae.dv$temp, ae.dv$rate, model_name = 'pawar_2018')

d_fits <- nest(ae.dv, data = c(temp, rate)) %>%
           mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                                                       data = .x,
                                                                       iter = c(3,3,3,3),
                                                                       start_lower = start_vals - 10,
                                                                       start_upper = start_vals + 10,
                                                                       supp_errors = 'Y',
                                                                       convergence_count = FALSE)))
#____________________


d_preds <- mutate(d_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot panel

ggplot(d_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, ae.dv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM1 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                               data = ae.dv,
                               start = coef(d_fits$pawar[[1]]),
                               lower = get_lower_lims(ae.dv$temp, ae.dv$rate, model_name = 'pawar_2018'),
                               upper = get_upper_lims(ae.dv$temp, ae.dv$rate, model_name = 'pawar_2018'),
                               weights = rep(1, times = nrow(ae.dv)))

# bootstrap using case resampling
boot1 <- Boot(fit_nlsLM1, method = 'residual')

boot1_preds <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(ae.dv$temp), max(ae.dv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot1_conf_preds <- group_by(boot1_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs
  
  ggplot() +
  geom_line(aes(temp, .fitted), d_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), ae.dv, size = 2, alpha = 0.5) +
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

#- Topt 

topt1 <- as_tibble(ci_extra_params1[2,])
topt1$species <- as.character("Aedes aegypti")

ggplot(topt1, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())


#____________________________________________________

# Juvenile mortality rate

ae.zj <- ae.ae %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

ae.zj <- ae.zj %>% select(standardisedtraitvalue,temp)
ae.zj <- ae.zj %>% rename(rate = standardisedtraitvalue)

ae.zj$rate <- 1/ae.zj$rate


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(ae.zj$temp, ae.zj$rate, model_name = 'pawar_2018')

a_fits <- nest(ae.zj, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


a_preds <- mutate(a_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot panel

ggplot(a_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, ae.zj)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM2 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = ae.zj,
                                start = coef(a_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(ae.zj)))

# bootstrap using case resampling
boot2 <- Boot(fit_nlsLM2, method = 'residual')

boot2_preds <- boot2$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(ae.zj$temp), max(ae.zj$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot2_conf_preds <- group_by(boot2_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

  ggplot() +
  geom_line(aes(temp, .fitted), a_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot2_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), ae.zj, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = '1/zj')

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

#- Topt 

topt2 <- as_tibble(ci_extra_params2[2,])
topt2$species <- as.character("Aedes aegypti")

ggplot(topt2, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())


#______________________________________

# Female mortality rate

ae.z <- ae.ae %>% filter(standardisedtraitname == "Female Mortality Rate")

ae.z <- ae.z %>% select(standardisedtraitvalue,temp)
ae.z <- ae.z %>% rename(rate = standardisedtraitvalue)

ae.z$rate <- 1/ae.z$rate


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(ae.z$temp, ae.z$rate, model_name = 'pawar_2018')

b_fits <- nest(ae.z, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


b_preds <- mutate(b_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot panel

ggplot(b_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, ae.z)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM3 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = ae.z,
                                start = coef(b_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(ae.z)))

# bootstrap using case resampling
boot3 <- Boot(fit_nlsLM3, method = 'residual')

boot3_preds <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(ae.z$temp), max(ae.z$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot3_conf_preds <- group_by(boot3_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

p3 <- ggplot() +
  geom_line(aes(temp, .fitted), b_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot3_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), ae.z, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = '1/z');p3

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

#- Topt 

topt3 <- as_tibble(ci_extra_params3[2,])
topt3$species <- as.character("Aedes aegypti")

ggplot(topt3, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'z TPC; using residual resampling')+
  theme(axis.title.y = element_blank())

ae.aeT <- rbind(topt1,topt2,topt3)

ae.aeT$param <- as.character(c("1/development time","juvenile mortality rate","adult mortality rate"))

ae.aeT$euler_param <- as.character(c("1/alpha","zj","z"))

ae.aeT$stage <- as.factor(c("juvenile","juvenile","adult"))
#ae.aeT$stage <- ordered(c("juvenile","juvenile", "adult"))

ae.aeT$stage <- relevel(ae.aeT$stage, "juvenile")

ggplot(ae.aeT, aes(estimate, species, shape = euler_param,colour=stage)) +
  geom_point(size = 2.5) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())+
  scale_color_brewer(type = 'qual', palette = 2)


#_________________________________________________

##               Aedes albopictus

ae.al <- df %>% filter(species == "Aedes albopictus", standardisedtraitname != "NA")

#- Development Rate

al.dv <- ae.al %>% filter(originaltraitname == "Juvenile Development Rate")

al.dv <- al.dv %>% select(originaltraitvalue,temp)
al.dv <- al.dv %>% rename(rate = originaltraitvalue)

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(al.dv$temp, al.dv$rate, model_name = 'pawar_2018')

f_fits <- nest(al.dv, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


f_preds <- mutate(f_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot 

ggplot(f_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, al.dv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM4 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = al.dv,
                                start = coef(f_fits$pawar[[1]]),
                                lower = get_lower_lims(al.dv$temp, al.dv$rate, model_name = 'pawar_2018'),
                                upper = get_upper_lims(al.dv$temp, al.dv$rate, model_name = 'pawar_2018'),
                                weights = rep(1, times = nrow(al.dv)))

# bootstrap using case resampling
boot4 <- Boot(fit_nlsLM4, method = 'residual')

boot4_preds <- boot4$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(al.dv$temp), max(al.dv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot4_conf_preds <- group_by(boot4_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

p4 <- ggplot() +
  geom_line(aes(temp, .fitted), f_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot4_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), al.dv, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate'); p4

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

#- Topt 

topt4 <- as_tibble(ci_extra_params4[2,])
topt4$species <- as.character("Aedes albopictus")

ggplot(topt4, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())


#____________________________________________________

# Juvenile mortality rate

al.zj <- ae.al %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

al.zj <- al.zj %>% select(standardisedtraitvalue,temp)
al.zj <- al.zj %>% rename(rate = standardisedtraitvalue)

al.zj$rate <- 1/al.zj$rate

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(al.zj$temp, al.zj$rate, model_name = 'pawar_2018')

g_fits <- nest(al.zj, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


g_preds <- mutate(g_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot 

ggplot(g_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, al.zj)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM5 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = al.zj,
                                start = coef(g_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(al.zj)))

# bootstrap using case resampling
boot5 <- Boot(fit_nlsLM5, method = 'residual')

boot5_preds <- boot5$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(al.zj$temp), max(al.zj$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot5_conf_preds <- group_by(boot5_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), g_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot5_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), al.zj, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = '1/zj')

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

#- Topt 

topt5 <- as_tibble(ci_extra_params2[2,])
topt5$species <- as.character("Aedes albopictus")

ggplot(topt5, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())


#______________________________________

# Fecundity rate

al.bpk <- ae.al %>% filter(originaltraitname == "Gross Reproductive Rate ")

al.bpk <- al.bpk %>% select(originaltraitvalue,temp)
al.bpk <- al.bpk %>% rename(rate = originaltraitvalue)

al.bpk$rate <- al.bpk$rate + 10^-6


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(al.bpk$temp, al.bpk$rate, model_name = 'pawar_2018')

t_fits <- nest(al.bpk, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


t_preds <- mutate(t_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot panel

ggplot(t_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, al.bpk)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'Rate')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM6 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = al.bpk,
                                start = coef(t_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(al.bpk)))

# bootstrap using case resampling
boot6 <- Boot(fit_nlsLM6, method = 'residual')

boot6_preds <- boot6$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(al.bpk$temp), max(al.bpk$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot6_conf_preds <- group_by(boot6_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), t_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot6_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), al.bpk, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'bpk')

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

#- Topt 

topt6 <- as_tibble(ci_extra_params6[2,])
topt6$species <- as.character("Aedes albopictus")

ggplot(topt6, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())

ae.alT <- rbind(topt4,topt5,topt6)

ae.alT$param <- as.character(c("1/development time","juvenile mortality rate","fecundity rate"))

ae.alT$euler_param <- as.character(c("1/alpha","zj","bpk"))

ae.alT$stage <- as.factor(c("juvenile","juvenile","adult"))

ae.alT$stage <- relevel(ae.alT$stage, "juvenile")

aedes <- rbind(ae.aeT,ae.alT)

ggplot(aedes, aes(estimate, species, shape = euler_param,colour=stage)) +
  geom_point(size = 2.5) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())+
  scale_color_brewer(type = 'qual', palette = 2)


#_________________________________________________

#_________________________________________________

##               Culex annulirostris            ##

cl.al <- df %>% filter(species == "Culex annulirostris", standardisedtraitname != "NA")

#- Development Rate

cl.dv <- cl.al %>% filter(standardisedtraitname == "Juvenile Development Rate")

cl.dv <- cl.dv %>% select(standardisedtraitvalue,temp)
cl.dv <- cl.dv %>% rename(rate = standardisedtraitvalue)

cl.dv$rate <- cl.dv$rate + 10^-6

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(cl.dv$temp, cl.dv$rate, model_name = 'pawar_2018')

z_fits <- nest(cl.dv, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


z_preds <- mutate(z_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot 

ggplot(z_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, cl.dv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM7 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = cl.dv,
                                start = coef(z_fits$pawar[[1]]),
                                lower = get_lower_lims(cl.dv$temp, cl.dv$rate, model_name = 'pawar_2018'),
                                upper = get_upper_lims(cl.dv$temp, cl.dv$rate, model_name = 'pawar_2018'),
                                weights = rep(1, times = nrow(cl.dv)))

# bootstrap using case resampling
boot7 <- Boot(fit_nlsLM7, method = 'residual')

boot7_preds <- boot7$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(cl.dv$temp), max(cl.dv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot7_conf_preds <- group_by(boot7_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

# plot bootstrapped CIs

  ggplot() +
  geom_line(aes(temp, .fitted), z_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot7_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), cl.dv, size = 2, alpha = 0.5) +
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

#- Topt 

topt7 <- as_tibble(ci_extra_params7[2,])
topt7$species <- as.character("Culex annulirostris")

ggplot(topt7, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())


#____________________________________________________

# Juvenile mortality rate

cl.zj <- cl.al %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

cl.zj <- cl.zj %>% select(standardisedtraitvalue,temp)
cl.zj <- cl.zj %>% rename(rate = standardisedtraitvalue)

cl.zj$rate <- 1/cl.zj$rate

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(cl.zj$temp, cl.zj$rate, model_name = 'pawar_2018')

s_fits <- nest(cl.zj, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


s_preds <- mutate(s_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot 

ggplot(s_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, cl.zj)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM8 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = cl.zj,
                                start = coef(s_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(cl.zj)))

# bootstrap using case resampling
boot8 <- Boot(fit_nlsLM8, method = 'residual')

boot8_preds <- boot8$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(cl.zj$temp), max(cl.zj$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot8_conf_preds <- group_by(boot8_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

  ggplot() +
  geom_line(aes(temp, .fitted), s_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot8_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), cl.zj, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = '1/zj')

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

#- Topt 

topt8 <- as_tibble(ci_extra_params8[2,])
topt8$species <- as.character("Culex annulirostris")

ggplot(topt8, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  labs(title = 'calculation of Topt with CIs',
       subtitle = 'dev rate TPC; using residual resampling')+
  theme(axis.title.y = element_blank())


#______________________________________

# Fecundity rate

cl.bpk <- cl.al %>% filter(standardisedtraitname == "Average Fecundity Rate")

cl.bpk <- cl.bpk %>% select(originaltraitvalue,temp)
cl.bpk <- cl.bpk %>% rename(rate = originaltraitvalue)

cl.bpk$rate <- cl.bpk$rate + 10^-6


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(cl.bpk$temp, cl.bpk$rate, model_name = 'pawar_2018')

p_fits <- nest(cl.bpk, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


p_preds <- mutate(p_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot panel

ggplot(p_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, cl.bpk)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'Rate')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM9 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = cl.bpk,
                                start = coef(p_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(cl.bpk)))

# bootstrap using case resampling
boot9 <- Boot(fit_nlsLM9, method = 'residual')

boot9_preds <- boot9$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(cl.bpk$temp), max(cl.bpk$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot9_conf_preds <- group_by(boot9_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), p_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot9_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), cl.bpk, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'bpk')

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

#- Topt 

topt9 <- as_tibble(ci_extra_params9[2,])
topt9$species <- as.character("Culex annulirostris")

ggplot(topt9, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())


cl.alT <- rbind(topt7,topt8,topt9)

cl.alT$param <- as.character(c("1/development time","juvenile mortality rate","fecundity rate"))

cl.alT$euler_param <- as.character(c("1/alpha","zj","bpk"))

cl.alT$stage <- as.factor(c("juvenile","juvenile","adult"))

cl.alT$stage <- relevel(cl.alT$stage, "juvenile")

moz <- rbind(aedes,cl.alT)

ggplot(moz, aes(estimate, species, shape = euler_param,colour=stage)) +
  geom_point(size = 2.5) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())+
  scale_color_brewer(type = 'qual', palette = 2)


#_________________________________________________

#_________________________________________________

##            Anopheles gambiae s.s.

an.g <- df %>% filter(species == "Anopheles gambiae s.s.", standardisedtraitname != "NA")

#- Development Rate

an.dv <- an.g %>% filter(standardisedtraitname == "Juvenile Development Rate")

an.dv <- an.dv %>% select(standardisedtraitvalue,temp)
an.dv <- an.dv %>% rename(rate = standardisedtraitvalue)

an.dv$rate <- an.dv$rate + 10^-6

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(an.dv$temp, an.dv$rate, model_name = 'pawar_2018')

e_fits <- nest(an.dv, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


e_preds <- mutate(e_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot 

ggplot(e_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, an.dv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM10 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = an.dv,
                                start = coef(e_fits$pawar[[1]]),
                                lower = get_lower_lims(an.dv$temp, an.dv$rate, model_name = 'pawar_2018'),
                                upper = get_upper_lims(an.dv$temp, an.dv$rate, model_name = 'pawar_2018'),
                                weights = rep(1, times = nrow(an.dv)))

# bootstrap using case resampling
boot10 <- Boot(fit_nlsLM10, method = 'residual')

boot10_preds <- boot10$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(an.dv$temp), max(an.dv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot10_conf_preds <- group_by(boot10_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), e_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot10_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), an.dv, size = 2, alpha = 0.5) +
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

#- Topt 

topt10 <- as_tibble(ci_extra_params10[2,])
topt10$species <- as.character("Anopheles gambiae s.s.")

ggplot(topt10, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())






#____________________________________________________

# Juvenile mortality rate

an.zj <- an.g %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

an.zj <- an.zj %>% select(standardisedtraitvalue,temp)
an.zj <- an.zj %>% rename(rate = standardisedtraitvalue)

an.zj$rate <- 1/an.zj$rate

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(an.zj$temp, an.zj$rate, model_name = 'pawar_2018')

m_fits <- nest(an.zj, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


m_preds <- mutate(m_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot 

ggplot(m_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, an.zj)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM11 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = an.zj,
                                start = coef(m_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(an.zj)))

# bootstrap using case resampling
boot11 <- Boot(fit_nlsLM11, method = 'residual')

boot11_preds <- boot11$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(an.zj$temp), max(an.zj$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot11_conf_preds <- group_by(boot11_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), m_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot11_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), an.zj, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = '1/zj')

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

#- Topt 

topt11 <- as_tibble(ci_extra_params11[2,])
topt11$species <- as.character("Anopheles gambiae s.s.")

ggplot(topt11, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())

#--

an.gT <- rbind(topt10,topt11)

an.gT$param <- as.character(c("1/development time","juvenile mortality rate"))

an.gT$euler_param <- as.character(c("1/alpha","zj"))

an.gT$stage <- as.factor(c("juvenile","juvenile"))

an.gT$stage <- relevel(an.gT$stage, "juvenile")

mozz <- rbind(moz,an.gT)

ggplot(mozz, aes(estimate, species, shape = euler_param,colour=stage)) +
  geom_point(size = 2.5) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())+
  scale_color_brewer(type = 'qual', palette = 2)

#________________________________________________________

#                   Aphis gossypii                      #

#________________________________________________________


ap.g <- df %>% filter(species == "Aphis gossypii", standardisedtraitname != "NA")

#- Development Rate

ap.dv <- ap.g %>% filter(standardisedtraitname == "Juvenile Development Rate")

ap.dv <- ap.dv %>% select(standardisedtraitvalue,temp)
ap.dv <- ap.dv %>% rename(rate = standardisedtraitvalue)

ap.dv$rate <- ap.dv$rate + 10^-6

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(ap.dv$temp, ap.dv$rate, model_name = 'pawar_2018')

v_fits <- nest(ap.dv, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


v_preds <- mutate(v_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot 

ggplot(v_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, ap.dv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM12 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = ap.dv,
                                start = coef(v_fits$pawar[[1]]),
                                lower = get_lower_lims(ap.dv$temp, ap.dv$rate, model_name = 'pawar_2018'),
                                upper = get_upper_lims(ap.dv$temp, ap.dv$rate, model_name = 'pawar_2018'),
                                weights = rep(1, times = nrow(ap.dv)))

# bootstrap using case resampling
boot12 <- Boot(fit_nlsLM12, method = 'residual')

boot12_preds <- boot12$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(ap.dv$temp), max(ap.dv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot12_conf_preds <- group_by(boot12_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), v_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot12_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), ap.dv, size = 2, alpha = 0.5) +
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

#- Topt 

topt12 <- as_tibble(ci_extra_params12[2,])
topt12$species <- as.character("Aphis gossypii")

ggplot(topt12, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())


#____________________________________________________

# Juvenile mortality rate

ap.zj <- ap.g %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

ap.zj <- ap.zj %>% select(standardisedtraitvalue,temp)
ap.zj <- ap.zj %>% rename(rate = standardisedtraitvalue)

ap.zj$rate <- 1/ap.zj$rate

ap.zj<-ap.zj[-14,]

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(ap.zj$temp, ap.zj$rate, model_name = 'pawar_2018')

h_fits <- nest(ap.zj, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


h_preds <- mutate(h_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot 

ggplot(h_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, ap.zj)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM13 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = ap.zj,
                                start = coef(h_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(ap.zj)))

# bootstrap using case resampling
boot13 <- Boot(fit_nlsLM13, method = 'residual')

boot13_preds <- boot13$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(ap.zj$temp), max(ap.zj$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot13_conf_preds <- group_by(boot13_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), h_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot13_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), ap.zj, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = '1/zj')

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

#- Topt 

topt13 <- as_tibble(ci_extra_params13[2,])
topt13$species <- as.character("Aphis gossypii")

ggplot(topt13, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())


#______________________________________

# Fecundity rate

ap.bpk <- ap.g %>% filter(standardisedtraitname == "Fecundity Rate")

ap.bpk <- ap.bpk %>% select(standardisedtraitvalue,temp)
ap.bpk <- ap.bpk %>% rename(rate = standardisedtraitvalue)

ap.bpk$rate <- ap.bpk$rate + 10^-6


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(ap.bpk$temp, ap.bpk$rate, model_name = 'pawar_2018')

j_fits <- nest(ap.bpk, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


j_preds <- mutate(j_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot panel

ggplot(j_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, ap.bpk)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'Rate')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM14 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = ap.bpk,
                                start = coef(j_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(ap.bpk)))

# bootstrap using case resampling
boot14 <- Boot(fit_nlsLM14, method = 'residual')

boot14_preds <- boot14$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(ap.bpk$temp), max(ap.bpk$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot14_conf_preds <- group_by(boot14_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), j_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot14_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), ap.bpk, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'bpk')

# calculate params with CIs

extra_params14 <- calc_params(fit_nlsLM14) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params14 <- Boot(fit_nlsLM14, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM14)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params14 <- left_join(ci_extra_params14, extra_params14)

#- Topt 

topt14 <- as_tibble(ci_extra_params14[2,])
topt14$species <- as.character("Aphis gossypii")

ggplot(topt14, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())


ap.gT <- rbind(topt12,topt13,topt14)

ap.gT$param <- as.character(c("1/development time","juvenile mortality rate","fecundity rate"))

ap.gT$euler_param <- as.character(c("1/alpha","zj","bpk"))

ap.gT$stage <- as.factor(c("juvenile","juvenile","adult"))

ap.gT$stage <- relevel(ap.gT$stage, "juvenile")

mozzi <- rbind(mozz,ap.gT)

ggplot(mozzi, aes(estimate, species, shape = euler_param,colour=stage)) +
  geom_point(size = 2.5) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())+
  scale_color_brewer(type = 'qual', palette = 2)


#_________________________________________________________#
#                     Telenomus isis                      #
#_________________________________________________________#

te.i <- df %>% filter(species == "Telenomus isis", standardisedtraitname != "NA")

#- Development Rate

te.dv <- te.i %>% filter(standardisedtraitname == "Juvenile Development Rate")

te.dv <- te.dv %>% select(standardisedtraitvalue,temp)
te.dv <- te.dv %>% rename(rate = standardisedtraitvalue)

te.dv$rate <- te.dv$rate + 10^-6

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(te.dv$temp, te.dv$rate, model_name = 'pawar_2018')

o_fits <- nest(te.dv, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


o_preds <- mutate(o_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot 

ggplot(o_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, te.dv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM15 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = te.dv,
                                 start = coef(o_fits$pawar[[1]]),
                                 lower = get_lower_lims(te.dv$temp, te.dv$rate, model_name = 'pawar_2018'),
                                 upper = get_upper_lims(te.dv$temp, te.dv$rate, model_name = 'pawar_2018'),
                                 weights = rep(1, times = nrow(te.dv)))

# bootstrap using case resampling
boot15 <- Boot(fit_nlsLM15, method = 'residual')

boot15_preds <- boot15$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(te.dv$temp), max(te.dv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot15_conf_preds <- group_by(boot15_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), o_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot15_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), te.dv, size = 2, alpha = 0.5) +
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

#- Topt 

topt15 <- as_tibble(ci_extra_params15[2,])
topt15$species <- as.character("Telenomus isis")

ggplot(topt15, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())

#____________________________________________________

# Juvenile mortality rate

te.zj <- te.i %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

te.zj <- te.zj %>% select(standardisedtraitvalue,temp)
te.zj <- te.zj %>% rename(rate = standardisedtraitvalue)

te.zj$rate <- 1/te.zj$rate

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(te.zj$temp, te.zj$rate, model_name = 'pawar_2018')

i_fits <- nest(te.zj, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


i_preds <- mutate(i_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot 

ggplot(i_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, te.zj)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM16 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = te.zj,
                                 start = coef(i_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(te.zj)))

# bootstrap using case resampling
boot16 <- Boot(fit_nlsLM16, method = 'residual')

boot16_preds <- boot16$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(te.zj$temp), max(te.zj$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot16_conf_preds <- group_by(boot16_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), i_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot16_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), te.zj, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = '1/zj')

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

#- Topt 

topt16 <- as_tibble(ci_extra_params16[2,])
topt16$species <- as.character("Telenomus isis")

ggplot(topt16, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())


te.iT <- rbind(topt15,topt16)

te.iT$param <- as.character(c("1/development time","juvenile mortality rate"))

te.iT$euler_param <- as.character(c("1/alpha","zj"))

te.iT$stage <- as.factor(c("juvenile","juvenile"))

te.iT$stage <- relevel(te.iT$stage, "juvenile")

#_____________________________________

mozzie <- rbind(mozzi,te.iT)



#_____________________________________




#_________________________________________________________#
#                    Anthonomus grandis                   #
#_________________________________________________________#

ag.i <- df %>% filter(species == "Anthonomus grandis", standardisedtraitname != "NA")

#- Development Rate

ag.dv <- ag.i %>% filter(standardisedtraitname == "Juvenile Development Rate ")

ag.dv <- ag.dv %>% select(standardisedtraitvalue,temp)
ag.dv <- ag.dv %>% rename(rate = standardisedtraitvalue)

ag.dv$rate <- ag.dv$rate + 10^-6

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(ag.dv$temp, ag.dv$rate, model_name = 'pawar_2018')

pp_fits <- nest(ag.dv, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


pp_preds <- mutate(pp_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot 

ggplot(pp_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, ag.dv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM17 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = ag.dv,
                                 start = coef(pp_fits$pawar[[1]]),
                                 lower = get_lower_lims(ag.dv$temp, ag.dv$rate, model_name = 'pawar_2018'),
                                 upper = get_upper_lims(ag.dv$temp, ag.dv$rate, model_name = 'pawar_2018'),
                                 weights = rep(1, times = nrow(ag.dv)))

# bootstrap using case resampling
boot17 <- Boot(fit_nlsLM17, method = 'residual')

boot17_preds <- boot17$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(ag.dv$temp), max(ag.dv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot17_conf_preds <- group_by(boot17_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), pp_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot17_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), ag.dv, size = 2, alpha = 0.5) +
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

#- Topt 

topt17 <- as_tibble(ci_extra_params17[2,])
topt17$species <- as.character("Anthonomus grandis")

ggplot(topt17, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())

#____________________________________________________

# Juvenile mortality rate

ag.zj <- ag.i %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

ag.zj <- ag.zj %>% select(standardisedtraitvalue,temp)
ag.zj <- ag.zj %>% rename(rate = standardisedtraitvalue)

ag.zj$rate <- 1/ag.zj$rate

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(ag.zj$temp, ag.zj$rate, model_name = 'pawar_2018')

ppp_fits <- nest(ag.zj, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


ppp_preds <- mutate(ppp_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot 

ggplot(ppp_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, te.zj)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM18 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = ag.zj,
                                 start = coef(ppp_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(ag.zj)))

# bootstrap using case resampling
boot18 <- Boot(fit_nlsLM18, method = 'residual')

boot18_preds <- boot18$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(ag.zj$temp), max(ag.zj$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot18_conf_preds <- group_by(boot18_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), ppp_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot18_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), ag.zj, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = '1/zj')

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

#- Topt 

topt18 <- as_tibble(ci_extra_params18[2,])
topt18$species <- as.character("Anthonomus grandis")

ggplot(topt18, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())

# female mortality rate
####################################################################
#### Come back to this
# ag.z <- ag.i %>% filter(standardisedtraitname == "Female Mortality Rate")
# 
# ag.z <- ag.z %>% select(standardisedtraitvalue,temp)
# ag.z <- ag.z %>% rename(rate = standardisedtraitvalue)
# 
# ag.z$rate <- 1/ag.z$rate
# 
# # fit chosen model formulation in rTPC
# 
# start_vals <- get_start_vals(ag.z$temp, ag.z$rate, model_name = 'pawar_2018')
# 
# pppp_fits <- nest(ag.z, data = c(temp, rate)) %>%
#   mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
#                                           data = .x,
#                                           iter = c(3,3,3,3),
#                                           start_lower = start_vals - 10,
#                                           start_upper = start_vals + 10,
#                                           supp_errors = 'Y',
#                                           convergence_count = FALSE)))
# #____________________
# 
# 
# pppp_preds <- mutate(pppp_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
#   # get rid of original data column
#   select(., -data) %>%
#   # stack models into a single column, with an id column for model_name
#   pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
#   # create new list column containing the predictions
#   # this uses both fit and new_data list columns
#   mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
#   # select only the columns we want to keep
#   select(model_name, preds) %>%
#   # unlist the preds list column
#   unnest(preds)
# 
# 
# # plot 
# 
# ggplot(pppp_preds) +
#   geom_line(aes(temp, .fitted, col = model_name)) +
#   geom_point(aes(temp, rate),size=0.2,alpha=0.5, te.zj)+
#   theme_bw() +
#   theme(legend.position = 'none') +
#   scale_color_brewer(type = 'qual', palette = 2) +
#   labs(x = 'Temperature (?C)',
#        y = 'rate',
#        title = 'Dev rate thermal performance curves')+
#   theme_bw(base_size = 10)+
#   theme(legend.position = "none")
# 
# #____________________
# 
# fit_nlsLM18 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
#                                  data = ag.z,
#                                  start = coef(pppp_fits$pawar[[1]]),
#                                  weights = rep(1, times = nrow(ag.z)))
# 
# # bootstrap using case resampling
# boot18 <- Boot(fit_nlsLM18, method = 'residual')
# 
# boot18_preds <- boot18$t %>%
#   as.data.frame() %>%
#   drop_na() %>%
#   mutate(iter = 1:n()) %>%
#   group_by_all() %>%
#   do(data.frame(temp = seq(min(ag.z$temp), max(ag.z$temp), length.out = 100))) %>%
#   ungroup() %>%
#   mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))
# 
# # calculate bootstrapped confidence intervals
# boot18_conf_preds <- group_by(boot18_preds, temp) %>%
#   summarise(conf_lower = quantile(pred, 0.025),
#             conf_upper = quantile(pred, 0.975)) %>%
#   ungroup()
# 
# 
# # plot bootstrapped CIs
# 
# ggplot() +
#   geom_line(aes(temp, .fitted), pppp_preds, col = 'blue') +
#   geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot18_conf_preds, fill = 'blue', alpha = 0.3) +
#   geom_point(aes(temp, rate), ag.z, size = 2, alpha = 0.5) +
#   theme_bw(base_size = 12) +
#   labs(x = 'Temperature (?C)',
#        y = '1/zj')
# 
# # calculate params with CIs
# 
# extra_params18 <- calc_params(fit_nlsLM18) %>%
#   pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
# 
# ci_extra_params18 <- Boot(fit_nlsLM18, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM18)), R = 200, method = 'residual') %>%
#   confint(., method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'residual bootstrap')
# 
# ci_extra_params18 <- left_join(ci_extra_params18, extra_params18)
# 
# #- Topt 
# 
# topt18 <- as_tibble(ci_extra_params18[2,])
# topt18$species <- as.character("Anthonomus grandis")
# 
# ggplot(topt18, aes(estimate, species)) +
#   geom_point(size = 4) + 
#   geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
#   theme_bw(base_size = 12) +
#   scale_x_continuous('') +
#   theme(axis.title.y = element_blank())



ag.iT <- rbind(topt17,topt18)

ag.iT$param <- as.character(c("1/development time","juvenile mortality rate"))

ag.iT$euler_param <- as.character(c("1/alpha","zj"))

ag.iT$stage <- as.factor(c("juvenile","juvenile"))

ag.iT$stage <- relevel(ag.iT$stage, "juvenile")

mozzie <- rbind(mozzi,ag.iT)
#_____________________________________


#_________________________________________________________#
#                    Corythucha ciliata                   #
#_________________________________________________________#

acc.i <- df %>% filter(species == "Corythucha ciliata", standardisedtraitname != "NA")

#- Development Rate

acc.dv <- acc.i %>% filter(standardisedtraitname == "Juvenile Development Rate")

acc.dv <- acc.dv %>% select(standardisedtraitvalue,temp)
acc.dv <- acc.dv %>% rename(rate = standardisedtraitvalue)

acc.dv$rate <- acc.dv$rate + 10^-6

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(acc.dv$temp, acc.dv$rate, model_name = 'pawar_2018')

qq_fits <- nest(acc.dv, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


qq_preds <- mutate(qq_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot 

ggplot(qq_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, acc.dv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM19 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = acc.dv,
                                 start = coef(qq_fits$pawar[[1]]),
                                 lower = get_lower_lims(acc.dv$temp, acc.dv$rate, model_name = 'pawar_2018'),
                                 upper = get_upper_lims(acc.dv$temp, acc.dv$rate, model_name = 'pawar_2018'),
                                 weights = rep(1, times = nrow(acc.dv)))

# bootstrap using case resampling
boot19 <- Boot(fit_nlsLM19, method = 'residual')

boot19_preds <- boot19$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(acc.dv$temp), max(acc.dv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot19_conf_preds <- group_by(boot19_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), qq_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot19_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), acc.dv, size = 2, alpha = 0.5) +
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

#- Topt 

topt19 <- as_tibble(ci_extra_params19[2,])
topt19$species <- as.character("Corythucha ciliata")

ggplot(topt19, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())

#____________________________________________________

# Juvenile mortality rate

acc.zj <- acc.i %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

acc.zj <- acc.zj %>% select(standardisedtraitvalue,temp)
acc.zj <- acc.zj %>% rename(rate = standardisedtraitvalue)

acc.zj$rate <- 1/acc.zj$rate

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(acc.zj$temp, acc.zj$rate, model_name = 'pawar_2018')

qqq_fits <- nest(acc.zj, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


qqq_preds <- mutate(qqq_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot 

ggplot(qqq_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, te.zj)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM20 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = acc.zj,
                                 start = coef(qqq_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(acc.zj)))

# bootstrap using case resampling
boot20 <- Boot(fit_nlsLM20, method = 'residual')

boot20_preds <- boot20$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(acc.zj$temp), max(acc.zj$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot20_conf_preds <- group_by(boot20_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), qqq_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot20_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), acc.zj, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = '1/zj')

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

#- Topt 

topt20 <- as_tibble(ci_extra_params20[2,])
topt20$species <- as.character("Corythucha ciliata")

ggplot(topt20, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())

#______________________________________

# Fecundity rate

acc.bpk <- acc.i %>% filter(standardisedtraitname == "Fecundity Rate")

acc.bpk <- acc.bpk %>% select(standardisedtraitvalue,temp)
acc.bpk <- acc.bpk %>% rename(rate = standardisedtraitvalue)

acc.bpk$rate <- acc.bpk$rate + 10^-6


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(acc.bpk$temp, acc.bpk$rate, model_name = 'pawar_2018')

jq_fits <- nest(acc.bpk, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


jq_preds <- mutate(jq_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot panel

ggplot(jq_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, acc.bpk)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'Rate')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM21 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = acc.bpk,
                                 start = coef(jq_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(acc.bpk)))

# bootstrap using case resampling
boot21 <- Boot(fit_nlsLM21, method = 'residual')

boot21_preds <- boot21$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(acc.bpk$temp), max(acc.bpk$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot21_conf_preds <- group_by(boot21_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), jq_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot21_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), acc.bpk, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'bpk')

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

#- Topt 

topt21 <- as_tibble(ci_extra_params21[2,])
topt21$species <- as.character("Corythucha ciliata")

ggplot(topt21, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())

acc.accT <- rbind(topt19,topt20,topt21)

acc.accT$param <- as.character(c("1/development time","juvenile mortality rate","fecundity rate"))

acc.accT$euler_param <- as.character(c("1/alpha","zj","bpk"))

acc.accT$stage <- as.factor(c("juvenile","juvenile","adult"))


acc.accT$stage <- relevel(acc.accT$stage, "juvenile")


mozzies <- rbind(mozzie,acc.accT)

#_____________________________________________________________________________
# 
#_________________________________________________________#
#                    Corythucha ciliata                   #
#_________________________________________________________#

acc.i <- df %>% filter(species == "Corythucha ciliata", standardisedtraitname != "NA")

#- Development Rate

acc.dv <- acc.i %>% filter(standardisedtraitname == "Juvenile Development Rate")

acc.dv <- acc.dv %>% select(standardisedtraitvalue,temp)
acc.dv <- acc.dv %>% rename(rate = standardisedtraitvalue)

acc.dv$rate <- acc.dv$rate + 10^-6

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(acc.dv$temp, acc.dv$rate, model_name = 'pawar_2018')

qq_fits <- nest(acc.dv, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


qq_preds <- mutate(qq_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot 

ggplot(qq_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, acc.dv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM19 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = acc.dv,
                                 start = coef(qq_fits$pawar[[1]]),
                                 lower = get_lower_lims(acc.dv$temp, acc.dv$rate, model_name = 'pawar_2018'),
                                 upper = get_upper_lims(acc.dv$temp, acc.dv$rate, model_name = 'pawar_2018'),
                                 weights = rep(1, times = nrow(acc.dv)))

# bootstrap using case resampling
boot19 <- Boot(fit_nlsLM19, method = 'residual')

boot19_preds <- boot19$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(acc.dv$temp), max(acc.dv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot19_conf_preds <- group_by(boot19_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), qq_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot19_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), acc.dv, size = 2, alpha = 0.5) +
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

#- Topt 

topt19 <- as_tibble(ci_extra_params19[2,])
topt19$species <- as.character("Corythucha ciliata")

ggplot(topt19, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())

#____________________________________________________

# Juvenile mortality rate

acc.zj <- acc.i %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

acc.zj <- acc.zj %>% select(standardisedtraitvalue,temp)
acc.zj <- acc.zj %>% rename(rate = standardisedtraitvalue)

acc.zj$rate <- 1/acc.zj$rate

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(acc.zj$temp, acc.zj$rate, model_name = 'pawar_2018')

qqq_fits <- nest(acc.zj, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


qqq_preds <- mutate(qqq_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot 

ggplot(qqq_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, te.zj)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM20 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = acc.zj,
                                 start = coef(qqq_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(acc.zj)))

# bootstrap using case resampling
boot20 <- Boot(fit_nlsLM20, method = 'residual')

boot20_preds <- boot20$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(acc.zj$temp), max(acc.zj$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot20_conf_preds <- group_by(boot20_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), qqq_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot20_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), acc.zj, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = '1/zj')

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

#- Topt 

topt20 <- as_tibble(ci_extra_params20[2,])
topt20$species <- as.character("Corythucha ciliata")

ggplot(topt20, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())

#______________________________________

# Fecundity rate

acc.bpk <- acc.i %>% filter(standardisedtraitname == "Fecundity Rate")

acc.bpk <- acc.bpk %>% select(standardisedtraitvalue,temp)
acc.bpk <- acc.bpk %>% rename(rate = standardisedtraitvalue)

acc.bpk$rate <- acc.bpk$rate + 10^-6


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(acc.bpk$temp, acc.bpk$rate, model_name = 'pawar_2018')

jq_fits <- nest(acc.bpk, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


jq_preds <- mutate(jq_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot panel

ggplot(jq_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, acc.bpk)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'Rate')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM21 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = acc.bpk,
                                 start = coef(jq_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(acc.bpk)))

# bootstrap using case resampling
boot21 <- Boot(fit_nlsLM21, method = 'residual')

boot21_preds <- boot21$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(acc.bpk$temp), max(acc.bpk$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot21_conf_preds <- group_by(boot21_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), jq_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot21_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), acc.bpk, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'bpk')

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

#- Topt 

topt21 <- as_tibble(ci_extra_params21[2,])
topt21$species <- as.character("Corythucha ciliata")

ggplot(topt21, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())

acc.accT <- rbind(topt19,topt20,topt21)

acc.accT$param <- as.character(c("1/development time","juvenile mortality rate","fecundity rate"))

acc.accT$euler_param <- as.character(c("1/alpha","zj","bpk"))

acc.accT$stage <- as.factor(c("juvenile","juvenile","adult"))


acc.accT$stage <- relevel(acc.accT$stage, "juvenile")


mozziess <- rbind(mozzies,acc.accT)

#_____________________________________________________________________________

#_________________________________________________________#
#                    Cydia pomonella                   #
#_________________________________________________________#

acp.i <- df %>% filter(species == "Cydia pomonella", standardisedtraitname != "NA")

#- Development Rate

acp.dv <- acp.i %>% filter(standardisedtraitname == "Juvenile Development Rate")

acp.dv <- acp.dv %>% select(standardisedtraitvalue,temp)
acp.dv <- acp.dv %>% rename(rate = standardisedtraitvalue)

acp.dv$rate <- acp.dv$rate + 10^-6

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(acp.dv$temp, acp.dv$rate, model_name = 'pawar_2018')

acp_fits <- nest(acp.dv, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


acp_preds <- mutate(acp_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot 

ggplot(acp_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, acp.dv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM22 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = acp.dv,
                                 start = coef(acp_fits$pawar[[1]]),
                                 lower = get_lower_lims(acp.dv$temp, acp.dv$rate, model_name = 'pawar_2018'),
                                 upper = get_upper_lims(acp.dv$temp, acp.dv$rate, model_name = 'pawar_2018'),
                                 weights = rep(1, times = nrow(acp.dv)))

# bootstrap using case resampling
boot22 <- Boot(fit_nlsLM22, method = 'residual')

boot22_preds <- boot22$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(acp.dv$temp), max(acp.dv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot22_conf_preds <- group_by(boot22_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), acp_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot22_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), acp.dv, size = 2, alpha = 0.5) +
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

#- Topt 

topt22 <- as_tibble(ci_extra_params22[2,])
topt22$species <- as.character("Cydia pomonella")

ggplot(topt22, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())


acp.acpT <- topt22

acp.acpT$param <- as.character(c("1/development time"))

acp.acpT$euler_param <- as.character(c("1/alpha"))

acp.acpT$stage <- as.factor(c("juvenile"))


acp.acpT$stage <- relevel(acp.acpT$stage, "juvenile")


mozziesss <- rbind(mozziess,acp.acpT)

#_____________________________________________________________________________

#_________________________________________________________#
#                    Helicoverpa armigera                   #
#_________________________________________________________#

ha.i <- df %>% filter(species == "Helicoverpa armigera", standardisedtraitname != "NA")

#- Development Rate

ha.dv <- ha.i %>% filter(standardisedtraitname == "Juvenile Development Time")

ha.dv <- ha.dv %>% select(standardisedtraitvalue,temp)
ha.dv <- ha.dv %>% rename(rate = standardisedtraitvalue)

ha.dv = head(ha.dv,-1)
################################################################################
#convert back to dev rate?
ha.dv$rate <- ha.dv$rate + 10^-6
ha.dv$rate <- 1/ha.dv$rate
ha.dv$rate[is.na(ha.dv$rate)] <- 0
# fit chosen model formulation in rTPC

start_vals <- get_start_vals(ha.dv$temp, ha.dv$rate, model_name = 'pawar_2018')

ha_fits <- nest(ha.dv, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


ha_preds <- mutate(ha_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot 

ggplot(ha_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, ha.dv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM23 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = ha.dv,
                                 start = coef(ha_fits$pawar[[1]]),
                                 lower = get_lower_lims(ha.dv$temp, ha.dv$rate, model_name = 'pawar_2018'),
                                 upper = get_upper_lims(ha.dv$temp, ha.dv$rate, model_name = 'pawar_2018'),
                                 weights = rep(1, times = nrow(ha.dv)))

# bootstrap using case resampling
boot23 <- Boot(fit_nlsLM23, method = 'residual')

boot23_preds <- boot23$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(ha.dv$temp), max(ha.dv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot23_conf_preds <- group_by(boot23_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), ha_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot23_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), ha.dv, size = 2, alpha = 0.5) +
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

#- Topt 

topt23 <- as_tibble(ci_extra_params23[2,])
topt23$species <- as.character("Helicoverpa armigera")

ggplot(topt23, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())

#____________________________________________________

# Juvenile mortality rate

ha.zj <- ha.i %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

ha.zj <- ha.zj %>% select(standardisedtraitvalue,temp)
ha.zj <- ha.zj %>% rename(rate = standardisedtraitvalue)

ha.zj$rate <- 1/ha.zj$rate

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(ha.zj$temp, ha.zj$rate, model_name = 'pawar_2018')

haq_fits <- nest(ha.zj, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


haq_preds <- mutate(haq_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot 

ggplot(haq_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, ha.zj)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM24 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = ha.zj,
                                 start = coef(haq_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(ha.zj)))

# bootstrap using case resampling
boot24 <- Boot(fit_nlsLM24, method = 'residual')

boot24_preds <- boot24$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(ha.zj$temp), max(ha.zj$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot24_conf_preds <- group_by(boot24_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), haq_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot24_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), ha.zj, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = '1/zj')

# calculate params with CIs

extra_params24 <- calc_params(fit_nlsLM24) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params24 <- Boot(fit_nlsLM24, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM20)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params24 <- left_join(ci_extra_params24, extra_params24)

#- Topt 

topt24 <- as_tibble(ci_extra_params24[2,])
topt24$species <- as.character("Helicoverpa armigera")

ggplot(topt24, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())

#______________________________________

# Fecundity rate

ha.bpk <- ha.i %>% filter(standardisedtraitname == "Fecundity Rate")

ha.bpk <- ha.bpk %>% select(standardisedtraitvalue,temp)
ha.bpk <- ha.bpk %>% rename(rate = standardisedtraitvalue)

ha.bpk$rate <- ha.bpk$rate + 10^-6


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(ha.bpk$temp, ha.bpk$rate, model_name = 'pawar_2018')

haf_fits <- nest(ha.bpk, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


haf_preds <- mutate(haf_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot panel

ggplot(haf_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, ha.bpk)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'Rate')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM24 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = ha.bpk,
                                 start = coef(haf_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(ha.bpk)))

# bootstrap using case resampling
boot24 <- Boot(fit_nlsLM24, method = 'residual')

boot24_preds <- boot24$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(ha.bpk$temp), max(ha.bpk$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot24_conf_preds <- group_by(boot24_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), haf_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot24_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), ha.bpk, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'bpk')

# calculate params with CIs

extra_params24 <- calc_params(fit_nlsLM24) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params24 <- Boot(fit_nlsLM24, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM24)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params24 <- left_join(ci_extra_params24, extra_params24)

#- Topt 

topt24 <- as_tibble(ci_extra_params24[2,])
topt24$species <- as.character("Helicoverpa armigera")

ggplot(topt24, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())

ha.accT <- rbind(topt22,topt23,topt24)

ha.accT$param <- as.character(c("1/development time","juvenile mortality rate","fecundity rate"))

ha.accT$euler_param <- as.character(c("1/alpha","zj","bpk"))

ha.accT$stage <- as.factor(c("juvenile","juvenile","adult"))


ha.accT$stage <- relevel(ha.accT$stage, "juvenile")


mozziessss <- rbind(mozziesss,ha.accT)

#_____________________________________________________________________________

#_________________________________________________________#
#                    Lepinotus reticulatus               #
#_________________________________________________________#

lr.i <- df %>% filter(species == "Lepinotus reticulatus", standardisedtraitname != "NA")

#- Development Rate

lr.dv <- lr.i %>% filter(standardisedtraitname == "Juvenile Development Rate")

lr.dv <- lr.dv %>% select(standardisedtraitvalue,temp)
lr.dv <- lr.dv %>% rename(rate = standardisedtraitvalue)


lr.dv$rate <- lr.dv$rate + 10^-6
# fit chosen model formulation in rTPC

start_vals <- get_start_vals(lr.dv$temp, lr.dv$rate, model_name = 'pawar_2018')

lr_fits <- nest(lr.dv, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


lr_preds <- mutate(lr_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot 

ggplot(lr_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, lr.dv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM25 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = lr.dv,
                                 start = coef(lr_fits$pawar[[1]]),
                                 lower = get_lower_lims(lr.dv$temp, lr.dv$rate, model_name = 'pawar_2018'),
                                 upper = get_upper_lims(lr.dv$temp, lr.dv$rate, model_name = 'pawar_2018'),
                                 weights = rep(1, times = nrow(lr.dv)))

# bootstrap using case resampling
boot25 <- Boot(fit_nlsLM25, method = 'residual')

boot25_preds <- boot25$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(lr.dv$temp), max(lr.dv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot25_conf_preds <- group_by(boot25_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), lr_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot25_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), lr.dv, size = 2, alpha = 0.5) +
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

#- Topt 

topt25 <- as_tibble(ci_extra_params25[2,])
topt25$species <- as.character("Lepinotus reticulatus")

ggplot(topt25, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())

#____________________________________________________

# Juvenile mortality rate
###############################################################################
#Didn't work

# lr.zj <- lr.i %>% filter(standardisedtraitname == "Juvenile Mortality Rate")
# 
# lr.zj <- lr.zj %>% select(standardisedtraitvalue,temp)
# lr.zj <- lr.zj %>% rename(rate = standardisedtraitvalue)
# 
# lr.zj$rate <- 1/lr.zj$rate
# 
# # fit chosen model formulation in rTPC
# 
# start_vals <- get_start_vals(lr.zj$temp, lr.zj$rate, model_name = 'pawar_2018')
# 
# lrq_fits <- nest(lr.zj, data = c(temp, rate)) %>%
#   mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
#                                           data = .x,
#                                           iter = c(3,3,3,3),
#                                           start_lower = start_vals - 10,
#                                           start_upper = start_vals + 10,
#                                           supp_errors = 'Y',
#                                           convergence_count = FALSE)))
# #____________________
# 
# 
# lrq_preds <- mutate(lrq_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
#   # get rid of original data column
#   select(., -data) %>%
#   # stack models into a single column, with an id column for model_name
#   pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
#   # create new list column containing the predictions
#   # this uses both fit and new_data list columns
#   mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
#   # select only the columns we want to keep
#   select(model_name, preds) %>%
#   # unlist the preds list column
#   unnest(preds)
# 
# 
# # plot 
# 
# ggplot(lrq_preds) +
#   geom_line(aes(temp, .fitted, col = model_name)) +
#   geom_point(aes(temp, rate),size=0.2,alpha=0.5, lr.zj)+
#   theme_bw() +
#   theme(legend.position = 'none') +
#   scale_color_brewer(type = 'qual', palette = 2) +
#   labs(x = 'Temperature (?C)',
#        y = 'rate',
#        title = 'Dev rate thermal performance curves')+
#   theme_bw(base_size = 10)+
#   theme(legend.position = "none")
# 
# #____________________
# 
# fit_nlsLM24 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
#                                  data = lr.zj,
#                                  start = coef(lrq_fits$pawar[[1]]),
#                                  weights = rep(1, times = nrow(lr.zj)))
# 
# # bootstrap using case resampling
# boot24 <- Boot(fit_nlsLM24, method = 'residual')
# 
# boot24_preds <- boot24$t %>%
#   as.data.frame() %>%
#   drop_na() %>%
#   mutate(iter = 1:n()) %>%
#   group_by_all() %>%
#   do(data.frame(temp = seq(min(lr.zj$temp), max(lr.zj$temp), length.out = 100))) %>%
#   ungroup() %>%
#   mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))
# 
# # calculate bootstrapped confidence intervals
# boot24_conf_preds <- group_by(boot24_preds, temp) %>%
#   summarise(conf_lower = quantile(pred, 0.025),
#             conf_upper = quantile(pred, 0.975)) %>%
#   ungroup()
# 
# 
# # plot bootstrapped CIs
# 
# ggplot() +
#   geom_line(aes(temp, .fitted), haq_preds, col = 'blue') +
#   geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot24_conf_preds, fill = 'blue', alpha = 0.3) +
#   geom_point(aes(temp, rate), lr.zj, size = 2, alpha = 0.5) +
#   theme_bw(base_size = 12) +
#   labs(x = 'Temperature (?C)',
#        y = '1/zj')
# 
# # calculate params with CIs
# 
# extra_params24 <- calc_params(fit_nlsLM24) %>%
#   pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
# 
# ci_extra_params24 <- Boot(fit_nlsLM24, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM20)), R = 200, method = 'residual') %>%
#   confint(., method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'residual bootstrap')
# 
# ci_extra_params24 <- left_join(ci_extra_params24, extra_params24)
# 
# #- Topt 
# 
# topt24 <- as_tibble(ci_extra_params24[2,])
# topt24$species <- as.character("Helicoverpa armigera")
# 
# ggplot(topt24, aes(estimate, species)) +
#   geom_point(size = 4) + 
#   geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
#   theme_bw(base_size = 12) +
#   scale_x_continuous('') +
#   theme(axis.title.y = element_blank())


jc.accT <- topt25

jc.accT$param <- as.character(c("1/development time"))

jc.accT$euler_param <- as.character(c("1/alpha"))

jc.accT$stage <- as.factor(c("juvenile"))


jc.accT$stage <- relevel(jc.accT$stage, "juvenile")


mozziesssss <- rbind(mozziessss,jc.accT)

#______________________________________
#_________________________________________________________#
#                  Macrocentrus iridescens               #
#_________________________________________________________#

mi.i <- df %>% filter(species == "Macrocentrus iridescens", standardisedtraitname != "NA")

mi.zj <- mi.i %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

mi.zj <- mi.zj %>% select(standardisedtraitvalue,temp)
mi.zj <- mi.zj %>% rename(rate = standardisedtraitvalue)

mi.zj$rate <- 1/mi.zj$rate

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(mi.zj$temp, mi.zj$rate, model_name = 'pawar_2018')

miq_fits <- nest(mi.zj, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


miq_preds <- mutate(miq_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot

ggplot(miq_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, mi.zj)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM26 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = mi.zj,
                                 start = coef(miq_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(mi.zj)))

# bootstrap using case resampling
boot26 <- Boot(fit_nlsLM26, method = 'residual')

boot26_preds <- boot26$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(mi.zj$temp), max(mi.zj$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot26_conf_preds <- group_by(boot26_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), miq_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot26_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), mi.zj, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = '1/zj')

# calculate params with CIs

extra_params26 <- calc_params(fit_nlsLM26) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params26 <- Boot(fit_nlsLM26, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM20)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params26 <- left_join(ci_extra_params26, extra_params26)

#- Topt

topt26 <- as_tibble(ci_extra_params26[2,])
topt26$species <- as.character("Macrocentrus iridescens")

ggplot(topt26, aes(estimate, species)) +
  geom_point(size = 4) +
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())


mi.accT <- topt26

mi.accT$param <- as.character(c("juvenile mortality rate"))

mi.accT$euler_param <- as.character(c("zj"))

mi.accT$stage <- as.factor(c("juvenile"))


mi.accT$stage <- relevel(mi.accT$stage, "juvenile")


mozziessssss <- rbind(mozziesssss,mi.accT)
#_____________________________________________________________________________


#_________________________________________________________#
#                    Planococcus citri                    #
#_________________________________________________________#

pc.i <- df %>% filter(species == "Planococcus citri", standardisedtraitname != "NA")

#- Development Rate

pc.dv <- pc.i %>% filter(standardisedtraitname == "Juvenile Development Rate")

pc.dv <- pc.dv %>% select(standardisedtraitvalue,temp)
pc.dv <- pc.dv %>% rename(rate = standardisedtraitvalue)


pc.dv$rate <- pc.dv$rate + 10^-6
# fit chosen model formulation in rTPC

start_vals <- get_start_vals(pc.dv$temp, pc.dv$rate, model_name = 'pawar_2018')

pc_fits <- nest(pc.dv, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


pc_preds <- mutate(pc_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot 

ggplot(pc_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, pc.dv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM27 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = pc.dv,
                                 start = coef(pc_fits$pawar[[1]]),
                                 # lower = get_lower_lims(pc.dv$temp, pc.dv$rate, model_name = 'pawar_2018'),
                                 # upper = get_upper_lims(pc.dv$temp, pc.dv$rate, model_name = 'pawar_2018'),
                                 weights = rep(1, times = nrow(pc.dv)))

# bootstrap using case resampling
boot27 <- Boot(fit_nlsLM27, method = 'residual')

boot27_preds <- boot27$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(pc.dv$temp), max(pc.dv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot27_conf_preds <- group_by(boot27_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.027),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), pc_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot27_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), pc.dv, size = 2, alpha = 0.5) +
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

#- Topt 

topt27 <- as_tibble(ci_extra_params27[2,])
topt27$species <- as.character("Planococcus citri")

ggplot(topt27, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())

#____________________________________________________

###############################################################################
#didn't work
# # female mortality rate
# 
# 
# pc.z <- pc.i %>% filter(standardisedtraitname == "Female Mortality Rate")
# 
# pc.z <- pc.z %>% select(standardisedtraitvalue,temp)
# pc.z <- pc.z %>% rename(rate = standardisedtraitvalue)
# 
# pc.z$rate <- 1/pc.z$rate
# 
# # fit chosen model formulation in rTPC
# 
# start_vals <- get_start_vals(pc.z$temp, pc.z$rate, model_name = 'pawar_2018')
# 
# pcq_fits <- nest(pc.z, data = c(temp, rate)) %>%
#   mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
#                                           data = .x,
#                                           iter = c(3,3,3,3),
#                                           start_lower = start_vals - 10,
#                                           start_upper = start_vals + 10,
#                                           supp_errors = 'Y',
#                                           convergence_count = FALSE)))
# #____________________
# 
# 
# pcq_preds <- mutate(pcq_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>%
#   # get rid of original data column
#   select(., -data) %>%
#   # stack models into a single column, with an id column for model_name
#   pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
#   # create new list column containing the predictions
#   # this uses both fit and new_data list columns
#   mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
#   # select only the columns we want to keep
#   select(model_name, preds) %>%
#   # unlist the preds list column
#   unnest(preds)
# 
# 
# # plot
# 
# ggplot(pcq_preds) +
#   geom_line(aes(temp, .fitted, col = model_name)) +
#   geom_point(aes(temp, rate),size=0.2,alpha=0.5, pc.z)+
#   theme_bw() +
#   theme(legend.position = 'none') +
#   scale_color_brewer(type = 'qual', palette = 2) +
#   labs(x = 'Temperature (?C)',
#        y = 'rate',
#        title = 'Dev rate thermal performance curves')+
#   theme_bw(base_size = 10)+
#   theme(legend.position = "none")
# 
# #____________________
# 
# fit_nlsLM24 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
#                                  data = pc.z,
#                                  start = coef(pcq_fits$pawar[[1]]),
#                                  weights = rep(1, times = nrow(pc.z)))
# 
# # bootstrap using case resampling
# boot24 <- Boot(fit_nlsLM24, method = 'residual')
# 
# boot24_preds <- boot24$t %>%
#   as.data.frame() %>%
#   drop_na() %>%
#   mutate(iter = 1:n()) %>%
#   group_by_all() %>%
#   do(data.frame(temp = seq(min(pc.z$temp), max(pc.z$temp), length.out = 100))) %>%
#   ungroup() %>%
#   mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))
# 
# # calculate bootstrapped confidence intervals
# boot24_conf_preds <- group_by(boot24_preds, temp) %>%
#   summarise(conf_lower = quantile(pred, 0.025),
#             conf_upper = quantile(pred, 0.975)) %>%
#   ungroup()
# 
# 
# # plot bootstrapped CIs
# 
# ggplot() +
#   geom_line(aes(temp, .fitted), haq_preds, col = 'blue') +
#   geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot24_conf_preds, fill = 'blue', alpha = 0.3) +
#   geom_point(aes(temp, rate), pc.z, size = 2, alpha = 0.5) +
#   theme_bw(base_size = 12) +
#   labs(x = 'Temperature (?C)',
#        y = '1/z')
# 
# # calculate params with CIs
# 
# extra_params24 <- calc_params(fit_nlsLM24) %>%
#   pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
# 
# ci_extra_params24 <- Boot(fit_nlsLM24, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM20)), R = 200, method = 'residual') %>%
#   confint(., method = 'bca') %>%
#   as.data.frame() %>%
#   rename(conf_lower = 1, conf_upper = 2) %>%
#   rownames_to_column(., var = 'param') %>%
#   mutate(method = 'residual bootstrap')
# 
# ci_extra_params24 <- left_join(ci_extra_params24, extra_params24)
# 
# #- Topt
# 
# topt24 <- as_tibble(ci_extra_params24[2,])
# topt24$species <- as.character("Planococcus citri")
# 
# ggplot(topt24, aes(estimate, species)) +
#   geom_point(size = 4) +
#   geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
#   theme_bw(base_size = 12) +
#   scale_x_continuous('') +
#   theme(axis.title.y = element_blank())
#______________________________________

# Fecundity rate

pc.bpk <- pc.i %>% filter(standardisedtraitname == "Fecundity Rate")

pc.bpk <- pc.bpk %>% select(originaltraitvalue,temp)
pc.bpk <- pc.bpk %>% rename(rate = originaltraitvalue)

pc.bpk$rate <- pc.bpk$rate + 10^-6


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(pc.bpk$temp, pc.bpk$rate, model_name = 'pawar_2018')

pc_fits <- nest(pc.bpk, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


pc_preds <- mutate(pc_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)


# plot panel

ggplot(pc_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, pc.bpk)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'Rate')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM28 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = pc.bpk,
                                start = coef(pc_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(pc.bpk)))

# bootstrap using case resampling
boot28 <- Boot(fit_nlsLM28, method = 'residual')

boot28_preds <- boot28$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(pc.bpk$temp), max(pc.bpk$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot28_conf_preds <- group_by(boot28_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), pc_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot28_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), pc.bpk, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'bpk')

# calculate params with CIs

extra_params28 <- calc_params(fit_nlsLM28) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params28 <- Boot(fit_nlsLM28, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM9)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params28 <- left_join(ci_extra_params28, extra_params28)

#- Topt 

topt28 <- as_tibble(ci_extra_params28[2,])
topt28$species <- as.character("Planococcus citri")

ggplot(topt28, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())


pc.accT <- rbind(topt27, topt28)

pc.accT$param <- as.character(c("1/development time","fecundity rate"))

pc.accT$euler_param <- as.character(c("1/alpha","bpk"))

pc.accT$stage <- as.factor(c("juvenile","adult"))


pc.accT$stage <- relevel(pc.accT$stage, "juvenile")


mozziesssssss <- rbind(mozziessssss,pc.accT)

#______________________________________
#_________________________________________________________#
#                   Macrocentrus iridescens               #
#_________________________________________________________#

mi.i

#- Development Rate

mi.dv <- mi.i %>% filter(standardisedtraitname == "Juvenile Development Rate")

mi.dv <- mi.dv %>% select(originaltraitvalue,temp)
mi.dv <- mi.dv %>% rename(rate = originaltraitvalue)


mi.dv$rate <- mi.dv$rate + 10^-6
# fit chosen model formulation in rTPC

start_vals <- get_start_vals(mi.dv$temp, mi.dv$rate, model_name = 'pawar_2018')

mi_fits <- nest(mi.dv, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


mi_preds <- mutate(mi_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot 

ggplot(mi_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, mi.dv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  labs(x = 'Temperature (?C)',
       y = 'rate',
       title = 'Dev rate thermal performance curves')+
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM29 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = mi.dv,
                                 start = coef(mi_fits$pawar[[1]]),
                                 lower = get_lower_lims(mi.dv$temp, mi.dv$rate, model_name = 'pawar_2018'),
                                 upper = get_upper_lims(mi.dv$temp, mi.dv$rate, model_name = 'pawar_2018'),
                                 weights = rep(1, times = nrow(mi.dv)))

# bootstrap using case resampling
boot29 <- Boot(fit_nlsLM29, method = 'residual')

boot29_preds <- boot29$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(mi.dv$temp), max(mi.dv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot29_conf_preds <- group_by(boot29_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.029),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), mi_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot29_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), mi.dv, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (?C)',
       y = 'dev rate')

# calculate params with CIs

extra_params29 <- calc_params(fit_nlsLM29) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params29 <- Boot(fit_nlsLM29, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM29)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params29 <- left_join(ci_extra_params29, extra_params29)

#- Topt 

topt29 <- as_tibble(ci_extra_params29[2,])
topt29$species <- as.character("Macrocentrus iridescens")

ggplot(topt29, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())


mi2.iT <- topt29

mi2.iT$param <- as.character(c("1/development time"))

mi2.iT$euler_param <- as.character(c("1/alpha"))

mi2.iT$stage <- as.factor(c("juvenile"))

mi2.iT$stage <- relevel(mi2.iT$stage, "juvenile")

mnfinalmoz <- rbind(mozziesssssss,mi2.iT)
#____________________________________________________

mnuniquefinalmoz <- unique(mnfinalmoz)

mnuniquefinalmoz$species <- as.factor(mnuniquefinalmoz$species)
  
ggplot(mnuniquefinalmoz, aes(estimate, species, shape = euler_param,colour=stage)) +
  geom_point(size = 3.5) + 
  geom_errorbar(aes(xmin = conf_lower, xmax = conf_upper),lwd=0.5,width=0.2) +
  theme_bw(base_size = 16) +
  scale_x_continuous('') +
  theme(axis.title.y = element_blank())+
  scale_color_brewer(type = 'qual', palette = 6)

#write.csv(mozziesssssss, '../data/eulerparammn.csv')
