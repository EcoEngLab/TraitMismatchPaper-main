# load packages

require('tidyverse')
require('nls.multstart')
require('broom')
require('rTPC')
require('data.table')
require('car')
require('patchwork')
require('minpack.lm')
require('boot')

rm(list=ls())
graphics.off()

setwd("~/Dropbox/ph_thesis/Topt_paper/data")

#read in the trait data

df <- as_tibble(read.csv('mismatch_dat.csv'))

df<- df %>%
  rename(temp = ambienttemp,species = interactor1)

df$species <- as.factor(df$species)

#_____________________________________________________#
#         Tetraneura nigriabdominalis                 #
#_____________________________________________________#

        
te.ni <- df %>% filter(species == "Tetraneura nigriabdominalis")

#- Development Rate

te.nidv <- te.ni %>% filter(standardisedtraitname == "Juvenile Development Rate")

te.nidv <- te.nidv %>% select(standardisedtraitvalue,temp)
te.nidv <- te.nidv %>% rename(rate = standardisedtraitvalue)

te.nidv$rate <- te.nidv$rate + 10^-6

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(te.nidv$temp, te.nidv$rate, model_name = 'pawar_2018')

d_fits <- nest(te.nidv, data = c(temp, rate)) %>%
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
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, te.nidv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  theme(legend.position = "none")

#____________________

fit_nlsLM1 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = te.nidv,
                                start = coef(d_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(te.nidv)))

# bootstrap using residual resampling
boot1 <- Boot(fit_nlsLM1, method = 'residual')

boot1_preds <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(te.nidv$temp), max(te.nidv$temp), length.out = 100))) %>%
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
  geom_point(aes(temp, rate), te.nidv, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12)

# calculate params with CIs

extra_params1 <- calc_params(fit_nlsLM1) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params1 <- Boot(fit_nlsLM1, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM1)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')


# params
ci_extra_params1 <- left_join(ci_extra_params1, extra_params1)

# dataset for preds with 95% CIs
teni <- cbind(d_preds,boot1_conf_preds) 
teni <- teni[-1:-2]
teni$.fitted <- 1/teni$.fitted
teni$conf_lower <- 1/teni$conf_lower
teni$conf_upper <- 1/teni$conf_upper

teni <- teni %>% rename(alpha = .fitted, alpha_lwr = conf_lower, alpha_upr = conf_upper)

teni$species <- 'Tetraneura nigriabdominalis'

teni <- teni %>% select(temp, alpha, alpha_lwr,alpha_upr, species)

#____________________________________________________

# Juvenile mortality rate

te.nizj <- te.ni %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

te.nizj <- te.nizj %>% select(standardisedtraitvalue,temp)
te.nizj <- te.nizj %>% rename(rate = standardisedtraitvalue)

te.nizj$rate <- 1/te.nizj$rate

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(te.nizj$temp, te.nizj$rate, model_name = 'pawar_2018')

a_fits <- nest(te.nizj, data = c(temp, rate)) %>%
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
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, te.nizj)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) 

#____________________

fit_nlsLM2 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = te.nizj,
                                start = coef(a_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(te.nizj)))

# bootstrap using residual resampling
boot2 <- Boot(fit_nlsLM2, method = 'residual')

boot2_preds <- boot2$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(te.nizj$temp), max(te.nizj$temp), length.out = 100))) %>%
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
  geom_point(aes(temp, rate), te.nizj, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/zj')

# calculate params with CIs

extra_params2 <- calc_params(fit_nlsLM2) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params2 <- Boot(fit_nlsLM2, f = function(x){unlist(calc_params(x))}, 
                    labels = names(calc_params(fit_nlsLM2)), R = 200, method = 'residual') %>%
                    confint(., method = 'bca') %>%
                    as.data.frame() %>%
                    rename(conf_lower = 1, conf_upper = 2) %>%
                    rownames_to_column(., var = 'param') %>%
                    mutate(method = 'residual bootstrap')

# params
ci_extra_params2 <- left_join(ci_extra_params2, extra_params2)

# dataset for preds with 95% CIs

tenzj <- cbind(a_preds,boot2_conf_preds) 
tenzj <- tenzj[-1:-2]
tenzj$.fitted <- 1/tenzj$.fitted
tenzj$conf_lower <- 1/tenzj$conf_lower
tenzj$conf_upper <- 1/tenzj$conf_upper
tenzj$species <- 'Tetraneura nigriabdominalis'

tenzj <- tenzj %>% rename(zj = .fitted, zj_lwr = conf_lower, zj_upr = conf_upper)
tenzj <- tenzj %>% select(zj, zj_lwr,zj_upr)
tenzj <- tenzj %>% mutate(zj = replace(zj, zj > 1, 1))
tenzj <- tenzj %>% mutate(zj_lwr = replace(zj_lwr, zj_lwr > 1, 1))

teni <- cbind(teni,tenzj)

#______________________________________

# Female mortality rate

te.z <- te.ni %>% filter(standardisedtraitname == "Female Mortality Rate")

te.z <- te.z %>% select(standardisedtraitvalue,temp)
te.z <- te.z %>% rename(rate = standardisedtraitvalue)

te.z$rate <- 1/te.z$rate

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(te.z$temp, te.z$rate, model_name = 'pawar_2018')

b_fits <- nest(te.z, data = c(temp, rate)) %>%
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
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, te.z)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM3 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = te.z,
                                start = coef(b_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(te.z)))

# bootstrap using residual resampling
boot3 <- Boot(fit_nlsLM3, method = 'residual')

boot3_preds <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(te.z$temp), max(te.z$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot3_conf_preds <- group_by(boot3_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

  ggplot() +
  geom_line(aes(temp, .fitted), b_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot3_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), te.z, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/z')

# calculate params with CIs
extra_params3 <- calc_params(fit_nlsLM3) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params3 <- Boot(fit_nlsLM3, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM3)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params3 <- left_join(ci_extra_params3, extra_params3)

# dataset for preds with 95% CIs
tenz <- cbind(b_preds,boot3_conf_preds) 
tenz <- tenz[-1:-2]
tenz$.fitted <- 1/tenz$.fitted
tenz$conf_lower <- 1/tenz$conf_lower
tenz$conf_upper <- 1/tenz$conf_upper

tenz <- tenz %>% rename(z = .fitted, z_lwr = conf_lower, z_upr = conf_upper)

tenz$species <- 'Tetraneura nigriabdominalis'
tenz <- tenz %>% select(z, z_lwr,z_upr)
teni <- cbind(teni,tenz)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Fecundity rate
te.f <- te.ni %>% filter(standardisedtraitname == "Fecundity Rate")

te.f <- te.f %>% select(standardisedtraitvalue,temp)
te.f <- te.f %>% rename(rate = standardisedtraitvalue)

# fit chosen model formulation in rTPC
start_vals <- get_start_vals(te.f$temp, te.f$rate, model_name = 'pawar_2018')

f_fits <- nest(te.f, data = c(temp, rate)) %>%
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


# plot panel

ggplot(f_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, te.f)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM4 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = te.f,
                                start = coef(f_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(te.f)))

# bootstrap using residual resampling
boot4 <- Boot(fit_nlsLM4, method = 'residual')

boot4_preds <- boot4$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(te.f$temp), max(te.f$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot4_conf_preds <- group_by(boot4_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs
ggplot() +
  geom_line(aes(temp, .fitted), f_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot4_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), te.f, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = 'f')

# calculate params with CIs
extra_params4 <- calc_params(fit_nlsLM4) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params4 <- Boot(fit_nlsLM4, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM4)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params4 <- left_join(ci_extra_params4, extra_params4)

# dataset for preds with 95% CIs

tenf <- cbind(f_preds,boot4_conf_preds) 
tenf <- tenf[-1:-2]
tenf$trait <- 'bpk'
tenf$species <- 'Tetraneura nigriabdominalis'
tenf <- tenf %>% rename(bpk = .fitted, bpk_lwr = conf_lower, bpk_upr = conf_upper)
tenf <- tenf %>% select(bpk, bpk_lwr,bpk_upr)

tenaz <- cbind(teni,tenf)
tenaz$kappa <- 10^-6

#---------------------------------- Define parameters ----------------------------------#

zj    <- tenaz$zj
alpha <- tenaz$alpha
z     <- tenaz$z
bpk   <- tenaz$bpk
k     <- tenaz$kappa


# Calculate rmax

tenaz$rmax <- (((k+z)*((log(bpk/(k+z)))-(alpha*zj)))/(alpha*(k+z)+1))

# lower

zj_lwr    <- tenaz$zj_lwr
alpha_lwr <- tenaz$alpha_lwr
z_lwr     <- tenaz$z_lwr
bpk_lwr   <- tenaz$bpk_lwr
k         <- tenaz$kappa

tenaz$rmax_lwr <- (((k+z_lwr)*((log(bpk_lwr/(k+z_lwr)))-(alpha_lwr*zj_lwr)))/(alpha_lwr*(k+z_lwr)+1))

# upper 

zj_upr    <- tenaz$zj_upr
alpha_upr <- tenaz$alpha_upr
z_upr     <- tenaz$z_upr
bpk_upr   <- tenaz$bpk_upr
k         <- tenaz$kappa

tenaz$rmax_upr <- (((k+z_upr)*((log(bpk_upr/(k+z_upr)))-(alpha_upr*zj_upr)))/(alpha_upr*(k+z_upr)+1))

tenaz <- tenaz %>%                              
  mutate(rmax = replace(rmax, rmax < -0.1, -0.1), 
         rmax_lwr = replace(rmax_lwr, rmax_lwr < -0.1, -0.1),
         rmax_upr = replace(rmax_upr, rmax_upr < -0.1, -0.1))

ter_max <- select(tenaz, temp, species, rmax)
ter_max <- ter_max %>% rename(estimate=rmax)
ter_max <- filter(ter_max,temp > 15)
ter_max <- filter(ter_max,temp < 34)
                       
p1 <- ggplot(tenaz, aes(temp, rmax,fill=species,colour=species,shape=species))+
  ggtitle("a")+
  geom_ribbon(aes(temp, ymin = rmax_lwr, ymax = rmax_upr), alpha=0.4,col = NA)+
  theme_bw(base_size = 12)+
  geom_line(aes(temp,estimate),ter_max)



#_____________________________________________________#
#                Anthonomus grandis                   #
#_____________________________________________________#


ag <- df %>% filter(species == "Anthonomus grandis")

#- Development Rate

agdv <- ag %>% filter(standardisedtraitname == "Juvenile Development Rate")

agdv <- agdv %>% select(standardisedtraitvalue,temp)
agdv <- agdv %>% rename(rate = standardisedtraitvalue)

agdv$rate <- agdv$rate + 10^-6

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(agdv$temp, agdv$rate, model_name = 'pawar_2018')

v_fits <- nest(agdv, data = c(temp, rate)) %>%
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

# plot panel

ggplot(v_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, agdv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  theme(legend.position = "none")

#____________________

fit_nlsLM5 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = agdv,
                                start = coef(v_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(agdv)))

# bootstrap using residual resampling
boot5 <- Boot(fit_nlsLM5, method = 'residual')

boot5_preds <- boot5$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(agdv$temp), max(agdv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot5_conf_preds <- group_by(boot5_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), v_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot5_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), agdv, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12)

# calculate params with CIs

extra_params5 <- calc_params(fit_nlsLM5) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params5 <- Boot(fit_nlsLM5, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM5)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')


# params
ci_extra_params5 <- left_join(ci_extra_params5, extra_params5)

# dataset for preds with 95% CIs
agi <- cbind(v_preds,boot5_conf_preds) 
agi <- agi[-1:-2]
agi$.fitted <- 1/agi$.fitted
agi$conf_lower <- 1/agi$conf_lower
agi$conf_upper <- 1/agi$conf_upper
agi <- agi %>% rename(alpha = .fitted, alpha_lwr = conf_lower, alpha_upr = conf_upper)
agi$species <- 'Anthonomus grandis'
agi <- agi %>% select(temp, alpha, alpha_lwr,alpha_upr, species)

#____________________________________________________

# Juvenile mortality rate

agzj <- ag %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

agzj <- agzj %>% select(standardisedtraitvalue,temp)
agzj <- agzj %>% rename(rate = standardisedtraitvalue)

agzj$rate <- 1/agzj$rate

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(agzj$temp, agzj$rate, model_name = 'pawar_2018')

ag_fits <- nest(agzj, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


ag_preds <- mutate(ag_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
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
ggplot(ag_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, agzj)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) 

#____________________

fit_nlsLM6 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = agzj,
                                start = coef(ag_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(agzj)))

# bootstrap using residual resampling
boot6 <- Boot(fit_nlsLM6, method = 'residual')

boot6_preds <- boot6$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(agzj$temp), max(agzj$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot6_conf_preds <- group_by(boot6_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), ag_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot6_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), agzj, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/zj')

# calculate params with CIs

extra_params6 <- calc_params(fit_nlsLM6) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params6 <- Boot(fit_nlsLM6, f = function(x){unlist(calc_params(x))}, 
                         labels = names(calc_params(fit_nlsLM6)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params6 <- left_join(ci_extra_params6, extra_params6)

# dataset for preds with 95% CIs
angzj <- cbind(ag_preds,boot6_conf_preds) 
angzj <- angzj[-1:-2]
angzj$.fitted <- 1/angzj$.fitted
angzj$conf_lower <- 1/angzj$conf_lower
angzj$conf_upper <- 1/angzj$conf_upper
angzj$species <- 'Anthonomus grandis'

angzj <- angzj %>% rename(zj = .fitted, zj_lwr = conf_lower, zj_upr = conf_upper)
angzj <- angzj %>% select(zj, zj_lwr,zj_upr)
angzj <- angzj %>% mutate(zj = replace(zj, zj > 1, 1))
angzj <- angzj %>% mutate(zj_lwr = replace(zj_lwr, zj_lwr > 1, 1))

agi <- cbind(agi,angzj)

#______________________________________

# Female mortality rate

agz <- ag %>% filter(standardisedtraitname == "Female Mortality Rate")

agz <- agz %>% select(standardisedtraitvalue,temp)
agz <- agz %>% rename(rate = standardisedtraitvalue)

agz$rate <- 1/agz$rate


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(agz$temp, agz$rate, model_name = 'pawar_2018')

bg_fits <- nest(agz, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


bg_preds <- mutate(bg_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
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

ggplot(bg_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, agz)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM7 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = agz,
                                start = coef(bg_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(agz)))

# bootstrap using residual resampling
boot7 <- Boot(fit_nlsLM7, method = 'residual')

boot7_preds <- boot7$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(agz$temp), max(agz$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot7_conf_preds <- group_by(boot7_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), bg_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot7_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), agz, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/z')

# calculate params with CIs

extra_params7 <- calc_params(fit_nlsLM7) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params7 <- Boot(fit_nlsLM7, f = function(x){unlist(calc_params(x))}, 
                         labels = names(calc_params(fit_nlsLM7)), R = 200, method = 'residual') %>%
                              confint(., method = 'bca') %>%
                              as.data.frame() %>%
                              rename(conf_lower = 1, conf_upper = 2) %>%
                              rownames_to_column(., var = 'param') %>%
                              mutate(method = 'residual bootstrap')

# params
ci_extra_params7 <- left_join(ci_extra_params7, extra_params7)

# dataset for preds with 95% CIs

agnz <- cbind(bg_preds,boot7_conf_preds) 
agnz <- agnz[-1:-2]
agnz$.fitted <- 1/agnz$.fitted
agnz$conf_lower <- 1/agnz$conf_lower
agnz$conf_upper <- 1/agnz$conf_upper

agnz <- agnz %>% rename(z = .fitted, z_lwr = conf_lower, z_upr = conf_upper)

agnz$species <- 'Anthonomus grandis'

agnz <- agnz %>% select(z, z_lwr,z_upr)

agnz <- agnz %>%                              
  mutate(z = replace(z, z > 1, 1))

agnz <- agnz %>%                              
  mutate(z_lwr = replace(z_lwr, z_lwr > 1, 1))

agi <- cbind(agi,agnz)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4

# Fecundity rate

agf <- ag %>% filter(standardisedtraitname == "Fecundity Rate")

agf <- agf %>% select(standardisedtraitvalue,temp)
agf <- agf %>% rename(rate = standardisedtraitvalue)


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(agf$temp, agf$rate, model_name = 'pawar_2018')

af_fits <- nest(agf, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


af_preds <- mutate(af_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
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

ggplot(af_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, agf)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM8 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = agf,
                                start = coef(af_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(agf)))

# bootstrap using residual resampling
boot8 <- Boot(fit_nlsLM8, method = 'residual')

boot8_preds <- boot8$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(agf$temp), max(agf$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot8_conf_preds <- group_by(boot8_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs
ggplot() +
  geom_line(aes(temp, .fitted), af_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot8_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), agf, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12)

# calculate params with CIs
extra_params8 <- calc_params(fit_nlsLM8) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params8 <- Boot(fit_nlsLM8, f = function(x){unlist(calc_params(x))}, 
                         labels = names(calc_params(fit_nlsLM8)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params8 <- left_join(ci_extra_params8, extra_params8)

# dataset for preds with 95% CIs

agnf <- cbind(af_preds,boot8_conf_preds) 
agnf <- agnf[-1:-2]
agnf$trait <- 'bpk'
agnf$species <- 'Anthonomus grandis'
agnf <- agnf %>% rename(bpk = .fitted, bpk_lwr = conf_lower, bpk_upr = conf_upper)
agnf <- agnf %>% select(bpk, bpk_lwr,bpk_upr)
agaz <- cbind(agi,agnf)

agaz$kappa <- 10^-6

#---------------------------------- Define parameters ----------------------------------#

zj    <- agaz$zj
alpha <- agaz$alpha
z     <- agaz$z
bpk   <- agaz$bpk
k     <- agaz$kappa


# Calculate rmax

agaz$rmax <- (((k+z)*((log(bpk/(k+z)))-(alpha*zj)))/(alpha*(k+z)+1))

# lower

zj_lwr    <- agaz$zj_lwr
alpha_lwr <- agaz$alpha_lwr
z_lwr     <- agaz$z_lwr
bpk_lwr   <- agaz$bpk_lwr
k         <- agaz$kappa

agaz$rmax_lwr <- (((k+z_lwr)*((log(bpk_lwr/(k+z_lwr)))-(alpha_lwr*zj_lwr)))/(alpha_lwr*(k+z_lwr)+1))

# upper 

zj_upr    <- agaz$zj_upr
alpha_upr <- agaz$alpha_upr
z_upr     <- agaz$z_upr
bpk_upr   <- agaz$bpk_upr
k         <- agaz$kappa

agaz$rmax_upr <- (((k+z_upr)*((log(bpk_upr/(k+z_upr)))-(alpha_upr*zj_upr)))/(alpha_upr*(k+z_upr)+1))

agaz <- agaz %>%                              
  mutate(rmax = replace(rmax, rmax < -0.1, -0.1), 
         rmax_lwr = replace(rmax_lwr, rmax_lwr < -0.1, -0.1),
         rmax_upr = replace(rmax_upr, rmax_upr < -0.1, -0.1))

ag_max <- select(agaz, temp, species, rmax)
ag_max <- ag_max %>% rename(estimate=rmax)
ag_max <- filter(ag_max,temp > 12.8)
ag_max <- filter(ag_max,temp < 37)


p2 <- ggplot(agaz, aes(temp, rmax,fill=species,colour=species,shape=species))+
  ggtitle("b")+
  geom_ribbon(aes(temp, ymin = rmax_lwr, ymax = rmax_upr), alpha=0.4,col = NA)+
  theme_bw(base_size = 12)+
  geom_line(aes(temp,estimate), ag_max)

p1+p2


#_____________________________________________________#
#               Rhopalosiphum maidis                  #
#_____________________________________________________#


rm <- df %>% filter(species == "Rhopalosiphum maidis")

#- Development Rate

rmdv <- rm %>% filter(standardisedtraitname == "Juvenile Development Rate")

rmdv <- rmdv %>% select(standardisedtraitvalue,temp)
rmdv <- rmdv %>% rename(rate = standardisedtraitvalue)

rmdv$rate <- rmdv$rate + 10^-6

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(rmdv$temp, rmdv$rate, model_name = 'pawar_2018')

rm_fits <- nest(rmdv, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


rm_preds <- mutate(rm_fits, new_data = map(data, ~tibble(temp = seq(min(0), max(.x$temp), length.out = 100))))%>% 
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

ggplot(rm_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, rmdv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  theme(legend.position = "none")

#____________________

fit_nlsLM9 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = rmdv,
                                start = coef(rm_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(rmdv)))

# bootstrap using residual resampling
boot9 <- Boot(fit_nlsLM9, method = 'residual')

boot9_preds <- boot9$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(0), max(rmdv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot9_conf_preds <- group_by(boot9_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs
ggplot() +
  geom_line(aes(temp, .fitted), rm_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot9_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), rmdv, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12)

# calculate params with CIs
extra_params9 <- calc_params(fit_nlsLM9) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params9 <- Boot(fit_nlsLM9, f = function(x){unlist(calc_params(x))}, 
                         labels = names(calc_params(fit_nlsLM9)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')


# params
ci_extra_params9 <- left_join(ci_extra_params9, extra_params9)

# dataset for preds with 95% CIs
rmi <- cbind(rm_preds,boot9_conf_preds) 
rmi <- rmi[-1:-2]
rmi$.fitted <- 1/rmi$.fitted
rmi$conf_lower <- 1/rmi$conf_lower
rmi$conf_upper <- 1/rmi$conf_upper

rmi <- rmi %>% rename(alpha = .fitted, alpha_lwr = conf_lower, alpha_upr = conf_upper)

rmi$species <- 'Rhopalosiphum maidis'

rmi <- rmi %>% select(temp, alpha, alpha_lwr,alpha_upr, species)

#____________________________________________________

# Juvenile mortality rate

rmzj <- rm %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

rmzj <- rmzj %>% select(standardisedtraitvalue,temp)
rmzj <- rmzj %>% rename(rate = standardisedtraitvalue)

rmzj$rate <- 1/rmzj$rate

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(rmzj$temp, rmzj$rate, model_name = 'pawar_2018')

rzj_fits <- nest(rmzj, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


rzj_preds <- mutate(rzj_fits, new_data = map(data, ~tibble(temp = seq(min(0), max(.x$temp), length.out = 100))))%>% 
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

ggplot(rzj_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, rmzj)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) 

#____________________

fit_nlsLM10 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = rmzj,
                                 start = coef(rzj_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(rmzj)))

# bootstrap using residual resampling
boot10 <- Boot(fit_nlsLM10, method = 'residual')

boot10_preds <- boot10$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(0), max(rmzj$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot10_conf_preds <- group_by(boot10_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), rzj_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot10_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), rmzj, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/zj')

# calculate params with CIs

extra_params10 <- calc_params(fit_nlsLM10) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params10 <- Boot(fit_nlsLM10, f = function(x){unlist(calc_params(x))}, 
                          labels = names(calc_params(fit_nlsLM10)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params10 <- left_join(ci_extra_params10, extra_params10)

# dataset for preds with 95% CIs

rmazj <- cbind(rzj_preds,boot10_conf_preds) 
rmazj <- rmazj[-1:-2]
rmazj$.fitted <- 1/rmazj$.fitted
rmazj$conf_lower <- 1/rmazj$conf_lower
rmazj$conf_upper <- 1/rmazj$conf_upper
rmazj$species <- 'Rhopalosiphum maidis'

rmazj <- rmazj %>% rename(zj = .fitted, zj_lwr = conf_lower, zj_upr = conf_upper)
rmazj <- rmazj %>% select(zj, zj_lwr,zj_upr)
rmazj <- rmazj %>% mutate(zj = replace(zj, zj > 1, 1))
rmazj <- rmazj %>% mutate(zj_lwr = replace(zj_lwr, zj_lwr > 1, 1))

rmi <- cbind(rmi,rmazj)

#______________________________________

# Female mortality rate

rmz <- rm %>% filter(standardisedtraitname == "Female Mortality Rate")

rmz <- rmz %>% select(standardisedtraitvalue,temp)
rmz <- rmz %>% rename(rate = standardisedtraitvalue)

rmz$rate <- 1/rmz$rate


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(rmz$temp, rmz$rate, model_name = 'pawar_2018')

rmz_fits <- nest(rmz, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


rmz_preds <- mutate(rmz_fits, new_data = map(data, ~tibble(temp = seq(min(0), max(.x$temp), length.out = 100))))%>% 
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

ggplot(rmz_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, rmz)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM12 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = rmz,
                                 start = coef(rmz_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(rmz)))

# bootstrap using residual resampling
boot12 <- Boot(fit_nlsLM12, method = 'residual')

boot12_preds <- boot12$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(0), max(rmz$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot12_conf_preds <- group_by(boot12_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs
ggplot() +
  geom_line(aes(temp, .fitted), rmz_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot12_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), rmz, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/z')

# calculate params with CIs

extra_params12 <- calc_params(fit_nlsLM12) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params12 <- Boot(fit_nlsLM12, f = function(x){unlist(calc_params(x))}, 
                          labels = names(calc_params(fit_nlsLM12)), R = 200, method = 'residual') %>%
                          confint(., method = 'bca') %>%
                          as.data.frame() %>%
                          rename(conf_lower = 1, conf_upper = 2) %>%
                          rownames_to_column(., var = 'param') %>%
                          mutate(method = 'residual bootstrap')

# params
ci_extra_params12 <- left_join(ci_extra_params12, extra_params12)

# dataset for preds with 95% CIs

rmaz <- cbind(rmz_preds,boot12_conf_preds) 
rmaz <- rmaz[-1:-2]
rmaz$.fitted <- 1/rmaz$.fitted
rmaz$conf_lower <- 1/rmaz$conf_lower
rmaz$conf_upper <- 1/rmaz$conf_upper

rmaz <- rmaz %>% rename(z = .fitted, z_lwr = conf_lower, z_upr = conf_upper)

rmaz$species <- 'Rhopalosiphum maidis'

rmaz <- rmaz %>% select(z, z_lwr,z_upr)
rmaz <- rmaz %>% mutate(z = replace(z, z > 1, 1))
rmaz <- rmaz %>% mutate(z_lwr = replace(z_lwr, z_lwr > 1, 1))

rmi <- cbind(rmi,rmaz)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4

# Fecundity rate

rmf <- rm %>% filter(standardisedtraitname == "Fecundity Rate")
rmf <- rmf %>% select(standardisedtraitvalue,temp)
rmf <- rmf %>% rename(rate = standardisedtraitvalue)


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(rmf$temp, rmf$rate, model_name = 'pawar_2018')

rmf_fits <- nest(rmf, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


rmf_preds <- mutate(rmf_fits, new_data = map(data, ~tibble(temp = seq(min(0), max(.x$temp), length.out = 100))))%>% 
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

ggplot(rmf_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, rmf)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM13 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = rmf,
                                 start = coef(rmf_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(rmf)))

# bootstrap using residual resampling
boot13 <- Boot(fit_nlsLM13, method = 'residual')

boot13_preds <- boot13$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(0), max(rmf$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot13_conf_preds <- group_by(boot13_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), rmf_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot13_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), rmf, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12)

# calculate params with CIs

extra_params13 <- calc_params(fit_nlsLM13) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params13 <- Boot(fit_nlsLM13, f = function(x){unlist(calc_params(x))}, 
                          labels = names(calc_params(fit_nlsLM13)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params13 <- left_join(ci_extra_params13, extra_params13)

# dataset for preds with 95% CIs

rhmf <- cbind(rmf_preds,boot13_conf_preds) 
rhmf <- rhmf[-1:-2]
rhmf$trait <- 'bpk'
rhmf$species <- 'Rhopalosiphum maidis'
rhmf <- rhmf %>% rename(bpk = .fitted, bpk_lwr = conf_lower, bpk_upr = conf_upper)
rhmf <- rhmf %>% select(bpk, bpk_lwr,bpk_upr)
rhaz <- cbind(rmi,rhmf)

rhaz$kappa <- 10^-6

#---------------------------------- Define parameters ----------------------------------#

zj    <- rhaz$zj
alpha <- rhaz$alpha
z     <- rhaz$z
bpk   <- rhaz$bpk
k     <- rhaz$kappa


# Calculate rmax

rhaz$rmax <- (((k+z)*((log(bpk/(k+z)))-(alpha*zj)))/(alpha*(k+z)+1))

# lower

zj_lwr    <- rhaz$zj_lwr
alpha_lwr <- rhaz$alpha_lwr
z_lwr     <- rhaz$z_lwr
bpk_lwr   <- rhaz$bpk_lwr
k         <- rhaz$kappa

rhaz$rmax_lwr <- (((k+z_lwr)*((log(bpk_lwr/(k+z_lwr)))-(alpha_lwr*zj_lwr)))/(alpha_lwr*(k+z_lwr)+1))

# upper 

zj_upr    <- rhaz$zj_upr
alpha_upr <- rhaz$alpha_upr
z_upr     <- rhaz$z_upr
bpk_upr   <- rhaz$bpk_upr
k         <- rhaz$kappa

rhaz$rmax_upr <- (((k+z_upr)*((log(bpk_upr/(k+z_upr)))-(alpha_upr*zj_upr)))/(alpha_upr*(k+z_upr)+1))

rhaz <- rhaz %>%                              
  mutate(rmax = replace(rmax, rmax < -0.1, -0.1), 
         rmax_lwr = replace(rmax_lwr, rmax_lwr < -0.1, -0.1),
         rmax_upr = replace(rmax_upr, rmax_upr < -0.1, -0.1))

rm_max <- select(rhaz, temp, species, rmax)
rm_max <- rm_max %>% rename(estimate=rmax)
rm_max <- filter(rm_max,temp < 32.5)


p3 <- ggplot(rhaz, aes(temp, rmax,fill=species,colour=species,shape=species))+
  ggtitle("a")+
  geom_ribbon(aes(temp, ymin = rmax_lwr, ymax = rmax_upr), alpha=0.4,col = NA)+
  theme_bw(base_size = 12)+
  geom_line(aes(temp,estimate), rm_max)


p2+p1+p3

#_____________________________________________________#
#               Muscidifurax zaraptor                  #
#_____________________________________________________#

mz <- df %>% filter(species == "Muscidifurax zaraptor")

#- Development Rate

mzdv <- mz %>% filter(standardisedtraitname == "Juvenile Development Rate")

mzdv <- mzdv %>% select(standardisedtraitvalue,temp)
mzdv <- mzdv %>% rename(rate = standardisedtraitvalue)

mzdv$rate <- mzdv$rate + 10^-6

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(mzdv$temp, mzdv$rate, model_name = 'pawar_2018')

mz_fits <- nest(mzdv, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


mz_preds <- mutate(mz_fits, new_data = map(data, ~tibble(temp = seq(min(8), max(40), length.out = 100))))%>% 
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

ggplot(mz_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, mzdv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  theme(legend.position = "none")

#____________________

fit_nlsLM14 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = mzdv,
                                 start = coef(mz_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(mzdv)))

# bootstrap using residual resampling
boot14 <- Boot(fit_nlsLM14, method = 'residual')

boot14_preds <- boot14$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(8), max(40), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot14_conf_preds <- group_by(boot14_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), mz_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot14_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), mzdv, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12)

# calculate params with CIs

extra_params14 <- calc_params(fit_nlsLM14) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params14 <- Boot(fit_nlsLM14, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM14)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')


# params
ci_extra_params14 <- left_join(ci_extra_params14, extra_params14)

# dataset for preds with 95% CIs
mzi <- cbind(mz_preds,boot14_conf_preds) 
mzi <- mzi[-1:-2]
mzi$.fitted <- 1/mzi$.fitted
mzi$conf_lower <- 1/mzi$conf_lower
mzi$conf_upper <- 1/mzi$conf_upper

mzi <- mzi %>% rename(alpha = .fitted, alpha_lwr = conf_lower, alpha_upr = conf_upper)

mzi$species <- 'Muscidifurax zaraptor'

mzi <- mzi %>% select(temp, alpha, alpha_lwr,alpha_upr, species)

#____________________________________________________

# Juvenile mortality rate

mzzj <- mz %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

mzzj <- mzzj %>% select(standardisedtraitvalue,temp)
mzzj <- mzzj %>% rename(rate = standardisedtraitvalue)

mzzj$rate <- 1/mzzj$rate

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(mzzj$temp, mzzj$rate, model_name = 'pawar_2018')

mzj_fits <- nest(mzzj, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


mzj_preds <- mutate(mzj_fits, new_data = map(data, ~tibble(temp = seq(min(8), max(40), length.out = 100))))%>% 
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

ggplot(mzj_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, mzzj)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) 

#____________________

fit_nlsLM15 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = mzzj,
                                 start = coef(mzj_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(mzzj)))

# bootstrap using residual resampling
boot15 <- Boot(fit_nlsLM15, method = 'residual')

boot15_preds <- boot15$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(8), max(40), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot15_conf_preds <- group_by(boot15_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), mzj_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot15_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), mzzj, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/zj')

# calculate params with CIs

extra_params15 <- calc_params(fit_nlsLM15) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params15 <- Boot(fit_nlsLM15, f = function(x){unlist(calc_params(x))}, 
                          labels = names(calc_params(fit_nlsLM15)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params15 <- left_join(ci_extra_params15, extra_params15)

# dataset for preds with 95% CIs

mzaji <- cbind(mzj_preds,boot15_conf_preds) 
mzaji <- mzaji[-1:-2]
mzaji$.fitted <- 1/mzaji$.fitted
mzaji$conf_lower <- 1/mzaji$conf_lower
mzaji$conf_upper <- 1/mzaji$conf_upper
mzaji$species <- 'Muscidifurax zaraptor'

mzaji <- mzaji %>% rename(zj = .fitted, zj_lwr = conf_lower, zj_upr = conf_upper)
mzaji <- mzaji %>% select(zj, zj_lwr,zj_upr)
mzaji <- mzaji %>% mutate(zj = replace(zj, zj > 1, 1))
mzaji <- mzaji %>% mutate(zj_lwr = replace(zj_lwr, zj_lwr > 1, 1))

mzi <- cbind(mzi,mzaji)

#______________________________________

# Female mortality rate

mzz <- mz %>% filter(standardisedtraitname == "Female Mortality Rate")

mzz <- mzz %>% select(standardisedtraitvalue,temp)
mzz <- mzz %>% rename(rate = standardisedtraitvalue)

mzz$rate <- 1/mzz$rate


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(mzz$temp, mzz$rate, model_name = 'pawar_2018')

mzz_fits <- nest(mzz, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


mzz_preds <- mutate(mzz_fits, new_data = map(data, ~tibble(temp = seq(min(8), max(40), length.out = 100))))%>% 
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

ggplot(mzz_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, mzz)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM16 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = mzz,
                                 start = coef(mzz_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(mzz)))

# bootstrap using residual resampling
boot16 <- Boot(fit_nlsLM16, method = 'residual')

boot16_preds <- boot16$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(8), max(40), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot16_conf_preds <- group_by(boot16_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), mzz_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot16_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), mzz, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/z')

# calculate params with CIs

extra_params16 <- calc_params(fit_nlsLM16) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params16 <- Boot(fit_nlsLM16, f = function(x){unlist(calc_params(x))}, 
                          labels = names(calc_params(fit_nlsLM16)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params16 <- left_join(ci_extra_params16, extra_params16)

# dataset for preds with 95% CIs

mzaz <- cbind(mzz_preds,boot16_conf_preds) 
mzaz <- mzaz[-1:-2]
mzaz$.fitted <- 1/mzaz$.fitted
mzaz$conf_lower <- 1/mzaz$conf_lower
mzaz$conf_upper <- 1/mzaz$conf_upper

mzaz <- mzaz %>% rename(z = .fitted, z_lwr = conf_lower, z_upr = conf_upper)
mzaz$species <- 'Muscidifurax zaraptor'
mzaz <- mzaz %>% select(z, z_lwr,z_upr)
mzaz <- mzaz %>% mutate(z = replace(z, z > 1, 1))
mzaz <- mzaz %>% mutate(z_lwr = replace(z_lwr, z_lwr > 1, 1))

mzi <- cbind(mzi,mzaz)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4

# Fecundity rate

mzf <- mz %>% filter(standardisedtraitname == "Fecundity Rate")

mzf <- mzf %>% select(standardisedtraitvalue,temp)
mzf <- mzf %>% rename(rate = standardisedtraitvalue)


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(mzf$temp, mzf$rate, model_name = 'pawar_2018')

mzf_fits <- nest(mzf, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


mzf_preds <- mutate(mzf_fits, new_data = map(data, ~tibble(temp = seq(min(8), max(40), length.out = 100))))%>% 
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

ggplot(mzf_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, mzf)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM17 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                 data = mzf,
                                 start = coef(mzf_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(mzf)))

# bootstrap using residual resampling
boot17 <- Boot(fit_nlsLM17, method = 'residual')

boot17_preds <- boot17$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(8), max(40), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot17_conf_preds <- group_by(boot17_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), mzf_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot17_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), mzf, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12)

# calculate params with CIs

extra_params17 <- calc_params(fit_nlsLM17) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params17 <- Boot(fit_nlsLM17, f = function(x){unlist(calc_params(x))}, 
                          labels = names(calc_params(fit_nlsLM17)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params17 <- left_join(ci_extra_params17, extra_params17)

# dataset for preds with 95% CIs

mzaf <- cbind(mzf_preds,boot17_conf_preds) 
mzaf <- mzaf[-1:-2]
mzaf$trait <- 'bpk'
mzaf$species <- 'Muscidifurax zaraptor'
mzaf <- mzaf %>% rename(bpk = .fitted, bpk_lwr = conf_lower, bpk_upr = conf_upper)
mzaf <- mzaf %>% select(bpk, bpk_lwr,bpk_upr)
mhaz <- cbind(mzi,mzaf)

mhaz$kappa <- 10^-6

#---------------------------------- Define parameters ----------------------------------#

zj    <- mhaz$zj
alpha <- mhaz$alpha
z     <- mhaz$z
bpk   <- mhaz$bpk
k     <- mhaz$kappa


# Calculate rmax

mhaz$rmax <- (((k+z)*((log(bpk/(k+z)))-(alpha*zj)))/(alpha*(k+z)+1))

# lower

zj_lwr    <- mhaz$zj_lwr
alpha_lwr <- mhaz$alpha_lwr
z_lwr     <- mhaz$z_lwr
bpk_lwr   <- mhaz$bpk_lwr
k         <- mhaz$kappa

mhaz$rmax_lwr <- (((k+z_lwr)*((log(bpk_lwr/(k+z_lwr)))-(alpha_lwr*zj_lwr)))/(alpha_lwr*(k+z_lwr)+1))

# upper 

zj_upr    <- mhaz$zj_upr
alpha_upr <- mhaz$alpha_upr
z_upr     <- mhaz$z_upr
bpk_upr   <- mhaz$bpk_upr
k         <- mhaz$kappa

mhaz$rmax_upr <- (((k+z_upr)*((log(bpk_upr/(k+z_upr)))-(alpha_upr*zj_upr)))/(alpha_upr*(k+z_upr)+1))

mhaz <- mhaz %>%                              
  mutate(rmax = replace(rmax, rmax < -0.1, -0.1), 
         rmax_lwr = replace(rmax_lwr, rmax_lwr < -0.1, -0.1),
         rmax_upr = replace(rmax_upr, rmax_upr < -0.1, -0.1))

mz_max <- select(mhaz, temp, species, rmax)
mz_max <- mz_max %>% rename(estimate=rmax)
mz_max <- mz_max %>% filter(temp < 37.7)



p4 <- ggplot(mhaz, aes(temp, rmax,fill=species,colour=species,shape=species))+
  ggtitle("a")+
  geom_ribbon(aes(temp, ymin = rmax_lwr, ymax = rmax_upr), alpha=0.4,col = NA)+
  theme_bw(base_size = 12)+
  geom_line(aes(temp,estimate), mz_max)


((p2|p1)/(p3|p4))


#_____________________________________________________#
#               Aphis nasturtii                       #
#_____________________________________________________#


an <- df %>% filter(species == "Aphis nasturtii")

#- Development Rate

andv <- an %>% filter(standardisedtraitname == "Juvenile Development Rate")

andv <- andv %>% select(standardisedtraitvalue,temp)
andv <- andv %>% rename(rate = standardisedtraitvalue)

andv$rate <- andv$rate + 10^-6

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(andv$temp, andv$rate, model_name = 'pawar_2018')

an_fits <- nest(andv, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 19),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


an_preds <- mutate(an_fits, new_data = map(data, ~tibble(temp = seq(min(0), max(.x$temp), length.out = 100))))%>% 
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

ggplot(an_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, andv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  theme(legend.position = "none")

#____________________

fit_nlsLM18 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 19),
                                 data = andv,
                                 start = coef(an_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(andv)))

# bootstrap using residual resampling
boot18 <- Boot(fit_nlsLM18, method = 'residual')

boot18_preds <- boot18$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(0), max(andv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 19))

# calculate bootstrapped confidence intervals
boot18_conf_preds <- group_by(boot18_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), an_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot18_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), andv, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12)

# calculate params with CIs

extra_params18 <- calc_params(fit_nlsLM18) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params18 <- Boot(fit_nlsLM18, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM18)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')


# params
ci_extra_params18 <- left_join(ci_extra_params18, extra_params18)

# dataset for preds with 95% CIs
ani <- cbind(an_preds,boot18_conf_preds) 
ani <- ani[-1:-2]
ani$.fitted <- 1/ani$.fitted
ani$conf_lower <- 1/ani$conf_lower
ani$conf_upper <- 1/ani$conf_upper

ani <- ani %>% rename(alpha = .fitted, alpha_lwr = conf_lower, alpha_upr = conf_upper)

ani$species <- 'Aphis nasturtii'

ani <- ani %>% select(temp, alpha, alpha_lwr,alpha_upr, species)

#____________________________________________________

# Juvenile mortality rate

anzj <- an %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

anzj <- anzj %>% select(standardisedtraitvalue,temp)
anzj <- anzj %>% rename(rate = standardisedtraitvalue)

anzj$rate <- 1/anzj$rate

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(anzj$temp, anzj$rate, model_name = 'pawar_2018')

anzj_fits <- nest(anzj, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 19),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


anzj_preds <- mutate(anzj_fits, new_data = map(data, ~tibble(temp = seq(min(0), max(.x$temp), length.out = 100))))%>% 
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

ggplot(anzj_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, anzj)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) 

#____________________

fit_nlsLM19 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 19),
                                 data = anzj,
                                 start = coef(anzj_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(anzj)))

# bootstrap using residual resampling
boot19 <- Boot(fit_nlsLM19, method = 'residual')

boot19_preds <- boot19$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(0), max(anzj$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 19))

# calculate bootstrapped confidence intervals
boot19_conf_preds <- group_by(boot19_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), anzj_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot19_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), anzj, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/zj')

# calculate params with CIs

extra_params19 <- calc_params(fit_nlsLM19) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params19 <- Boot(fit_nlsLM19, f = function(x){unlist(calc_params(x))}, 
                          labels = names(calc_params(fit_nlsLM19)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params19 <- left_join(ci_extra_params19, extra_params19)

# dataset for preds with 95% CIs

azaji <- cbind(anzj_preds,boot19_conf_preds) 
azaji <- azaji[-1:-2]
azaji$.fitted <- 1/azaji$.fitted
azaji$conf_lower <- 1/azaji$conf_lower
azaji$conf_upper <- 1/azaji$conf_upper
azaji$species <- 'Aphis nasturtii'

azaji <- azaji %>% rename(zj = .fitted, zj_lwr = conf_lower, zj_upr = conf_upper)
azaji <- azaji %>% select(zj, zj_lwr,zj_upr)
azaji <- azaji %>% mutate(zj = replace(zj, zj > 1, 1))
azaji <- azaji %>% mutate(zj_lwr = replace(zj_lwr, zj_lwr > 1, 1))

ani <- cbind(ani,azaji)

#______________________________________

# Female mortality rate

anz <- an %>% filter(standardisedtraitname == "Female Mortality Rate")

anz <- anz %>% select(standardisedtraitvalue,temp)
anz <- anz %>% rename(rate = standardisedtraitvalue)

anz$rate <- 1/anz$rate


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(anz$temp, anz$rate, model_name = 'pawar_2018')

anz_fits <- nest(anz, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 19),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


anz_preds <- mutate(anz_fits, new_data = map(data, ~tibble(temp = seq(min(0), max(.x$temp), length.out = 100))))%>% 
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

ggplot(anz_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, anz)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM20 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 19),
                                 data = anz,
                                 start = coef(anz_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(anz)))

# bootstrap using residual resampling
boot20 <- Boot(fit_nlsLM20, method = 'residual')

boot20_preds <- boot20$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(0), max(anz$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 19))

# calculate bootstrapped confidence intervals
boot20_conf_preds <- group_by(boot20_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), anz_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot20_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), anz, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/z')

# calculate params with CIs

extra_params20 <- calc_params(fit_nlsLM20) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params20 <- Boot(fit_nlsLM20, f = function(x){unlist(calc_params(x))}, 
                          labels = names(calc_params(fit_nlsLM20)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params20 <- left_join(ci_extra_params20, extra_params20)

# dataset for preds with 95% CIs

azaz <- cbind(anz_preds,boot20_conf_preds) 
azaz <- azaz[-1:-2]
azaz$.fitted <- 1/azaz$.fitted
azaz$conf_lower <- 1/azaz$conf_lower
azaz$conf_upper <- 1/azaz$conf_upper

azaz <- azaz %>% rename(z = .fitted, z_lwr = conf_lower, z_upr = conf_upper)

azaz$species <- 'Aphis nasturtii'

azaz <- azaz %>% select(z, z_lwr,z_upr)

azaz <- azaz %>%                              
  mutate(z = replace(z, z > 1, 1))

azaz <- azaz %>%                              
  mutate(z_lwr = replace(z_lwr, z_lwr > 1, 1))

ani <- cbind(ani,azaz)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4

# Fecundity rate

anf <- an %>% filter(standardisedtraitname == "Fecundity Rate")

anf <- anf %>% select(standardisedtraitvalue,temp)
anf <- anf %>% rename(rate = standardisedtraitvalue)


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(anf$temp, anf$rate, model_name = 'pawar_2018')

anf_fits <- nest(anf, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 19),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


anf_preds <- mutate(anf_fits, new_data = map(data, ~tibble(temp = seq(min(0), max(.x$temp), length.out = 100))))%>% 
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

ggplot(anf_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, anf)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM21 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 19),
                                 data = anf,
                                 start = coef(anf_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(anf)))

# bootstrap using residual resampling
boot21 <- Boot(fit_nlsLM21, method = 'residual')

boot21_preds <- boot21$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(0), max(anf$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 19))

# calculate bootstrapped confidence intervals
boot21_conf_preds <- group_by(boot21_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), anf_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot21_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), anf, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12)

# calculate params with CIs

extra_params21 <- calc_params(fit_nlsLM21) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params21 <- Boot(fit_nlsLM21, f = function(x){unlist(calc_params(x))}, 
                          labels = names(calc_params(fit_nlsLM21)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params21 <- left_join(ci_extra_params21, extra_params21)

# dataset for preds with 95% CIs

anaf <- cbind(anf_preds,boot21_conf_preds) 
anaf <- anaf[-1:-2]
anaf$trait <- 'bpk'
anaf$species <- 'Aphis nasturtii'
anaf <- anaf %>% rename(bpk = .fitted, bpk_lwr = conf_lower, bpk_upr = conf_upper)
anaf <- anaf %>% select(bpk, bpk_lwr,bpk_upr)
anaz <- cbind(ani,anaf)

anaz$kappa <- 10^-6

#---------------------------------- Define parameters ----------------------------------#

zj    <- anaz$zj
alpha <- anaz$alpha
z     <- anaz$z
bpk   <- anaz$bpk
k     <- anaz$kappa


# Calculate rmax

anaz$rmax <- (((k+z)*((log(bpk/(k+z)))-(alpha*zj)))/(alpha*(k+z)+1))

# lower

zj_lwr    <- anaz$zj_lwr
alpha_lwr <- anaz$alpha_lwr
z_lwr     <- anaz$z_lwr
bpk_lwr   <- anaz$bpk_lwr
k         <- anaz$kappa

anaz$rmax_lwr <- (((k+z_lwr)*((log(bpk_lwr/(k+z_lwr)))-(alpha_lwr*zj_lwr)))/(alpha_lwr*(k+z_lwr)+1))

# upper 

zj_upr    <- anaz$zj_upr
alpha_upr <- anaz$alpha_upr
z_upr     <- anaz$z_upr
bpk_upr   <- anaz$bpk_upr
k         <- anaz$kappa

anaz$rmax_upr <- (((k+z_upr)*((log(bpk_upr/(k+z_upr)))-(alpha_upr*zj_upr)))/(alpha_upr*(k+z_upr)+1))

anaz <- anaz %>%                              
  mutate(rmax = replace(rmax, rmax < -0.1, -0.1), 
         rmax_lwr = replace(rmax_lwr, rmax_lwr < -0.1, -0.1),
         rmax_upr = replace(rmax_upr, rmax_upr < -0.1, -0.1))

an_max <- select(anaz, temp, species, rmax)
an_max <- an_max %>% rename(estimate=rmax)


p5 <- ggplot(anaz, aes(temp, rmax,fill=species,colour=species,shape=species))+
  ggtitle("a")+
  geom_ribbon(aes(temp, ymin = rmax_lwr, ymax = rmax_upr), alpha=0.4,col = NA)+
  theme_bw(base_size = 12)+
  geom_line(aes(temp,estimate), an_max)


((p2|p1)/(p3|p4)/(p5+plot_spacer()))


#_____________________________________________________#
#               Paracoccus marginatu                  #
#_____________________________________________________#


pm <- df %>% filter(species == "Paracoccus marginatu")

#- Development Rate

pmdv <- pm %>% filter(standardisedtraitname == "Juvenile Development Rate")

pmdv <- pmdv %>% select(standardisedtraitvalue,temp)
pmdv <- pmdv %>% rename(rate = standardisedtraitvalue)

pmdv$rate <- pmdv$rate + 10^-6

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(pmdv$temp, pmdv$rate, model_name = 'pawar_2018')

pm_fits <- nest(pmdv, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 19),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


pm_preds <- mutate(pm_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
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

ggplot(pm_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, pmdv)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) +
  theme(legend.position = "none")

#____________________

fit_nlsLM22 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 19),
                                 data = pmdv,
                                 start = coef(pm_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(pmdv)))

# bootstrap using residual resampling
boot22 <- Boot(fit_nlsLM22, method = 'residual')

boot22_preds <- boot22$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(pmdv$temp), max(pmdv$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 19))

# calculate bootstrapped confidence intervals
boot22_conf_preds <- group_by(boot22_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), pm_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot22_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), pmdv, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12)

# calculate params with CIs

extra_params22 <- calc_params(fit_nlsLM22) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params22 <- Boot(fit_nlsLM22, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM22)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')


# params
ci_extra_params22 <- left_join(ci_extra_params22, extra_params22)

# dataset for preds with 95% CIs
pmi <- cbind(pm_preds,boot22_conf_preds) 
pmi <- pmi[-1:-2]
pmi$.fitted <- 1/pmi$.fitted
pmi$conf_lower <- 1/pmi$conf_lower
pmi$conf_upper <- 1/pmi$conf_upper

pmi <- pmi %>% rename(alpha = .fitted, alpha_lwr = conf_lower, alpha_upr = conf_upper)

pmi$species <- 'Paracoccus marginatu'

pmi <- pmi %>% select(temp, alpha, alpha_lwr,alpha_upr, species)

#____________________________________________________

# Juvenile mortality rate

pmzj <- pm %>% filter(standardisedtraitname == "Juvenile Mortality Rate")

pmzj <- pmzj %>% select(standardisedtraitvalue,temp)
pmzj <- pmzj %>% rename(rate = standardisedtraitvalue)

pmzj$rate <- 1/pmzj$rate

# fit chosen model formulation in rTPC

start_vals <- get_start_vals(pmzj$temp, pmzj$rate, model_name = 'pawar_2018')

pmzj_fits <- nest(pmzj, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 19),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


pmzj_preds <- mutate(pmzj_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
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

ggplot(pmzj_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, pmzj)+
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_brewer(type = 'qual', palette = 2) 

#____________________

fit_nlsLM23 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 19),
                                 data = pmzj,
                                 start = coef(pmzj_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(pmzj)))

# bootstrap using residual resampling
boot23 <- Boot(fit_nlsLM23, method = 'residual')

boot23_preds <- boot23$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(pmzj$temp), max(pmzj$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 19))

# calculate bootstrapped confidence intervals
boot23_conf_preds <- group_by(boot23_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), pmzj_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot23_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), pmzj, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/zj')

# calculate params with CIs

extra_params23 <- calc_params(fit_nlsLM23) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params23 <- Boot(fit_nlsLM23, f = function(x){unlist(calc_params(x))}, 
                          labels = names(calc_params(fit_nlsLM23)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params23 <- left_join(ci_extra_params23, extra_params23)

# dataset for preds with 95% CIs

pzaji <- cbind(pmzj_preds,boot23_conf_preds) 
pzaji <- pzaji[-1:-2]
pzaji$.fitted <- 1/pzaji$.fitted
pzaji$conf_lower <- 1/pzaji$conf_lower
pzaji$conf_upper <- 1/pzaji$conf_upper
pzaji$species <- 'Paracoccus marginatu'

pzaji <- pzaji %>% rename(zj = .fitted, zj_lwr = conf_lower, zj_upr = conf_upper)
pzaji <- pzaji %>% select(zj, zj_lwr,zj_upr)
pzaji <- pzaji %>% mutate(zj = replace(zj, zj > 1, 1))
pzaji <- pzaji %>% mutate(zj_lwr = replace(zj_lwr, zj_lwr > 1, 1))

pmi <- cbind(pmi,pzaji)

#______________________________________

# Female mortality rate

pmz <- pm %>% filter(standardisedtraitname == "Female Mortality Rate")

pmz <- pmz %>% select(standardisedtraitvalue,temp)
pmz <- pmz %>% rename(rate = standardisedtraitvalue)

pmz$rate <- 1/pmz$rate


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(pmz$temp, pmz$rate, model_name = 'pawar_2018')

pmz_fits <- nest(pmz, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 19),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


pmz_preds <- mutate(pmz_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), 
                                                                      max(.x$temp), length.out = 100))))%>% 
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

ggplot(pmz_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, pmz)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM24 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 19),
                                 data = pmz,
                                 start = coef(pmz_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(pmz)))

# bootstrap using residual resampling
boot24 <- Boot(fit_nlsLM24, method = 'residual')

boot24_preds <- boot24$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(pmz$temp), max(pmz$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 19))

# calculate bootstrapped confidence intervals
boot24_conf_preds <- group_by(boot24_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), pmz_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot24_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), pmz, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/z')

# calculate params with CIs

extra_params24 <- calc_params(fit_nlsLM24) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params24 <- Boot(fit_nlsLM24, f = function(x){unlist(calc_params(x))}, 
                          labels = names(calc_params(fit_nlsLM24)), R = 240, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params24 <- left_join(ci_extra_params24, extra_params24)

# dataset for preds with 95% CIs

pzaz <- cbind(pmz_preds,boot24_conf_preds) 
pzaz <- pzaz[-1:-2]
pzaz$.fitted <- 1/pzaz$.fitted
pzaz$conf_lower <- 1/pzaz$conf_lower
pzaz$conf_upper <- 1/pzaz$conf_upper

pzaz <- pzaz %>% rename(z = .fitted, z_lwr = conf_lower, z_upr = conf_upper)

pzaz$species <- 'Paracoccus marginatu'

pzaz <- pzaz %>% select(z, z_lwr,z_upr)
pzaz <- pzaz %>% mutate(z = replace(z, z > 1, 1))
pzaz <- pzaz %>% mutate(z_lwr = replace(z_lwr, z_lwr > 1, 1))

pmi <- cbind(pmi,pzaz)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4

# Fecundity rate

pmf <- pm  %>% filter(standardisedtraitname == "Fecundity Rate")
pmf <- pmf %>% select(standardisedtraitvalue,temp)
pmf <- pmf %>% rename(rate = standardisedtraitvalue)


# fit chosen model formulation in rTPC

start_vals <- get_start_vals(pmf$temp, pmf$rate, model_name = 'pawar_2018')

pmf_fits <- nest(pmf, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 19),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________


pmf_preds <- mutate(pmf_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
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

ggplot(pmf_preds) +
  geom_line(aes(temp, .fitted, col = model_name)) +
  geom_point(aes(temp, rate),size=0.2,alpha=0.5, pmf)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________

fit_nlsLM25 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 19),
                                 data = pmf,
                                 start = coef(pmf_fits$pawar[[1]]),
                                 weights = rep(1, times = nrow(pmf)))

# bootstrap using residual resampling
boot25 <- Boot(fit_nlsLM25, method = 'residual')

boot25_preds <- boot25$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(pmf$temp), max(pmf$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 19))

# calculate bootstrapped confidence intervals
boot25_conf_preds <- group_by(boot25_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs

ggplot() +
  geom_line(aes(temp, .fitted), pmf_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot25_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, rate), pmf, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12)

# calculate params with CIs

extra_params25 <- calc_params(fit_nlsLM25) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params25 <- Boot(fit_nlsLM25, f = function(x){unlist(calc_params(x))}, 
                     labels = names(calc_params(fit_nlsLM25)), R = 200, method = 'residual') %>%
                     confint(., method = 'bca') %>%
                     as.data.frame() %>%
                     rename(conf_lower = 1, conf_upper = 2) %>%
                     rownames_to_column(., var = 'param') %>%
                     mutate(method = 'residual bootstrap')

# params
ci_extra_params25 <- left_join(ci_extra_params25, extra_params25)

# dataset for preds with 95% CIs

pmaf <- cbind(pmf_preds,boot25_conf_preds) 
pmaf <- pmaf[-1:-2]
pmaf$trait <- 'bpk'
pmaf$species <- 'Paracoccus marginatu'
pmaf <- pmaf %>% rename(bpk = .fitted, bpk_lwr = conf_lower, bpk_upr = conf_upper)
pmaf <- pmaf %>% select(bpk, bpk_lwr,bpk_upr)
pmaz <- cbind(pmi,pmaf)

pmaz$kappa <- 10^-6

#---------------------------------- Define parameters ----------------------------------#

zj    <- pmaz$zj
alpha <- pmaz$alpha
z     <- pmaz$z
bpk   <- pmaz$bpk
k     <- pmaz$kappa


# Calculate rmax

pmaz$rmax <- (((k+z)*((log(bpk/(k+z)))-(alpha*zj)))/(alpha*(k+z)+1))

# lower

zj_lwr    <- pmaz$zj_lwr
alpha_lwr <- pmaz$alpha_lwr
z_lwr     <- pmaz$z_lwr
bpk_lwr   <- pmaz$bpk_lwr
k         <- pmaz$kappa

pmaz$rmax_lwr <- (((k+z_lwr)*((log(bpk_lwr/(k+z_lwr)))-(alpha_lwr*zj_lwr)))/(alpha_lwr*(k+z_lwr)+1))

# upper 

zj_upr    <- pmaz$zj_upr
alpha_upr <- pmaz$alpha_upr
z_upr     <- pmaz$z_upr
bpk_upr   <- pmaz$bpk_upr
k         <- pmaz$kappa

pmaz$rmax_upr <- (((k+z_upr)*((log(bpk_upr/(k+z_upr)))-(alpha_upr*zj_upr)))/(alpha_upr*(k+z_upr)+1))

pmaz <- pmaz %>%                              
  mutate(rmax = replace(rmax, rmax < -0.1, -0.1), 
         rmax_lwr = replace(rmax_lwr, rmax_lwr < -0.1, -0.1),
         rmax_upr = replace(rmax_upr, rmax_upr < -0.1, -0.1))

pm_max <- select(pmaz, temp, species, rmax)
pm_max <- pm_max %>% rename(estimate=rmax)
pm_max <- pm_max %>% filter(temp > 15.3)


p6 <- ggplot(pmaz, aes(temp, rmax,fill=species,colour=species,shape=species))+
  ggtitle("a")+
  geom_ribbon(aes(temp, ymin = rmax_lwr, ymax = rmax_upr), alpha=0.4,col = NA)+
  theme_bw(base_size = 12)+
  geom_line(aes(temp,estimate), pm_max)


((p1|p2)/(p3|p4)/(p5|p6))


# dataset for rmax plot 

# preds 
tn_max <- ter_max %>% filter(temp < 33)
tn_max[72,3] <- 0
tn_max <- tn_max %>% filter(temp > 16.1)
tn_max[1,3] <- 0


ar_max <- ag_max %>% filter(temp > 13.2)
ar_max[1,3] <- 0
ar_max <- ar_max %>% filter(temp < 35.2)
ar_max[78,3] <- 0

rh_max <- rm_max
rh_max[92,3] <- 0
rh_max <- rh_max %>% filter(temp > 2.2)
rh_max[1,3] <- 0

mu_max <- mz_max %>% filter(temp < 37)
mu_max[90,3] <- 0
mu_max <- mu_max %>% filter(temp > 11.6)
mu_max[1,3] <- 0

at_max <- an_max %>% filter(temp < 34)
at_max[97,3] <- 0
at_max <- at_max %>% filter(temp > 9.2)
at_max[1,3] <- 0

pa_max <- pm_max %>% filter(temp > 16)
pa_max[1,3] <- 0
pa_max <- pa_max %>% filter(temp < 31.3)
pa_max[76,3] <- 0

rmpredictions <- as_tibble(rbind(tn_max,ar_max,rh_max,mu_max,at_max,pa_max))
rmpredictions$species <- as.factor(rmpredictions$species)

rmpredictions$species <- fct_relevel(rmpredictions$species, "Anthonomus grandis")
rmpredictions$species <- fct_relevel(rmpredictions$species, "Rhopalosiphum maidis", after = Inf)
rmpredictions$species <- fct_relevel(rmpredictions$species, "Aphis nasturtii", after = 1)
rmpredictions$species <- fct_relevel(rmpredictions$species, "Muscidifurax zaraptor", after = 3)



rmestimates <-  ggplot(rmpredictions, aes(temp,estimate,colour=species))+
  scale_x_continuous(expression(plain(paste(" Temperature, ",degree,"C"))))+
  scale_y_continuous(expression(plain(paste(" Population growth rate ("~italic(r)[m]~")"))),
                      limits=c(-0.001,0.255),
                      expand = c(0.01, 0),
                      breaks=seq(0,0.25, by=0.05))+
                      theme_bw(base_size = 12)+
  geom_line(size=0.85)+
  geom_hline(aes(yintercept = 0), linetype = 2, show.legend = FALSE)+
  scale_colour_manual(values = c("#d9d9d9","#bdbdbd","#969696","#737373","#525252","#252525"),
                    name=expression(bold("species")),
                    labels = c("A. grandis","A. nasturtii","P. marginatu",
                               "M. zaraptor","T. nigriabdominalis","R. maidis"),
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

ggsave("~/Dropbox/ph_thesis/Topt_paper/results/rmestimates.pdf",
       rmestimates, width = 10, height = 10, units = "cm",device = cairo_pdf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# code from here is unfinished. I prepared this to address Samraat's prior request for:

# "for each of the two mismatch types, 
# can you please plot T_opt of r_m vs mismatch level 
# as well as niche width vs mismatch 
# in addtion to  r_opt vs mismatch level, 
# thereby generating a two 4-panel figure 
# (in each, the first figure being of TPCs colored by mismatch level). 
# So basically, emulating the theory figure."


# dataset for ropt vs mismatch level analysis 

px  <- dplyr::arrange(pmaz, desc(rmax))
ppk <- data.frame(px[1,])

ax  <- dplyr::arrange(anaz, desc(rmax))
apk <- data.frame(ax[1,])

mx  <- dplyr::arrange(mhaz, desc(rmax))
mpk <- data.frame(mx[1,])

rx  <- dplyr::arrange(rhaz, desc(rmax))
rpk <- data.frame(rx[1,])

agx  <- dplyr::arrange(agaz, desc(rmax))
agpk <- data.frame(agx[1,])

tx  <- dplyr::arrange(tenaz, desc(rmax))
tpk <- data.frame(tx[1,])

ropt_mzmtch <- as_tibble(rbind(tpk,agpk,rpk,mpk,apk,ppk))

ropt_mzmtch <- ropt_mzmtch %>% select(temp,species,rmax,rmax_lwr,rmax_upr)

ropt_mzmtch <- arrange(ropt_mzmtch, desc(rmax))

# compile temp mismatch datasets

alpha  <- as_tibble(read.csv('alpha_Tpks.csv', header = TRUE))
zj     <- as_tibble(read.csv('zj_Tpks.csv', header = TRUE))
z      <- as_tibble(read.csv('z_Tpks.csv', header = TRUE))
bpk    <- as_tibble(read.csv('bpk_Tpks.csv', header = TRUE))


# alpha minus zj








# zj minus z

zj$species <- as.factor(zj$species)
zj$trait   <- as.factor(zj$trait)

z$species <- as.factor(z$species)
z$trait   <- as.factor(z$trait)

zjz <- zj  %>% filter(trait == "juvenile mortality rate")
zjz <- zjz %>% filter(param == "topt")
zjz <- zjz %>% filter(species!= "Stethorus punctillum" # remove species that didn't have all required traits for rm calculations
                      & species!= "Aphis gossypii" 
                      & species!= "Aedes krombeini"
                      & species!= "Telenomus isis"
                      & species!= "Aedes aegypti"
                      & species!= "Aedes albopictus"
                      & species!= "Bemisia tabaci"
                      & species!= "Anopheles gambiae s.s."
                      & species!= "Tetranychus mcdanieli")

zjz <- zjz %>% rename(zj=estimate,zj_lwr=conf_lower,zj_upr=conf_upper)

#%%

zz <- z  %>% filter(trait == "adult mortality rate")
zz <- zz %>% filter(param == "topt")
zz <- zz %>% filter(species!= "Stethorus punctillum"
                    & species!= "Aphis gossypii"
                    & species!= "Aedes krombeini"
                    & species!= "Telenomus isis"
                    & species!= "Aedes aegypti"
                    & species!= "Aedes albopictus"
                    & species!= "Bemisia tabaci"
                    & species!= "Anopheles gambiae s.s."
                    & species!= "Tetranychus mcdanieli")

zz <- zz %>% rename(z=estimate,
                    z_lwr=conf_lower,
                    z_upr=conf_upper,specz=species,trtz=trait)

# organise data for plotting 
mortz <- cbind(zjz,zz)
mortz <- select(mortz,trait, zj,zj_lwr,zj_upr,species,trtz,z,z_lwr,z_upr)

# subtract z from zj
mortz$zjminusz     <- mortz$zj-mortz$z
mortz$zjminusz_lwr <- mortz$zj_lwr-mortz$z_lwr
mortz$zjminusz_upr <- mortz$zj_upr-mortz$z_upr


mortz <- arrange(mortz, desc(zjminusz))










