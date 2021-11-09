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
require('doMC')
require('foreach')

rm(list=ls())
graphics.off() 

# setwd("~/Dropbox/ph_thesis/Topt_paper/data")

#read in the trait data 

df <- as_tibble(read.csv('mismatch_dat.csv'))

df<- df %>%
  rename(temp = ambienttemp,species = interactor1)

dv <- df %>% filter(standardisedtraitname == 'Fecundity Rate')

dv <- dv %>% select(curve_ID, species, temp, standardisedtraitname, standardisedtraitvalue, stage)

dv <- dv %>% rename(rate = standardisedtraitvalue)

dv$rate <- dv$rate + 0.001

dv <- dv %>% arrange(dv, curve_ID)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# fit TPC model for each species

start_vals <- get_start_vals(dv$temp, dv$rate, model_name = 'pawar_2018')
low_lims <- get_lower_lims(dv$temp, dv$rate, model_name = 'pawar_2018')
upper_lims <- get_upper_lims(dv$temp, dv$rate, model_name = 'pawar_2018')


dv_fits <- nest(dv, data = c(temp, rate)) %>%
           mutate(fit = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 10,
                                          start_upper = start_vals + 10,
                                          lower = low_lims,
                                          upper=upper_lims,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))


dv_preds <- mutate(dv_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100)))) %>%
  # get rid of original data column
  select(., -data) %>%
  # stack models into a single column, with an id column for model_name
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(fit)) %>%
  # create new list column containing the predictions
  # this uses both fit and new_data list columns
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  # select only the columns we want to keep
  select(curve_ID, species, standardisedtraitname, stage, preds) %>%
  # unlist the preds list column
  unnest(preds)

glimpse(dv_preds)

ggplot(dv_preds) +
  geom_line(aes(temp, .fitted, col = curve_ID)) +
  geom_point(aes(temp, rate), dv) +
  facet_wrap(~species, scales = 'free_y', ncol = 6) +
  theme_bw() +
  theme(legend.position = 'none')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#€€€€€ Anthonomus grandis

an.g <- dv %>% filter(curve_ID == 4)

An_g <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                         data = an.g,
                         start = coef(dv_fits$fit[[1]]),
                         lower = get_lower_lims(an.g$temp, an.g$rate, model_name = 'pawar_2018'),
                         upper = get_upper_lims(an.g$temp, an.g$rate, model_name = 'pawar_2018'),
                         weights = rep(1, times = nrow(an.g)))


extra_params5 <- calc_params(An_g) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params5 <- Boot(An_g, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(An_g)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params5 <- left_join(ci_extra_params5, extra_params5)

topt5 <- as_tibble(rbind(ci_extra_params5[1,],ci_extra_params5[2,]))
topt5$species <- as.character("Anthonomus grandis")

#€€€€€ Aphis gossypii

ap.g <- dv %>% filter(curve_ID == 5)

Ap_g <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                          data = ap.g,
                          start = coef(dv_fits$fit[[2]]),
                          weights = rep(1, times = nrow(ap.g)))


extra_params6 <- calc_params(Ap_g) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params6 <- Boot(Ap_g, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(Ap_g)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params6 <- left_join(ci_extra_params6, extra_params6)

topt6 <- as_tibble(rbind(ci_extra_params6[1,],ci_extra_params6[2,]))
topt6$species <- as.character("Aphis gossypii")

#€€€€€ Aphis nasturtii

ap.n <- dv %>% filter(curve_ID == 6)

Ap_n <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                          data = ap.n,
                          start = coef(dv_fits$fit[[3]]),
                          weights = rep(1, times = nrow(ap.n)))


extra_params7 <- calc_params(Ap_n) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params7 <- Boot(Ap_n, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(Ap_n)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params7 <- left_join(ci_extra_params7, extra_params7)

topt7 <- as_tibble(rbind(ci_extra_params7[1,],ci_extra_params7[2,]))
topt7$species <- as.character("Aphis nasturtii")

#€€€€€ Muscidifurax zaraptor

m.z <- dv %>% filter(curve_ID == 11)

M_z <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                         data = m.z,
                         start = coef(dv_fits$fit[[4]]),
                         lower = get_lower_lims(m.z$temp, m.z$rate, model_name = 'pawar_2018'),
                         upper = get_upper_lims(m.z$temp, m.z$rate, model_name = 'pawar_2018'),
                         weights = rep(1, times = nrow(m.z)))


extra_params12 <- calc_params(M_z) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params12 <- Boot(M_z, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(M_z)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params12 <- left_join(ci_extra_params12, extra_params12)

topt12 <- as_tibble(rbind(ci_extra_params12[1,],ci_extra_params12[2,]))
topt12$species <- as.character("Muscidifurax zaraptor")


#€€€€€ Paracoccus marginatu

p.m <- dv %>% filter(curve_ID == 12)

P_m <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                         data = p.m,
                         start = coef(dv_fits$fit[[5]]),
                         lower = get_lower_lims(p.m$temp, p.m$rate, model_name = 'pawar_2018'),
                         upper = get_upper_lims(p.m$temp, p.m$rate, model_name = 'pawar_2018'),
                         weights = rep(1, times = nrow(p.m)))


extra_params13 <- calc_params(P_m) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params13 <- Boot(P_m, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(P_m)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params13 <- left_join(ci_extra_params13, extra_params13)

topt13         <- as_tibble(rbind(ci_extra_params13[1,],ci_extra_params13[2,]))
topt13$species <- as.character("Paracoccus marginatu")


#€€€€€ Planococcus citri

p.c <- dv %>% filter(curve_ID == 13)
P_c <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                          data = p.c,
                          start = coef(dv_fits$fit[[6]]),
                          weights = rep(1, times = nrow(p.c)))


extra_params14 <- calc_params(P_c) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params14 <- Boot(P_c, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(P_c)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params14 <- left_join(ci_extra_params14, extra_params14)

topt14 <- as_tibble(rbind(ci_extra_params14[1,],ci_extra_params14[2,]))
topt14$species <- as.character("Planococcus citri")

#€€€€€ Rhopalosiphum maidis

r.m <- dv %>% filter(curve_ID == 14)

R_m <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                         data = r.m,
                         start = coef(dv_fits$fit[[7]]),
                         lower = get_lower_lims(r.m$temp, r.m$rate, model_name = 'pawar_2018'),
                         upper = get_upper_lims(r.m$temp, r.m$rate, model_name = 'pawar_2018'),
                         weights = rep(1, times = nrow(r.m)))


extra_params15 <- calc_params(R_m) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params15 <- Boot(R_m, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(R_m)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params15 <- left_join(ci_extra_params15, extra_params15)

topt15 <- as_tibble(rbind(ci_extra_params15[1,],ci_extra_params15[2,]))
topt15$species <- as.character("Rhopalosiphum maidis")


#€€€€€ Tetraneura nigriabdominalis

t.n <- dv %>% filter(curve_ID == 17)

T_n <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                         data = t.n,
                         start = coef(dv_fits$fit[[8]]),
                         lower = get_lower_lims(t.n$temp, t.n$rate, model_name = 'pawar_2018'),
                         upper = get_upper_lims(t.n$temp, t.n$rate, model_name = 'pawar_2018'),
                         weights = rep(1, times = nrow(t.n)))


extra_params18 <- calc_params(T_n) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params18 <- Boot(T_n, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(T_n)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params18 <- left_join(ci_extra_params18, extra_params18)

topt18 <- as_tibble(rbind(ci_extra_params18[1,],ci_extra_params18[2,]))
topt18$species <- as.character("Tetraneura nigriabdominalis")

#¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢¢

dv_tz <- rbind(topt5,topt6,topt7,topt12,topt13,topt14,topt15,topt18)

dv_tz$trait <- c("fecundity")

ggplot(dv_tz, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

write.csv(dv_tz, "bpk_Tpks.csv")



