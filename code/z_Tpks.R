## load packages

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

#read in and prepare the trait data 

df <- as_tibble(read.csv('mismatch_dat.csv'))

df<- df %>%
  rename(temp = ambienttemp,species = interactor1)

dv <- df %>% filter(standardisedtraitname == "Female Mortality Rate")

dv <- dv %>% select(curve_ID, species, temp, standardisedtraitname, standardisedtraitvalue, stage)

dv <- dv %>% rename(rate = standardisedtraitvalue)

dv$rate <- 1/dv$rate

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

# Bootstrap by species to estimate CIs for Tpks and Bpks

#€€€€€ Aedes aegypti

a.e <- dv %>% filter(curve_ID == 1)

A_e <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                               data = a.e,
                               start = coef(dv_fits$fit[[1]]),
                               weights = rep(1, times = nrow(a.e)))


extra_params1 <- calc_params(A_e) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params1 <- Boot(A_e, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(A_e)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params1 <- left_join(ci_extra_params1, extra_params1)

topt1 <- as_tibble(rbind(ci_extra_params1[1,],ci_extra_params1[2,]))
topt1$species <- as.character("Aedes aegypti")

#€€€€€ Aedes albopictus

a.b <- dv %>% filter(curve_ID == 2)

A_b <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                         data = a.b,
                         start = coef(dv_fits$fit[[2]]),
                         weights = rep(1, times = nrow(a.b)))


extra_params2 <- calc_params(A_b) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params2 <- Boot(A_b, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(A_b)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params2 <- left_join(ci_extra_params2, extra_params2)

topt2 <- as_tibble(rbind(ci_extra_params2[1,],ci_extra_params2[2,]))
topt2$species <- as.character("Aedes albopictus")


#€€€€€ Aedes krobeini

a.k <- dv %>% filter(curve_ID == 3)

A_k <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                         data = a.k,
                         start = coef(dv_fits$fit[[3]]),
                         lower = get_lower_lims(a.k$temp, a.k$rate, model_name = 'pawar_2018'),
                         upper = get_upper_lims(a.k$temp, a.k$rate, model_name = 'pawar_2018'),
                         weights = rep(1, times = nrow(a.k)))


extra_params3 <- calc_params(A_k) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params3 <- Boot(A_k, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(A_k)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params3 <- left_join(ci_extra_params3, extra_params3)

topt3 <- as_tibble(rbind(ci_extra_params3[1,],ci_extra_params3[2,]))
topt3$species <- as.character("Aedes krombeini")


#€€€€€ Anopheles gambiae

a.g <- dv %>% filter(curve_ID == 4)

A_g <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                         data = a.g,
                         start = coef(dv_fits$fit[[4]]),
                         weights = rep(1, times = nrow(a.g)))


extra_params4 <- calc_params(A_g) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params4 <- Boot(A_g, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(A_g)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params4 <- left_join(ci_extra_params4, extra_params4)

topt4 <- as_tibble(rbind(ci_extra_params4[1,],ci_extra_params4[2,]))
topt4$species <- as.character("Anopheles gambiae s.s.")

#€€€€€ Anthonomus grandis

an.g <- dv %>% filter(curve_ID == 5)

An_g <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                         data = an.g,
                         start = coef(dv_fits$fit[[5]]),
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

#€€€€€ Aphis nasturtii    

ap.n <- dv %>% filter(curve_ID == 7)

Ap_n <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                          data = ap.n,
                          start = coef(dv_fits$fit[[6]]),
                          lower = get_lower_lims(ap.n$temp, ap.n$rate, model_name = 'pawar_2018'),
                          upper = get_upper_lims(ap.n$temp, ap.n$rate, model_name = 'pawar_2018'),
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

#€€€€€ Bemisia tabaci

b.t <- dv %>% filter(curve_ID == 8)

B_t <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                          data = b.t,
                          start = coef(dv_fits$fit[[7]]),
                          weights = rep(1, times = nrow(b.t)))


extra_params8 <- calc_params(B_t) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params8 <- Boot(B_t, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(B_t)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params8 <- left_join(ci_extra_params8, extra_params8)

topt8 <- as_tibble(rbind(ci_extra_params8[1,],ci_extra_params8[2,]))
topt8$species <- as.character("Bemisia tabaci")

#€€€€€ Muscidifurax zaraptor

m.z <- dv %>% filter(curve_ID == 12)

M_z <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                         data = m.z,
                         start = coef(dv_fits$fit[[8]]),
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

p.m <- dv %>% filter(curve_ID == 13)

P_m <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                         data = p.m,
                         start = coef(dv_fits$fit[[9]]),
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

#€€€€€ Rhopalosiphum maidis

r.m <- dv %>% filter(curve_ID == 15)

R_m <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                         data = r.m,
                         start = coef(dv_fits$fit[[10]]),
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

#€€€€€ Tetraneura nigri abdominalis

t.n <- dv %>% filter(curve_ID == 16)

T_n <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                         data = t.n,
                         start = coef(dv_fits$fit[[11]]),
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

#€€€€€ Tetranychus mcdanieli

t.m <- dv %>% filter(curve_ID == 17)

T_m <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                         data = t.m,
                         start = coef(dv_fits$fit[[12]]),
                         weights = rep(1, times = nrow(t.m)))


extra_params19 <- calc_params(T_m) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')

ci_extra_params19 <- Boot(T_m, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(T_m)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

ci_extra_params19 <- left_join(ci_extra_params19, extra_params19)

topt19 <- as_tibble(rbind(ci_extra_params19[1,],ci_extra_params19[2,]))
topt19$species <- as.character("Tetranychus mcdanieli")

#%%%%%%%%%%% bind z Tpks and Bpks

dv_tz <- rbind(topt1,topt2,topt3,topt4,
               topt5,topt7,topt8,topt12,
               topt13,topt15,topt18,topt19)

dv_tz$trait <- c("adult mortality rate")

dv_topt          <- dv_tz %>% filter(param =="topt")
dv_rm            <- dv_tz %>% filter(param =="rmax")
dv_rm$temp       <- dv_topt$estimate
dv_rm$estimate   <- 1/dv_rm$estimate # invert for rate
dv_rm$conf_lower <- 1/dv_rm$conf_lower
dv_rm$conf_upper <- 1/dv_rm$conf_upper

#%%%%%%%%% plot Tpks

ggplot(dv_topt, aes(estimate, species)) +
  geom_point(size = 4) + 
  geom_linerange(aes(xmin = conf_lower, xmax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') 

#%%%%%%%%% plot Bpks

ggplot(dv_rm, aes(temp, estimate)) +
  geom_point(size = 2.5) + 
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw(base_size = 12) +
  scale_x_continuous('') 


#%%%%%%%%%%

write.csv(dv_tz, "z_Tpks.csv")



