



#______________________________________
#  Tetraneura nigriabdominalis

Tni <- z %>% filter(species == 'Tetraneura nigriabdominalis')

# fit chosen model formulation in rTPC
start_vals <- get_start_vals(Tni$temp, Tni$rate, model_name = 'pawar_2018')

b_fits <- nest(Tni, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 1,
                                          start_upper = start_vals + 1,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________
b_preds <- mutate(b_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  select(., -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  select(model_name, preds) %>%
  unnest(preds)

# plot panel
ggplot(b_preds) +
  geom_line(aes(temp, 1/.fitted, col = model_name)) +
  geom_point(aes(temp, 1/rate),size=0.2,alpha=0.5, Tni)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________
fit_nlsLM3 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = Tni,
                                start = coef(b_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(Tni)))

# bootstrap using residual resampling
boot3 <- Boot(fit_nlsLM3, method = 'residual')

boot3_preds <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(Tni$temp), max(Tni$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot3_conf_preds <- group_by(boot3_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs
ggplot() +
  geom_line(aes(temp, 1/.fitted), b_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = 1/conf_lower, ymax = 1/conf_upper), 
              boot3_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, 1/rate), Tni, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/z')

# calculate params with CIs
extra_params3 <- calc_params(fit_nlsLM3) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')%>%
  mutate(trait="z")

ci_extra_params3 <- Boot(fit_nlsLM3, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM3)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params3 <- left_join(ci_extra_params3, extra_params3)

##Get other parameters
params <- broom::tidy(fit_nlsLM3) %>% select(param = term, estimate)
BootOut <- Boot(fit_nlsLM3, method = 'residual')
## Get the param Names that has multiple values:
paramName <- colnames(BootOut[[2]])[which(apply(BootOut[[2]], 2, function(x) length(unique(x))>1))]

params_cis <- BootOut %>%
  confint(.,parm = paramName, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

params_cis <- bind_rows(params_cis) %>%
  left_join(., params) %>% mutate(trait = 'z') %>% filter(param != 'topt')


TniParams     <- as_tibble(rbind(ci_extra_params3,params_cis)) %>%
  mutate(species = 'Tetraneura nigriabdominalis')

TniZetaFits <- b_preds %>% select(.fitted) %>% bind_cols(boot3_conf_preds) %>% 
  rename(z = .fitted, zLwr = conf_lower, zUpr = conf_upper) %>%
  mutate(z = 1/z, zLwr = 1/zLwr, zUpr = 1/zUpr)

#______________________________________
#  Acyrthosiphon pisum

Api <- z %>% filter(species == 'Acyrthosiphon pisum')

# fit chosen model formulation in rTPC
start_vals <- get_start_vals(Api$temp, Api$rate, model_name = 'pawar_2018')

b_fits <- nest(Api, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 1,
                                          start_upper = start_vals + 1,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________
b_preds <- mutate(b_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  select(., -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  select(model_name, preds) %>%
  unnest(preds)

# plot panel
ggplot(b_preds) +
  geom_line(aes(temp, 1/.fitted, col = model_name)) +
  geom_point(aes(temp, 1/rate),size=0.2,alpha=0.5, Api)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________
fit_nlsLM3 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = Api,
                                start = coef(b_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(Api)))

# bootstrap using residual resampling
boot3 <- Boot(fit_nlsLM3, method = 'residual')

boot3_preds <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(Api$temp), max(Api$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot3_conf_preds <- group_by(boot3_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs
ggplot() +
  geom_line(aes(temp, 1/.fitted), b_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = 1/conf_lower, ymax = 1/conf_upper), 
              boot3_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, 1/rate), Api, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/z')

# calculate params with CIs
extra_params3 <- calc_params(fit_nlsLM3) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')%>%
  mutate(trait="z")

ci_extra_params3 <- Boot(fit_nlsLM3, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM3)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params3 <- left_join(ci_extra_params3, extra_params3)

## Get other parameters
params <- broom::tidy(fit_nlsLM3) %>% select(param = term, estimate)
BootOut <- Boot(fit_nlsLM3, method = 'residual')

## Get the param Names that has multiple values:
paramName <- colnames(BootOut[[2]])[which(apply(BootOut[[2]], 2, function(x) length(unique(x))>1))]

params_cis <- BootOut %>%
  confint(.,parm = paramName, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

params_cis <- bind_rows(params_cis) %>%
  left_join(., params) %>% mutate(trait = 'z') %>% filter(param != 'topt')

ApiParams   <- as_tibble(rbind(ci_extra_params3,params_cis)) %>%
  mutate(species = 'Acyrthosiphon pisum')

ApiZetaFits <-  b_preds %>% select(.fitted) %>% bind_cols(boot3_conf_preds) %>% 
  rename(z = .fitted, zLwr = conf_lower, zUpr = conf_upper) %>%
  mutate(z = 1/z, zLwr = 1/zLwr, zUpr = 1/zUpr)

#______________________________________
#  Stethorus punctillum

Spu <- z %>% filter(species == 'Stethorus punctillum') # %>% 
# mutate_at(vars(c(rate)), 
#~ifelse(rate == 10, 50, .)) %>% mutate(rate = as.numeric(rate))


# fit chosen model formulation in rTPC
start_vals <- get_start_vals(Spu$temp, Spu$rate, model_name = 'pawar_2018')

b_fits <- nest(Spu, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 1,
                                          start_upper = start_vals + 1,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________
b_preds <- mutate(b_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  select(., -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  select(model_name, preds) %>%
  unnest(preds)

# plot panel
ggplot(b_preds) +
  geom_line(aes(temp, 1/.fitted, col = model_name)) +
  geom_point(aes(temp, 1/rate),size=0.2,alpha=0.5, Spu)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________
fit_nlsLM3 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = Spu,
                                start = coef(b_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(Spu)))

# bootstrap using residual resampling
boot3 <- Boot(fit_nlsLM3, method = 'residual')

boot3_preds <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(Spu$temp), max(Spu$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot3_conf_preds <- group_by(boot3_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


boot3_conf_preds <- boot3_conf_preds %>% 
  mutate(conf_lower = 1/conf_lower, conf_upper = 1/conf_upper) #%>%
#filter(temp > 14.32324) %>% 
#mutate_at(vars(c(conf_lower)), 
#          ~ifelse(conf_lower > 0.1, 0.1, .))


# plot bootstrapped CIs
ggplot() +
  geom_line(aes(temp, 1/.fitted), b_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), 
              boot3_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, 1/rate), Spu, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/z')

# calculate params with CIs
extra_params3 <- calc_params(fit_nlsLM3) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')%>%
  mutate(trait="z")

ci_extra_params3 <- Boot(fit_nlsLM3, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM3)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params3 <- left_join(ci_extra_params3, extra_params3)

##Get other parameters
params <- broom::tidy(fit_nlsLM3) %>% select(param = term, estimate)
BootOut <- Boot(fit_nlsLM3, method = 'residual')
## Get the param Names that has multiple values:
paramName <- colnames(BootOut[[2]])[which(apply(BootOut[[2]], 2, function(x) length(unique(x))>1))]

params_cis <- BootOut %>%
  confint(.,parm = paramName, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

params_cis <- bind_rows(params_cis) %>%
  left_join(., params) %>% mutate(trait = 'z') %>% filter(param != 'topt')

SpuParams   <- as_tibble(rbind(ci_extra_params3,params_cis)) %>%
  mutate(species = 'Stethorus punctillum')

SpuZetaFits <-  b_preds %>% select(.fitted) %>% bind_cols(boot3_conf_preds) %>% 
  rename(z = .fitted, zLwr = conf_lower, zUpr = conf_upper) %>%
  mutate(z = 1/z, zLwr = 1/zLwr, zUpr = 1/zUpr)

#______________________________________
#  Planococcus citri

Pci <- z %>% filter(species == 'Planococcus citri') %>% filter(temp != 35)

# fit chosen model formulation in rTPC
start_vals <- get_start_vals(Pci$temp, Pci$rate, model_name = 'pawar_2018')

b_fits <- nest(Pci, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 1,
                                          start_upper = start_vals + 1,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________
b_preds <- mutate(b_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  select(., -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot panel
ggplot(b_preds) +
  geom_line(aes(temp, 1/.fitted, col = model_name)) +
  geom_point(aes(temp, 1/rate),size=0.2,alpha=0.5, Pci)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________
fit_nlsLM3 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = Pci,
                                start = coef(b_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(Pci)))

# bootstrap using residual resampling
boot3 <- Boot(fit_nlsLM3, method = 'residual')

boot3_preds <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(Pci$temp), max(Pci$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot3_conf_preds <- group_by(boot3_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs
ggplot() +
  geom_line(aes(temp, 1/.fitted), b_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = 1/conf_lower, ymax = 1/conf_upper), 
              boot3_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, 1/rate), Pci, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/z')

# calculate params with CIs
extra_params3 <- calc_params(fit_nlsLM3) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')%>%
  mutate(trait="z")

ci_extra_params3 <- Boot(fit_nlsLM3, f = function(x){unlist(calc_params(x))}, 
                         labels = names(calc_params(fit_nlsLM3)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params3 <- left_join(ci_extra_params3, extra_params3)

##Get other parameters
params <- broom::tidy(fit_nlsLM3) %>% select(param = term, estimate)
BootOut <- Boot(fit_nlsLM3, method = 'residual')
## Get the param Names that has multiple values:
paramName <- colnames(BootOut[[2]])[which(apply(BootOut[[2]], 2, function(x) length(unique(x))>1))]

params_cis <- BootOut %>%
  confint(.,parm = paramName, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

params_cis <- bind_rows(params_cis) %>%
  left_join(., params) %>% mutate(trait = 'z') %>% filter(param != 'topt')

PciParams   <- as_tibble(rbind(ci_extra_params3,params_cis)) %>%
  mutate(species = 'Planococcus citri')

PciZetaFits <-  b_preds %>% select(.fitted) %>% bind_cols(boot3_conf_preds) %>% 
  rename(z = .fitted, zLwr = conf_lower, zUpr = conf_upper) %>%
  mutate(z = 1/z, zLwr = 1/zLwr, zUpr = 1/zUpr)

#______________________________________
#  Anthonomus grandis

Agr <- z %>% filter(species == 'Anthonomus grandis')

# fit chosen model formulation in rTPC
start_vals <- get_start_vals(Agr$temp, Agr$rate, model_name = 'pawar_2018')

b_fits <- nest(Agr, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 1,
                                          start_upper = start_vals + 1,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________
b_preds <- mutate(b_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  select(., -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot panel
ggplot(b_preds) +
  geom_line(aes(temp, 1/.fitted, col = model_name)) +
  geom_point(aes(temp, 1/rate),size=0.2,alpha=0.5, Agr)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________
fit_nlsLM3 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = Agr,
                                start = coef(b_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(Agr)))

# bootstrap using residual resampling
boot3 <- Boot(fit_nlsLM3, method = 'residual')

boot3_preds <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(Agr$temp), max(Agr$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot3_conf_preds <- group_by(boot3_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs
ggplot() +
  geom_line(aes(temp, 1/.fitted), b_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = 1/conf_lower, ymax = 1/conf_upper), 
              boot3_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, 1/rate), Agr, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/z')

# calculate params with CIs
extra_params3 <- calc_params(fit_nlsLM3) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')%>%
  mutate(trait="z")

ci_extra_params3 <- Boot(fit_nlsLM3, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM3)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params3 <- left_join(ci_extra_params3, extra_params3)

##Get other parameters
params <- broom::tidy(fit_nlsLM3) %>% select(param = term, estimate)
BootOut <- Boot(fit_nlsLM3, method = 'residual')
## Get the param Names that has multiple values:
paramName <- colnames(BootOut[[2]])[which(apply(BootOut[[2]], 2, function(x) length(unique(x))>1))]

params_cis <- BootOut %>%
  confint(.,parm = paramName, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

params_cis <- bind_rows(params_cis) %>%
  left_join(., params) %>% mutate(trait = 'z') %>% filter(param != 'topt')

AgrParams   <- as_tibble(rbind(ci_extra_params3,params_cis)) %>%
  mutate(species = 'Anthonomous grandis')

AgrZetaFits <-  b_preds %>% select(.fitted) %>% bind_cols(boot3_conf_preds) %>% 
  rename(z = .fitted, zLwr = conf_lower, zUpr = conf_upper) %>%
  mutate(z = 1/z, zLwr = 1/zLwr, zUpr = 1/zUpr)


#______________________________________
#  Tetranychus mcdanieli

Tmc <- z %>% filter(species == 'Tetranychus mcdanieli')

# fit chosen model formulation in rTPC
start_vals <- get_start_vals(Tmc$temp, Tmc$rate, model_name = 'pawar_2018')

b_fits <- nest(Tmc, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 1,
                                          start_upper = start_vals + 1,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________
b_preds <- mutate(b_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  select(., -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot panel
ggplot(b_preds) +
  geom_line(aes(temp, 1/.fitted, col = model_name)) +
  geom_point(aes(temp, 1/rate),size=0.2,alpha=0.5, Tmc)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________
fit_nlsLM3 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = Tmc,
                                start = coef(b_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(Tmc)))

# bootstrap using residual resampling
boot3 <- Boot(fit_nlsLM3, method = 'residual')

boot3_preds <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(Tmc$temp), max(Tmc$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot3_conf_preds <- group_by(boot3_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()


# plot bootstrapped CIs
ggplot() +
  geom_line(aes(temp, 1/.fitted), b_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = 1/conf_lower, ymax = 1/conf_upper), 
              boot3_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, 1/rate), Tmc, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/z')

# calculate params with CIs
extra_params3 <- calc_params(fit_nlsLM3) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')%>%
  mutate(trait="z")

ci_extra_params3 <- Boot(fit_nlsLM3, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM3)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params3 <- left_join(ci_extra_params3, extra_params3)

##Get other parameters
params <- broom::tidy(fit_nlsLM3) %>% select(param = term, estimate)
BootOut <- Boot(fit_nlsLM3, method = 'residual')
## Get the param Names that has multiple values:
paramName <- colnames(BootOut[[2]])[which(apply(BootOut[[2]], 2, function(x) length(unique(x))>1))]

params_cis <- BootOut %>%
  confint(.,parm = paramName, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

params_cis <- bind_rows(params_cis) %>%
  left_join(., params) %>% mutate(trait = 'z') %>% filter(param != 'topt')

TmcParams   <- as_tibble(rbind(ci_extra_params3,params_cis)) %>%
  mutate(species = 'Tetranychus mcdanieli')

TmcZetaFits <-  b_preds %>% select(.fitted) %>% bind_cols(boot3_conf_preds) %>% 
  rename(z = .fitted, zLwr = conf_lower, zUpr = conf_upper) %>%
  mutate(z = 1/z, zLwr = 1/zLwr, zUpr = 1/zUpr)

#______________________________________
#  Aphis nasturtii

Ana <- z %>% filter(species == 'Aphis nasturtii')

# fit chosen model formulation in rTPC
start_vals <- get_start_vals(Ana$temp, Ana$rate, model_name = 'pawar_2018')

b_fits <- nest(Ana, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 1,
                                          start_upper = start_vals + 1,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________
b_preds <- mutate(b_fits, new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))))%>% 
  select(., -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  select(model_name, preds) %>%
  # unlist the preds list column
  unnest(preds)

# plot panel
ggplot(b_preds) +
  geom_line(aes(temp, 1/.fitted, col = model_name)) +
  geom_point(aes(temp, 1/rate),size=0.2,alpha=0.5, Ana)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________
fit_nlsLM3 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = Ana,
                                start = coef(b_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(Ana)))

# bootstrap using residual resampling
boot3 <- Boot(fit_nlsLM3, method = 'residual')

boot3_preds <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(Ana$temp), max(Ana$temp), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = pawar_2018(temp, r_tref, e, eh, topt, tref = 15))

# calculate bootstrapped confidence intervals
boot3_conf_preds <- group_by(boot3_preds, temp) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup()

# plot bootstrapped CIs
ggplot() +
  geom_line(aes(temp, 1/.fitted), b_preds, col = 'blue') +
  geom_ribbon(aes(temp, ymin = 1/conf_lower, ymax = 1/conf_upper), 
              boot3_conf_preds, fill = 'blue', alpha = 0.3) +
  geom_point(aes(temp, 1/rate), Ana, size = 2, alpha = 0.5) +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (C)',
       y = '1/z')

# calculate params with CIs
extra_params3 <- calc_params(fit_nlsLM3) %>%
  pivot_longer(everything(), names_to =  'param', values_to = 'estimate')%>%
  mutate(trait="z")

ci_extra_params3 <- Boot(fit_nlsLM3, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM3)), R = 200, method = 'residual') %>%
  confint(., method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

# params
ci_extra_params3 <- left_join(ci_extra_params3, extra_params3)

##Get other parameters
params <- broom::tidy(fit_nlsLM3) %>% select(param = term, estimate)
BootOut <- Boot(fit_nlsLM3, method = 'residual')
## Get the param Names that has multiple values:
paramName <- colnames(BootOut[[2]])[which(apply(BootOut[[2]], 2, function(x) length(unique(x))>1))]

params_cis <- BootOut %>%
  confint(.,parm = paramName, method = 'bca') %>%
  as.data.frame() %>%
  rename(conf_lower = 1, conf_upper = 2) %>%
  rownames_to_column(., var = 'param') %>%
  mutate(method = 'residual bootstrap')

params_cis <- bind_rows(params_cis) %>%
  left_join(., params) %>% mutate(trait = 'z') %>% filter(param != 'topt')

AnaParams   <- as_tibble(rbind(ci_extra_params3,params_cis)) %>%
  mutate(species = 'Tetranychus mcdanieli')

AnaZetaFits <-  b_preds %>% select(.fitted) %>% bind_cols(boot3_conf_preds) %>% 
  rename(z = .fitted, zLwr = conf_lower, zUpr = conf_upper) %>%
  mutate(z = 1/z, zLwr = 1/zLwr, zUpr = 1/zUpr)