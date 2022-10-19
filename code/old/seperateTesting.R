


Cto <- as_tibble(read.csv('TraitData.csv')) 

Cto <- Cto %>% rename(temp = interactor1temp, species = interactor1, rate = standardisedtraitvalue) %>%
  select(species, temp, standardisedtraitname, rate) %>%
  filter(standardisedtraitname == 'z', rate != 'NA', 
         species == 'Anoplophora glabripennis') %>%
  mutate(curve_ID = case_when(species == 'Anoplophora glabripennis' ~ '21')) %>% 
  mutate(temp = as.numeric(temp), rate = 1/rate)


# fit chosen model formulation in rTPC
start_vals <- get_start_vals(Cto$temp, Cto$rate, model_name = 'pawar_2018')

b_fits <- nest(Cto, data = c(temp, rate)) %>%
  mutate(pawar = map(data, ~nls_multstart(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                          data = .x,
                                          iter = c(3,3,3,3),
                                          start_lower = start_vals - 1,
                                          start_upper = start_vals + 1,
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)))
#____________________
b_preds <- mutate(b_fits, new_data = map(data, ~tibble(temp = seq(min(0), max(45), length.out = 1000))))%>% 
  select(., -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', c(pawar)) %>%
  mutate(preds = map2(fit, new_data, ~augment(.x, newdata = .y))) %>%
  select(model_name, preds) %>%
  unnest(preds)

# plot panel
ggplot(b_preds) +
  geom_line(aes(temp, 1/.fitted, col = model_name)) +
  geom_point(aes(temp, 1/rate),size=0.2,alpha=0.5, Cto)+
  scale_color_brewer(type = 'qual', palette = 2) +
  theme_bw(base_size = 10)+
  theme(legend.position = "none")

#____________________
fit_nlsLM3 <- minpack.lm::nlsLM(rate~pawar_2018(temp = temp, r_tref,e,eh,topt, tref = 15),
                                data = Cto,
                                start = coef(b_fits$pawar[[1]]),
                                weights = rep(1, times = nrow(Cto)))

# bootstrap using residual resampling
boot3 <- Boot(fit_nlsLM3, method = 'residual')

boot3_preds <- boot3$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(temp = seq(min(0), max(45), length.out = 1000))) %>%
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
  geom_point(aes(temp, 1/rate), Cto, size = 2, alpha = 0.5) +
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

CtoParams   <- as_tibble(rbind(ci_extra_params3,params_cis)) %>%
  mutate(species = 'Clavigralla tomentosicollis')

write_csv(CtoParams, 'CtoParams.csv')

CtoZeta <- boot3_conf_preds %>% select(conf_lower,conf_upper)
CtoZetaFits <- b_preds %>% select(temp,.fitted) %>% bind_cols(CtoZeta) %>%
  mutate(standardisedtraitname = 'z', species = 'Clavigralla tomentosicollis', curve_ID = '20') %>%
  select(curve_ID, species, standardisedtraitname, temp, .fitted, conf_lower, conf_upper)
