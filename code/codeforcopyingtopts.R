#_______________________________________________
d_fits %>% print(n = Inf)

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

