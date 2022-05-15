# -*- coding: utf-8 -*-
# ±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±±

# get the exponent for how rm_opt scales with mass or assume 
# rm_opt scales to the quarter power of body mass

# r_m,opt = r_0 * m^-y * e-E/(k T_opt) 

# Therefore, taking logs of both sides,

# log(r_m,opt) = log(r_0) - y log (m) - E (1/kT_opt) 

# r_m,opt / m^-y  = r_0 * e-E/(kT_opt)

# load in data 

Tc <- as_tibble(read.csv('dataForTempCorrections.csv'))

#plot uncorrected data in log-log scale
Tc %>%
  ggplot(aes(x=log(mass),y = log(rm_opt)))+
  geom_point()+
  geom_smooth(method = 'lm')

# extract parameter values using model coefficients
r0_log <- coef(model)[1]
y      <- coef(model)[2]


# linear model (note the allometry is linear in log-log scale)
Tcmodel <- lm(log(rm_opt) ~ log(mass) * 1/(k*T_opt), data = Tc)
summary(Tcmodel); 
coef(Tcmodel); 
confint(Tcmodel,level = .95)
Anova(Tcmodel)

# define mass range to plot over
# M_range <- range(log(rmMass$mass))
# M_log   <- seq(M_range[1],M_range[2],length.out = 100)


# define function to get log(r0) at a given temperature 

# temp_norm_r0 <- function(Temp,Temp_r){
#  R0_norm    <- r0_log - ( (E/k) * ((1/Temp) - (1/Temp_r)) ) 
#  return(R0_norm) 
# }

# ref temperatures
# T_20 <- 293.15; T_22 <- 295.15; 
# T_26 <- 299.15
# T_32 <- 305.15

# calculate  Temp corrections
r0_tc_log <- temp_norm_r0(Temp = T_opt+273.15 , Temp_r = T_26)


# get log rm_opt 
# r_20_log <- r0_20_log + (y * M_log) 
# r_22_log <- r0_22_log + (y * M_log)
# r_26_log <- r0_26_log + (y * M_log)
# r_32_log <- r0_32_log + (y * M_log)


# plot results
