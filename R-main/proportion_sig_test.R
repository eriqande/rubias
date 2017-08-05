library(ggplot2)
library(dplyr)

load("~/R/rubias/cjfas_data.RData")

# Exploring the relationship between RU bias and the difference of the true rho and Nc/P
### STASTICS IN THE PAPER APPEAR AT THE BOTTOM OF THE DOCUMENT ####

# First, in coalescent data; calculate residuals both with (prop_diff) and without (diff)
# expressing as a proportion of the true rho
coal_rho_data$prop_diff <- (coal_rho_data$Estimate - coal_rho_data$true_rho) / coal_rho_data$true_rho
coal_rho_data$diff <- (coal_rho_data$Estimate - coal_rho_data$true_rho)
coal_rho_data$Np_C <- rep(c(2/17, 3/17, 12/17), 100)
coal_rho_data$Np_diff <- coal_rho_data$Np_C - coal_rho_data$true_rho
# Model based on proportional residuals
coal_mod_prop <- lm(prop_diff ~ Np_diff * method, data = coal_rho_data)
plot(coal_mod_prop)
# Appears non-normal and heteroscedastic, so simple regression probably inappropriate.
summary(coal_mod_prop)
# Significant relationship between the residual as a proportion of the true value
# and the difference between the true rho and the Np/C expectation (P = 5.17e-10).
# Also a signficant effect of method (P = .0234, PB reduces slope of relationship).
# Interaction of the two is signficant (P = 5.34e-05), i.e. PB changes the slope
# of the relationship in the two by shrinking residuals towards 0.

# Try again with plain residuals
coal_mod <- lm(diff ~ Np_diff * method, data = coal_rho_data)
plot(coal_mod)
# Appears to meet assumptions of linear models
summary(coal_mod)
# Highly significant relationship between the plain residual and the
# difference between the true rho and the Np/C expectation (P = <2e-16)
# But no effect of method (expected, since PB has a net zero effect when averaged
# across positively and negatively biased populations). Highly signficant interaction
# (P = <2e-16) confirms this.

#Graphs of the above relationship, split by reporting unit and method
cp <- ggplot2::ggplot(coal_rho_data, aes(x = Np_diff, y = prop_diff, colour = repunit)) +
  geom_point() +
  facet_grid(repunit ~ method) +
  labs(x='Np/C - True Rho', y = 'Proportional Residual') +
  scale_color_brewer(palette = "Set1")
# Clearly not linear, so the proportional bias statistics must be disregarded.

# Now with just plain residuals
c <- ggplot2::ggplot(coal_rho_data, aes(x = Np_diff, y = diff, colour = repunit)) +
  geom_point() +
  facet_grid(repunit ~ method) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm") +
  labs(x='Np/C - True Rho', y = 'Residual') +
  scale_color_brewer(palette = "Set1")

# Similar graph, but not split by reporting unit
cbp <- ggplot2::ggplot(coal_rho_data, aes(x = Np_diff, y = diff)) +
  geom_point(aes(colour = repunit)) +
  facet_grid(~ method) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", colour = "darkgrey", level = 0) +
  labs(x='Np/C - True Rho', y = 'Residual') +
  scale_color_brewer(palette = "Set1")

# Now the same, but with Hasselman alewife data
Hass_rho_data$prop_diff <- (Hass_rho_data$Estimate - Hass_rho_data$true_rho) / Hass_rho_data$true_rho
Hass_rho_data$diff <- (Hass_rho_data$Estimate - Hass_rho_data$true_rho)
Hass_rho_data$Np_C <- rep(c(6/21, 3/21, 12/21), 100)
Hass_rho_data$Np_diff <- Hass_rho_data$Np_C - Hass_rho_data$true_rho

hass_mod_prop <- lm(prop_diff ~ Np_diff * method, data = Hass_rho_data)
plot(hass_mod_prop)
# very non-normal, assumptions violated
summary(hass_mod_prop)
# No signficant relationships, should ignore anyways

hass_mod <- lm(diff ~ Np_diff * method, data = Hass_rho_data)
plot(hass_mod) # Looks better, although not as normal as coalescent
summary(hass_mod)
# Significant relationship between the plain residual
# and the difference between the true rho and the Np/C expectation (P = 0.0417)
# No signficant interaction (consistent with bias correction less successful than
# in coalescent simulation)

hp <- ggplot2::ggplot(Hass_rho_data, aes(x = Np_diff, y = prop_diff, colour = repunit)) +
  geom_point() +
  facet_grid(repunit ~ method) +
  labs(x='Np/C - True Rho', y = 'Proportional Residual') +
  scale_color_brewer(palette = "Set1")
# Yep, again, nonlinear, and so regression is inappropriate

h <- ggplot2::ggplot(Hass_rho_data, aes(x = Np_diff, y = diff, colour = repunit)) +
  geom_point() +
  facet_grid(repunit ~ method) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm") +
  labs(x='Np/C - True Rho', y = 'Residual') +
  scale_color_brewer(palette = "Set1")


# Similar graph, but not split by reporting unit
hbp <- ggplot2::ggplot(Hass_rho_data, aes(x = Np_diff, y = diff)) +
  geom_point(aes(colour = repunit)) +
  facet_grid(~ method) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", colour = "darkgrey", level = 0) +
  labs(x='Np/C - True Rho', y = 'Residual') +
  scale_color_brewer(palette = "Set1")

# Separate regressions for coalescent data of each method
##### ACTUAL STATISTICS CITED IN PAPER  #####
coal_pb <- coal_rho_data %>% filter(method == 'PB')
coal_pb_mod <- lm(diff ~ Np_diff, data = coal_pb)
summary(coal_pb_mod)
coal_mcmc <- coal_rho_data %>% filter(method == 'MCMC')
coal_mcmc_mod <- lm(diff ~ Np_diff, data = coal_mcmc)
summary(coal_mcmc_mod)
