library(ggplot2)
library(dplyr)

load("~/R/rubias/cjfas_data.RData")

# Remove the data from the failed BH correction
coal_data <- coal_data %>% filter(method != "BH")
hass_data <- Hass[[1]] %>% filter(method != "bh_rho")

# testing the relationship between RU bias and the difference of the true rho and Nc/P
# first, in coalescent data; calculate residuals both with (prop_diff) and without (diff)
# expressing as a proportion of the true rho
coal_data$prop_diff <- (coal_data$Estimate - coal_data$true_rho) / coal_data$true_rho
coal_data$diff <- (coal_data$Estimate - coal_data$true_rho)
coal_data$Np_C <- rep(c(2/17, 3/17, 12/17), 100)
coal_data$Np_diff <- coal_data$Np_C - coal_data$true_rho
coal_mod_prop <- lm(prop_diff ~ Np_diff * method, data = coal_data)
summary(coal_mod_prop)
# Significant relationship between the residual as a proportion of the true value
# and the difference between the true rho and the Np/C expectation (P = 0.00552)
coal_mod <- lm(diff ~ Np_diff * method, data = coal_data)
summary(coal_mod)
# Marginally significant relationship between the plain residual and the
# difference between the true rho and the Np/C expectation (P = 0.0689)

# Now the same, but with Hasselman alewife data
hass_data$prop_diff <- (hass_data$rho_est - hass_data$true_rho) / hass_data$true_rho
hass_data$diff <- (hass_data$rho_est - hass_data$true_rho)
hass_data$Np_C <- rep(c(6/21, 3/21, 12/21), 100)
hass_data$Np_diff <- hass_data$Np_C - hass_data$true_rho

hass_mod_prop <- lm(prop_diff ~ Np_diff * method, data = hass_data)
summary(hass_mod_prop)
# Significant relationship between the residual as a proportion of the true value
# and the difference between the true rho and the Np/C expectation (P < 2.2e-16)

hass_mod <- lm(diff ~ Np_diff * method, data = hass_data)
summary(hass_mod)
# Significant relationship between the plain residual
# and the difference between the true rho and the Np/C expectation (P < 2.2e-16)

#Graphs of the above relationship, split by reporting unit and method
cp <- ggplot2::ggplot(coal_data, aes(x = Np_diff, y = prop_diff, colour = repunit)) +
  geom_point() +
  facet_grid(repunit ~ method) +
  labs(x='Np/C - True Rho', y = 'Proportional Residual') +
  scale_color_brewer(palette = "Set1")
print(cp)

c <- ggplot2::ggplot(coal_data, aes(x = Np_diff, y = diff, colour = repunit)) +
  geom_point() +
  facet_grid(repunit ~ method) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm") +
  labs(x='Np/C - True Rho', y = 'Residual') +
  scale_color_brewer(palette = "Set1")
print(c)

hp <- ggplot2::ggplot(hass_data, aes(x = Np_diff, y = prop_diff, colour = repunit)) +
  geom_point() +
  facet_grid(repunit ~ method) +
  labs(x='Np/C - True Rho', y = 'Proportional Residual') +
  scale_color_brewer(palette = "Set1")
print(hp)


h <- ggplot2::ggplot(hass_data, aes(x = Np_diff, y = diff, colour = repunit)) +
  geom_point() +
  facet_grid(repunit ~ method) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm") +
  labs(x='Np/C - True Rho', y = 'Residual') +
  scale_color_brewer(palette = "Set1")
print(h)

# Similar graphs, but not split by reporting unit
cbp <- ggplot2::ggplot(coal_data, aes(x = Np_diff, y = diff)) +
  geom_point(aes(colour = repunit)) +
  facet_grid(~ method) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", colour = "darkgrey", level = 0) +
  labs(x='Np/C - True Rho', y = 'Residual') +
  scale_color_brewer(palette = "Set1")
print(cbp)

hbp <- ggplot2::ggplot(hass_data, aes(x = Np_diff, y = diff)) +
  geom_point(aes(colour = repunit)) +
  facet_grid(~ method) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", colour = "darkgrey", level = 0) +
  labs(x='Np/C - True Rho', y = 'Residual') +
  scale_color_brewer(palette = "Set1")
print(hbp)

coal_pb <- coal_data %>% filter(method == 'PB')
coal_pb_mod <- lm(diff ~ Np_diff, data = coal_pb)
summary(coal_pb_mod)
coal_mcmc <- coal_data %>% filter(method == 'MCMC')
coal_mcmc_mod <- lm(diff ~ Np_diff, data = coal_mcmc)
summary(coal_mcmc_mod)
