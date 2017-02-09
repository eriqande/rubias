library(ggplot2)

load("~/R/rubias/cjfas_data.RData")

# testing the relationship between RU bias and the difference of the true rho and Np/C
# first, in coalescent data; calculate residuals both with (prop_diff) and without (diff)
# expressing as a proportion of the true rho
coal_data$prop_diff <- (coal_data$Estimate - coal_data$true_rho) / coal_data$true_rho
coal_data$diff <- (coal_data$Estimate - coal_data$true_rho)
coal_data$Np_C <- rep(c(2/17, 3/17, 12/17), 150)
coal_data$Np_diff <- coal_data$Np_C - coal_data$true_rho
coal_mod_prop <- lm(prop_diff ~ Np_diff, data = coal_data)
summary(coal_mod_prop)
# Significant relationship between the residual as a proportion of the true value
# and the difference between the true rho and the Np/C expectation (P = 0.00552)
coal_mod <- lm(diff ~ Np_diff, data = coal_data)
summary(coal_mod)
# Marginally significant relationship between the plain residual and the
# difference between the true rho and the Np/C expectation (P = 0.0689)

# Now the same, but with Hasselman alewife data
hass_data <- Hass[[1]]
hass_data$prop_diff <- (hass_data$rho_est - hass_data$true_rho) / hass_data$true_rho
hass_data$diff <- (hass_data$rho_est - hass_data$true_rho)
hass_data$Np_C <- rep(c(6/21, 3/21, 12/21), 150)
hass_data$Np_diff <- hass_data$Np_C - hass_data$true_rho
hass_mod_prop <- lm(prop_diff ~ Np_diff, data = hass_data)
summary(hass_mod_prop)
# Significant relationship between the residual as a proportion of the true value
# and the difference between the true rho and the Np/C expectation (P < 2.2e-16)
hass_mod <- lm(diff ~ Np_diff, data = hass_data)
summary(hass_mod)
# Significant relationship between the plain residual
# and the difference between the true rho and the Np/C expectation (P < 2.2e-16)

cp <- ggplot2::ggplot(coal_data, aes(x = Np_diff, y = prop_diff, colour = repunit)) +
  geom_point() +
  facet_grid(~ repunit) +
  labs(x='Np/C - True Rho', y = 'Proportional Residual') +
  scale_color_brewer(palette = "Set1")
print(cp)

c <- ggplot2::ggplot(coal_data, aes(x = Np_diff, y = diff, colour = repunit)) +
  geom_point() +
  #facet_grid(~ repunit) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x='Np/C - True Rho', y = 'Residual') +
  scale_color_brewer(palette = "Set1")
print(c)

hp <- ggplot2::ggplot(hass_data, aes(x = Np_diff, y = prop_diff, colour = repunit)) +
  geom_point() +
  facet_grid(~ repunit) +
  labs(x='Np/C - True Rho', y = 'Proportional Residual') +
  scale_color_brewer(palette = "Set1")
print(hp)

h <- ggplot2::ggplot(hass_data, aes(x = Np_diff, y = diff, colour = repunit)) +
  geom_point() +
  #facet_grid(~ repunit) +
  #geom_abline(intercept = coal_mod$coefficients[1], slope = coal_mod$coefficients[2]) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x='Np/C - True Rho', y = 'Residual') +
  scale_color_brewer(palette = "Set1")
print(h)
