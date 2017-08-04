

#### Hasselman recreation through Leave-One-Out Cross-Validation ####
Hass_new <- assess_reference_loo(alewife, 17, mixsize = 1000)

D <- rbind(alewife, alewife[1:1000,])
D$sample_type[1071:2070] <- "mixture"
Hass_new$repunit <- factor(Hass_new$repunit, levels = unique(Hass_new$repunit))

rho_mcmc <- Hass_new %>%
  dplyr::group_by(iter, repunit) %>%
  dplyr::summarise(true_rho = sum(true_pi),
                   post_mean_rho = sum(post_mean_pi)) %>%
  dplyr::ungroup()

boot_hass <- lapply(1:50, function(rep) {
  print(rep)
  pi_est <- Hass_new$post_mean_pi[Hass_new$iter == rep]
  rho_est <- rho_mcmc$post_mean_rho[rho_mcmc$iter == rep]
  bootstrap_rho(rho_est, pi_est, D, 17)
  }) %>% unlist()

Hass_rho_data <- cbind(rho_mcmc, boot_hass) %>%
  tidyr::gather(key = "method", value = "Estimate", post_mean_rho:boot_hass)
Hass_rho_data$repunit <- rep(c("N New England", "S New England", "Mid-Atlantic"), 100)
Hass_rho_data$method <- rep(c("MCMC", "Bootstrap"), each = 150)
Hass_rho_data$method <- factor(Hass_rho_data$method, levels = unique(Hass_rho_data$method))

g <- ggplot(Hass_rho_data, aes(x = true_rho, y = Estimate, colour = repunit)) +
  geom_point() +
  facet_grid(method ~ repunit) +
  scale_color_brewer(palette = "Set1") +
  geom_abline(intercept = 0, slope = 1)
print(g)

Hass_rho_data$method <- rep(c("MCMC", "PB"), each = 150)
Hass_rho_data$method <- as.factor(Hass_rho_data$method)

Hass_rho_dev <- Hass_rho_data %>%
  mutate(dev = (true_rho - Estimate)^2) %>%
  mutate(prop_bias = (Estimate-true_rho) / true_rho) %>%
  mutate(bias = Estimate-true_rho) %>%
  group_by(repunit, method) %>%
  summarise(MSE = mean(dev), mean_prop_bias = mean(prop_bias), mean_bias = mean(bias))

# the standard mean error of each method and reporting unit
d <- ggplot(data = Hass_rho_dev, aes(x = method, y = MSE, fill = repunit)) +
  geom_bar(stat = "identity") +
  guides(fill = F) +
  theme(axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~repunit)
print(d)

# proportional bias makes error of equal size have larger effects
# on small values than large values
prop.bias <- ggplot(data = Hass_rho_dev, aes(x = method, y = mean_prop_bias, fill = repunit)) +
  geom_bar(stat = "identity") +
  facet_wrap(~repunit) +
  scale_fill_brewer(palette = "Set1") +
  guides(fill = F)
print(prop.bias)

# the unscaled mean bias
m.bias <- ggplot(data = Hass_rho_dev, aes(x = method, y = mean_bias, fill = repunit)) +
  geom_bar(stat = "identity") +
  guides(fill = F) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~repunit)
print(m.bias)

### everything checks out ###

#### Monte Carlo Cross-Validation ####

mc_new <- assess_pb_bias_correction(alewife, 17)

mc_rho_data <- mc_new %>%
  tidyr::gather(key = "method", value = "Estimate", rho_mcmc:rho_pb)
mc_rho_data$repunit <- rep(c("N New England", "S New England", "Mid-Atlantic"), 100)
mc_rho_data$method <- rep(c("MCMC", "Bootstrap"), each = 150)
mc_rho_data$method <- factor(mc_rho_data$method, levels = unique(mc_rho_data$method))

g <- ggplot(mc_rho_data, aes(x = true_rho, y = Estimate, colour = repunit)) +
  geom_point() +
  facet_grid(method ~ repunit) +
  scale_color_brewer(palette = "Set1") +
  geom_abline(intercept = 0, slope = 1)
print(g)

mc_rho_data$method <- rep(c("MCMC", "PB"), each = 150)
mc_rho_data$method <- as.factor(mc_rho_data$method)

mc_rho_dev <- mc_rho_data %>%
  mutate(dev = (true_rho - Estimate)^2) %>%
  mutate(prop_bias = (Estimate-true_rho) / true_rho) %>%
  mutate(bias = Estimate-true_rho) %>%
  group_by(repunit, method) %>%
  summarise(MSE = mean(dev), mean_prop_bias = mean(prop_bias), mean_bias = mean(bias))

# the standard mean error of each method and reporting unit
d <- ggplot(data = mc_rho_dev, aes(x = method, y = MSE, fill = repunit)) +
  geom_bar(stat = "identity") +
  guides(fill = F) +
  theme(axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~repunit)
print(d)

# proportional bias makes error of equal size have larger effects
# on small values than large values
prop.bias <- ggplot(data = mc_rho_dev, aes(x = method, y = mean_prop_bias, fill = repunit)) +
  geom_bar(stat = "identity") +
  facet_wrap(~repunit) +
  scale_fill_brewer(palette = "Set1") +
  guides(fill = F)
print(prop.bias)

# the unscaled mean bias
m.bias <- ggplot(data = mc_rho_dev, aes(x = method, y = mean_bias, fill = repunit)) +
  geom_bar(stat = "identity") +
  guides(fill = F) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~repunit)
print(m.bias)


### Slight increase in MSE, reduction in bias
