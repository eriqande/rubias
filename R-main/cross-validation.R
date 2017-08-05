

#### Hasselman recreation through Leave-One-Out Cross-Validation ####

# Perform Leave-One-Out CV-ML with same parameters & dataset as Hasselman et al. 2015
Hass_new <- assess_reference_loo(alewife, 17, mixsize = 1000)

# Reformat, calculate pis, and bootstrap
D <- rbind(alewife, alewife[1:1000,])       # mixture added for bootstrap_rho,
D$sample_type[1071:2070] <- "mixture"       # only num. of mix indiv. is important
Hass_new$repunit <- factor(Hass_new$repunit, levels = unique(Hass_new$repunit))

rho_mcmc <- Hass_new %>%
  dplyr::group_by(iter, repunit) %>%
  dplyr::summarise(true_rho = sum(true_pi),
                   post_mean_rho = sum(post_mean_pi)) %>%
  dplyr::ungroup()

boot_hass <- lapply(1:50, function(rep) {
  message(paste("Hasselman bootstrap repetition", rep))
  pi_est <- Hass_new$post_mean_pi[Hass_new$iter == rep]
  rho_est <- rho_mcmc$post_mean_rho[rho_mcmc$iter == rep]
  bootstrap_rho(rho_est, pi_est, D, 17)
  }) %>% unlist()

# Reformatting for graphing
Hass_rho_data <- cbind(rho_mcmc, boot_hass) %>%
  tidyr::gather(key = "method", value = "Estimate", post_mean_rho:boot_hass)
Hass_rho_data$repunit <- rep(c("N New England", "S New England", "Mid-Atlantic"), 100)
Hass_rho_data$repunit <- factor(Hass_rho_data$repunit, levels = unique(Hass_rho_data$repunit))
Hass_rho_data$method <- rep(c("MCMC", "Bootstrap"), each = 150)
Hass_rho_data$method <- factor(Hass_rho_data$method, levels = unique(Hass_rho_data$method))

# Graph of individual MCMC estimates vs. true 1:1 line
h.lg <- ggplot(Hass_rho_data, aes(x = true_rho, y = Estimate, colour = repunit)) +
  geom_point() +
  facet_grid(method ~ repunit) +
  scale_color_brewer(palette = "Set1") +
  geom_abline(intercept = 0, slope = 1)

Hass_rho_data$method <- rep(c("MCMC", "PB"), each = 150)
Hass_rho_data$method <- as.factor(Hass_rho_data$method)

# Calculating summary statistics
Hass_rho_dev <- Hass_rho_data %>%
  mutate(dev = (true_rho - Estimate)^2) %>%
  mutate(prop_bias = (Estimate-true_rho) / true_rho) %>%
  mutate(bias = Estimate-true_rho) %>%
  group_by(repunit, method) %>%
  summarise(MSE = mean(dev), mean_prop_bias = mean(prop_bias), mean_bias = mean(bias)) %>%
  ungroup()

# Graph the standard mean error of each method and reporting unit
h.mse <- ggplot(data = Hass_rho_dev, aes(x = method, y = MSE, fill = repunit)) +
  geom_bar(stat = "identity") +
  guides(fill = F) +
  theme(axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~repunit)

# Graph proportional bias (makes error of equal size have
# larger effects on small values than large values)
h.prop.bias <- ggplot(data = Hass_rho_dev, aes(x = method, y = mean_prop_bias, fill = repunit)) +
  geom_bar(stat = "identity") +
  facet_wrap(~repunit) +
  scale_fill_brewer(palette = "Set1") +
  guides(fill = F)

# Graph the unscaled mean bias, or mean residual
h.m.bias <- ggplot(data = Hass_rho_dev, aes(x = method, y = mean_bias, fill = repunit)) +
  geom_bar(stat = "identity") +
  guides(fill = F) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~repunit)

### Printing plots shows a decrease in MSE, mean proportional bias, ###
### and mean bias (mean residual) across all reporting units        ###


#### Monte Carlo Cross-Validation, aka sampling mix. from ref. without replacement ####

# Perform Monte Carlo CV with default parameters
mc_new <- assess_pb_bias_correction(alewife, 17)

# Reformatting for graphing
mc_rho_data <- mc_new %>%
  tidyr::gather(key = "method", value = "Estimate", rho_mcmc:rho_pb)
mc_rho_data$repunit <- rep(c("N New England", "S New England", "Mid-Atlantic"), 100)
mc_rho_data$repunit <- factor(mc_rho_data$repunit, levels = unique(mc_rho_data$repunit))
mc_rho_data$method <- rep(c("MCMC", "Bootstrap"), each = 150)
mc_rho_data$method <- factor(mc_rho_data$method, levels = unique(mc_rho_data$method))

# Graph of individual MCMC estimates vs. true 1:1 line
mc.lg <- ggplot(mc_rho_data, aes(x = true_rho, y = Estimate, colour = repunit)) +
  geom_point() +
  facet_grid(method ~ repunit) +
  scale_color_brewer(palette = "Set1") +
  geom_abline(intercept = 0, slope = 1)

mc_rho_data$method <- rep(c("MCMC", "PB"), each = 150)
mc_rho_data$method <- as.factor(mc_rho_data$method)

# Calculating summary statistics
mc_rho_dev <- mc_rho_data %>%
  mutate(dev = (true_rho - Estimate)^2) %>%
  mutate(prop_bias = (Estimate-true_rho) / true_rho) %>%
  mutate(bias = Estimate-true_rho) %>%
  group_by(repunit, method) %>%
  summarise(MSE = mean(dev), mean_prop_bias = mean(prop_bias), mean_bias = mean(bias)) %>%
  ungroup()

# Graph the standard mean error of each method and reporting unit
mc.mse <- ggplot(data = mc_rho_dev, aes(x = method, y = MSE, fill = repunit)) +
  geom_bar(stat = "identity") +
  guides(fill = F) +
  theme(axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~repunit)

# Graph proportional bias (makes error of equal size have
# larger effects on small values than large values)
mc.prop.bias <- ggplot(data = mc_rho_dev, aes(x = method, y = mean_prop_bias, fill = repunit)) +
  geom_bar(stat = "identity") +
  facet_wrap(~repunit) +
  scale_fill_brewer(palette = "Set1") +
  guides(fill = F)

# Graph the unscaled mean bias, or mean residual
mc.m.bias <- ggplot(data = mc_rho_dev, aes(x = method, y = mean_bias, fill = repunit)) +
  geom_bar(stat = "identity") +
  guides(fill = F) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~repunit)


### Printing bar plots shows slight increases in MSE, variable    ###
### decreases in mean proportional bias, and consistent decreases ###
### in mean bias (mean residual) across all reporting units       ###
