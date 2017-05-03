library(ggplot2)
library(tidyr)
library(dplyr)
load("~/R/rubias/cjfas_data.RData")
################################# Hasselman Recreation Graphs


rho_data <- Hass[[1]]


rho_data$Estimate <- rho_data$rho_est
rho_data$Repunit <- rho_data$repunit
rho_data <- rho_data[151:450,]
rho_data$Method <- rep(c("MCMC", "Bootstrap"), each = 150)
rho_data$Method <- as.factor(rho_data$Method)
rho_data$Method <- factor(rho_data$Method, levels = levels(rho_data$Method)[2:1])

g <- ggplot2::ggplot(rho_data, ggplot2::aes(x = true_rho, y = Estimate, colour = Repunit)) +
  ggplot2::geom_point() +
  ggplot2::facet_grid(Method ~ Repunit) +
  ggplot2::labs(x='True (Simulated) Mixing Proportion', y = 'Estimated Mixing Proportion') +
  ggplot2::guides(color=F) +
  scale_color_brewer(palette = "Set1") +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed")
print(g)

rho_dev <- Hass[[2]]
rho_dev$Method <- rep(c("BH", "MCMC", "PB"), 3) %>%
  as.factor()
rho_dev$Method <- factor(rho_dev$Method, levels = levels(rho_dev$Method)[c(2,3,1)])
rho_dev$repunit <- rep(c("NNE", "SNE", "MAT"), each = 3) %>%
  as.factor()
rho_dev$repunit <- factor(rho_dev$repunit, levels = levels(rho_dev$repunit)[c(2,3,1)])
rho_dev <- rho_dev %>% dplyr::filter(Method != "BH")

d <- ggplot2::ggplot(data = rho_dev, ggplot2::aes(x = Method, y = MSE, fill = repunit)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::guides(fill = F) +
  ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
  scale_fill_brewer(palette = "Set1") +
  ggplot2::facet_wrap(~repunit)
print(d)

b <- ggplot2::ggplot(data = rho_dev, ggplot2::aes(x = Method, y = mean_bias, fill = repunit)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::guides(fill = F) +
  scale_fill_brewer(palette = "Set1") +
  ggplot2::facet_wrap(~repunit)
print(b)


###################################### Cross-Validation Graphs

rho50 <- boot$rho_iterations
names(rho50) <- 1:50
rho50x <- rho50 %>% dplyr::bind_rows(.id = "iter")
rho50x$repunit <- rep(unique(alewife$repunit), 50)

rho_data <- rho50x %>%
  tidyr::gather(key = "Method", value = "Estimate", rho_mcmc:rho_pb)
rho_data$Method <- rep(c("MCMC", "BH", "PB"), each = 150)
rho_data <- rho_data[-(151:300),]
g <- ggplot2::ggplot(rho_data, ggplot2::aes(x = true_rho, y = Estimate, colour = repunit)) +
  ggplot2::geom_point() +
  ggplot2::facet_grid(Method ~ repunit) +
  ggplot2::labs(x='True (Simulated) Mixing Proportion', y = 'Estimated Mixing Proportion') +
  ggplot2::guides(color=F) +
  scale_color_brewer(palette = "Set1") +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed")
print(g)


rho_dev <- boot$rho_dev
rho_dev$Method <- rep(c("BH", "MCMC", "PB"), 3)
rho_dev <- rho_dev %>% dplyr::filter(Method != "BH")
rho_dev$MSE <- rho_dev$mse
d <- ggplot2::ggplot(data = rho_dev, ggplot2::aes(x = Method, y = MSE, fill = repunit)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::guides(fill = F) +
  scale_fill_brewer(palette = "Set1") +
  ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
  ggplot2::facet_wrap(~repunit)
print(d)

b <- ggplot2::ggplot(data = rho_dev, ggplot2::aes(x = Method, y = mean_bias, fill = repunit)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::guides(fill = F) +
  scale_fill_brewer(palette = "Set1") +
  ggplot2::facet_wrap(~repunit)
print(b)

#####################     Coalescent simulation Graphs
rho_data <- coal_data

rho_data$repunit <- rep(c("RU 1", "RU 2", "RU 3"), 150)
rho_data$method <- rep(c("MCMC", "Bayes. Hierarch.", "Bootstrap"), each = 150)
rho_data <- rho_data[ -(151:300),]
rho_data$method <- as.factor(rho_data$method)
rho_data$method <- factor(rho_data$method, levels = levels(rho_data$method)[3:1])


g <- ggplot(rho_data, aes(x = true_rho, y = Estimate, colour = repunit)) +
  geom_point() +
  facet_grid(repunit ~ method) +
  ggplot2::labs(x='True (Simulated) Mixing Proportion', y = 'Estimated Mixing Proportion') +
  ggplot2::guides(color=F) +
  scale_color_brewer(palette = "Set1") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")
print(g)

rho_data$method <- rep(c("MCMC", "PB"), each = 150)
rho_data$method <- as.factor(rho_data$method)

rho_dev <- rho_data %>%
  mutate(dev = (true_rho - Estimate)^2) %>%
  mutate(prop_bias = (Estimate-true_rho) / true_rho) %>%
  mutate(bias = Estimate-true_rho) %>%
  group_by(repunit, method) %>%
  summarise(MSE = mean(dev), mean_prop_bias = mean(prop_bias), mean_bias = mean(bias))
rho_dev$method <- rep(c("MCMC", "PB"), 3)

# the standard mean error of each method and reporting unit
d <- ggplot(data = rho_dev, aes(x = method, y = MSE, fill = repunit)) +
  geom_bar(stat = "identity") +
  guides(fill = F) +
  theme(axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~repunit)
print(d)


# the unscaled mean bias
m.bias <- ggplot(data = rho_dev, aes(x = method, y = mean_bias, fill = repunit)) +
  geom_bar(stat = "identity") +
  guides(fill = F) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~repunit)
print(m.bias)

######################## Graph depicting effects of Nc/P on the residuals
coal_data <- coal_data %>% filter(method != "BH")
coal_data$prop_diff <- (coal_data$Estimate - coal_data$true_rho) / coal_data$true_rho
coal_data$diff <- (coal_data$Estimate - coal_data$true_rho)
coal_data$Np_C <- rep(c(2/17, 3/17, 12/17), 100)
coal_data$Np_diff <- coal_data$Np_C - coal_data$true_rho

cbp <- ggplot2::ggplot(coal_data, aes(x = Np_diff, y = diff)) +
  geom_point() +
  facet_grid(method ~ .) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", colour = "red", level = 0) +
  labs(x= expression(over(N[C],P)~~-~rho[r]^{sim}), y = expression(Residual~(hat(rho)[r]-rho[r]^{sim} ))) +
  scale_color_brewer(palette = "Set1")
print(cbp)
