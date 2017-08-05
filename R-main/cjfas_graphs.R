library(ggplot2)
library(tidyr)
library(dplyr)
load("~/R/rubias/cjfas_data.RData")
################################# Figure 2 (Coalescent Simulation Cross-Validation)

rho_data <- coal_rho_data

rho_data$repunit <- rep(c("RU~1~~(N[C]==2)", "RU~2~~(N[C]==3)","RU~3~~(N[C]==12)"), 100)
rho_data$method <- rep(c("MCMC", "Bootstrap"), each = 150)
rho_data$method <- as.factor(rho_data$method)
rho_data$method <- factor(rho_data$method, levels = levels(rho_data$method)[2:1])

fig.2 <- ggplot2::ggplot(rho_data, ggplot2::aes(x = true_rho, y = Estimate, colour = repunit)) +
  ggplot2::geom_point() +
  ggplot2::facet_grid(repunit ~ method, labeller = label_parsed) +
  ggplot2::labs(x='True (Simulated) Mixing Proportion', y = 'Estimated Mixing Proportion') +
  ggplot2::guides(color=F) +
  scale_color_brewer(palette = "Set1") +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed")
print(fig.2)

# exported as PDF, dimensions 4.5 x 6 inches


######################## Figure 3 (Effects of Nc/P on the residuals)
coal_bias_data <- coal_rho_data
coal_bias_data$prop_diff <- (coal_bias_data$Estimate - coal_bias_data$true_rho) / coal_bias_data$true_rho
coal_bias_data$diff <- (coal_bias_data$Estimate - coal_bias_data$true_rho)
coal_bias_data$Np_C <- rep(c(2/17, 3/17, 12/17), 100)
coal_bias_data$Np_diff <- coal_bias_data$Np_C - coal_bias_data$true_rho

fig.3 <- ggplot2::ggplot(coal_bias_data, aes(x = Np_diff, y = diff)) +
  geom_point() +
  theme(axis.title.x = element_text(margin = margin(t = 2, b = 0)),
        axis.title.y = element_text(margin = margin(t = 0, b = 2)),
        plot.margin = margin(t=6, r=6, b = 0, l=6)) +
  facet_grid(method ~ .) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", colour = "red", level = 0) +
  labs(x= expression(over(N[C],P)~~-~rho[r]^{sim}), y = expression(Residual~(tilde(rho)[r]-rho[r]^{sim}))) +

  scale_color_brewer(palette = "Set1")
print(fig.3)

## Exported as PDF, dimensions 4.14 x 4.53 inches
