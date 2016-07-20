library(rubias)
out <- Hasselman_simulation_pipeline(alewife, 15)

out$MSE <- out$mean_dev
out$method <- rep(c("BH", "EM", "MCMC"), 3)
prop.bias <- ggplot2::ggplot(dplyr::filter(out, method %in% c("BH", "MCMC")), ggplot2::aes(x = method, y = MSE, fill = repunit)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::guides(fill = FALSE) +
  ggplot2::facet_wrap( ~ repunit, nrow = 3)
print(prop.bias)


boot <- bias_comparison(alewife, 15)
rho50 <- boot$rho_iterations
names(rho50) <- 1:50
rho50x <- rho50 %>% dplyr::bind_rows(.id = "iter")
rho50x$repunit <- rep(unique(alewife$repunit), 50)

rho_data <- rho50x %>%
  tidyr::gather(key = "Method", value = "Estimate", rho_mcmc:rho_pb)
rho_data$Method <- rep(c("MCMC", "BH", "PB"), each = 150)
g <- ggplot2::ggplot(rho_data, ggplot2::aes(x = true_rho, y = Estimate, colour = repunit)) +
  ggplot2::geom_point() +
  ggplot2::facet_grid(repunit ~ Method) +
  ggplot2::geom_abline(intercept = 0, slope = 1)
print(g)


rho_dev <- boot$rho_dev
rho_dev$Method <- rep(c("BH", "MCMC", "PB"), 3)
rho_dev$MSE <- rho_dev$mse
d <- ggplot2::ggplot(data = rho_dev, ggplot2::aes(x = Method, y = MSE, fill = repunit)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::guides(fill = F) +
  ggplot2::facet_wrap(~repunit)
print(d)

b <- ggplot2::ggplot(data = rho_dev, ggplot2::aes(x = Method, y = mean_bias, fill = repunit)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::guides(fill = F) +
  ggplot2::facet_wrap(~repunit)
print(b)

