library(GSImulator)
library(dplyr)
library(tidyr)
library(ggplot2)

Nmix <- 1000
Dpar <- 1.5
ruDpar <- 1.5
Npop <- 17
rs <- rep(144, Npop)  # refsizes
ru <- c(2, 3, 12)

reppy_frame <- data_frame(repunit = paste("repu", rep(1:3, ru), sep = "_"),
                                 collection = paste("Pop", 1:Npop, sep = "_"))

rho50 <- lapply(1:50, function(rep) {
  print(rep)
  ms <- {
    rut <- rgamma(length(ru), shape = ruDpar, scale = 10)
    rut <- rut / sum(rut)
    pvec <- lapply(1:length(rut), function(x) {
      tmp <- rgamma(ru[x], shape = Dpar, scale = 10)
      tmp <- tmp / sum(tmp)
      tmp <- tmp * rut[x]
    }) %>%
      unlist

    rmultinom(1, size = Nmix, prob = pvec)[,1]

  }
  TruePi <- pvec
  TrueRho <- rut


  suppressMessages(GSImulate(refsizes = rs,
                                         mixsizes = ms,
                                         num_loci = 20,
                                         num_sets = 1,
                                         variability_string = " -t 0.3 ",
                                         marker_pars = " -u .15 3 ",
                                         M_matrix = sh_island_mig_mat(ru = c(2, 3, 12), Min = 70, Mbt = 7.0))
  )

  indat <- gsim2bhru(1, reppy_frame)
  ref <- filter(indat, sample_type == "reference")
  mix <- filter(indat, sample_type == "mixture")

  pi_mcmc <- ref_and_mix_pipeline(ref, mix, 5, method = "MCMC")$mean$pi
  rho_mcmc <- lapply(levels(ref$repunit), function(ru){
    reppy <- reppy_frame %>%
      mutate(coll_int = as.integer(factor(reppy_frame$collection, levels = unique(reppy_frame$collection)))) %>%
      filter(repunit == ru)
    out <- sum(pi_mcmc[reppy$coll_int])
  }) %>% unlist()
  rho_bh <- ref_and_mix_pipeline(ref, mix, 5, method = "BH")$mean$rho
  rho_pb <- bootstrap_rho(rho_est = rho_mcmc,
                          pi_est = pi_mcmc,
                          D = indat,
                          gen_start_col = 5)

  out <- list("true_rho" = rut, "rho_mcmc" = rho_mcmc, "rho_bh" = rho_bh, "rho_pb" = rho_pb)
})

names(rho50) <- 1:50
rho50x <- rho50 %>% bind_rows(.id = "iter")
rho50x$repunit <- rep(unique(reppy_frame$repunit), 50)

rho_data <- rho50x %>%
  gather(key = "method", value = "Estimate", rho_mcmc:rho_pb)
rho_data$repunit <- rep(c("Reporting Unit 1", "Reporting Unit 2", "Reporting Unit 3"), 150)
rho_data$method <- rep(c("MCMC", "Bayes. Hierarch.", "Bootstrap"), each = 150)
rho_data$method <- as.factor(rho_data$method)
rho_data$method <- factor(rho_data$method, levels = levels(rho_data$method)[3:1])


g <- ggplot(rho_data, aes(x = true_rho, y = Estimate, colour = repunit)) +
  geom_point() +
  facet_grid(method ~ repunit) +
  scale_color_brewer(palette = "Set1") +
  geom_abline(intercept = 0, slope = 1)
print(g)

rho_data$method <- rep(c("MCMC", "BH", "PB"), each = 150)
rho_data$method <- as.factor(rho_data$method)
rho_data$method <- factor(rho_data$method, levels = levels(rho_data$method)[c(2,3,1)])

rho_dev <- rho_data %>%
  mutate(dev = (true_rho - Estimate)^2) %>%
  mutate(prop_bias = (Estimate-true_rho) / true_rho) %>%
  mutate(bias = Estimate-true_rho) %>%
  group_by(repunit, method) %>%
  summarise(MSE = mean(dev), mean_prop_bias = mean(prop_bias), mean_bias = mean(bias))

# the standard mean error of each method and reporting unit
d <- ggplot(data = rho_dev, aes(x = method, y = MSE, fill = repunit)) +
  geom_bar(stat = "identity") +
  guides(fill = F) +
  theme(axis.title.x = element_blank()) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~repunit)
print(d)

# proportional bias makes error of equal size have larger effects
# on small values than large values
prop.bias <- ggplot(data = rho_dev, aes(x = method, y = mean_prop_bias, fill = repunit)) +
  geom_bar(stat = "identity") +
  facet_wrap(~repunit) +
  scale_fill_brewer(palette = "Set1") +
  guides(fill = F)
print(prop.bias)

# the unscaled mean bias
m.bias <- ggplot(data = rho_dev, aes(x = method, y = mean_bias, fill = repunit)) +
  geom_bar(stat = "identity") +
  guides(fill = F) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~repunit)
print(m.bias)
