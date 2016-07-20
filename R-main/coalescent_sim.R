
Nmix <- 1000
Dpar <- 1.5
ruDpar <- 1.5
Npop <- 17
rs <- rep(144, Npop)  # refsizes
ru <- c(2, 3, 12)

reppy_frame <- dplyr::data_frame(repunit = paste("repu", rep(1:3, ru), sep = "_"),
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


  suppressMessages(GSImulator::GSImulate(refsizes = rs,
                                         mixsizes = ms,
                                         num_loci = 20,
                                         num_sets = 1,
                                         variability_string = " -t 0.3 ",
                                         marker_pars = " -u .15 3 ",
                                         M_matrix = GSImulator::sh_island_mig_mat(ru = c(2, 3, 12), Min = 50, Mbt = 5.0))
  )

  indat <- GSImulator::gsim2bhru(1, reppy_frame)
  ref <- dplyr::filter(indat, sample_type == "reference")
  mix <- dplyr::filter(indat, sample_type == "mixture")

  pi_mcmc <- ref_and_mix_pipeline(ref,mix, 5, method = "MCMC")$mean$pi
  rho_mcmc <- lapply(levels(ref$repunit), function(ru){
    reppy <- reppy_frame %>%
      dplyr::mutate(coll_int = as.integer(factor(reppy_frame$collection, levels = unique(reppy_frame$collection)))) %>%
      dplyr::filter(repunit == ru)
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
rho50x <- rho50 %>% dplyr::bind_rows(.id = "iter")
rho50x$repunit <- rep(unique(reppy_frame$repunit), 50)

rho_data <- rho50x %>%
  tidyr::gather(key = "method", value = "rho_est", rho_mcmc:rho_pb)

rho_dev <- rho_data %>%
  dplyr::mutate(dev = (true_rho - rho_est)^2) %>%
  dplyr::mutate(prop_bias = (rho_est-true_rho) / true_rho) %>%
  dplyr::mutate(bias = rho_est-true_rho) %>%
  dplyr::group_by(repunit, method) %>%
  dplyr::summarise(mse = mean(dev), mean_prop_bias = mean(prop_bias), mean_bias = mean(bias))


g <- ggplot2::ggplot(rho_data, ggplot2::aes(x = true_rho, y = rho_est, colour = repunit)) +
  ggplot2::geom_point() +
  ggplot2::facet_grid(repunit ~ method) +
  ggplot2::geom_abline(intercept = 0, slope = 1)
print(g)

# the standard mean error of each method and reporting unit
d <- ggplot2::ggplot(data = rho_dev, ggplot2::aes(x = method, y = mse, fill = repunit)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::facet_wrap(~repunit)
print(d)

# proportional bias makes error of equal size have larger effects
# on small values than large values
prop.bias <- ggplot2::ggplot(data = rho_dev, ggplot2::aes(x = method, y = mean_prop_bias, fill = repunit)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::facet_wrap(~repunit)
print(prop.bias)

# the unscaled mean bias
m.bias <- ggplot2::ggplot(data = rho_dev, ggplot2::aes(x = method, y = mean_bias, fill = repunit)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::facet_wrap(~repunit)
print(m.bias)
