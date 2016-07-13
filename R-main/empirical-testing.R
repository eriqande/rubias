library(ggplot2)

repidxs <- alewife %>%
  dplyr::mutate(coll_int = as.integer(collection)) %>%
  dplyr::select(repunit, coll_int) %>%
  dplyr::group_by(repunit, coll_int) %>%
  dplyr::tally()

ref_params <- tcf2param_list(alewife, 15, summ = F)
set.seed(5)
#Let's get a random rho centered around the true values in the dataset

rho50 <- lapply(1:50, function(rr) {
  rho <- gtools::rdirichlet(1, c(1.5, 1.5, 1.5))
  drawn <- mixture_draw(alewife, rhos = rho, N = 100, min_remaining = .005)
  # get estimates of rho from vanilla mcmc
  pi_mcmc <- ref_and_mix_pipeline(drawn$reference, drawn$mixture, 15, method = "MCMC")$mean$pi
  rho_mcmc <- lapply(levels(alewife$repunit), function(ru){
    out <- sum(pi_mcmc[repidxs$coll_int[repidxs$repunit == ru]])
  }) %>% unlist()
  # and from finagled mcmc
  rho_bh <- ref_and_mix_pipeline(drawn$reference, drawn$mixture, 15, method = "BH")$mean$rho

  delin <- rbind(drawn$reference, drawn$mixture)
  ref_star_params <- tcf2param_list(delin, 15, samp_type = "reference")
  rho_mean <- lapply(1:100, function(rep) {
    sim_ns <- rmultinom(n = 1, size = length(pi_mcmc), prob = pi_mcmc)
    sim_colls <- lapply(1:length(sim_ns), function(coll){
      rep(coll, sim_ns[coll])
      }) %>%
      unlist()
    sim_inds <- gprob_sim_ind(ref_star_params, sim_colls)
    SL <- apply(exp(sim_inds), 2, function(x) x/sum(x))
    pi_pb <- gsi_mcmc_1(SL = SL,
                        Pi_init = rep(1 / ref_star_params$C, ref_star_params$C),
                        lambda = rep(1 / ref_star_params$C, ref_star_params$C),
                        reps = 2000,
                        burn_in = 100,
                        sample_int_Pi = 0,
                        sample_int_PofZ = 0)
    rho_pb <- lapply(levels(alewife$repunit), function(ru){
      out <- sum(pi_pb$mean$pi[repidxs$coll_int[repidxs$repunit == ru]])
    }) %>% unlist()
  }) %>%
    simplify2array() %>%
    rowMeans()

  rho_pb <- rho_mcmc - (rho_mean - rho_mcmc)

  out <- list("true_rho" = rho, "rho_mcmc" = rho_mcmc, "rho_bh" = rho_bh, "rho_pb" = rho_pb)
  })

names(rho50) <- 1:50
rho50 <- rho50 %>% dplyr::bind_rows(.id = "iter")

rho50x <- rho50 #just so I can mess with the format without having to rerun the loop
rho50x$repunit <- rep(c("NNE", "SNE", "MAT"), 50)

rho_data <- rho50x %>%
  tidyr::gather(key = "method", value = "rho_est", rho_mcmc:rho_pb)

rho_dev <- rho_data %>%
  dplyr::mutate(dev = (true_rho - rho_est)^2) %>%
  dplyr::mutate(prop_bias = (rho_est-true_rho) / true_rho) %>%
  dplyr::mutate(bias = rho_est-true_rho) %>%
  dplyr::group_by(repunit, method) %>%
  dplyr::summarise(mean_dev = mean(dev), mean_prop_bias = mean(prop_bias), mean_bias = mean(bias))

g <- ggplot2::ggplot(rho_data, aes(x = true_rho, y = rho_est, colour = repunit)) +
  geom_point() +
  facet_grid(repunit ~ method) +
  geom_abline(intercept = 0, slope = 1)
print(g)

d <- ggplot(data = rho_dev, aes(x = method, y = mean_dev, fill = repunit)) +
  geom_bar(stat = "identity") +
  facet_wrap(~repunit)
print(d)

prop.bias <- ggplot(data = rho_dev, aes(x = method, y = mean_prop_bias, fill = repunit)) +
  geom_bar(stat = "identity") +
  facet_wrap(~repunit)
print(prop.bias)

m.bias <- ggplot(data = rho_dev, aes(x = method, y = mean_bias, fill = repunit)) +
  geom_bar(stat = "identity") +
  facet_wrap(~repunit)
print(m.bias)

rho_dev
