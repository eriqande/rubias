



library(tidyverse)
library(rubias)

baseline <- read_table(
  "development/ppns_less_than_1_over_n/data/from_e_rondeau/bco_SNP_Skeena_v5.0.1_2025-03-18_CUs-Bulk.txt.gz",
  col_types = cols(.default = col_character())
)
mixture24 <- read_table(
  "development/ppns_less_than_1_over_n/data/from_e_rondeau/mco_2024.txt",
  col_types = cols(.default = col_character())
)

# get the sample size to pretend that it is the total catch
SS <- mixture24 %>%
  count(collection) %>%
  pull(n)

tct <- tibble(collection = "Skeena_TF(24)", tot_catch = list(SS))

set.seed(5)
mix_res <- infer_mixture(baseline, mixture24, 5, total_catch_tib = tct, reps = 2e4, burn_in = 1e3)


plotty <- plot_PofZ_vs_pi(mix_res, mix_coll = "Skeena_TF(24)", left_edge_fract = 0.18)

ggsave(plot = plotty$plot, filename = "development/ppns_less_than_1_over_n/outputs/SkeenaTF24_normal.pdf", height = 5, width = 5)


plotty2 <- plot_PofZ_vs_pi(
  mix_res,
  mix_coll = "Skeena_TF(24)",
  left_edge_fract = 0.18,
  plot_flavor = "bars_fully_expanded",
  arrange_scaled_likelihood_ascending = TRUE
)


ggsave(plot = plotty2$plot, filename = "development/ppns_less_than_1_over_n/outputs/SkeenaTF24_bars_expanded_sc_like_sorted.pdf", height = 5, width = 5)



# now, do the same for 23
mixture23 <- read_table(
  "development/ppns_less_than_1_over_n/data/from_e_rondeau/mco_2023.txt",
  col_types = cols(.default = col_character())
)

SS23 <- mixture23 %>%
  count(collection) %>%
  pull(n)

tct23 <- tibble(collection = "Skeena_TF(23)", tot_catch = list(SS23))

set.seed(15)
mix_res23 <- infer_mixture(baseline, mixture23, 5, total_catch_tib = tct23, reps = 2e4, burn_in = 1e3)


# and plot them
plotty23 <- plot_PofZ_vs_pi(mix_res23, mix_coll = "Skeena_TF(23)", left_edge_fract = 0.18)

ggsave(plot = plotty23$plot, filename = "development/ppns_less_than_1_over_n/outputs/SkeenaTF23_normal.pdf", height = 5, width = 5)


plotty23_2 <- plot_PofZ_vs_pi(
  mix_res23,
  mix_coll = "Skeena_TF(23)",
  left_edge_fract = 0.18,
  plot_flavor = "bars_fully_expanded",
  arrange_scaled_likelihood_ascending = TRUE
)


ggsave(plot = plotty23_2$plot, filename = "development/ppns_less_than_1_over_n/outputs/SkeenaTF23_bars_expanded_sc_like_sorted.pdf", height = 5, width = 5)



