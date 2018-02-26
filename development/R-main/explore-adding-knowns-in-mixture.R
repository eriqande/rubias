

# Some code to get us up and ready to dive into infer_mixture
library(rubias)
library(tidyverse)

c2 <- chinook %>%
  mutate(known_collection = collection) %>%
  select(sample_type:collection, known_collection, indiv, everything())

c2_mix <- chinook_mix %>%
  mutate(known_collection = NA) %>%
  select(sample_type:collection, known_collection, indiv, everything())

# here we make a mixture data set where the first 4 fish in each mixture collection
# is known to be from some population.  We will sort fish first etc.
set.seed(5)
c2_kmix <- c2_mix %>%
  arrange(collection) %>%
  group_by(collection) %>%
  mutate(known_collection = c(sample(unique(c2$collection), size = 4, replace = FALSE), rep(NA, n() - 4))) %>%
  ungroup()

reference <- c2
mixture <- c2_kmix
gen_start_col <- 6
method = "MCMC"
alle_freq_prior = list("const_scaled" = 1)
reps = 2000
burn_in = 100
pb_iter = 100
sample_int_Pi = 1

# now we can start running through the infer_mixture function


# here we can just run the function
result <- infer_mixture(reference = c2, mixture = c2_kmix, gen_start_col = 6)

# now screw it up a bit
c2_kmix_mess <- c2_kmix
c2_kmix_mess[7, 4] <- "boingo"
c2_kmix_mess[10, 4] <- "Cheez-whizz"

result <- infer_mixture(reference = c2, mixture = c2_kmix_mess, gen_start_col = 6)
