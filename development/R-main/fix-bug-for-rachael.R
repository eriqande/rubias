
library(tidyverse)
library(rubias)

load("development/data/rubiasTrial.RData")


breedFrame <- as_tibble(breedFrame)  # reference frame
unkFrame <- as_tibble(unkFrame) # "real" frame of unknowns
fake <- as_tibble(fake)  # fake frame of unknowns
realunit     # actual populations of unknowns in "fake"

# make a data frame of the realunits
ru <- tibble(indiv = fake$indiv, true_unit = realunit)


# first, see what this looks like
inf_m1 <- infer_mixture(breedFrame, fake, gen_start_col = 6)$indiv_posteriors

# and get the MAP assignments and compare
inf_m1 %>%
  group_by(indiv) %>%
  filter(near(PofZ, max(PofZ))) %>%
  left_join(ru) %>%
  group_by(repunit, true_unit) %>%
  tally()

# so things are scrambled, no doubt.

# now, let's see if we get the same thing after we have converted all factors to character...
bf <- lapply(breedFrame, as.character) %>%
  as_data_frame()
ff <- lapply(fake, as.character) %>%
  as_data_frame()

inf_m2 <- infer_mixture(bf, ff, gen_start_col = 6)$indiv_posteriors


inf_m2 %>%
  group_by(indiv) %>%
  filter(near(PofZ, max(PofZ))) %>%
  left_join(ru) %>%
  group_by(repunit, true_unit) %>%
  tally()

# aha! that looks OK...

# how about if we coerce repunit and collection to integers?
bfi <- bf
bfi$repunit <- as.integer(bfi$repunit)
bfi$collection <- as.integer(bfi$collection)

ffi <- ff
ffi$repunit <- as.integer(ffi$repunit)

inf_m3 <- infer_mixture(bfi, ffi, gen_start_col = 6)$indiv_posteriors

inf_m3 %>%
  group_by(indiv) %>%
  filter(near(PofZ, max(PofZ))) %>%
  left_join(ru) %>%
  group_by(repunit, true_unit) %>%
  tally()

# that looks OK, too...


# let's see if it still bombs when things are factors but not integers...
bf2 <- breedFrame %>%
  mutate(repunit = factor(paste0("g", repunit), levels = paste0("g", 1:4))) %>%
  mutate(collection = factor(paste0("g", collection), levels = paste0("g", 1:4)))

inf_m4 <- infer_mixture(bf2, fake, gen_start_col = 6)$indiv_posteriors

inf_m4 %>%
  group_by(indiv) %>%
  filter(near(PofZ, max(PofZ))) %>%
  left_join(ru) %>%
  group_by(repunit, true_unit) %>%
  tally()


# OK, yep.  Problem. Somehow, when things are getting re-factorized it is getting screwed up.
