# should be run at the top level of the repository

#### let's say we wanted to re-do the simulations for alewife.  ####
# we would do it like this:
library(dplyr)
library(rubias)
library(ggplot2)
library(readr)
library(stringr)

ale <- assess_reference_loo(reference = alewife, gen_start_col = 15, reps = 50, seed = 55)

# then plot them by collection
ggplot(ale, aes(x = omega, y = mle, colour = repunit)) +
  geom_point() +
  facet_wrap(~ collection) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")


# have a look by reporting unit
ale_rep <- ale %>%
  group_by(iter, repunit) %>%
  summarise(rep_pi_true = sum(omega),
            rep_n_true = sum(n),
            rep_post_mean = sum(post_mean),
            rep_mle = sum(mle))

# then plot that
ggplot(ale_rep, aes(x = rep_pi_true, y = rep_mle, colour = repunit)) +
  geom_point() +
  facet_wrap(~ repunit) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")



#### Now, do it with Chinook we just need to get it in the right format and do the same thing ####

### Get Dan's data and deal with the junky loci
chin <- read.csv("development/data/Ots-baseline-numerical.csv", stringsAsFactors = FALSE) %>%
  tbl_df()



# but the file from Dan has a bunch of crap in it.  Specifically, loci with invalid characters
# like "X" and "-", but it seems like just a few loci, and the next few steps dump them
col_classes <- sapply(chin, class)
chars <- which(col_classes == "character")
dump_these <- chars[chars >= 14]

chin2 <- chin[, -dump_these]

# those loci were:
names(dump_these)[c(T,F)]

# now, get locus names + ".1" on the second column of each locus
loccols <- seq(14, ncol(chin2), by = 2)
names(chin2)[loccols + 1] <- paste(names(chin2)[loccols], "1", sep = ".")

# finally, dump the columns that we really don't need and name the columns as might be better
chin3 <- chin2[, -c(2:10, 12, 13)] %>%
  rename(pop_number = Pop.order, indiv = Drainage.code)

# turn the 0's into NAs
chin3[chin3 == 0] <- NA


### Figure out what reporting unit everything is in
grps <- read_delim("development/data/Ots-19-grps.txt", delim = "\t") %>%
  rename(rep_num = RepGroupNum,
         repunit = RepGroup,
         collection = Pop) %>%
  mutate(pop_number = as.integer(str_replace(collection, "OTS", "")),
         collection = factor(collection, levels = unique(collection)),
         repunit = factor(repunit, levels = unique(repunit)),
         sample_type = "reference")

### Then stick that all into single data frame that we can use.
chinook <- left_join(grps, chin3) %>%
  select(sample_type, collection, repunit, indiv, everything())


### Finally....run it
chin_sim <- assess_reference_loo(reference = chinook, gen_start_col = 7, reps = 300, seed = 55, alpha_repunit = 0.5, alpha_collection = 1.5)


# have a look at the individuals pops.
popsplot <- ggplot(chin_sim, aes(x = omega, y = mle, colour = repunit)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ collection) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

ggsave(popsplot, filename = "dan_chinook_pops.pdf", width = 15, height = 10)

# have a look by reporting unit
chin_rep <- chin_sim %>%
  group_by(iter, repunit) %>%
  summarise(rep_pi_true = sum(omega),
            rep_n_true = sum(n),
            rep_post_mean = sum(post_mean),
            rep_mle = sum(mle))


# then plot that
repsplot <- ggplot(chin_rep, aes(x = rep_pi_true, y = rep_mle, colour = repunit)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ repunit) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

ggsave(repsplot, filename = "dan_chinook_repunits.pdf", width = 15, height = 10)
