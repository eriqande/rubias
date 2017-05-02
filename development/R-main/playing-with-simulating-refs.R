


# working on a few things to with simulating lots of mixtures at randomly
# chosen mixing proportions.

# ultimately I think we will want to be able to specify the
# alpha-parameters so that we can get the proportions of some reporting units
# much larger, on average.

# one way that would be nice to be able to parameterize things would be to
# have an expected fraction from each repunit and then have a weight applied to
# that.

# at any rate, this is just some code to play with the chinook in this context.

simmed <- assess_reference_loo(chinook, 5, reps = 50, mixsize = 2000)
whop <- simmed %>% group_by(iter, repunit) %>% summarise(rho = sum(omega), rho_post = sum(post_mean))

ggplot(whop, aes(x = rho, y = rho_post)) +
  geom_point() +
  facet_wrap(~repunit) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "red")
