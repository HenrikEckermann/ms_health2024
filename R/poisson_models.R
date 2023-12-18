# H1a mat sens in 1st year of life --> fewer health complaints during first 14 years of life
# H1b higher attachment security in 1st year of life --> fewer health complains during first 14 years of life
# H1b EXPLORE: secure vs insecure and 4 attachment categories
# H2 higher mat sens across first 14 years of life --> fewer health complaints

library(tidyverse)
library(brms)
library(glue)
library(patchwork)
library(HDInterval)

load(here::here("data/rdata/d.Rds"))
# we need an indicator of whether the child is 1 year or older for the models
# because the data have slighty different data generating processes
dlong$firstyear <- ifelse(dlong$year == 1, 1, 0)


# cov names I need repeatedly
covs <- c(
  "sex",
  "edu",
  "sib",
  "cc",
  "bf"
)
# for labelling in plots
symp_names <- tibble(
  symptom = c(
    "resp", "dig", "skin", "gen", "totalhealth"
  ),
  name = c(
    "Respiratory Symptoms", "Digestive Symptoms", "Skin Symptoms", "General Symptoms", "Total Symptoms"
  )
)


# get an impression of data distribution
integers <- c("resp", "skin", "dig", "gen", "totalhealth")
dists <- map(integers, function(symptom) {
  dlong[dlong[[symptom]] < 500, ] %>%
    mutate(across(all_of(symptom), function(x) as.numeric(x))) %>%
    ggplot(aes_string(x = symptom)) +
    geom_histogram()
})
dists[[3]]


###############################################################################
#######################     Prior Predictive checks     #######################
###############################################################################


# prior selection with poisson models is less straight forward, therefore I
# do prior predictive checks rather than sticking to default priors of brms.
# a wide uninformative prior assumes expected values in the range of 8e17
lb <- rlnorm(n = 1e4, mean = 0, sd = 3)
plot(density(lb))
# what about a narrower prior
lb <- rlnorm(n = 1e4, mean = 2, sd = 0.5)
plot(density(lb))
# this prior works already but we would like more likelihood at 0 and wider tail
lb <- rlnorm(n = 1e4, mean = 2, sd = 0.75)
plot(density(lb))
# the expected value would be ~10 complains per child in 12 months
exp(2 + 0.75^2 / 2)
# this would be ok for a model with only an intercept, now we add a prior for year
# we assume that there are more health complaints in the first few years due to
# higher (skin) sensitivity of infants as well as entrance to childcare in NL



# number of curves 
n <- 100
p <- tibble(
  i = 1:n,
  a = rnorm(n, mean = 2, sd = 0.75)
) %>%
  mutate(
    `beta%~%Normal(-0.15*', '*0.15)`  = rnorm(n, mean = -0.15 , sd = 0.15),
    `beta%~%Normal(0*', '*0.2)`  = rnorm(n, mean = 0 , sd = 0.2)
  ) %>%
  pivot_longer(contains("beta"), values_to = "b", names_to = "prior") %>%
  expand_grid(x = seq(from = 1, to = 14, length.out = 100)) %>%
  ggplot(aes(x = x, y = exp(a + b * x), group = i)) +
  geom_line(linewidth = 1/4, alpha = 2/3) +
  labs(x = "years", y = "symptoms") +
  facet_wrap(~prior, scales = "free_y")
p 
# we can  see that these priors are barely informative allowing for relationships
# in other directions than we would assume. At the same time they dont lead to 
# totally unrealistic prior predicitons. They will be easily overwhelmed
# by the data.

# and finally for all other coefs where we dont have prior information except 
# that they should not allow unrealistic number of symptoms
p <- tibble(
  i = 1:n,
  a = rnorm(n, mean = 2, sd = 0.75)
) %>%
  mutate(
    `beta%~%Normal(0*', '*0.2)`  = rnorm(n, mean = 0 , sd = 0.2),
    `beta%~%Normal(0*', '*0.5)`  = rnorm(n, mean = 0 , sd = 0.5)
  ) %>%
  pivot_longer(contains("beta"), values_to = "b", names_to = "prior") %>%
  expand_grid(x = seq(from = -2, to = 2, length.out = 100)) %>%
  ggplot(aes(x = x, y = exp(a + b * x), group = i)) +
  geom_line(linewidth = 1/4, alpha = 2/3) +
  labs(x = "x", y = "symptoms") +
  facet_wrap(~prior, scales = "free_y")


colnames(dlong)
# the N(0, 0.2) prior works well. Lets confirm this with a prior predictive check
# in brms
get_prior(resp ~ year + sensitivity_firstyear, data = dlong)
m0_prior <- brm(
  family = poisson(),
  formula = resp ~ year,
  data = dlong,
  prior = c(
    prior(normal(-0.15, 0.15), coef = "year"),
    #prior(normal(0, 0.5), class = "b"),
    prior(normal(2, 0.75), class = "Intercept")
    
  ),
  sample_prior = "only"
)

prior <- as_draws_df(m0_prior)
prior %>%
  slice_sample(n = 50) %>%
  rownames_to_column("draw") %>%
  expand_grid(year = 1:15) %>%
  mutate(lambda = exp(b_Intercept + b_year * year)) %>%
  ggplot(aes(x = year, y = lambda)) +
  geom_line(aes(group = draw), color = "firebrick", alpha = 0.4)

# looks good. Now for any other predictor with unknown effect.
m1_prior <- brm(
  family = poisson(),
  formula = resp ~ sensitivity_firstyear,
  data = dlong,
  prior = c(
    #prior(normal(-0.2, 0.2), coef = "year"),
    prior(normal(0, 0.2), class = "b"),
    prior(normal(2, 0.75), class = "Intercept")
    
  ),
  sample_prior = "only"
)

prior <- as_draws_df(m1_prior)
prior %>%
  slice_sample(n = 50) %>%
  rownames_to_column("draw") %>%
  expand_grid(sensitivity_firstyear = seq(-2, 2, length.out = 100)) %>%
  mutate(lambda = exp(b_Intercept + b_sensitivity_firstyear * sensitivity_firstyear)) %>%
  ggplot(aes(x = sensitivity_firstyear, y = lambda)) +
  geom_line(aes(group = draw), color = "firebrick", alpha = 0.4)


# also looks good. Now we are ready to start building the models
# note that in the following I work with the variable resp only
# to explore optimal model structure. This I will use to select candidate models
# then later, I will let these candidate models compete using loo for each 
# variable and imputed dataset. Then we report the statistics for the winning 
# model.



###############################################################################
#############################     Model Building   ############################
###############################################################################


# the base model
m0 <- brm(
  family = poisson(),
  formula = resp ~ year,
  data = dlong,
  prior = c(
    prior(normal(-0.15, 0.15), coef = "year"),
    prior(normal(2, 0.75), class = "Intercept")
  ),
  file = here::here("data/rdata/m0_resp.Rds")
)

post <- as_draws_df(m0)
post %>%
  slice_sample(n = 50) %>%
  rownames_to_column("draw") %>%
  expand_grid(year = 0:15) %>%
  mutate(lambda = exp(b_Intercept + b_year * year)) %>%
  ggplot(aes(x = year, y = lambda)) +
  geom_line(aes(group = draw), color = "firebrick", alpha = 0.4) +
  scale_x_continuous(breaks = 0:15)

filter(dlong, year == 1) %>%
  summarise(m = mean(resp, na.rm = TRUE))

# average symptoms count at age 1 is ok but 8 lower compared to sample average
# lets see if we get different results if I leave out prior specification
m0_2 <- brm(
  family = poisson(),
  formula = resp ~ year,
  data = dlong,
  file = here::here("data/rdata/m0_noprior_resp.Rds")
)
# results are the same

post <- as_draws_df(m0_2)
post %>%
  slice_sample(n = 50) %>%
  rownames_to_column("draw") %>%
  expand_grid(year = 0:15) %>%
  mutate(lambda = exp(b_Intercept + b_year * year)) %>%
  ggplot(aes(x = year, y = lambda)) +
  geom_line(aes(group = draw), color = "firebrick", alpha = 0.4) +
  scale_x_continuous(breaks = 0:15)



m1 <- brm(
  family = poisson(),
  formula = resp ~ year  + I(year^2),
  data = dlong,
  prior = c(
    prior(normal(-0.15, 0.15), coef = "year"),
    prior(normal(0, 0.5), class = "b"),
    prior(normal(2, 0.75), class = "Intercept")
  ),
  file = here::here("data/rdata/m1_resp.Rds")
)
post <- as_draws_df(m1)
head(post)
post %>%
  slice_sample(n = 50) %>%
  rownames_to_column("draw") %>%
  expand_grid(year = 0:15) %>%
  mutate(lambda = exp(b_Intercept + b_year * year + b_IyearE2 * year)) %>%
  ggplot(aes(x = year, y = lambda)) +
  geom_line(aes(group = draw), color = "firebrick", alpha = 0.4) +
  scale_x_continuous(breaks = 0:15)

# this model fits much better. lets confirm this with loo
loo_m0 <- add_criterion(
  m0,
  "loo",
  file = here::here("data/rdata/loo_m0"),
  moment_match = FALSE
)
loo_m1 <- add_criterion(
  m1,
  "loo",
  file = here::here("data/rdata/loo_m1"),
  moment_match = FALSE
)




# now lets add the random intercepts and see if that improves model fit further
# as expected

m2 <- brm(
  family = poisson(),
  formula = resp ~ year + I(year^2) +  (1 | id),
  data = dlong,
  prior = c(
    prior(normal(-0.15, 0.15), coef = "year"),
    prior(normal(2, 0.75), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  file = here::here("data/rdata/m2_refit_resp.Rds")
)
summary(m2)
loo_m2 <- add_criterion(
  m2,
  "loo",
  file = here::here("data/rdata/loo_m2"),
  moment_match = FALSE
)

# again clearly improving model fit as expected. The random intercept also
# should nicely capture differences in reporting tendencies independent of 
# maternal sensitivity. I will do posterior predictive 
# checks in order to find out where the models fails

pred <- posterior_predict(m2) %>%
  as.data.frame()
cnames <- filter(dlong, !is.na(resp)) %>% 
  mutate(idnew = glue("id_{id}_{year}")) %>%
  .$idnew
colnames(pred) <- cnames
colnames(pred)
pred_sum <- pred %>%
  pivot_longer(
    everything(),
    names_to = "idfull",
    values_to = "pred"
  ) %>%
  group_by(idfull) %>%
  summarise(
    m = median(pred),
    lower = quantile(pred, 0.025),
    upper = quantile(pred, 0.975)
  ) %>%
  mutate(
    id = as.integer(str_extract(idfull, "(\\d{3})")),
    year = as.numeric(str_extract(idfull, "\\d\\d?\\.?\\d?$"))
  ) %>%
  left_join(dlong, by = c("id", "year"))

all_preds1 <- map(unique(dlong$id), function(idn) {
  filter(pred_sum, id == idn) %>%
    ggplot(aes(x = year, y = resp)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), color = "skyblue") +
    geom_point(aes(y = m), color = "skyblue", size = 2) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = 1:14) +
    scale_y_continuous(breaks = seq(0, max(dlong$resp, na.rm = TRUE), 5), limits = c(0, max(dlong$resp, na.rm = TRUE))) +
    theme_bw(base_size = 20)
})
# here, for each individual, we can check posterior predictions
all_preds1[[16]]


# posterior checks show that in the first year, the model typically
# underestimates number of symptoms and overestimates in the next timepoint.
# this is what we expect from the data generation process as we did not yet 
# implement the difference of first year vs rest of data, which we do in the next
# step using a variable called firstyear which basically indicates if it is the first year or 
# not, which we possibly let interact with cc and bf as those factors should have 
# strongest effects in the first year (cc causes more symptoms as you go but then
# no more or even less over time (immunity))
# Other things I tried but deleted are cubic effects and random effects
# for the time vars

# m2_cubic <- brm(
#   family = poisson(),
#   formula = resp ~ year + I(year^2) + I(year^3) +  (1 | id),
#   data = dlong,
#   prior = c(
#     prior(normal(-0.15, 0.15), coef = "year"),
#     prior(normal(2, 0.75), class = "Intercept"),
#     prior(exponential(1), class = "sd")
#   ),
#   file = here::here("data/rdata/m2_cubic_resp.Rds")
# )
# summary(m2_cubic)
# loo_m2_cubic <- add_criterion(
#   m2_cubic,
#   "loo",
#   file = here::here("data/rdata/loo_m2_cubic"),
#   moment_match = FALSE
# )
# # cubic does not fit well at all
# 
# # maybe random slope for year will help
# m2_rs <- brm(
#   family = poisson(),
#   formula = resp ~ year + I(year^2) +  (1 + year | id),
#   data = dlong,
#   prior = c(
#     prior(normal(-0.15, 0.15), coef = "year"),
#     prior(normal(2, 0.75), class = "Intercept"),
#     prior(exponential(1), class = "sd")
#   ),
#   file = here::here("data/rdata/m2_rs_resp.Rds")
# )
# summary(m2_rs)
# loo_m2_rs <- add_criterion(
#   m2_rs,
#   "loo",
#   file = here::here("data/rdata/loo_m2_rs"),
#   moment_match = FALSE
# )
# 
# # does not fit well either
# m2_rs2 <- brm(
#   family = poisson(),
#   formula = resp ~ year + I(year^2) +  (1 + year + I(year^2) | id),
#   data = dlong,
#   prior = c(
#     prior(normal(-0.15, 0.15), coef = "year"),
#     prior(normal(2, 0.75), class = "Intercept"),
#     prior(exponential(1), class = "sd")
#   ),
#   file = here::here("data/rdata/m2_rs2_resp.Rds")
# )
# summary(m2_rs2)
# loo_m2_rs2 <- add_criterion(
#   m2_rs2,
#   "loo",
#   file = here::here("data/rdata/loo_m2_rs2"),
#   moment_match = FALSE
# )

# finally, we want to decide if exact age works better (in order to compare
# i must impute missing age with the year)
m2_age <- brm(
  family = poisson(),
  formula = resp ~ age + I(age^2) +  (1 | id),
  data = mutate(dlong, age = ifelse(is.na(age), year, age)),
  prior = c(
    prior(normal(-0.3, 0.2), coef = "age"),
    prior(normal(0, 0.5), class = "b"),
    prior(normal(2, 0.75), class = "Intercept")
  ),
  file = here::here("data/rdata/m2_age2_resp.Rds")
)
summary(m2_age)
loo_m2_age <- add_criterion(
  m2_age,
  "loo",
  file = here::here("data/rdata/loo_m2_age_resp"),
  moment_match = FALSE
)

# surprisingly the age model performs worse, indicating that somewhere noise
# is introduced. For sake of completeness I will still try a model with age diff
m2_age_diff <- brm(
  family = poisson(),
  formula = resp ~ year + I(year^2) + age_diff +  (1 | id),
  data = mutate(dlong, age_diff = ifelse(is.na(age_diff), 0, age_diff)),
  prior = c(
    prior(normal(-0.15, 0.15), coef = "year"),
    prior(normal(0, 0.5), class = "b"),
    prior(normal(2, 0.75), class = "Intercept")
  ),
  file = here::here("data/rdata/m2_age_diff_resp.Rds")
)

loo_m2_age_diff <- add_criterion(
  m2_age_diff,
  "loo",
  file = here::here("data/rdata/loo_m2_age_diff"),
  moment_match = FALSE
)
lcomp <- loo_compare(loo_m0, loo_m1, loo_m2, loo_m2_age, loo_m2_age_diff)
lcomp
# there is zero added benefit of adding the little deviations from the targeted 
# age of collection, therefore we continue with the year variable alone
# lets add the first year var and check if post pred checks are better

m3 <- brm(
  family = poisson(),
  formula = resp ~ year + I(year^2) + firstyear + (1 | id),
  data = dlong,
  prior = c(
    prior(normal(-0.15, 0.15), coef = "year"),
    prior(normal(2, 0.75), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  file = here::here("data/rdata/m3_resp.Rds")
)
summary(m3)
loo_m3 <- add_criterion(
  m3,
  "loo",
  file = here::here("data/rdata/loo_m3"),
  moment_match = FALSE
)

pred <- posterior_predict(m3) %>%
  as.data.frame()
cnames <- filter(dlong, !is.na(resp)) %>% 
  mutate(idnew = glue("id_{id}_{year}")) %>%
  .$idnew
colnames(pred) <- cnames
colnames(pred)
pred_sum <- pred %>%
  pivot_longer(
    everything(),
    names_to = "idfull",
    values_to = "pred"
  ) %>%
  group_by(idfull) %>%
  summarise(
    m = median(pred),
    lower = quantile(pred, 0.025),
    upper = quantile(pred, 0.975)
  ) %>%
  mutate(
    id = as.integer(str_extract(idfull, "(\\d{3})")),
    year = as.numeric(str_extract(idfull, "\\d\\d?\\.?\\d?$"))
  ) %>%
  left_join(dlong, by = c("id", "year"))

all_preds <- map(unique(dlong$id), function(idn) {
  filter(pred_sum, id == idn) %>%
    ggplot(aes(x = year, y = resp)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), color = "skyblue") +
    geom_point(aes(y = m), color = "skyblue", size = 2) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = 1:14) +
    theme_bw(base_size = 20)
})
all_preds[[11]]

# clearly improved fit, now with CC interaction


m4 <- brm(
  family = poisson(),
  formula = resp ~ year + I(year^2) + firstyear * cc + (1 | id),
  data = dlongimp[[1]],
  prior = c(
    prior(normal(-0.15, 0.15), coef = "year"),
    prior(normal(2, 0.75), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  file = here::here("data/rdata/m4_resp.Rds")
)
summary(m4)
loo_m4 <- add_criterion(
  m4,
  "loo",
  file = here::here("data/rdata/loo_m4"),
  moment_match = FALSE
)

m4_2 <- brm(
  family = poisson(),
  formula = resp ~ year + I(year^2) + firstyear + firstyear:cc + (1 | id),
  data = dlongimp[[1]],
  prior = c(
    prior(normal(-0.15, 0.15), coef = "year"),
    prior(normal(2, 0.75), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  file = here::here("data/rdata/m4_2_resp.Rds")
)
summary(m4_2)
loo_m4_2 <- add_criterion(
  m4_2,
  "loo",
  file = here::here("data/rdata/loo_m4_2"),
  moment_match = FALSE
)

pred <- posterior_predict(m4) %>%
  as.data.frame()
cnames <- filter(dlong, !is.na(resp)) %>% 
  mutate(idnew = glue("id_{id}_{year}")) %>%
  .$idnew
colnames(pred) <- cnames
colnames(pred)
pred_sum <- pred %>%
  pivot_longer(
    everything(),
    names_to = "idfull",
    values_to = "pred"
  ) %>%
  group_by(idfull) %>%
  summarise(
    m = median(pred),
    lower = quantile(pred, 0.025),
    upper = quantile(pred, 0.975)
  ) %>%
  mutate(
    id = as.integer(str_extract(idfull, "(\\d{3})")),
    year = as.numeric(str_extract(idfull, "\\d\\d?\\.?\\d?$"))
  ) %>%
  left_join(dlong, by = c("id", "year"))

all_preds2 <- map(unique(dlong$id), function(idn) {
  filter(pred_sum, id == idn) %>%
    ggplot(aes(x = year, y = resp)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), color = "skyblue") +
    geom_point(aes(y = m), color = "skyblue", size = 2) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = 1:14) +
    scale_y_continuous(breaks = seq(0, max(dlong$resp, na.rm = TRUE), 5), limits = c(0, max(dlong$resp, na.rm = TRUE))) +
    theme_bw(base_size = 20)
})

all_preds <- map(1:length(all_preds), function(i) {
  all_preds1[[i]] + all_preds2[[i]] + theme(axis.text.y = element_blank(),
                                            axis.ticks.y = element_blank(),
                                            axis.title.y = element_blank())
})


all_preds[[12]]

# even better! lastly, a similar model would make sense for bf. 
# Bf protects somewhat from respiratory diseases but mainly in the first year.
m5 <- brm(
  family = poisson(),
  formula = resp ~ year + I(year^2) + firstyear * cc + bf * firstyear + (1 | id),
  data = dlongimp[[1]],
  prior = c(
    prior(normal(-0.15, 0.15), coef = "year"),
    prior(normal(2, 0.75), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  file = here::here("data/rdata/m5_resp.Rds")
)
summary(m5)
loo_m5 <- add_criterion(
  m5,
  "loo",
  file = here::here("data/rdata/loo_m5"),
  moment_match = FALSE
)

# this model is a tiny bit better, but not significantly. As the effect size
# is in the expected direction and it makes theoretically sense, I choose this
# model as the new base model. Now I will add edu and then in a next step sib

m6 <- brm(
  family = poisson(),
  formula = resp ~ year + I(year^2) + firstyear * cc + bf * firstyear + edu + (1 | id),
  data = dlongimp[[1]],
  prior = c(
    prior(normal(-0.15, 0.15), coef = "year"),
    prior(normal(2, 0.75), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  file = here::here("data/rdata/m6_resp.Rds")
)
summary(m6)
loo_m6 <- add_criterion(
  m6,
  "loo",
  file = here::here("data/rdata/loo_m6"),
  moment_match = FALSE
)

m7 <- brm(
  family = poisson(),
  formula = resp ~ year + I(year^2) + firstyear * cc + bf * firstyear + edu + sib + (1 | id),
  data = dlongimp[[1]],
  prior = c(
    prior(normal(-0.15, 0.15), coef = "year"),
    prior(normal(2, 0.75), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  file = here::here("data/rdata/m7_resp.Rds")
)
summary(m7)
loo_m7 <- add_criterion(
  m7,
  "loo",
  file = here::here("data/rdata/loo_m7"),
  moment_match = FALSE
)


loo_compare(loo_m1, loo_m2, loo_m3, loo_m4, loo_m4_2, loo_m5, loo_m6, loo_m7)

# so m4 is the best candidate when also trying to keep complexity lowest but
# m5 and m6 are also okay. m7 is worse fit. Now we are ready to test our
# hypothesis One more thing: It makes sense to fit a random slope for the variable
# of interest BUT this cannot be tested with the base models as they do not 
# include those variables. Therefore, I will test whether or not such a model
# fits better or not using the candidate models. 



###############################################################################
##############################     Final Model Fit  ###########################
###############################################################################

# in the following code, we fit several candidate models to the imputed data,
# then determine which model fits best. We work with the posterior distributions
# of the best fitting model to report test statistics and visualize the data.
# we repeat these steps for each outcome and the different hypotheses.


integers <- c("resp", "skin", "dig", "gen", "totalhealth")


# 1st year sensitivity
firstyear_sens <- map(integers, function(symptom) {
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + sensitivity_firstyear +  (1 | id)"))
  m8_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.)))
    m8 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m8_imp_{m}_{symptom}.Rds"))
    )
    loo_m8 <- add_criterion(
      m8,
      "loo",
      ##file = here::here(glue("data/rdata/loo_m8_imp_{m}_{symptom}")),
      moment_match = FALSE
    )
    list(
      m8,
      loo_m8
    )
  })
  m8_post <- m8_post <- map_dfr(m8_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  m8_h <- map(m8_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "sensitivity_firstyear < 0")
  })
  
  
  
  # with random slopes 
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + sensitivity_firstyear +  (1 | id) + (sensitivity_firstyear | yearf)"))
  m11_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.))) %>%
              mutate(yearf = as.factor(year))
    m11 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m11_imp_{m}_{symptom}.Rds"))
    )
    loo_m11 <- add_criterion(
      m11,
      "loo",
      ##file = here::here(glue("data/rdata/loo_m11_imp_{m}_{symptom}")),
      moment_match = FALSE
    )
    list(
      m11,
      loo_m11
    )
  })
  m11_post <- map_dfr(m11_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  m11_h <- map(m11_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "sensitivity_firstyear < 0")
  })
  
  
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + bf * firstyear + edu + sib + sensitivity_firstyear + (1 | id)"))
  m10_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.)))
    m10 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m10_imp_{m}_{symptom}.Rds"))
    )
    loo_m10 <- add_criterion(
      m10,
      "loo",
      ##file = here::here(glue("data/rdata/loo_m10_imp_{m}_{symptom}")),
      moment_match = FALSE
    )
    list(
      m10,
      loo_m10
    )
    
  })
  m10_post <- map_dfr(m10_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  m10_h <- map(m10_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "sensitivity_firstyear < 0")
  })
  
  
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + bf * firstyear + edu + sib + sensitivity_firstyear + (1 | id) + (sensitivity_firstyear | yearf)"))
  m12_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.))) %>%
              mutate(yearf = as.factor(year))
    m12 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m12_imp_{m}_{symptom}.Rds"))
    )
    loo_m12 <- add_criterion(
      m12,
      "loo",
      ##file = here::here(glue("data/rdata/loo_m12_imp_{m}_{symptom}")),
      moment_match = FALSE
    )
    list(
      m12,
      loo_m12
    )
    
  })
  m12_post <- m12_post <- map_dfr(m12_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  m12_h <- map(m12_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "sensitivity_firstyear < 0")
  })
  
  
  
  loo <- loo_compare(m8_models[[1]][[2]], m10_models[[1]][[2]], m11_models[[1]][[2]], m12_models[[1]][[2]])
  winner <- rownames(loo)[1]
  list(
    m8 = m8_models,
    m10 = m10_models,
    loo = loo,
    winner = winner,
    m8_h = m8_h,
    m10_h = m10_h,
    m8_post = m8_post,
    m10_post = m10_post,
    
    m11 = m11_models,
    m12 = m12_models,
    m11_h = m11_h,
    m12_h = m12_h,
    m11_post = m11_post,
    m12_post = m12_post
  )
})

# lets check which model we should interpret
map(firstyear_sens, ~.x$winner)
map2(1:5, firstyear_sens, function(no, mlist) {
  if (no != 3) {
    mlist$m12_h
  } else {
    mlist$m11_h
  }
})


# results indicate that mat sens in the first year mainly effects resp symptoms
# a similar trend can be observed across other symptoms such that also 3 out of
# 5 models are rejecting the null for total health, but the other single symptoms
# cannot reject the null.

# to report effect size i use combined posteriors. For directional hypothesis
b <- map2_dfr(firstyear_sens, 1:5, function(mlist, no) {
  if (no != 3) {
    posttemp <- mlist$m12_post
  } else {
    posttemp <- mlist$m11_post
  }
  select(posttemp, b_sensitivity_firstyear) %>%
    summarise(
      m = median(b_sensitivity_firstyear),
      lower = hdi(b_sensitivity_firstyear)[1],
      upper = hdi(b_sensitivity_firstyear)[2],
      lower_d = quantile(b_sensitivity_firstyear, 0),
      upper_d = quantile(b_sensitivity_firstyear, 0.95)
    ) %>%
    mutate(y = integers[no])
})
p <- map2_dfr(1:5, firstyear_sens, function(no, mlist) {
  if (no != 3) {
    htemp <- mlist$m12_h
  } else {
    htemp <- mlist$m11_h
  }
  p <- map_dbl(htemp, function(h) {
    h$hypothesis$Post.Prob
  })
  tibble(y = integers[no], p = median(p))
})
b <- full_join(b, p, by = "y") %>%
  mutate(across(where(is.numeric), function(x) format(round(x, 3), digits = 3)))
save(b, file = here::here("data/rdata/b_sensitivity_firstyear.Rds"))



# plot
year <- unique(dlong$year)
post_year <- map_dfr(year, function(yr) {
  # generate outcomes for effect size prediction
  nd <- expand_grid(
    cc = c(0, 1), 
    bf = median(dlongimp[[1]]$bf),
    year = yr,
    edu = c(1, 2),
    sib = c(0, 1),
    #sensitivity_all = seq(-2, 2, 1),
    sensitivity_firstyear = seq(-2, 2, 1),
    #attach_dummy = c(0, 1),
    #attach_cat = 1:4
    firstyear = ifelse(yr == 1, 1, 0)
  )
  nd$yearf <- factor(nd$year)
  
  post_pred <- map_dfr(1:5, function(no) {
      models <- "m11"
    if (b[b$y == integers[no], "p"] >= 0.3) {
      post_pred <- map_dfr(firstyear_sens[[no]][[models]], function(model) {
        posterior_predict(model[[1]], newdata = nd, re.form =  "~(sensitivity_firstyear|yearf)") %>%
          as.data.frame() %>%
          pivot_longer(
            everything(),
            names_to = "var",
            values_to = "values"
          ) %>%
          mutate(
            #values = values - median(values),
            y = integers[no]
          )
      })
    }
  })
  
  match_names <- tibble(
    var = unique(post_pred$var), 
    sensitivity_firstyear = nd$sensitivity_firstyear,
    firstyear = 0,
    cc = nd$cc, 
    bf = nd$bf,
    year = nd$year,
    edu = nd$edu,
    sib = nd$sib
  )
  post_pred %>%
    left_join(match_names, by = "var") %>%
    select(-var) 
})


summarised <- post_year %>%
  group_by(sensitivity_firstyear, year, y) %>%
  summarise(
    m = median(values),
    lower = hdi(values)[1],
    upper = hdi(values)[2],
    lower_50 = hdi(values, credMass = 0.68)[1],
    upper_50 = hdi(values, credMass = 0.68)[2]
  )  

post_nested <- group_by(post_year, y) %>% nest()
post_nested[[2]][[1]]
p_sens_firstyear <- map2(post_nested$y, post_nested$data, function(ystring, data) {
  p <- filter(summarised, y == ystring) %>%
      mutate(year = as.factor(year)) %>%
    ggplot(aes(sensitivity_firstyear, m, group = year, color = year)) +
      geom_point() +
      geom_line() +
      xlab(glue("Early Sensitivity")) + 
      ylab(filter(symp_names, symptom == ystring) %>% .$name) +
      labs(color = "Year") +
      theme_bw(base_size = 20) 
  p
})

save(p_sens_firstyear, file = here::here("data/rdata/firstyear_sens_newplots.Rds"))



# after all, Carolina wants to split it up into first year and rest:
p_sens_firstyear <- map2(post_nested$y, post_nested$data, function(ystring, data) {
  map(firstyear, function(fy) {
    if (fy) {
      p <- temp <- filter(summarised, y == ystring, year == 1)
        ggplot(temp, aes(sensitivity_firstyear, m)) +
          geom_point(color = "black") +
          geom_line(color = "black") +
          scale_y_continuous(breaks = seq(min(temp$m), max(temp$m), 1)) +
          xlab(glue("Early Sensitivity")) + 
          ylab(filter(symp_names, symptom == ystring) %>% .$name) +
          theme_bw(base_size = 20)
    } else {
      p <- temp <- filter(summarised, y == ystring, year != 1)
        mutate(temp, year = factor(year, levels = c("1", "2.5", "4", "5", "6", "7", "8", "10", "11", "12.5", "14"))) %>%
        ggplot(aes(sensitivity_firstyear, m, group = year, color = year)) +
          geom_point() +
          geom_line() +
          scale_y_continuous(breaks = seq(min(temp$m), max(temp$m), 1)) +
          scale_color_manual(values = c(
            "1" = "black", "2.5" = "#a6cee3", "4" = "#b2df8a", "5" = "#33a02c", "6" = "#fb9a99", "7" = "#e31a1c", "8" = "#fdbf6f", "10" = "#ff7f00", "11" = "#cab2d6", "12.5" = "#6a3d9a", "14" = "#1f78b4"), drop = FALSE) +
          xlab(glue("Early Sensitivity")) + 
          ylab(filter(symp_names, symptom == ystring) %>% .$name) +
          labs(color = "Year") +
          theme_bw(base_size = 20)
    }
  })
})


p_sens_firstyear[[5]][[2]] + p_sens_firstyear[[5]][[1]] + ylab("")
save(p_sens_firstyear, file = here::here("data/rdata/firstyear_sens_newplots.Rds"))







# all time sensitivity
sens_all <- map(integers, function(symptom) {
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + sensitivity_all +  (1 | id)"))
  m8_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.)))
    m8 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m8_imp_{m}_{symptom}_sensall.Rds"))
    )
    loo_m8 <- add_criterion(
      m8,
      "loo",
      ##file = here::here(glue("data/rdata/loo_m8_imp_{m}_{symptom}_sensall")),
      moment_match = FALSE
    )
    list(
      m8,
      loo_m8
    )
  })
  
  m8_post <- m8_post <- map_dfr(m8_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m8_h <- map(m8_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "sensitivity_all < 0")
  })
  
  
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + sensitivity_all +  (1 | id) + (sensitivity_all | yearf)"))
  m11_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.))) %>%
      mutate(yearf = as.factor(year))
    m11 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m11_imp_{m}_{symptom}_sensall.Rds"))
    )
    loo_m11 <- add_criterion(
      m11,
      "loo",
      ##file = here::here(glue("data/rdata/loo_m11_imp_{m}_{symptom}_sensall")),
      moment_match = FALSE
    )
    list(
      m11,
      loo_m11
    )
  })
  
  m11_post <- m11_post <- map_dfr(m11_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m11_h <- map(m11_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "sensitivity_all < 0")
  })
  
  
  
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + bf * firstyear + edu + sib + sensitivity_all + (1 | id)"))
  m10_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.)))
    m10 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m10_imp_{m}_{symptom}_sensall.Rds"))
    )
    loo_m10 <- add_criterion(
      m10,
      "loo",
      ##file = here::here(glue("data/rdata/loo_m10_imp_{m}_{symptom}_sensall")),
      moment_match = FALSE
    )
    list(
      m10,
      loo_m10
    )
  })
  
  m10_post <- m10_post <- map_dfr(m10_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  
  m10_h <- map(m10_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "sensitivity_all < 0")
  })
  
  
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + bf * firstyear + edu + sib + sensitivity_all + (1 | id) + (sensitivity_all | yearf)"))
  m12_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.))) %>%
      mutate(yearf = as.factor(year))
    m12 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m12_imp_{m}_{symptom}_sensall.Rds"))
    )
    loo_m12 <- add_criterion(
      m12,
      "loo",
      ##file = here::here(glue("data/rdata/loo_m12_imp_{m}_{symptom}_sensall")),
      moment_match = FALSE
    )
    list(
      m12,
      loo_m12
    )
  })
  
  m12_post <- m12_post <- map_dfr(m12_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  
  m12_h <- map(m12_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "sensitivity_all < 0")
  })
  
  
  
  loo <- loo_compare(m8_models[[1]][[2]], m10_models[[1]][[2]], m11_models[[1]][[2]], m12_models[[1]][[2]])
  winner <- rownames(loo)[1]
  list(
    m8 = m8_models,
    m10 = m10_models,
    loo = loo,
    winner = winner,
    m8_h = m8_h,
    m10_h = m10_h,
    m8_post = m8_post,
    m10_post = m10_post,
    
    m11 = m11_models,
    m12 = m12_models,
    loo = loo,
    m11_h = m11_h,
    m12_h = m12_h,
    m11_post = m11_post,
    m12_post = m12_post
  )
})

# lets first check which model we should interpret
map(sens_all, ~.x$winner)
map2(1:5, sens_all, function(no, mlist) {
  mlist$m11_h
})



# to report effect size i use combined posteriors. For directional hypothesis
b <- map2_dfr(sens_all, 1:5, function(mlist, no) {
  posttemp <- mlist$m11_post
  select(posttemp, b_sensitivity_all) %>%
    summarise(
      m = median(b_sensitivity_all),
      lower = hdi(b_sensitivity_all)[1],
      upper = hdi(b_sensitivity_all)[2],
      lower_d = quantile(b_sensitivity_all, 0),
      upper_d = quantile(b_sensitivity_all, 0.95)
    ) %>%
    mutate(y = integers[no])
})
p <- map2_dfr(1:5, sens_all, function(no, mlist) {
  htemp <- mlist$m11_h
  p <- map_dbl(htemp, function(h) {
    h$hypothesis$Post.Prob
  })
  tibble(y = integers[no], p = median(p))
})
b <- full_join(b, p, by = "y") %>%
  mutate(across(where(is.numeric), function(x) format(round(x, 3), digits = 3)))
b
save(b, file = here::here("data/rdata/b_sensitivity_all.Rds"))



# plot
year <- unique(dlong$year)
post_year <- map_dfr(year, function(yr) {
  # generate outcomes for effect size prediction
  nd <- expand_grid(
    cc = c(0, 1), 
    bf = median(dlongimp[[1]]$bf),
    year = yr,
    edu = c(1, 2),
    sib = c(0, 1),
    sensitivity_all = seq(-2, 2, 1),
    firstyear = ifelse(yr == 1, 1, 0)
  )
  nd$yearf <- factor(nd$year)
  
  post_pred <- map_dfr(1:5, function(no) {
    models <- "m11"
    if (b[b$y == integers[no], "p"] >= 0.3) {
      post_pred <- map_dfr(sens_all[[no]][[models]], function(model) {
        posterior_predict(model[[1]], newdata = nd, re.form = ~(sensitivity_all | yearf)) %>%
          as.data.frame() %>%
          pivot_longer(
            everything(),
            names_to = "var",
            values_to = "values"
          ) %>%
          mutate(
            y = integers[no]
          )
      })
    }
  })
  
  match_names <- tibble(
    var = unique(post_pred$var), 
    sensitivity_all = nd$sensitivity_all,
    firstyear = 0,
    cc = nd$cc, 
    bf = nd$bf,
    year = nd$year,
    yearf = nd$yearf,
    edu = nd$edu,
    sib = nd$sib
  )
  post_pred %>%
    left_join(match_names, by = "var") %>%
    select(-var) 
})


summarised <- post_year %>%
  group_by(sensitivity_all, year, y) %>%
  summarise(
    m = median(values),
    lower = hdi(values)[1],
    upper = hdi(values)[2],
    lower_50 = hdi(values, credMass = 0.68)[1],
    upper_50 = hdi(values, credMass = 0.68)[2]
  )  




# alternative with patchwork
post_nested <- group_by(post_year, y) %>% nest()
p_sens_all <- map2(post_nested$y, post_nested$data, function(ystring, data) {
  p <- ggplot(filter(summarised, y == ystring), aes(sensitivity_all, m, group = year, color = as.factor(year))) +
    geom_point() +
    geom_line() +
    xlab(glue("Sensitivity Throughout Childhood")) + 
    ylab(filter(symp_names, symptom == ystring) %>% .$name) +
    labs(color = "Year") +
    theme_bw(base_size = 20) 
  p
})


save(p_sens_all, file = here::here("data/rdata/sensall_newplots.Rds"))


# after all, Carolina wants to split it up into first year and rest:
p_sens_all <- map2(post_nested$y, post_nested$data, function(ystring, data) {
  map(firstyear, function(fy) {
    if (fy) {
      p <- temp <- filter(summarised, y == ystring, year == 1)
        ggplot(temp, aes(sensitivity_all, m)) +
          geom_point(color = "black") +
          geom_line(color = "black") +
          scale_y_continuous(breaks = seq(min(temp$m), max(temp$m), 1)) +
          xlab(glue("Sensitivity Throughout Childhood")) + 
          ylab(filter(symp_names, symptom == ystring) %>% .$name) +
          theme_bw(base_size = 20)
    } else {
      p <- temp <- filter(summarised, y == ystring, year != 1)
        mutate(temp, year = factor(year, levels = c("1", "2.5", "4", "5", "6", "7", "8", "10", "11", "12.5", "14"))) %>%
        ggplot(aes(sensitivity_all, m, group = year, color = year)) +
          geom_point() +
          geom_line() +
          scale_y_continuous(breaks = seq(min(temp$m), max(temp$m), 1)) +
          scale_color_manual(values = c(
            "1" = "black", "2.5" = "#a6cee3", "4" = "#b2df8a", "5" = "#33a02c", "6" = "#fb9a99", "7" = "#e31a1c", "8" = "#fdbf6f", "10" = "#ff7f00", "11" = "#cab2d6", "12.5" = "#6a3d9a", "14" = "#1f78b4"), drop = FALSE) +
          xlab(glue("Sensitivity Throughout Childhood")) +  
          ylab(filter(symp_names, symptom == ystring) %>% .$name) +
          labs(color = "Year") +
          theme_bw(base_size = 20)
    }
  })
})


p_sens_all[[5]][[2]] + p_sens_all[[5]][[1]] + ylab("")
save(p_sens_all, file = here::here("data/rdata/sensall_newplots.Rds"))




attachment_sec_all <- map(integers, function(symptom) {
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + attachment_security +  (1 | id)"))
  m8_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.)))
    m8 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m8_imp_{m}_{symptom}_attachment_sec_all.Rds"))
    )
    loo_m8 <- add_criterion(
      m8,
      "loo",
      #file = here::here(glue("data/rdata/loo_m8_imp_{m}_{symptom}_attachment_sec__all")),
      moment_match = FALSE
    )
    list(
      m8,
      loo_m8
    )
  })
  
  m8_post <- m8_post <- map_dfr(m8_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  
  m8_h <- map(m8_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "attachment_security > 0")
  })



  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + attachment_security +  (1 | id) + (attachment_security | yearf)"))
  m11_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.))) %>%
              mutate(yearf = as.factor(year))
    m11 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m11_imp_{m}_{symptom}_attachment_sec_all.Rds"))
    )
    loo_m11 <- add_criterion(
      m11,
      "loo",
      #file = here::here(glue("data/rdata/loo_m11_imp_{m}_{symptom}_attachment_sec__all")),
      moment_match = FALSE
    )
    list(
      m11,
      loo_m11
    )
  })
  
  m11_post <- m11_post <- map_dfr(m11_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  
  m11_h <- map(m11_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "attachment_security > 0")
  })

  
  
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + bf * firstyear + edu + sib + attachment_security + (1 | id)"))
  m10_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.)))
    m10 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m10_imp_{m}_{symptom}_attachment_sec_all.Rds"))
    )
    loo_m10 <- add_criterion(
      m10,
      "loo",
      ##file = here::here(glue("data/rdata/loo_m10_imp_{m}_{symptom}_attachment_sec_all")),
      moment_match = FALSE
    )
    list(
      m10,
      loo_m10
    )
    
  })
  
  m10_post <- m10_post <- map_dfr(m10_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m10_h <- map(m10_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "attachment_security > 0")
  })


form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + bf * firstyear + edu + sib + attachment_security + (1 | id) + (attachment_security | yearf)"))
  m12_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.))) %>%
              mutate(yearf = as.factor(year))
    m12 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m12_imp_{m}_{symptom}_attachment_sec_all.Rds"))
    )
    loo_m12 <- add_criterion(
      m12,
      "loo",
      ##file = here::here(glue("data/rdata/loo_m12_imp_{m}_{symptom}_attachment_sec_all")),
      moment_match = FALSE
    )
    list(
      m12,
      loo_m12
    )
    
  })
  
  m12_post <- m12_post <- map_dfr(m12_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m12_h <- map(m12_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "attachment_security > 0")
  })



  
  loo <- loo_compare(m8_models[[1]][[2]], m10_models[[1]][[2]], m11_models[[1]][[2]], m12_models[[1]][[2]])
  winner <- rownames(loo)[1]
  list(
    m8 = m8_models,
    m10 = m10_models,
    loo = loo,
    winner = winner,
    m8_h = m8_h,
    m10_h = m10_h,
    m8_post = m8_post,
    m10_post = m10_post,
    
    m11 = m11_models,
    m12 = m12_models,
    m11_h = m11_h,
    m12_h = m12_h,
    m11_post = m11_post,
    m12_post = m12_post
  )
})

# lets first check which model we should interpret
map(attachment_sec_all, ~.x$winner)
map2(1:5, attachment_sec_all, function(no, mlist) {
  if (no <= 2) {
    mlist$m10_h
  } else {
    mlist$m8_h
  }
})

# to report effect size i use combined posteriors. For directional hypothesis
b <- map2_dfr(attachment_sec_all, 1:5, function(mlist, no) {
  if (no <= 2) {
    posttemp <- mlist$m10_post
  } else {
    posttemp <- mlist$m8_post
  }
  select(posttemp, b_attachment_security) %>%
    summarise(
      m = median(b_attachment_security),
      lower = hdi(b_attachment_security)[1],
      upper = hdi(b_attachment_security)[2],
      lower_d = quantile(b_attachment_security, 0),
      upper_d = quantile(b_attachment_security, 0.95)
    ) %>%
    mutate(y = integers[no])
})
p <- map2_dfr(1:5, attachment_sec_all, function(no, mlist) {
  if (no <= 2) {
    htemp <- mlist$m10_h
  } else {
    htemp <- mlist$m8_h
  }
  p <- map_dbl(htemp, function(h) {
    h$hypothesis$Post.Prob
  })
  tibble(y = integers[no], p = median(p))
})
b <- full_join(b, p, by = "y") %>%
  mutate(across(where(is.numeric), function(x) format(round(x, 3), digits = 3)))
b
save(b, file = here::here("data/rdata/b_attachments_security.Rds"))

# wasnt significant, therefore no plot





# attachment dummy
ad_all <- map(integers, function(symptom) {
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + attach_dummy +  (1 | id)"))
  m8_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.)))
    m8 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m8_imp_{m}_{symptom}_ad_all.Rds"))
    )
    loo_m8 <- add_criterion(
      m8,
      "loo",
      #file = here::here(glue("data/rdata/loo_m8_imp_{m}_{symptom}_ad_all")),
      moment_match = FALSE
    )
    list(
      m8,
      loo_m8
    )
  })
  
  m8_post <- m8_post <- map_dfr(m8_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m8_h <- map(m8_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "Intercept + attach_dummy1 < Intercept")
  })




  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + attach_dummy +  (1 | id) + (attach_dummy | yearf)"))
  m11_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.))) %>%
              mutate(yearf = as.factor(year))
    m11 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m11_imp_{m}_{symptom}_ad_all.Rds"))
    )
    loo_m11 <- add_criterion(
      m11,
      "loo",
      #file = here::here(glue("data/rdata/loo_m11_imp_{m}_{symptom}_ad_all")),
      moment_match = FALSE
    )
    list(
      m11,
      loo_m11
    )
  })
  
  m11_post <- m11_post <- map_dfr(m11_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m11_h <- map(m11_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "Intercept + attach_dummy1 < Intercept")
  })
  
  
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + bf * firstyear + edu + sib + attach_dummy + (1 | id)"))
  m10_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.)))
    m10 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m10_imp_{m}_{symptom}_ad_all.Rds"))
    )
    loo_m10 <- add_criterion(
      m10,
      "loo",
      ##file = here::here(glue("data/rdata/loo_m10_imp_{m}_{symptom}_ad_all")),
      moment_match = FALSE
    )
    list(
      m10,
      loo_m10
    )
    
  })
  
  m10_post <- m10_post <- map_dfr(m10_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m10_h <- map(m10_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "Intercept + attach_dummy1 < Intercept")
  })


  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + bf * firstyear + edu + sib + attach_dummy + (1 | id) + (attach_dummy | yearf)"))
  m12_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.))) %>% 
              mutate(yearf = as.factor(year))
    m12 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m12_imp_{m}_{symptom}_ad_all.Rds"))
    )
    loo_m12 <- add_criterion(
      m12,
      "loo",
      ##file = here::here(glue("data/rdata/loo_m12_imp_{m}_{symptom}_ad_all")),
      moment_match = FALSE
    )
    list(
      m12,
      loo_m12
    )
    
  })
  
  m12_post <- m12_post <- map_dfr(m12_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m12_h <- map(m12_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "Intercept + attach_dummy1 < Intercept")
  })

  
  loo <- loo_compare(m8_models[[1]][[2]], m10_models[[1]][[2]], m11_models[[1]][[2]], m12_models[[1]][[2]])
  winner <- rownames(loo)[1]
  list(
    m8 = m8_models,
    m10 = m10_models,
    loo = loo,
    winner = winner,
    m8_h = m8_h,
    m10_h = m10_h,
    m8_post = m8_post,
    m10_post = m10_post,
    
    m11 = m11_models,
    m12 = m12_models,
    m11_h = m11_h,
    m12_h = m12_h,
    m11_post = m11_post,
    m12_post = m12_post
  )
})
# lets first check which model we should interpret
map(ad_all, ~.x$winner)

map2(1:5, ad_all, function(no, mlist) {
  if (no %in% c(1, 3, 5)) {
    mlist$m12_h
  } else if  (no == 2) {
    mlist$m11_h
  } else {
    mlist$m8_h
  }
})

# to report effect size i use combined posteriors. For directional hypothesis
b <- map2_dfr(ad_all, 1:5, function(mlist, no) {
  if (no %in% c(1, 3, 5)) {
    posttemp <- mlist$m12_post
  } else if  (no == 2) {
    posttemp <- mlist$m11_post
  } else {
    posttemp <- mlist$m8_post
  }
  
  select(posttemp, b_attach_dummy1) %>%
    summarise(
      m = median(b_attach_dummy1),
      lower = hdi(b_attach_dummy1)[1],
      upper = hdi(b_attach_dummy1)[2],
      lower_d = quantile(b_attach_dummy1, 0),
      upper_d = quantile(b_attach_dummy1, 0.95)
    ) %>%
    mutate(y = integers[no])
})

p <- map2_dfr(1:5, ad_all, function(no, mlist) {
  if (no %in% c(1, 3, 5)) {
    htemp <- mlist$m12_h
  } else if  (no == 2) {
    htemp <- mlist$m11_h
  } else {
    htemp <- mlist$m8_h
  }

  p <- map_dbl(htemp, function(h) {
    h$hypothesis$Post.Prob
  })
  tibble(y = integers[no], p = median(p))
})
b <- full_join(b, p, by = "y") %>%
  mutate(across(where(is.numeric), function(x) format(round(x, 3), digits = 3)))
b
save(b, file = here::here("data/rdata/b_attach_dummy.Rds"))


# Attachment has even stronger effects and across all symptoms except skin symps


year <- unique(dlong$year)
post_year <- map_dfr(year, function(yr) {
  # generate outcomes for effect size prediction
  nd <- expand_grid(
    cc = c(0, 1), 
    bf = median(dlongimp[[1]]$bf),
    year = yr,
    edu = c(1, 2),
    sib = c(0, 1),
    # sensitivity_all = seq(-2, 2, 1),
    #sensitivity_firstyear = seq(-2, 2, 1),
    attach_dummy = c(0, 1),
    #attach_cat = 1:4
    firstyear = ifelse(yr == 1, 1, 0)
  )
  nd$yearf <- factor(nd$year)
  
  post_pred <- map_dfr(1:5, function(no) {
    if (no %in% c(1, 3, 5)) {
      models <- "m12"
    } else if  (no == 2) {
      models <- "m11"
    } else {
      models <- "m8"
    }
    if (b[b$y == integers[no], "p"] >= 0.3) {
      post_pred <- map_dfr(ad_all[[no]][[models]], function(model) {
        posterior_predict(model[[1]], newdata = nd, re.form = ~(attach_dummy | yearf)) %>%
          as.data.frame() %>%
          pivot_longer(
            everything(),
            names_to = "var",
            values_to = "values"
          ) %>%
          mutate(
            #values = values - median(values),
            y = integers[no]
          )
      })
    }
  })
  
  match_names <- tibble(
    var = unique(post_pred$var), 
    #sensitivity_all = nd$sensitivity_all,
    firstyear = 0,
    cc = nd$cc, 
    bf = nd$bf,
    year = nd$year,
    yearf = nd$yearf,
    edu = nd$edu,
    sib = nd$sib,
    #sens_factor = factor(nd$sensitivity_firstyear)
    #sensitivity_all = seq(-2, 2, 1),
    attach_dummy = nd$attach_dummy,
    #attach_cat = 1:4
  )
  post_pred %>%
    left_join(match_names, by = "var") %>%
    select(-var) 
})


summarised <- post_year %>%
  group_by(attach_dummy, year, y) %>%
  summarise(
    m = median(values),
    lower = hdi(values)[1],
    upper = hdi(values)[2],
    lower_50 = hdi(values, credMass = 0.68)[1],
    upper_50 = hdi(values, credMass = 0.68)[2]
  )  



post_nested <- group_by(post_year, y) %>% nest()

# alternative with patchwork
p_attach_dummy <- map2(post_nested$y, post_nested$data, function(ystring, data) {
  p <- ggplot(filter(summarised, y == ystring), aes(as.factor(attach_dummy), m, group = year, color = as.factor(year))) +
    geom_point() +
    geom_line() +
    # geom_hline(aes(yintercept = 0), linetype = "dashed") +
    #geom_jitter(alpha = 0.005) +
    #geom_linerange(data = filter(summarised, y == ystring), aes(y = m, ymin = lower, ymax = upper), linewidth = 1.25, color = "#bdbdbd") +
    #geom_pointrange(data = filter(summarised, y == ystring), aes(x = sensitivity_firstyear, y = m, ymin = lower_50, ymax = upper_50), size = 1.5, linewidth = 1.8, color = "#636363") +
    #scale_y_continuous(breaks = c(seq(min(data$values), -3, 3), 0, seq(3, max(data$values), 3))) +
    #scale_y_continuous(breaks = seq(-100, 100, 5)) +
    xlab(glue("Attachment")) + 
    scale_x_discrete(labels = c("Secure", "Insecure")) +
    ylab(filter(symp_names, symptom == ystring) %>% .$name) +
    labs(color = "Year") +
    theme_bw(base_size = 20) 
  p
})

# save(p_attach_dummy, file = here::here("data/rdata/attach_dummy_newplots.Rds"))
# p_attach_dummy[[5]]
#1, 2, 3, 5 are rs



# after all, Carolina wants to split it up into first year and rest:
p_attach_dummy <- map2(post_nested$y, post_nested$data, function(ystring, data) {
  map(firstyear, function(fy) {
    if (fy) {
      p <- temp <- filter(summarised, y == ystring, year == 1)
        ggplot(temp, aes(factor(attach_dummy), m, group = 1)) +
          geom_point(color = "black") +
          geom_line(color = "black") +
          xlab(glue("Attachment")) +  
          scale_x_discrete(labels = c("Secure", "Insecure")) +
          ylab(filter(symp_names, symptom == ystring) %>% .$name) +
          theme_bw(base_size = 20)
    } else {
      p <- temp <- filter(summarised, y == ystring, year != 1)
        mutate(temp, year = factor(year, levels = c("1", "2.5", "4", "5", "6", "7", "8", "10", "11", "12.5", "14"))) %>%
        ggplot(aes(as.factor(attach_dummy), m, group = year, color = year)) +
          geom_point() +
          geom_line() +
          scale_y_continuous(breaks = seq(min(temp$m), max(temp$m), 1)) +
          scale_color_manual(values = c(
            "1" = "black", "2.5" = "#a6cee3", "4" = "#b2df8a", "5" = "#33a02c", "6" = "#fb9a99", "7" = "#e31a1c", "8" = "#fdbf6f", "10" = "#ff7f00", "11" = "#cab2d6", "12.5" = "#6a3d9a", "14" = "#1f78b4"), drop = FALSE) +
          xlab(glue("Attachment")) +   
          scale_x_discrete(labels = c("Secure", "Insecure")) +
          ylab(filter(symp_names, symptom == ystring) %>% .$name) +
          labs(color = "Year") +
          theme_bw(base_size = 20)
    }
  })
})


p_attach_dummy[[5]][[2]] + p_attach_dummy[[5]][[1]] + ylab("")
save(p_attach_dummy, file = here::here("data/rdata/attach_dummy_newplots.Rds"))



# attachment cat
# Stefania asked me to calculate contrasts for all comparisons so i use a grid
# for the epred fun

nd <- expand_grid(
  attach_cat = 1:4,
  firstyear = 1,
  cc = 1, 
  bf = median(dlongimp[[1]]$bf),
  year = 1,
  edu = 1,
  sib = 1,
  sex = 1
)

ac_all <- map(integers, function(symptom) {
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + attach_cat +  (1 | id)"))
  m8_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.)))
    m8 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m8_imp_{m}_{symptom}_ac_all.Rds"))
    )
    loo_m8 <- add_criterion(
      m8,
      "loo",
      #file = here::here(glue("data/rdata/loo_m8_imp_{m}_{symptom}_ac_all")),
      moment_match = FALSE
    )
    list(
      m8,
      loo_m8
    )
  })
  
  m8_post <- m8_post <- map_dfr(m8_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m8_post2 <- map_dfr(m8_models, function(mlist) {
    posterior_epred(mlist[[1]], newdata = nd, re_formula = NA) %>%
      as.data.frame()
  })
  
  m8_h <- m8_post2 %>%
    mutate(
      a1_vs_a2 = V1 - V2,
      a1_vs_a3 = V1 - V3,
      a1_vs_a4 = V1 - V3,
      a2_vs_a3 = V2 - V3,
      a2_vs_a4 = V2 - V4,
      a3_vs_a4 = V3 - V4
    ) %>%
    select(
      a1_vs_a2,
      a1_vs_a3,
      a1_vs_a4,
      a2_vs_a3,
      a2_vs_a4,
      a3_vs_a4
    ) %>%
    pivot_longer(
      everything(),
      names_to = "contrast",
      values_to = "value"
    ) %>%
    group_by(contrast) %>%
    summarise(
      m = median(value),
      lower = hdi(value)[1],
      upper = hdi(value)[2]
    )
  
  
  
  
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + bf * firstyear + edu + sib + attach_cat + (1 | id)"))
  m10_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.)))
    m10 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m10_imp_{m}_{symptom}_ac_all.Rds"))
    )
    loo_m10 <- add_criterion(
      m10,
      "loo",
      ##file = here::here(glue("data/rdata/loo_m10_imp_{m}_{symptom}_ac_all")),
      moment_match = FALSE
    )
    list(
      m10,
      loo_m10
    )
    
  })
  
  m10_post <- m10_post <- map_dfr(m10_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m10_post2 <- map_dfr(m10_models, function(mlist) {
    posterior_epred(mlist[[1]], newdata = nd, re_formula = NA) %>%
      as.data.frame()
  })
  m10_h <- mutate(m10_post2, 
                  a1_vs_a2 = V1 - V2,
                  a1_vs_a3 = V1 - V3,
                  a1_vs_a4 = V1 - V3,
                  a2_vs_a3 = V2 - V3,
                  a2_vs_a4 = V2 - V4,
                  a3_vs_a4 = V3 - V4
  ) %>%
    select(
      a1_vs_a2,
      a1_vs_a3,
      a1_vs_a4,
      a2_vs_a3,
      a2_vs_a4,
      a3_vs_a4
    ) %>%
    pivot_longer(
      everything(),
      names_to = "contrast",
      values_to = "value"
    ) %>%
    group_by(contrast) %>%
    summarise(
      m = median(value),
      lower = hdi(value)[1],
      upper = hdi(value)[2]
    )
  
  loo <- loo_compare(m8_models[[1]][[2]], m10_models[[1]][[2]])
  winner <- rownames(loo)[1]
  list(
    m8 = m8_models,
    m10 = m10_models,
    loo = loo,
    winner = winner,
    m8_h = m8_h,
    m10_h = m10_h,
    m8_post = m8_post,
    m10_post = m10_post,
    m8_post2 = m8_post2,
    m10_post2 = m10_post2
  )
})



# lets first check which model we should interpret
map(ac_all, ~.x$winner)
map2(1:5, ac_all, function(no, mlist) {
  mlist$m8_h
})


cat_plots <- map(ac_all, function(mlist) {
  summarised <- mlist[["m8_post2"]] %>%
    pivot_longer(
      everything(),
      names_to = "category",
      values_to = "values"
    ) %>% 
    group_by(category) %>%
    summarise(
      m = median(values),
      lower = hdi(values)[1],
      upper = hdi(values)[2],
      lower_50 = hdi(values, 0.68)[1],
      upper_50 = hdi(values, 0.68)[2]
    )
  mlist[["m8_post2"]] %>%
    pivot_longer(
      everything(),
      names_to = "category",
      values_to = "values"
    ) %>%
    ggplot(aes(category, values)) +
    # geom_jitter(width = 0.2, alpha = 0.01) +
    #geom_hline(aes(yintercept = 0), linetype = "dashed") +
    geom_linerange(data = summarised, aes(y = m, ymin = lower, ymax = upper), linewidth = 1.25, color = "#bdbdbd") +
    geom_pointrange(data = summarised, aes(y = m, ymin = lower_50, ymax = upper_50), size = 1.5, linewidth = 1.8, color = "#636363") +
    
    #geom_errorbar(data = summarised, aes(ymin = lower, ymax = upper, y = m, color = "red")) +
    #geom_point(data = summarised, aes(y = m), color = "red") +
    theme_bw(base_size = 20)
})
cat_plots[[4]]


firstyear <- c(0, 1)
ps <- map(firstyear, function(fy) {
  # generate outcomes for effect size prediction
  nd <- expand_grid(
    cc = c(0, 1), 
    #bf = median(dlongimp[[1]]$bf),
    year = if (fy == 0) c(1, 2.5, 4, 5, 6, 7, 8, 10, 11, 12.5, 14) else 1,
    #edu = c(1, 2),
    #sib = c(0, 1),
    #sensitivity_all = seq(-2, 2, 1),
    #sensitivity_firstyear = seq(-2, 2, 1),
    #attach_dummy = c(0, 1),
    attach_cat = 1:4,
    firstyear = fy
  )
  
  
  post_pred <- map_dfr(1:5, function(no) {
      models <- "m8"
      post_pred <- map_dfr(ac_all[[no]][[models]], function(model) {
        posterior_predict(model[[1]], newdata = nd, re.form = NA) %>%
          as.data.frame() %>%
          pivot_longer(
            everything(),
            names_to = "var",
            values_to = "values"
          ) %>%
          mutate(
            #values = values - median(values),
            y = integers[no]
          )
      })
  })
  
  match_names <- tibble(
    var = unique(post_pred$var), 
    #sensitivity_firstyear = nd$sensitivity_firstyear,
    firstyear = nd$firstyear,
    cc = nd$cc, 
    year = nd$year,
    #sens_factor = factor(nd$sensitivity_firstyear)
    #sensitivity_all = seq(-2, 2, 1),
    #attach_dummy = nd$attach_dummy
    attach_cat = nd$attach_cat
  )
  
  summarised <- post_pred %>%
    left_join(match_names, by = "var") %>%
    group_by(attach_cat, y) %>%
    summarise(
      m = median(values),
      lower = hdi(values)[1],
      upper = hdi(values)[2],
      lower_50 = hdi(values, credMass = 0.68)[1],
      upper_50 = hdi(values, credMass = 0.68)[2]
    )  
  
  
  post_pred_centered <- post_pred %>%
    left_join(match_names, by = "var") %>%
    select(-var) 
  
  # alternative with patchwork
  post_nested <- group_by(post_pred_centered, y) %>% nest()
  p_attach_cat <- map2(post_nested$y, post_nested$data, function(ystring, data) {
    p <- ggplot(data, aes(as.factor(attach_cat), values)) +
      # geom_hline(aes(yintercept = 0), linetype = "dashed") +
      #geom_jitter(alpha = 0.005) +
      geom_linerange(data = filter(summarised, y == ystring), aes(y = m, ymin = lower, ymax = upper), linewidth = 1.25, color = "#bdbdbd") +
      geom_pointrange(data = filter(summarised, y == ystring), aes(x = as.factor(attach_cat), y = m, ymin = lower_50, ymax = upper_50), size = 1.5, linewidth = 1.8, color = "#636363") +
      #scale_y_continuous(breaks = c(seq(min(data$values), -3, 3), 0, seq(3, max(data$values), 3))) +
      scale_y_continuous(breaks = seq(-100, 100, 5)) +
      scale_x_discrete(labels = c("Avoidant", "Secure", "Resistant", "Disorganized")) +
      xlab("Attachment") +
      ylab(filter(symp_names, symptom == ystring) %>% .$name) +
      theme_bw(base_size = 20) 
    p
  })
  
  p_attach_cat
})
ps








###############################################################################
######################    Explore Interaction Models   ########################
###############################################################################

# finally we want to see if the relationship depends on whether mother are sens
# or whether infants are insecure, respectively


# 1st year sensitivity
firstyear_sens_attach <- map(integers, function(symptom) {
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + sensitivity_firstyear * attach_dummy +  (1 | id)"))
  m8_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.)))
    m8 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m8_imp_{m}_{symptom}_interaction.Rds"))
    )
    loo_m8 <- add_criterion(
      m8,
      "loo",
      #file = here::here(glue("data/rdata/loo_m8_imp_{m}_{symptom}_interaction")),
      moment_match = FALSE
    )
    list(
      m8,
      loo_m8
    )
  })
  
  m8_post <- m8_post <- map_dfr(m8_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m8_h <- map(m8_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "sensitivity_firstyear:attach_dummy1 > 0")
  })



  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + sensitivity_firstyear * attach_dummy +  (1 | id) + (sensitivity_firstyear * attach_dummy | yearf)"))
  m11_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.))) %>%
              mutate(yearf = as.factor(year))
    m11 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m11_imp_{m}_{symptom}_interaction.Rds"))
    )
    loo_m11 <- add_criterion(
      m11,
      "loo",
      #file = here::here(glue("data/rdata/loo_m11_imp_{m}_{symptom}_interaction")),
      moment_match = FALSE
    )
    list(
      m11,
      loo_m11
    )
  })
  
  m11_post <- m11_post <- map_dfr(m11_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m11_h <- map(m11_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "sensitivity_firstyear:attach_dummy1 > 0")
  })
  
  
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + bf * firstyear + edu + sib + sensitivity_firstyear * attach_dummy + (1 | id)"))
  m10_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.)))
    m10 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m10_imp_{m}_{symptom}_interaction.Rds"))
    )
    loo_m10 <- add_criterion(
      m10,
      "loo",
      ##file = here::here(glue("data/rdata/loo_m10_imp_{m}_{symptom}_interaction")),
      moment_match = FALSE
    )
    list(
      m10,
      loo_m10
    )
  })
  
  m10_post <- m10_post <- map_dfr(m10_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m10_h <- map(m10_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "sensitivity_firstyear:attach_dummy1 > 0")
  })




  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + bf * firstyear + edu + sib + sensitivity_firstyear * attach_dummy + (1 | id) + (sensitivity_firstyear * attach_dummy | yearf)"))
  m12_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.))) %>% 
              mutate(yearf = as.factor(year))
    m12 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m12_imp_{m}_{symptom}_interaction.Rds"))
    )
    loo_m12 <- add_criterion(
      m12,
      "loo",
      ##file = here::here(glue("data/rdata/loo_m12_imp_{m}_{symptom}_interaction")),
      moment_match = FALSE
    )
    list(
      m12,
      loo_m12
    )
  })
  
  m12_post <- m12_post <- map_dfr(m12_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m12_h <- map(m12_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "sensitivity_firstyear:attach_dummy1 > 0")
  })


  
  loo <- loo_compare(m8_models[[1]][[2]], m10_models[[1]][[2]], m11_models[[1]][[2]], m12_models[[1]][[2]])
  winner <- rownames(loo)[1]
  list(
    m8 = m8_models,
    m10 = m10_models,
    loo = loo,
    winner = winner,
    m8_h = m8_h,
    m10_h = m10_h,
    m8_post = m8_post,
    m10_post = m10_post,
    
    m11 = m11_models,
    m12 = m12_models,
    m11_h = m11_h,
    m12_h = m12_h,
    m11_post = m11_post,
    m12_post = m12_post
  )
})

# lets first check which model we should interpret
map(firstyear_sens_attach, ~.x$winner)
map2(1:5, firstyear_sens_attach, function(no, mlist) {
  if (no == 1) {
    mlist$m12_h
  } else {
    mlist$m11_h
  }
})




# to report effect size i use combined posteriors. For directional hypothesis
b <- map2_dfr(firstyear_sens_attach, 1:5, function(mlist, no) {
  if (no == 1) {
    posttemp <- mlist$m12_post
  } else {
    posttemp <- mlist$m11_post
  }
  select(posttemp, `b_sensitivity_firstyear:attach_dummy1`) %>%
    summarise(
      m = median(`b_sensitivity_firstyear:attach_dummy1`),
      lower = hdi(`b_sensitivity_firstyear:attach_dummy1`)[1],
      upper = hdi(`b_sensitivity_firstyear:attach_dummy1`)[2],
      lower_d = quantile(`b_sensitivity_firstyear:attach_dummy1`, 0.025),
      upper_d = quantile(`b_sensitivity_firstyear:attach_dummy1`, 0.975),
      pp = mean(`b_sensitivity_firstyear:attach_dummy1` > 0)
    ) %>%
    mutate(y = integers[no])
})
p <- map2_dfr(1:5, firstyear_sens_attach, function(no, mlist) {
  if (no == 1) {
    htemp <- mlist$m12_h
  } else {
    htemp <- mlist$m11_h
  }
  p <- map_dbl(htemp, function(h) {
    h$hypothesis$Post.Prob
  })
  tibble(y = integers[no], p = median(p))
})
b <- full_join(b, p, by = "y") %>%
  mutate(across(where(is.numeric), function(x) format(round(x, 3), digits = 3)))
b
#load(here::here("data/rdata/b_firstyear_sens_attach.Rds"))
save(b, file = here::here("data/rdata/b_firstyear_sens_attach.Rds"))





year <- unique(dlong$year)
post_year <- map_dfr(year, function(yr) {
  # generate outcomes for effect size prediction
  nd <- expand_grid(
    cc = c(0, 1), 
    bf = median(dlongimp[[1]]$bf),
    year = yr,
    edu = c(1, 2),
    sib = c(0, 1),
    # sensitivity_all = seq(-2, 2, 1),
    sensitivity_firstyear = seq(-2, 2, 1),
    attach_dummy = c(0, 1),
    #attach_cat = 1:4
    firstyear = ifelse(yr == 1, 1, 0)
  )
  nd$yearf <- factor(nd$year)
  
  post_pred <- map_dfr(1:5, function(no) {
    if (no ==1) {
      models <- "m12"
    } else {
      models <- "m11"
    }
    if (b[b$y == integers[no], "p"] >= 0.3) {
      post_pred <- map_dfr(firstyear_sens_attach[[no]][[models]], function(model) {
        posterior_predict(model[[1]], newdata = nd, re.form = ~(sensitivity_firstyear * attach_dummy | yearf)) %>%
          as.data.frame() %>%
          pivot_longer(
            everything(),
            names_to = "var",
            values_to = "values"
          ) %>%
          mutate(
            #values = values - median(values),
            y = integers[no]
          )
      })
    }
  })
  
  match_names <- tibble(
    var = unique(post_pred$var), 
    sensitivity_firstyear = nd$sensitivity_firstyear,
    firstyear = 0,
    cc = nd$cc, 
    bf = nd$bf,
    year = nd$year,
    yearf = nd$yearf,
    edu = nd$edu,
    sib = nd$sib,
    attach_dummy = nd$attach_dummy
  )
  post_pred %>%
    left_join(match_names, by = "var") %>%
    select(-var) 
})


summarised <- post_year %>%
  group_by(attach_dummy, sensitivity_firstyear, year, y) %>%
  summarise(
    m = median(values),
    lower = hdi(values)[1],
    upper = hdi(values)[2],
    lower_50 = hdi(values, credMass = 0.68)[1],
    upper_50 = hdi(values, credMass = 0.68)[2]
  )  


# alternative with patchwork
attach_dummy_names <- c(
    `0` = "Early Sensitivity (Secure)",
    `1` = "Early Sensitivity (Insecure)"
  )
post_nested <- group_by(post_year, y) %>% nest()
p_attach_sens <- map2(post_nested$y, post_nested$data, function(ystring, data) {
  p <- ggplot(filter(summarised, y == ystring), aes(sensitivity_firstyear, m, group = year, color = as.factor(year))) +
    geom_point() +
    geom_line() +
    xlab("") + 
    facet_wrap(~attach_dummy, strip.position = "bottom", labeller = as_labeller(attach_dummy_names)) +
    ylab(filter(symp_names, symptom == ystring) %>% .$name) +
    labs(color = "Year") +
    theme_bw(base_size = 20) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.x = element_text(size = 14)
    ) 
  p
})


# save(p_attach_sens, file = here::here("data/rdata/attach_sens_newplots.Rds"))
# p_attach_sens[[5]]
#1, 2, 3, 5 are rs


# decide which plots to put together, make sure label includes first year, what to do about scaling y?

# after all, Carolina wants to split it up into first year and rest:
p_attach_sens <- map2(post_nested$y, post_nested$data, function(ystring, data) {
  map(firstyear, function(fy) {
    map(unique(summarised$attach_dummy), function(attachment_type) {
          if (fy) {
      temp <- filter(summarised, y == ystring, year == 1, attach_dummy == attachment_type)
      p <- mutate(temp, year = factor(year, levels = c("1", "2.5", "4", "5", "6", "7", "8", "10", "11", "12.5", "14"))) %>%
        ggplot(aes(sensitivity_firstyear, m)) +
        geom_point(color = "black") +
        geom_line(color = "black") +
        scale_y_continuous(breaks = seq(min(temp$m), max(temp$m), 1)) +
        xlab("Early Sensitivity") + 
        ylab(filter(symp_names, symptom == ystring) %>% .$name) +
        theme_bw(base_size = 20) +
          ggtitle(glue("{attachment_type}_{fy}"))
    } else {
        temp <- filter(summarised, y == ystring, year != 1, attach_dummy == attachment_type)
        p <- mutate(temp, year = factor(year, levels = c("1", "2.5", "4", "5", "6", "7", "8", "10", "11", "12.5", "14"))) %>%
          ggplot(aes(sensitivity_firstyear, m, group = year, color = as.factor(year))) +
          geom_point() +
          geom_line() +
          scale_y_continuous(breaks = seq(min(temp$m), max(temp$m), 1)) +
          scale_color_manual(values = c(
            "1" = "black", "2.5" = "#a6cee3", "4" = "#b2df8a", "5" = "#33a02c", "6" = "#fb9a99", "7" = "#e31a1c", "8" = "#fdbf6f", "10" = "#ff7f00", "11" = "#cab2d6", "12.5" = "#6a3d9a", "14" = "#1f78b4"), drop = FALSE) +
          xlab("Early Sensitivity") + 
          ylab(filter(symp_names, symptom == ystring) %>% .$name) +
          labs(color = "Year") +
          theme_bw(base_size = 20) +
          ggtitle(glue("{attachment_type}_{fy}"))
    }
    })
  })
})





p_attach_sens[[5]][[2]][[1]] + p_attach_sens[[5]][[2]][[2]] +   ylab("")
p_attach_sens[[5]][[1]][[1]] + theme(legend.position = "none") + p_attach_sens[[5]][[1]][[2]] +   ylab("")
save(p_attach_sens, file = here::here("data/rdata/attach_sens_newplots.Rds"))




# all time sensitivity
sens_all_attach <- map(integers, function(symptom) {
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + sensitivity_all * attach_dummy +  (1 | id)"))
  m8_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.)))
    m8 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m8_imp_{m}_{symptom}_sensall_attach.Rds"))
    )
    loo_m8 <- add_criterion(
      m8,
      "loo",
      # #file = here::here(glue("data/rdata/loo_m8_imp_{m}_{symptom}_sensall_attach")),
      moment_match = FALSE
    )
    list(
      m8,
      loo_m8
    )
  })
  
  m8_post <- m8_post <- map_dfr(m8_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m8_h <- map(m8_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "sensitivity_all:attach_dummy1 < 0")
  })



  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + sensitivity_all * attach_dummy +  (1 | id) + (sensitivity_all * attach_dummy | yearf)"))
  m11_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.))) %>% 
              mutate(yearf = as.factor(year))
    m11 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m11_imp_{m}_{symptom}_sensall_attach.Rds"))
    )
    loo_m11 <- add_criterion(
      m11,
      "loo",
      # #file = here::here(glue("data/rdata/loo_m11_imp_{m}_{symptom}_sensall_attach")),
      moment_match = FALSE
    )
    list(
      m11,
      loo_m11
    )
  })
  
  m11_post <- m11_post <- map_dfr(m11_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m11_h <- map(m11_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "sensitivity_all:attach_dummy1 < 0")
  })
  
  
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + bf * firstyear + edu + sib + sensitivity_all * attach_dummy + (1 | id)"))
  m10_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.)))
    m10 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m10_imp_{m}_{symptom}_sensall_attach.Rds"))
    )
    loo_m10 <- add_criterion(
      m10,
      "loo",
      # #file = here::here(glue("data/rdata/loo_m10_imp_{m}_{symptom}_sensall_attach")),
      moment_match = FALSE
    )
    list(
      m10,
      loo_m10
    )
    
  })
  
  m10_post <- m10_post <- map_dfr(m10_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m10_h <- map(m10_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "sensitivity_all:attach_dummy1 > 0")
  })
  


  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + bf * firstyear + edu + sib + sensitivity_all * attach_dummy + (1 | id) + (sensitivity_all * attach_dummy | yearf)"))
  m12_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.))) %>% 
              mutate(yearf = as.factor(year))
    m12 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m12_imp_{m}_{symptom}_sensall_attach.Rds"))
    )
    loo_m12 <- add_criterion(
      m12,
      "loo",
      # #file = here::here(glue("data/rdata/loo_m12_imp_{m}_{symptom}_sensall_attach")),
      moment_match = FALSE
    )
    list(
      m12,
      loo_m12
    )
    
  })
  
  m12_post <- m12_post <- map_dfr(m12_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m12_h <- map(m12_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "sensitivity_all:attach_dummy1 > 0")
  })



  loo <- loo_compare(m8_models[[1]][[2]], m10_models[[1]][[2]], m11_models[[1]][[2]], m12_models[[1]][[2]])
  winner <- rownames(loo)[1]
  list(
    m8 = m8_models,
    m10 = m10_models,
    loo = loo,
    winner = winner,
    m8_h = m8_h,
    m10_h = m10_h,
    m8_post = m8_post,
    m10_post = m10_post,
    
    m11 = m11_models,
    m12 = m12_models,
    m11_h = m11_h,
    m12_h = m12_h,
    m11_post = m11_post,
    m12_post = m12_post
  )
})

# lets first check which model we should interpret
map(sens_all_attach, ~.x$winner)
map2(1:5, sens_all_attach, function(no, mlist) {
  if (no != 5) {
    mlist$m11_h
  } else {
    mlist$m8_h
  }
})



# to report effect size i use combined posteriors. For directional hypothesis
b <- map2_dfr(sens_all_attach, 1:5, function(mlist, no) {
  if (no != 5) {
    posttemp <- mlist$m11_post
  } else {
    posttemp <- mlist$m8_post
  }
  select(posttemp, `b_sensitivity_all:attach_dummy1`) %>%
    summarise(
      m = median(`b_sensitivity_all:attach_dummy1`),
      lower = hdi(`b_sensitivity_all:attach_dummy1`)[1],
      upper = hdi(`b_sensitivity_all:attach_dummy1`)[2],
      lower_d = quantile(`b_sensitivity_all:attach_dummy1`, 0.025),
      upper_d = quantile(`b_sensitivity_all:attach_dummy1`, 0.975)
    ) %>%
    mutate(y = integers[no])
})
p <- map2_dfr(1:5, sens_all_attach, function(no, mlist) {
  if (no != 5) {
    htemp <- mlist$m11_h
  } else {
    htemp <- mlist$m8_h
  }
  p <- map_dbl(htemp, function(h) {
    h$hypothesis$Post.Prob
  })
  tibble(y = integers[no], p = median(p))
})
b <- full_join(b, p, by = "y") %>%
  mutate(across(where(is.numeric), function(x) format(round(x, 3), digits = 3)))
b

save(b, file = here::here("data/rdata/b_firstyear_sens_all_attach.Rds"))




year <- unique(dlong$year)
post_year <- map_dfr(year, function(yr) {
  # generate outcomes for effect size prediction
  nd <- expand_grid(
    cc = c(0, 1), 
    bf = median(dlongimp[[1]]$bf),
    year = yr,
    edu = c(1, 2),
    sib = c(0, 1),
    sensitivity_all = seq(-2, 2, 1),
    # sensitivity_firstyear = seq(-2, 2, 1),
    attach_dummy = c(0, 1),
    #attach_cat = 1:4
    firstyear = ifelse(yr == 1, 1, 0)
  )
  nd$yearf <- factor(nd$year)
  
  post_pred <- map_dfr(1:5, function(no) {
    if (no != 5) {
      models <- "m11"
    } else {
      models <- "m8"
    }
    if (b[b$y == integers[no], "p"] >= 0.3) {
      post_pred <- map_dfr(sens_all_attach[[no]][[models]], function(model) {
        posterior_predict(model[[1]], newdata = nd, re.form = ~(sensitivity_all * attach_dummy | yearf)) %>%
          as.data.frame() %>%
          pivot_longer(
            everything(),
            names_to = "var",
            values_to = "values"
          ) %>%
          mutate(
            #values = values - median(values),
            y = integers[no]
          )
      })
    }
  })
  
  match_names <- tibble(
    var = unique(post_pred$var), 
    # sensitivity_firstyear = nd$sensitivity_firstyear,
    firstyear = 0,
    cc = nd$cc, 
    bf = nd$bf,
    year = nd$year,
    yearf = nd$yearf,
    edu = nd$edu,
    sib = nd$sib,
    #sens_factor = factor(nd$sensitivity_firstyear)
    sensitivity_all = nd$sensitivity_all,
    attach_dummy = nd$attach_dummy,
    #attach_cat = 1:4
  )
  post_pred %>%
    left_join(match_names, by = "var") %>%
    select(-var) 
})


summarised <- post_year %>%
  group_by(attach_dummy, sensitivity_all, year, y) %>%
  summarise(
    m = median(values),
    lower = hdi(values)[1],
    upper = hdi(values)[2],
    lower_50 = hdi(values, credMass = 0.68)[1],
    upper_50 = hdi(values, credMass = 0.68)[2]
  )  


# alternative with patchwork
attach_dummy_names <- c(
  `0` = "Sensitivity Throughout Childhood (Secure)",
  `1` = "Sensitivity Throughout Childhood (Insecure)"
)
post_nested <- group_by(post_year, y) %>% nest()
p_attach_sensall <- map2(post_nested$y, post_nested$data, function(ystring, data) {
  p <- ggplot(filter(summarised, y == ystring), aes(sensitivity_all, m, group = year, color = as.factor(year))) +
    geom_point() +
    geom_line() +
    xlab("") + 
    facet_wrap(~attach_dummy, strip.position = "bottom", labeller = as_labeller(attach_dummy_names)) +
    ylab(filter(symp_names, symptom == ystring) %>% .$name) +
    labs(color = "Year") +
    theme_bw(base_size = 20) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.x = element_text(size = 14)
    ) 
  p
})


# save(p_attach_sensall, file = here::here("data/rdata/attach_sens_all_newplots.Rds"))
# p_attach_sensall[[3]]




# after all, Carolina wants to split it up into first year and rest:
p_attach_sens_all <- map2(post_nested$y, post_nested$data, function(ystring, data) {
  map(firstyear, function(fy) {
    map(unique(summarised$attach_dummy), function(attachment_type) {
          if (fy) {
      temp <- filter(summarised, y == ystring, year == 1, attach_dummy == attachment_type)
      p <- mutate(temp, year = factor(year, levels = c("1", "2.5", "4", "5", "6", "7", "8", "10", "11", "12.5", "14"))) %>%
        ggplot(aes(sensitivity_all, m)) +
        geom_point(color = "black") +
        geom_line(color = "black") +
        scale_y_continuous(breaks = seq(min(temp$m), max(temp$m), 1)) +
        xlab("Sensitivity Throughout Childhood") + 
        ylab(filter(symp_names, symptom == ystring) %>% .$name) +
        theme_bw(base_size = 20) +
          ggtitle(glue("{attachment_type}_{fy}"))
    } else {
        temp <- filter(summarised, y == ystring, year != 1, attach_dummy == attachment_type)
        p <- mutate(temp, year = factor(year, levels = c("1", "2.5", "4", "5", "6", "7", "8", "10", "11", "12.5", "14"))) %>%
          ggplot(aes(sensitivity_all, m, group = year, color = as.factor(year))) +
          geom_point() +
          geom_line() +
          scale_y_continuous(breaks = seq(min(temp$m), max(temp$m), 1)) +
          scale_color_manual(values = c(
            "1" = "black", "2.5" = "#a6cee3", "4" = "#b2df8a", "5" = "#33a02c", "6" = "#fb9a99", "7" = "#e31a1c", "8" = "#fdbf6f", "10" = "#ff7f00", "11" = "#cab2d6", "12.5" = "#6a3d9a", "14" = "#1f78b4"), drop = FALSE) +
          xlab("Sensitivity Throughout Childhood") + 
          ylab(filter(symp_names, symptom == ystring) %>% .$name) +
          labs(color = "Year") +
          theme_bw(base_size = 20) +
          ggtitle(glue("{attachment_type}_{fy}"))
    }
    })
  })
})




p_attach_sens_all[[3]][[2]][[1]] + p_attach_sens_all[[3]][[2]][[2]] +   ylab("")
p_attach_sens_all[[3]][[1]][[1]] + theme(legend.position = "none") + p_attach_sens_all[[3]][[1]][[2]] +   ylab("")
save(p_attach_sens_all, file = here::here("data/rdata/attach_sens_all_newplots.Rds"))








firstyear <- c(0, 1)
ps <- map(firstyear, function(fy) {
  # generate outcomes for effect size prediction
  nd <- expand_grid(
    cc = c(0, 1), 
    bf = median(dlongimp[[1]]$bf),
    year = if (fy == 0) c(1, 2.5, 4, 5, 6, 7, 8, 10, 11, 12.5, 14) else 1,
    edu = 1,
    sib = 0,
    sensitivity_all = seq(-2, 2, 1),
    #sensitivity_firstyear = seq(-2, 2, 1),
    attach_dummy = c(0, 1),
    #attach_cat = 1:4,
    firstyear = fy
  )
  
  
  post_pred <- map_dfr(1:5, function(no) {
    if (no == 1) {
      models <- "m10"
    } else {
      models <-"m8"
    }
    post_pred <- map_dfr(sens_all_attach[[no]][[models]], function(model) {
      posterior_predict(model[[1]], newdata = nd, re.form = NA) %>%
        as.data.frame() %>%
        pivot_longer(
          everything(),
          names_to = "var",
          values_to = "values"
        ) %>%
        mutate(
          #values = values - median(values),
          y = integers[no]
        )
    })
  })
  
  match_names <- tibble(
    var = unique(post_pred$var), 
    #sensitivity_firstyear = nd$sensitivity_firstyear,
    firstyear = nd$firstyear,
    cc = nd$cc, 
    year = nd$year,
    #sens_factor = factor(nd$sensitivity_firstyear)
    sensitivity_all = nd$sensitivity_all,
    attach_dummy = nd$attach_dummy
    #attach_cat = nd$attach_cat
  )
  
  summarised <- post_pred %>%
    left_join(match_names, by = "var") %>%
    group_by(sensitivity_all, attach_dummy, y) %>%
    summarise(
      m = median(values),
      lower = hdi(values)[1],
      upper = hdi(values)[2],
      lower_50 = hdi(values, credMass = 0.68)[1],
      upper_50 = hdi(values, credMass = 0.68)[2]
    )  
  
  
  post_pred_centered <- post_pred %>%
    left_join(match_names, by = "var") %>%
    select(-var) 
  
  attach_dummy_names <- c(
    `0` = "Maternal Sensitivity (Secure)",
    `1` = "Maternal Sensitivity (Insecure)"
  )
  
  # alternative with patchwork
  post_nested <- group_by(post_pred_centered, y) %>% nest()
  p_all_sens_attach <- map2(post_nested$y, post_nested$data, function(ystring, data) {
    p <- ggplot(data, aes(sensitivity_all, values)) +
      # geom_hline(aes(yintercept = 0), linetype = "dashed") +
      #geom_jitter(alpha = 0.005) +
      geom_linerange(data = filter(summarised, y == ystring), aes(y = m, ymin = lower, ymax = upper), linewidth = 1.25, color = "#bdbdbd") +
      geom_pointrange(data = filter(summarised, y == ystring), aes(x = sensitivity_all, y = m, ymin = lower_50, ymax = upper_50), size = 1.5, linewidth = 1.8, color = "#636363") +
      #scale_y_continuous(breaks = c(seq(min(data$values), -3, 3), 0, seq(3, max(data$values), 3))) +
      facet_wrap(~attach_dummy, strip.position = "bottom", labeller = as_labeller(attach_dummy_names)) +
      scale_y_continuous(breaks = seq(-100, 100, 5)) +
      #scale_x_discrete(labels = c("Avoidant", "Secure", "Resistant", "Disorganized")) +
      xlab("") +
      ylab(filter(symp_names, symptom == ystring) %>% .$name) +
      theme_bw(base_size = 20) +
      theme(
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20)
      ) 
    p
  })
  
  p_all_sens_attach
})

ps[[2]][[5]]

save(ps, file = here::here("data/rdata/sens_all_attach_plots.Rds"))





# attachment continuous
attachment_sec_all_matsens <- map(integers, function(symptom) {
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + attachment_security * sensitivity_all +  (1 | id)"))
  m8_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.)))
    m8 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m8_imp_{m}_{symptom}_attachment_sec_all_matsenall.Rds"))
    )
    loo_m8 <- add_criterion(
      m8,
      "loo",
      # #file = here::here(glue("data/rdata/loo_m8_imp_{m}_{symptom}_attachment_sec__all_matsenall")),
      moment_match = FALSE
    )
    list(
      m8,
      loo_m8
    )
  })
  
  m8_post <- m8_post <- map_dfr(m8_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m8_h <- map(m8_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "attachment_security:sensitivity_all > 0")
  })
  
  
  form <- bf(glue("{symptom} ~ year + I(year^2) + firstyear * cc + bf * firstyear + edu + sib + attachment_security*sensitivity_all + (1 | id)"))
  
  m10_models <- map(1:5, function(m) {
    dtemp <- filter(dlongimp[[m]], if_any(symptom, ~ !is.na(.)))
    m10 <- brm(
      family = poisson(),
      formula = form,
      data = dtemp,
      prior = c(
        prior(normal(-0.15, 0.15), coef = "year"),
        prior(normal(0, 0.5), class = "b"),
        prior(normal(2, 0.75), class = "Intercept"),
        prior(exponential(1), class = "sd")
      ),
      file = here::here(glue("data/rdata/m10_imp_{m}_{symptom}_attachment_sec_all_matsenall.Rds"))
    )
    loo_m10 <- add_criterion(
      m10,
      "loo",
      # #file = here::here(glue("data/rdata/loo_m10_imp_{m}_{symptom}_attachment_sec_all_matsenall")),
      moment_match = FALSE
    )
    list(
      m10,
      loo_m10
    )
    
  })
  
  m10_post <- m10_post <- map_dfr(m10_models, function(mlist) {
    as_draws_df(mlist[[1]])
  })
  
  m10_h <- map(m10_models, function(mlist) {
    m <- mlist[[1]]
    hypothesis(m, "attachment_security:sensitivity_all > 0")
  })
  
  loo <- loo_compare(m8_models[[1]][[2]], m10_models[[1]][[2]])
  winner <- rownames(loo)[1]
  list(
    m8 = m8_models,
    m10 = m10_models,
    loo = loo,
    winner = winner,
    m8_h = m8_h,
    m10_h = m10_h,
    m8_post = m8_post,
    m10_post = m10_post
  )
})

# lets first check which model we should interpret
map(attachment_sec_all_matsens, ~.x$winner)
map2(1:5, attachment_sec_all_matsens, function(no, mlist) {
  
  mlist$m10_h
  
})


# to report effect size i use combined posteriors. For directional hypothesis
b <- map2_dfr(attachment_sec_all_matsens, 1:5, function(mlist, no) {
  
  posttemp <- mlist$m10_post
  
  select(posttemp, `b_attachment_security:sensitivity_all`) %>%
    summarise(
      m = median(`b_attachment_security:sensitivity_all`),
      lower = hdi(`b_attachment_security:sensitivity_all`)[1],
      upper = hdi(`b_attachment_security:sensitivity_all`)[2],
      lower_d = quantile(`b_attachment_security:sensitivity_all`, 0.025),
      upper_d = quantile(`b_attachment_security:sensitivity_all`, 0.975)
    ) %>%
    mutate(y = integers[no])
})
p <- map2_dfr(1:5, attachment_sec_all_matsens, function(no, mlist) {
  
  htemp <- mlist$m10_h
  
  p <- map_dbl(htemp, function(h) {
    h$hypothesis$Post.Prob
  })
  tibble(y = integers[no], p = median(p))
})
b <- full_join(b, p, by = "y") %>%
  mutate(across(where(is.numeric), function(x) format(round(x, 3), digits = 3)))
b
attachment_sec_all_matsens[[5]]$m10[[1]]


nd <- expand_grid(
  firstyear = 1,
  cc = 1,
  bf = median(dlongimp[[1]]$bf),
  year = 1,
  edu = 1,
  sib = 1,
  sex = 1,
  sensitivity_all = seq(-2, 2, 1),
  #sensitivity_firstyear = seq(-2, 2, 1),
  attachment_security = c(-2, 2),
  #attach_dummy = c(0, 1),
  #attach_cat = 1:4
)


post_pred <- map_dfr(1:5, function(no) {
  
  models <- "m10"
  
  if (b[b$y == integers[no], "p"] >= 0.95) {
    post_pred <- map_dfr(attachment_sec_all_matsens[[no]][[models]], function(model) {
      posterior_predict(model[[1]], newdata = nd, re.form = NA) %>%
        as.data.frame() %>%
        pivot_longer(
          everything(),
          names_to = "var",
          values_to = "values"
        ) %>%
        mutate(
          values = values - median(values),
          y = integers[no]
        )
    })
  }
})

match_names <- tibble(
  var = unique(post_pred$var),
  attachment_security = nd$attachment_security,
  sensitivity_all = nd$sensitivity_all
)
head(post_pred)
summarised <- post_pred %>%
  group_by(var, y) %>%
  summarise(
    m = median(values),
    lower = hdi(values)[1],
    upper = hdi(values)[2],
    lower_50 = hdi(values, credMass = 0.68)[1],
    upper_50 = hdi(values, credMass = 0.68)[2]
  )  %>% left_join(match_names, by = "var") %>%
  select(-var)
summarised


post_pred_centered <- post_pred %>%
  left_join(match_names, by = "var") %>%
  select(-var)



# alternative with patchwork
post_nested <- group_by(post_pred_centered, y) %>% nest()

attachment_security_names <- c(
  `-2` = "-2SD Attchment Security",
  `2` = "+2SD Attchment Security"
)

p_sens_attach <- map2(post_nested$y, post_nested$data, function(ystring, data) {
  p <- ggplot(data, aes(sensitivity_all, values)) +
    #geom_jitter(alpha = 0.005) +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    geom_linerange(data = filter(summarised, y == ystring), aes(y = m, ymin = lower, ymax = upper), linewidth = 1.25, color = "#bdbdbd") +
    geom_pointrange(data = filter(summarised, y == ystring), aes(x = sensitivity_all, y = m, ymin = lower_50, ymax = upper_50), size = 1.5, linewidth = 1.8, color = "#636363") +
    #scale_y_continuous(breaks = c(seq(min(data$values), -3, 3), 0, seq(3, max(data$values), 3))) +
    scale_y_continuous(breaks = seq(-30, 30, 2)) +
    #scale_x_discrete(labels = c("secure", "insecure")) +
    facet_wrap(~as.factor(attachment_security), strip.position = "bottom", labeller = as_labeller(attachment_security_names)) +
    xlab("Average Sensitivity") +
    ylab(filter(symp_names, symptom == ystring) %>% .$name) +
    theme_bw(base_size = 20) +
    theme(
      legend.position = "none",
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.x = element_text(size = 10)
    )
  p
})

p_sens_attach[[2]]





