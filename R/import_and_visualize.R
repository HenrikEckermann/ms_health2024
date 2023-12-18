set.seed(1)
library(tidyverse)
library(see)

f <- here::here("data/all BIBO health_varofinterest.csv")
d <- read_csv2(
  f,
  na = c("99999", "")
)



d <- d %>%
  select(
    -contains("_months"), 
    -Sensitivity_5w,
    -ZSensitivity_5w,
    -Sensitivity12m_mean,
    -ZSensitivity12m_mean,
    -Sensitivity30m_mean,
    -ZSensitivity30m_mean,
    -Sensitivity10y_mean,
    -ZSensitivity10y_mean,
    -Sensitivity14y_mean,
    -ZSensitivity14y_mean
    ) %>%
  rename(
    bf = Borstvoeding12m_duration,
    cc = Childcare,
    sex = childsex,
    attach_dummy = Dummysecureversusinsecure,
    sib = Siblingsnumber,
    edu = Med_cat,
    attach_cat = Attach_allcateg,
    sensitivity_all = Sensitivity_history_mean5points
  ) 
colnames(d) <- str_to_lower(colnames(d))


# bring data to long format
colnames(d) <- str_replace(colnames(d), "12m", "1y")
colnames(d) <- str_replace(colnames(d), "(\\d+)y", "\\1")
colnames(d) <- str_replace(colnames(d), "_(\\d?\\d?\\.?\\d+)_years", "\\1")
colnames(d) <- str_replace(colnames(d), "(totalhealth)_(.*)", "\\1\\2")
colnames(d)

mlr::summarizeColumns(d)

factors <- c("sex", "edu", "attach_cat", "attach_dummy", "sib")
integers <- c("resp", "skin", "dig", "gen", "totalhealth")
continuous <- c("sensitivity_firstyear", "sensitivity_all")
dlong <- d %>%
  pivot_longer(
    matches("\\d+$"),
    names_to = c("variable", "year"),
    names_pattern = "([a-z_]+)(\\d+\\.?\\d?)$") %>%
  pivot_wider(
    names_from = variable,
    values_from = value) %>%
  mutate(
    year = as.numeric(year),
    id = as.integer(id),
    across(all_of(factors), function(x) as.factor(x)),
    across(all_of(integers), function(x) as.integer(x)),
    across(all_of(continuous), function(x) scale(x)[, 1]),
    age_diff = age - year
  ) 

# the data is highly skewed and not suited for a gaussian model. The residuals 
# of the model matter in the end but also from the data generation process, this
# is closer to count data and a poisson model seems more appropriate
vplots <- map(integers, function(var) {
  dlong %>%
    mutate(year = as.factor(year)) %>%
    ggplot(aes_string("year", var)) +
    geom_violinhalf() + 
    #scale_y_log10() +
    theme_bw()
})
vplots[[1]]

symp_names <- tibble(
  symptom = c(
    "resp", "dig", "skin", "gen", "totalhealth"
  ),
  name = c(
    "Respiratory Symptoms", "Digestive Symptoms", "Skin Symptoms", "General Symptoms", "Total Symptoms"
  )
)

# figure 1
bplots <- map(integers, function(var) {
  dlong %>%
    mutate(year = as.factor(year)) %>%
    ggplot(aes_string("year", var),) +
    geom_boxplot(outlier.alpha = 0, color = "#C00000", fill = "#C00000", alpha = 0.5) +
    geom_jitter(width = 0.2, alpha = 0.2,  color = "#C00000", fill = "#C00000") +
    #scale_y_log10() +
    geom_vline(xintercept = 1.5, linetype = "dashed", size = 0.8) +
    xlab("Age (years)") +
    ylab(symp_names[symp_names$symptom == var, "name"]) +
    theme_bw(base_size = 15)
})
bplots[[1]]
save(bplots, file = here::here("data/rdata/bplots.Rds"))

#ggsave(bplots[[1]], device = "png", file = "fig/dist2.png")

# imputation
dwide <- d %>%
    mutate(
      id = as.integer(id),
      across(all_of(factors), function(x) as.factor(x)),
      across(all_of(contains(integers)), function(x) as.integer(x)),
      across(all_of(continuous), function(x) scale(x)[, 1])
    )
colnames(dwide)
dimp <- mice::mice(dwide)
vars_to_impute <- c(
  "age1", "age2.5", "age4", "age5", "age6", "age7", "age8", "age10", "age11", "age12.5", 
  "age14", "attach_cat", "attach_dummy", "attachment_security", "cc", "sex", 
  "edu", "sib", "bf", "sensitivity_firstyear", "sensitivity_all"
)

dlongimp <- map(1:5, function(m) {
  complete(dimp, m) %>%
    select(id, all_of(vars_to_impute)) %>%
    right_join(select(dwide, -all_of(vars_to_impute)), by = "id") %>%
    pivot_longer(
      matches("\\d+$"),
      names_to = c("variable", "year"),
      names_pattern = "([a-z_]+)(\\d+\\.?\\d?)$") %>%
    pivot_wider(
      names_from = variable,
      values_from = value) %>%
    mutate(
      year = as.numeric(year),
      id = as.integer(id),
      across(all_of(factors), function(x) as.factor(x)),
      across(all_of(integers), function(x) as.integer(x)),
      across(all_of(continuous), function(x) scale(x)[, 1]),
      age_diff = age - year,
      firstyear = factor(ifelse(year == 1, 1, 0))
    ) 
})
# how much does age vary around year?
ggplot(dlongimp[[1]], aes(as.factor(year), age_diff)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.1, alpha = 0.3)
mlr::summarizeColumns(dlongimp[[1]])
save(dlong, d, dimp, dwide, dlongimp, file = "data/rdata/d.Rds")

missing_df <- mlr::summarizeColumns(dlong)
nobs <- dim(dlong)[1]
mutate(missing_df, na_perc = na/nobs * 100)
# double check
missing_df2 <- mlr::summarizeColumns(dwide)
nobs2 <- dim(dwide)[1]
mutate(missing_df2, na_perc = na/nobs2 * 100)
