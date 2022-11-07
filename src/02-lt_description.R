# General life-table analysis

# Init ------------------------------------------------------------

set.seed(1987)

library(here); library(glue)
library(tidyverse); library(yaml)
library(patchwork)
library(gt)

# Constants -------------------------------------------------------

wd <- here()
cnst <- list()
cnst <- within(cnst, {
  groups_for_analysis = c('nhw','nhb','latinx')
  path_out = glue('{wd}/out')
  path_tmp = glue('{wd}/tmp')
  # number of Poisson life-table replicates
  n_sim = 500
  years_ind = 2010:2019
})

dat <- list()
fig <- list()
tab <- list()

# Function --------------------------------------------------------

# simple piecewise-exponential life-table
CalculateLifeTable <-
  function (df, x, nx, Dx, Ex) {

    require(dplyr)

    df %>%
      transmute(
        x = {{x}},
        nx = {{nx}},
        mx = {{Dx}}/{{Ex}},
        px = exp(-mx*{{nx}}),
        qx = 1-px,
        lx = head(cumprod(c(1, px)), -1),
        dx = c(-diff(lx), tail(lx, 1)),
        Lx = ifelse(mx==0, lx*nx, dx/mx),
        Tx = rev(cumsum(rev(Lx))),
        ex = Tx/lx,
        edx = rev(cumsum(rev(c(c(.5*c(ex[-1L], 0) + (1-.5)*ex)[-length(x)],ex[length(x)])*dx)))/lx
      )

  }

# Data ------------------------------------------------------------

# input data for life-table calculation
dat$lt_input_100 <- readRDS(glue('{cnst$path_out}/mx_input_100.rds'))

# Subset to groups of interest ---------------------------------
dat$lt_input_100 <- dat$lt_input_100[reth %in% cnst$groups_for_analysis]

# Calculate annual life tables ------------------------------------

# life-tables by group, sex, and year
# open age_group 100+
dat$lt_100 <-
  dat$lt_input_100 %>%
  arrange(reth, sex, year, age) %>%
  group_by(reth, sex, year) %>%
  group_modify(~{
    CalculateLifeTable(.x, age, nx, deaths, pop)
  }) %>%
  ungroup()


# Create Poisson life-table replicates ----------------------------

# create life table replicates by group, sex, and year
# based on repeatedly sampling death counts from a Poisson
dat$lt_100_sim <-
  dat$lt_input_100%>%
  expand_grid(id_sim = 1:cnst$n_sim) %>%
  group_by(reth, sex, year, age) %>%
  mutate(death_total_sim = rpois(cnst$n_sim, deaths)) %>%
  arrange(age,reth, sex, year) %>%
  group_by(id_sim, reth, sex, year) %>%
  group_modify(~{
    CalculateLifeTable(.x, age, nx, death_total_sim, pop)
  }) %>%
  ungroup()

# Assemble table with ex statistics -------------------------------

# central estimates of life-expectancy, annual life-expectancy difference,
# and average annual life-expectancy difference 2015 to 2020

# 95% uncertainty intervals around the central estimates
dat$lt_ex_ci <-
  dat$lt_100_sim %>%
  filter(year %in% cnst$years_ind) %>%
  select(id_sim, reth, sex, year, x, mx, ex, edx) %>%
  arrange(id_sim, reth, sex, x, year) %>%
  group_by(reth, sex, x, year) %>%
  summarise(
    ex_q025 = quantile(ex, 0.025, na.rm = TRUE),
    ex_q975 = quantile(ex, 0.975, na.rm = TRUE),
    edx_q025 = quantile(edx, 0.025, na.rm = TRUE),
    edx_q975 = quantile(edx, 0.975, na.rm = TRUE)
  )


# assemble all the ex statistics in a single table
# for further computation
dat$lt_ex_long <-
  left_join(
    dat$lt_100,
    dat$lt_ex_ci
  )



# Export ----------------------------------------------------------

# save the regrouped life table input data
saveRDS(dat$lt_ex_long, file = glue('{wd}/out/lt_CI.rds'))
write.csv(dat$lt_ex_long, file = glue('{wd}/out/lt_CI.csv'))

