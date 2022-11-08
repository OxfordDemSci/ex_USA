
# Plot ex changes -------------------------------------------------

# General life-table analysis

# Init ------------------------------------------------------------
library(data.table)
library(here); library(glue)
library(tidyverse);library(ggplot2)
library(patchwork)
library(gt)

# Constants -------------------------------------------------------

wd <- here()
cnst <- list()
cnst <- within(cnst, {
  path_out = glue('{wd}/out')
  path_tmp = glue('{wd}/tmp')
  #years to show in graps
  years_fig = c(2010,2015,2019)
  cols_fig = c('#38907B','#E1C7CF','#4A232E')
})

dat <- list()
fig <- list()
tab <- list()



# Data ------------------------------------------------------------

# input data for life-table calculation
dat$lt_100 <- data.table(readRDS(glue('{cnst$path_out}/lt_CI.rds')))

# Plots -------------------------------------------------------------------


# plot  e0 and e60
# from 2019 to 2020 by sex and region and compare with average
# annual change over 2015 to 2019 period

fig$fig_e0 <- ggplot(data = dat$lt_100[x == 0 & year <= max(cnst$years_fig)]) +
  ggtitle('e0 by year, sex and subgroup')+
  geom_ribbon(
    aes(
      x = year,
      ymin = ex_q025,
      ymax = ex_q975,
      fill = reth,
      group = reth
    ),
    alpha = 0.3
  ) +
  geom_line(
    aes(x = year, y = ex, color = reth, group = reth)
  ) +
  geom_point(
    aes(x = year, y = ex, color = reth,group = reth),
    size = 1
  ) +
  labs(
    x = 'Year',
    y = 'Life expectancy at birth'
  ) +
  facet_grid(.~sex)+
  scale_colour_manual(values= cnst$cols_fig)+
  scale_fill_manual(values= cnst$cols_fig)+
  scale_x_continuous(breaks = c(2010,2012,2014,2016,2018,2019)) +
  theme_bw()

fig$fig_e0



fig$fig_e60 <- ggplot(data = dat$lt_100[x == 60 & year <= max(cnst$years_fig)]) +
  ggtitle('e60 by year, sex and subgroup')+
  geom_ribbon(
    aes(
      x = year,
      ymin = ex_q025,
      ymax = ex_q975,
      fill = reth,
      group = reth
    ),
    alpha = 0.3
  ) +
  geom_line(
    aes(x = year, y = ex, color = reth, group = reth)
  ) +
  geom_point(
    aes(x = year, y = ex, color = reth,group = reth),
    size = 1
  ) +
  labs(
    x = 'Year',
    y = 'Life expectancy at age 60'
  ) +
  facet_grid(.~sex)+
  scale_colour_manual(values= cnst$cols_fig)+
  scale_fill_manual(values= cnst$cols_fig)+
  scale_x_continuous(breaks = c(2010,2012,2014,2016,2018,2019)) +
  theme_bw()

fig$fig_e60


fig$fig_ex <- fig$fig_e0/fig$fig_e60


ggsave(glue('{cnst$path_out}/fig_ex.png'),plot = fig$fig_ex,width = 10,height = 10)


### life expectancy v lifespan inequality
fig$fig_ex_ed <- ggplot(data = dat$lt_100[x == 0]) +
  ggtitle('life expectancy vs lifespan inequality')+
  geom_point(
    aes(x = ex, y = edx, color = reth,shape = sex,group = reth),
    size = 1
  ) +
  labs(
    x = 'Life expectanacy',
    y = 'Lifespan inequality'
  )+
  theme_bw()
  
fig$fig_ex_ed

ggsave(glue('{cnst$path_out}/fig_ex_dx.png'),plot = fig$fig_ex_ed,width = 5,height = 5)
