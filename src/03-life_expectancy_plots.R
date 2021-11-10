
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

fig$fig_e0 <- ggplot(data = dat$lt_100[x == 0]) +
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
  facet_grid(.~sex)

fig$fig_e0



fig$fig_e60 <- ggplot(data = dat$lt_100[x == 60]) +
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
  facet_grid(.~sex)

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
  )
  
fig$fig_ex_ed

ggsave(glue('{cnst$path_out}/fig_ex_dx.png'),plot = fig$fig_ex_ed,width = 5,height = 5)


# Plot mx changes -------------------------------------------------

dat$mx_change <-
  dat$lt_85_sim %>%
  mutate(
    age_group = cut(x, c(0, 40, 60, 70, Inf), right = FALSE)
  ) %>%
  group_by(id_sim, region_iso, sex, year, age_group) %>%
  summarise(
    mx = sum(dx)/sum(Lx)
  ) %>%
  pivot_wider(names_from = year, values_from = mx) %>%
  mutate(diff = `2020`/`2019`) %>%
  group_by(region_iso, sex, age_group) %>%
  summarise(
    mean_diff = mean(diff),
    q025_diff = quantile(diff, 0.025),
    q975_diff = quantile(diff, 0.975),
    sig_elevated = ifelse(q025_diff > 1, TRUE, FALSE)
  )

fig$mx_change <-
  dat$mx_change %>%
  ggplot(aes(x = age_group, color = sex, group = sex)) +
  geom_hline(yintercept = 1, color = 'grey') +
  geom_pointrange(aes(
    y = mean_diff,ymin = q025_diff, ymax = q975_diff,
    alpha = sig_elevated
  ), fatten = 1, position = position_dodge(width = 0.5)) +
  scale_y_log10(breaks = c(0.8, 1, 1.2),
                labels = c('.8', '1', '1.2')) +
  scale_color_manual(values = fig_spec$sex_colors) +
  facet_wrap(~region_iso) +
  fig_spec$MyGGplotTheme(scaler = 0.8, panel_border = TRUE) +
  coord_flip() +
  scale_alpha_manual(values = c(0.2, 1)) +
  labs(
    title = 'Ratio of life-table death rates 2020 to 2019',
    subtitle = '95% CIs via Poisson simulation of raw death counts',
    x = 'Age group', y = 'Ratio 2020 to 2019'
  )
fig$mx_change

# Compare our ex estimates with wpp estimates ---------------------

walk(c(0, 60), ~{
  fig[[glue('e{.x}_consistency_check')]] <<-
    bind_rows(
      `85+` = filter(dat$lt_85, x == .x),
      `100+` = filter(dat$lt_100, x == .x),
      .id = 'open_age_group'
    ) %>%
    ggplot(aes(x = year, y = ex, color = sex)) +
    geom_segment(
      aes(x = 2015, xend = 2019, y = ex_wpp_estimate,
          yend = ex_wpp_estimate),
      size = 1, alpha = 0.5,
      data =
        dat$lt_input_85_sub %>%
        filter(age_start == .x, year == 2018)
    ) +
    geom_point(aes(shape = open_age_group)) +
    geom_line(
      aes(y = ex_hmd_estimate),
      data =
        dat$lt_input_85_sub %>%
        filter(age_start == .x)
    ) +
    facet_wrap(~region_iso, scales = 'free_y') +
    scale_x_continuous(
      breaks = 2015:2020,
      labels = c('', '2016', '', '2018', '', '2020')
    ) +
    scale_color_manual(values = fig_spec$sex_colors) +
    scale_shape_manual(values = c(`85+` = 1, `100+` = 3)) +
    fig_spec$MyGGplotTheme(panel_border = TRUE, grid = 'xy', scaler = 0.8) +
    labs(
      title = glue('Estimated yearly life expectancy at age {.x} compared with HMD (thin line) and WPP (bold line) 5 year average estimates'),
      y = glue('e{.x}'),
      shape = 'Open age group',
      color = 'Sex'
    )
})
fig$e0_consistency_check
fig$e60_consistency_check