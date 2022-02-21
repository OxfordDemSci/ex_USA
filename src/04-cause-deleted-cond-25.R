# Cause-deleted analysis

# Init ------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(patchwork)
library(here); library(glue)
library(data.table)
library(reshape2)


# Constants -------------------------------------------------------

wd <- here()
cnst <- list()
cnst <- within(cnst, {
  groups_for_analysis = c('nhw','nhb','latinx')
  path_out = glue('{wd}/out')
  path_tmp = glue('{wd}/tmp')
  # number of Poisson life-table replicates
  n_sim = 500
  years_fig = c(2010,2015,2020)
  condition_age = 25
})

dat <- list()
fig <- list()
tab <- list()


# Data ------------------------------------------------------------

# input data for life-table calculation
dat$cod_input_100 <- readRDS(glue('{cnst$path_out}/cod_input_100.rds'))

# Subset to groups of interest ---------------------------------
dat$cod_input_100 <- dat$cod_input_100[reth %in% cnst$groups_for_analysis]

dat$cod_input_100[, mx := deaths/pop]

dat$cod_input_100[,nx := 1]
dat$cod_input_100[age == 100,nx := Inf]

#check
#dat$cod_input_100[is.na(deaths_cause_prop)]

#view(dat$mx_input_decom)
# Functions -------------------------------------------------------

#General Life expectancy function
life.expectancy.fun  <-
  function (mx, x, nx, cond_age = 0) {

    px = exp(-mx*{{nx}})
    qx = 1-px
    lx = head(cumprod(c(1, px)), -1)
    dx = c(-diff(lx), tail(lx, 1))
    Lx = ifelse(mx==0, lx*nx, dx/mx)
    Tx = rev(cumsum(rev(Lx)))
    ex = Tx/lx

    ex[cond_age + 1]

  }

#source Preston et al 2002
cause.deleted.life.expectancy.fun  <-
  function (mx, x, nx, Ri, cond_age = 0) {
    
    #functions from master life table
    px = exp(-mx*{{nx}})
    lx = head(cumprod(c(1, px)), -1)
    dx = c(-diff(lx), tail(lx, 1))
    Lx = nx[-length(nx)] *lx[-1]+.5*dx[-length(nx)]
    Lx = c(Lx,(dx/mx)[length(nx)])
    #Lx = ifelse(mx==0, lx*nx, dx/mx)
    Tx = rev(cumsum(rev(Lx)))
    ex = Tx/lx
    
    #ASDLT 
    new.px = px^(1-Ri)
    new.lx = lx/lx
    new.lx = head(cumprod(c(1, new.px)), -1)
    new.dx = c(-diff(new.lx), tail(new.lx, 1))
    new.Lx = nx[-length(nx)] *new.lx[-1]+.5*new.dx[-length(nx)]
    new.Lx = c(new.Lx,new.lx[length(new.lx)]/(mx[length(mx)]*(1-Ri[length(Ri)])))
    new.Tx = rev(cumsum(rev(new.Lx)))
    new.ex = new.Tx/new.lx
    
    -(new.ex - ex)[cond_age + 1]
    
  }


# calculate life expectancy and cause deleted life expectancy ------------------------------------------

dat$cod_input_100[cause == 'rest']
unique(dat$cod_input_100$cause)

dat$cause_deleted_results <- dat$cod_input_100[, .(ex = life.expectancy.fun(mx = mx,x = age,nx = nx,cond_age = cnst$condition_age),
                                                          cause_contribution = cause.deleted.life.expectancy.fun(mx = mx,
                                                                                                  x = age_start,
                                                                                                  nx = nx,
                                                                                                  Ri= deaths_cause_prop ,
                                                                                                  cond_age = cnst$condition_age)), 
                                                      by = .(reth, scheme, sex, year, cause)]

dat$cause_deleted_results[, cause_contribution_months := cause_contribution*12]



# Obesity plot -------------------------------------------------------------

dat$fig_obesity_contributions <- dat$cause_deleted_results[cause %in% 'obesity' &
                                                             year %in% cnst$years_fig]

levels(dat$fig_obesity_contributions$reth) <- c('White', 'Black', 'Latino', 'Asian', 'Other')
levels(dat$fig_obesity_contributions$scheme) <- c('Acosta','Adair','GBD','Masters','Ucod')
levels(dat$fig_obesity_contributions$sex) <- c('Female','Male')
dat$fig_obesity_contributions$Scheme <- dat$fig_obesity_contributions$scheme
  

fig$fig_obesity_contributions <- ggplot(dat$fig_obesity_contributions) +
  #ggtitle('Contribution of obesity to life expectancy according to different schemes')+
  geom_point(
    aes(x = cause_contribution_months, y = reth, color = Scheme,group = Scheme,shape = Scheme),
    size = 2,
  ) +
  labs(
    x = 'Contribution (months)',
    y = ''
  ) +
  facet_grid(year~sex)+
  scale_colour_manual( values = c('#DAEDC2','#7ED5B8','#2BABC2','#6C6AB5','#80146E'))

fig$fig_obesity_contributions


ggsave(glue('{cnst$path_out}/obesity_contribution_25.png'),plot = fig$fig_obesity_contributions,width = 8,height = 5)






