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
    
    -(new.ex - ex)[1]
    
  }


# calculate life expectancy and cause deleted life expectancy ------------------------------------------

dat$cod_input_100[cause == 'rest']
unique(dat$cod_input_100$cause)

dat$cause_deleted_results <- dat$cod_input_100[, .(e0 = life.expectancy.fun(mx = mx,x = age,nx = nx,cond_age = 0),
                                                          cause_contribution = cause.deleted.life.expectancy.fun(mx = mx,
                                                                                                  x = age_start,
                                                                                                  nx = nx,
                                                                                                  Ri= deaths_cause_prop ,
                                                                                                  cond_age = 0)), 
                                                      by = .(reth, scheme, sex, year, cause)]


# Make a plot -------------------------------------------------------------

fig$fig.cause_deleted_females <- ggplot(data = dat$cause_deleted_results[cause != 'rest' & cause != 'male_cancer' & sex =='female'], aes(year, cause_contribution, col = reth, group = reth))+
  geom_line()+
  ggtitle('Females, cause contrib to e0')+
  facet_grid(cause~scheme,scales="free_y")

fig$fig.cause_deleted_females



fig$fig.cause_deleted_males <- ggplot(data = dat$cause_deleted_results[cause != 'rest'& cause != 'female_cancer' & sex =='male'], aes(year, cause_contribution, col = reth, group = reth))+
  geom_line()+
  ggtitle('Males, cause contrib to e0')+
  facet_grid(cause~scheme,scales="free_y")

fig$fig.cause_deleted_males

ggsave(glue('{cnst$path_out}/fig_cause_deleted_females.png'),plot = fig$fig.cause_deleted_females,width = 10,height = 10)


ggsave(glue('{cnst$path_out}/fig_cause_deleted_males.png'),plot = fig$fig.cause_deleted_males,width = 10,height = 10)


#save results
saveRDS(dat$cause_deleted_results,glue('{cnst$path_out}/cause-deleted-results.rds'))

write.csv(dat$cause_deleted_results,glue('{cnst$path_out}/cause-deleted-results.csv'))
