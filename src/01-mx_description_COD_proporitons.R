# General life-table analysis

# Init ------------------------------------------------------------
set.seed(1987)

library(data.table)
library(here); library(glue)
library(tidyverse); library(yaml)
library(patchwork)
library(gt); library(openxlsx)

# Constants -------------------------------------------------------

wd <- here()
cnst <- list()
config <- read_yaml(glue('{wd}/cfg/config.yaml'))
cnst <- within(cnst, {
  regions_for_analysis = config$regions_for_all_cause_analysis
  path_out = glue('{wd}/out')
  path_tmp = glue('{wd}/tmp')
  # number of Poisson life-table replicates
  n_sim = 500
  #years to show in graps
  years_fig = c(2010,2015,2019)
  cols_fig = c('#80146E',
                '#7A55AB',
                '#478EC1',
                '#34B8C0',
                '#7ED5B8',
                '#C6E8BC',
                '#F5F2D8',
               'lightgrey')
})

dat <- list()
fig <- list()
tab <- list()

#80146E
#7A55AB
#478EC1
#34B8C0
#7ED5B8
#C6E8BC
#F5F2D8

# '#FF9AA1',
# '#EC6068',
# '#AA4146',
# '#250F10',
# '#007575',
# '#00A6A4',
# '#64CBCF',

# functions ----------------------------------------------------------------

CI_s_mx <- function(dt = .SD, n_sim = cnst$n_sim){
  deaths_sim <- rpois(lambda = dt$deaths,n = n_sim)
  deaths_low <- quantile(deaths_sim, 0.025, na.rm = TRUE)
  deaths_high <- quantile(deaths_sim, 0.975, na.rm = TRUE)
  
  mx_low <- deaths_low/dt$pop
  mx_high <- deaths_high/dt$pop
  mx <-  dt$deaths/dt$pop
  
  return(data.table(mx,mx_low,mx_high))
}



# Data ------------------------------------------------------------

# input data for life-table calculation
# harmonized death counts and population exposures with open age group 100+
dat$deaths_input_100 <- readRDS(glue('{cnst$path_out}/deaths_input_100.rds'))

#create a dataset with total deaths only 
dat$mx_input_100 <- dat$deaths_input_100[scheme %in% 'Masters' & cause %in% 'Obesity',
                                         .(deaths = deaths, pop = pop), 
                                         by = .(reth,year,sex,age)]

# Description of mx -------------------------------------------------------

#calculate mx
dat$mx_input_100[,nx := 1]
dat$mx_input_100[age == 100,nx := Inf]

#recode some variables
#Reth	race/ethnic group variable, should have labels (1 "NHW" 2 "NHB" 3 "LatinX" 4 "AI/AN" 5 "Other")	
#Male	1=male, 0=female	
dat$mx_input_100[, reth:= factor(x = reth,labels = c('nhw','nhb','latinx','aian','other'))]
dat$mx_input_100[, sex:= factor(x = sex,labels = c('female','male'))]

dat$mx_input_100_to_save <- dat$mx_input_100

dat$mx_input_100 <- dat$mx_input_100[, CI_s_mx(dt = .SD,n_sim = cnst$n_sim),by = .(reth,year,sex,age)]

# Diagnostic plots of mx --------------------------------------------------

fig$fig_mx <- ggplot(data = dat$mx_input_100[year %in% cnst$years_fig]) +
  ggtitle('mx by year, sex and subgroup')+
  geom_ribbon(
    aes(
      x = age,
      ymin = log(mx_low),
      ymax = log(mx_high),
      fill = sex,
      group = sex
    ),
    alpha = 0.3
  ) +
  geom_line(
    aes(x = age, y = log(mx), color = sex, group = sex)
  ) +
  geom_point(
    aes(x = age, y = log(mx), color = sex,group = sex),
    size = 1
  ) +
  labs(
    x = 'Age',
    y = 'mx'
  ) +
  facet_grid(reth~year)+
  theme_bw()

fig$fig_mx
  
ggsave(glue('{cnst$path_out}/fig_mx.png'),plot = fig$fig_mx,width = 10,height = 10)

#NTS: need to potentially take care of aian and other using some modelling strategy


# Diagnostic plots of COD --------------------------------------------------

dat$deaths_input_100[, reth:= factor(x = reth,labels = c('nhw','nhb','latinx','aian','other'))]
dat$deaths_input_100[, sex:= factor(x = sex,labels = c('female','male'))]
dat$deaths_input_100[, cause:= as.factor(cause)]
dat$deaths_input_100[, scheme:= as.factor(scheme)]
dat$deaths_input_100[, deaths_cause_prop:= deaths_cause/deaths]


dat$fig_cod_female <- dat$deaths_input_100[sex == 'female' 
                                           & year == 2019 
                                           & reth %in% c('nhw','nhb','latinx')]

# dat$fig_cod_female <- dat$fig_cod_female[,cause_short := cause]
# dat$fig_cod_female[cause %in% c('drug','alcohol','suicide'),]$cause_short <- 'despair'
# dat$fig_cod_female[cause %in% c('female_cancer','male_cancer','rest_cancers'),]$cause_short <- 'cancers'
# dat$fig_cod_female <- dat$fig_cod_female[, .(deaths_cause = sum(deaths_cause)),
#                                          by = .(age,year,sex,reth,scheme,deaths,pop,cause_short)]
# dat$fig_cod_female[, deaths_cause_prop:= deaths_cause/deaths]
# dat$fig_cod_female[, cause_short:= as.character(cause_short)]
# 
# dat$fig_cod_female[, cause_short:=  factor(cause_short, c('obesity',
#                                                           'despair',
#                                                           'accidents',
#                                                           'cancers',
#                                                           'infectious',
#                                                           'respiratory',
#                                                           'rest'))]



#rename variables for plot
# levels(dat$fig_cod_female$cause_short) <- c('Obesity','Despair','Accidents','Cancers',
#                                             'Infectious','Respiratory','Rest')
levels(dat$fig_cod_female$reth) <- c('White', 'Black', 'Latino', 'Asian', 'Other')



fig$fig_COD_prop_females <- 
ggplot(data = dat$fig_cod_female[cause!='Covid'],
        aes(
          x = age,
          y = deaths_cause_prop,
          fill = cause)) +
  #ggtitle('COD structure over time, scheme and subgroup, females 2020')+
  geom_area(stat = 'identity', position = 'fill')+
  labs(
    x = 'Age',
    y = 'Proportion'
  ) +
  facet_grid(reth~scheme)+
   scale_fill_manual('Cause of Death',values = cnst$cols_fig)+
  theme_bw()


fig$fig_COD_prop_females

################################################################################

dat$fig_cod_male <- dat$deaths_input_100[sex == 'male' 
                                           & year == 2019 
                                           & reth %in% c('nhw','nhb','latinx')]

# dat$fig_cod_male <- dat$fig_cod_male[,cause_short := cause]
# dat$fig_cod_male[cause %in% c('drug','alcohol','suicide'),]$cause_short <- 'despair'
# dat$fig_cod_male[cause %in% c('female_cancer','male_cancer','rest_cancers'),]$cause_short <- 'cancers'
# dat$fig_cod_male <- dat$fig_cod_male[, .(deaths_cause = sum(deaths_cause)),
#                                          by = .(age,year,sex,reth,scheme,deaths,pop,cause_short)]
# dat$fig_cod_male[, deaths_cause_prop:= deaths_cause/deaths]
# dat$fig_cod_male[, cause_short:= as.character(cause_short)]
# dat$fig_cod_male[, cause_short:=  factor(cause_short, c('obesity',
#                                                           'despair',
#                                                           'accidents',
#                                                           'cancers',
#                                                           'infectious',
#                                                           'respiratory',
#                                                           'rest'))]

# levels(dat$fig_cod_male$cause_short) <- c('Obesity','Despair','Accidents','Cancers',
#                                             'Infectious','Respiratory','Rest')
levels(dat$fig_cod_male$reth) <- c('White', 'Black', 'Latino', 'Asian', 'Other')




fig$fig_COD_prop_males <- 
  ggplot(data = dat$fig_cod_male[cause!='Covid'],
         aes(
           x = age,
           y = deaths_cause_prop,
           fill = cause)) +
  #ggtitle('COD structure over time, scheme and subgroup, males 2020')+
  geom_area(stat = 'identity', position = 'fill')+
  labs(
    x = 'Age',
    y = 'Proportion'
  ) +
  facet_grid(reth~scheme)+
  scale_fill_manual('Cause of Death',values = cnst$cols_fig)+
  theme_bw()


fig$fig_COD_prop_males

ggsave(glue('{cnst$path_out}/fig_cod_prop_females.pdf'),plot = fig$fig_COD_prop_females,width = 10,height = 10, device = cairo_pdf)

ggsave(glue('{cnst$path_out}/fig_cod_prop_females.png'),plot = fig$fig_COD_prop_females,width = 10,height = 10)

ggsave(glue('{cnst$path_out}/fig_cod_prop_males.pdf'),plot = fig$fig_COD_prop_males,width = 10,height = 10, device = cairo_pdf)

ggsave(glue('{cnst$path_out}/fig_cod_prop_males.png'),plot = fig$fig_COD_prop_males,width = 10,height = 10)




#### save 

saveRDS(dat$mx_input_100_to_save,glue('{cnst$path_out}/mx_input_100.rds'))

saveRDS(dat$deaths_input_100,glue('{cnst$path_out}/cod_input_100.rds'))


#### Some results in the paper

# Under the UCOD scheme there are very few obesity-related deaths (less than XYZ %), 
dat$paper_results_1 <- dat$deaths_input_100[ year == 2019 
                     & reth %in% c('nhw','nhb','latinx'),]

dat$paper_results_1 <- dat$paper_results_1[, .(total_deaths = sum(deaths), deaths_cause = sum(deaths_cause)), 
                    by = .(year,sex,reth,scheme,cause)]

dat$paper_results_1[, prop_cause := deaths_cause/total_deaths*100]

#check
dat$paper_results_1[, sum(prop_cause), by =.(year,sex,reth,scheme)]

# results 
dat$paper_results_1[scheme %in% 'Underlying' & cause %in% 'Obesity', ]

dat$paper_results_1[scheme %in% 'Contributing' & cause %in% 'Obesity', ]




