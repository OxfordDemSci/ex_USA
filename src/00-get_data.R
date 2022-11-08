# Get data in functional form

# Init ------------------------------------------------------------
library(here); library(glue)
library(tidyverse); library(yaml)
library(patchwork); library(data.table)
library(gt)
library(haven)

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
  min_age = 15
  var_names = c('age','year','sex','reth','deaths','pop','Obesity',
                'Suicide','Substance','Cancers',
                'Respiratory','Infectious','Accidents','Covid','Other')
})

dat <- list()
fig <- list()
tab <- list()

# Read data --------------------------------------------------------

##Pop	population count in race/ethnic, sex, year group	
#Mort	total number of deaths  in race/ethnic, sex, year group	
#Reth	race/ethnic group variable, should have labels (1 "NHW" 2 "NHB" 3 "LatinX" 4 "AI/AN" 5 "Other")	
#Male	1=male, 0=female	


#some checks prior

dat$analytic <- data.table(read_dta(glue('{wd}/dat/analytic_obesity.dta')))

#view(dat$analytic)

dat$masters <- dat$analytic[,c('age','year','male','reth','mort','pop','m1',
                               'm2','m3','m4','m5','m6','m7','m8','m9')]

names(dat$masters) <- cnst$var_names

dat$masters$scheme <- 'Masters'

dat$acosta <- dat$analytic[,c('age','year','male','reth','mort','pop','ac1',
                               'ac2','ac3','ac4','ac5','ac6','ac7','ac8','ac9')]

names(dat$acosta) <- cnst$var_names

dat$acosta$scheme <- 'Acosta'

dat$adair <- dat$analytic[,c('age','year','male','reth','mort','pop','ad1',
                               'ad2','ad3','ad4','ad5','ad6','ad7','ad8','ad9')]

names(dat$adair) <- cnst$var_names

dat$adair$scheme <- 'Adair'

dat$gbd <- dat$analytic[,c('age','year','male','reth','mort','pop','cc1',
                           'cc2','cc3','cc4','cc5','cc6','cc7','cc8','cc9')]

names(dat$gbd) <- cnst$var_names

dat$gbd$scheme <- 'Contributing'

dat$ucod <- dat$analytic[,c('age','year','male','reth','mort','pop','o1',
                           'o2','o3','o4','o5','o6','o7','o8','o9')]

names(dat$ucod) <- cnst$var_names

dat$ucod$scheme <- 'Underlying'

## bind them to have only one dataset

dat$analytic_bind <- rbind(dat$masters,dat$acosta,dat$adair,dat$gbd,dat$ucod)
#view(dat$analytic_bind)

#check death rates, mort >  pop, consistency with totals, etc.

dat$analytic_bind[deaths > pop]


#move to long format, easier to calculate proportions

dat$deaths_input_100 <- melt.data.table(data = dat$analytic_bind,
                                          id.vars = c('age','year','sex','reth','scheme','deaths','pop'),
                                          variable.name = 'cause',value.name = 'deaths_cause')

dat$deaths_input_100$scheme <- factor(dat$deaths_input_100$scheme, levels =  c("Acosta", "Adair", 
                                                                               "Masters", "Contributing", "Underlying"))

#check for NA's

#view(check_NA)
dat$deaths_input_100[is.na(dat$deaths_input_100$deaths)]$deaths_cause <- 0
dat$deaths_input_100[is.na(dat$deaths_input_100$deaths)]$deaths <- 0

check_NA <- dat$deaths_input_100[is.na(dat$deaths_input_100$deaths)]
check_NA

#check (cause-specific) = total
sum(dat$deaths_input_100[,sum(deaths_cause) - deaths[1], by = .(age,year,sex,reth,scheme)]$V1)


#save dataset
saveRDS(dat$deaths_input_100,glue('{cnst$path_out}/deaths_input_100.rds'))


