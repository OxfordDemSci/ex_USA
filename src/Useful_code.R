#Useful code
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
        ex = Tx/lx
      )
    
  }


life.expectancy.cod.fun <-
  function(mx.cod, x, nx,cond_age = 0){
    dim(mx.cod) <- c(length(x),length(mx.cod)/length(x))
    mx          <- rowSums(mx.cod)
    life.expectancy.fun(mx,x,nx,cond_age)
  }


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


Decomp <-function (func, rates1, rates2, N, ...) {
  y1 <- func(rates1, ...)
  y2 <- func(rates2, ...)
  d <- rates2 - rates1
  n <- length(rates1)
  delta <- d/N
  x <- rates1 + d * matrix(rep(0.5:(N - 0.5)/N, length(rates1)), 
                           byrow = TRUE, ncol = N)
  cc <- matrix(0, nrow = n, ncol = N)
  for (j in 1:N) {
    for (i in 1:n) {
      z <- rep(0, n)
      z[i] <- delta[i]/2
      cc[i, j] <- func((x[, j] + z), ...) - func((x[, j] - 
                                                    z), ...)
    }
  }
  return(rowSums(cc))
}


models$gam.nb <-
  gam(
    deaths ~
      # special days
      first_week*sex*age.n.fct +
      last_week*sex*age.n.fct +
      week21*sex*age.n.fct +
      # log-linear long term trend
      time*sex*age.n.fct +
      # penalized spline for age effect
      s(age.n, bs = 'ps', k = 6, by = sex) +
      # penalized cyclic spline for seasonality
      s(week_ify, bs = 'cp', k = 8, by = sex) +
      # smooth interaction between age and seasonality
      ti(
        week_ify, age.n,
        bs = c('cp', 'ps'),
        k = c(8, 6),
        by = sex
      ) +
      offset(log(exposures)),
    data = dat$training.data,
    family = nb(link = 'log'),
    method = 'REML'
  )

# GAM Poisson -----------------------------------------------------

# Same as before, but with Poisson distribution.

models$gam.poisson <-
  gam(
    deaths ~
      # special days
      first_week*sex*age.n.fct +
      last_week*sex*age.n.fct +
      week21*sex*age.n.fct +
      # log-linear long term trend
      time*sex*age.n.fct +
      # penalized spline for age effect
      s(age.n, bs = 'ps', k = 6, by = sex) +
      # penalized cyclic spline for seasonality
      s(week_ify, bs = 'cp', k = 8, by = sex) +
      # smooth interaction between age and seasonality
      ti(
        week_ify, age.n,
        bs = c('cp', 'ps'),
        k = c(8, 6),
        by = sex
      ) +
      offset(log(exposures)),
    data = dat$training.data,
    family = poisson(link = 'log'),
    method = 'REML'
  )



AverageMortalityModel <-
  function (df, year, week, deaths, exposures, training.years, ...) {
    require(dplyr)
    .year = enquo(year); .strata = enquos(...)
    .week = enquo(week); .deaths = enquo(deaths);
    .exposures = enquo(exposures)
    
    avg_mx <-
      df %>%
      filter(!!.year %in% training.years) %>%
      group_by(!!!.strata, !!.week) %>%
      summarise(
        avg_mortality = mean(!!.deaths/!!.exposures),
        .groups = 'drop'
      )
    
    structure(list(avg = avg_mx), class = 'avgmx')
    
  }

predict.avgmx <- function (object, newdata, exposures, ...) {
  require(dplyr)
  left_join(newdata, object$avg) %>%
    transmute(deaths = !!enquo(exposures) * avg_mortality) %>%
    pull(deaths)
}

models$avg.mortality <-
  AverageMortalityModel(
    dat$training.data,
    year, week,
    deaths, exposures,
    training.years = cnst$avg.mortality.period,
    sex, age.n
  )

