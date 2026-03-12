
#' RR function to calculate the relative risk
#'
#' @param date # date (days in Date format)
#' @param temp # mean annual temperature
#' @param deaths # daily deaths
#' @param year_eval # years to be evaluated
#' @param temp_dec # decimal digits of the temperature data
#' @param ci_level # confidence level of the prediction
#'
#' @returns # MMT (minimum mortality temperature), RR (realtive risk)
#' @export
#'
#' @examples
RR_fun <- function(date, temp, deaths, year_eval, temp_dec = 0.1, ci_level = 0.95){

  # data preparation (ten years, time trend, and seasonality) ---------------

  dat <- tibble(date = date,
                temp = temp,
                deaths = deaths) |>
    mutate(year = year(date)) |>
    filter((year >= year_eval - 9) & (year <= year_eval)) |>
    arrange(date) |>
    mutate(
      time = as.numeric(date - min(date)) + 1,
      dow  = factor(weekdays(date),
                    levels = c("Monday","Tuesday","Wednesday","Thursday",
                               "Friday","Saturday","Sunday"))
    )


  # crossbasis setup --------------------------------------------------------


  # temperature: quadratic B-spline with nodes at the 75th percentile
  k_temp <- quantile(dat$temp, 0.75, na.rm = TRUE)

  # lags: 7 days at 2 knots on a log scale
  lag_knots <- logknots(7, nk = 2)

  # crossbasis
  cb_temp <- crossbasis(
    dat$temp,
    lag = 7,
    argvar = list(fun = "bs", degree = 2, knots = k_temp),
    arglag = list(fun = "ns", knots = lag_knots)
  )

  # dim(cb_temp)



  # Quasi-Poisson DLNM fitten -----------------------------------------------


  # including long-term trend and day of the week
  # standard in heat-mortality analyses
  # time trend (approximately 8 df/year, adjustable)

  model <- glm(
    deaths ~ cb_temp +
      ns(time, df = 8) +
      dow,
    family = quasipoisson(),
    data = dat
  )

  # summary(model)


  # minimum mortality temperature (MMT) -------------------------------------

  # prediction (with median temperature as centering value)

  pred0 <- crosspred(cb_temp,
                     model,
                     cen = median(dat$temp),
                     by = temp_dec,
                     ci.level = ci_level)


  # reference temperature (temperature with minimal mortality risk)
  MMT <- pred0$predvar[which.min(pred0$allRRfit)[1]]

  # MMT in output data
  dat$MMT <- MMT


  # prediction with MMT as centering value ----------------------------------

  pred <- crosspred(
    cb_temp,
    model,
    cen = MMT,
    by = temp_dec,
    ci.level = ci_level
  )


  # relative risk (RR) per day ----------------------------------------------

  # approximate the nearest forecast value for each day
  dat$RR <- approx(pred$predvar, pred$allRRfit, xout = dat$temp)$y
  dat$RR_low  <- approx(pred$predvar, pred$allRRlow,  xout = dat$temp)$y
  dat$RR_high <- approx(pred$predvar, pred$allRRhigh, xout = dat$temp)$y
  dat$RR_se <- approx(pred$predvar, pred$allse, xout = dat$temp)$y


  # output ------------------------------------------------------------------

  out <- dat |>
    mutate(period = paste0(min(dat$year), "-", max(dat$year)),
           year_pred = if_else(year == year_eval, 1, 0)) |>
    select(year, date, temp, deaths, MMT,
           RR, RR_low, RR_high, RR_se, period, year_pred)

  return(out)

}

