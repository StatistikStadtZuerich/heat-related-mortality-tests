# header ------------------------------------------------------------------

# heat-related mortality: first tests
# rok, March 2026



# preparation -------------------------------------------------------------

# packages
library(dlnm) # distributed lag non-linear models
library(lubridate) # date, time
library(mgcv) # flexible regressions
library(patchwork) # multiple graphics
library(renv) # reproducibility
library(splines) # non-linear models
library(styler) # code style
library(tidyverse) # tidyverse
library(truncnorm) # truncated normal distribution


# functions
source("R/RR_fun.R")


# colors
# https://github.com/StatistikStadtZuerich/zuericolors/blob/main/R/palettes.R

qual12 <- c("#3431DE", "#0A8DF6", "#23C3F1", "#7B4FB7", "#DB247D", "#FB737E",
            "#007C78", "#1F9E31", "#99C32E", "#9A5B01", "#FF720C", "#FBB900")




# paths -------------------------------------------------------------------

# deaths in Zurich
path_death <- "https://data.stadt-zuerich.ch/dataset/bev_tag_todesfaelle_quartier_geschl_ag_herkunft_od4211/download/BEV421OD4211.csv"

# daily temperature at SMA Zurich (harmonized time series)
path_temp <- "https://data.geo.admin.ch/ch.meteoschweiz.ogd-nbcn/sma/ogd-nbcn_sma_d_historical.csv"




# deaths ------------------------------------------------------------------

# import deaths
dea_imp <- read_csv(path_death) |>
  mutate(date = as.Date(GueltigAbDat),
         age = factor(if_else(AlterV20Kurz %in% c("80-99", "100 u. älter"),
                              "80+", "0-79"),
                      levels = c("0-79", "80+")),
         sex = factor(if_else(SexCd == 2, "female", "male"),
                      levels = c("female", "male"))) |>
  group_by(date, age, sex) |>
  summarise(deaths = sum(AnzSterWir, na.rm = TRUE), .groups = "drop")


# deaths for all days from May to September
dea <- expand_grid(date = seq(min(dea_imp$date), max(dea_imp$date), by = "day"),
                   age  = sort(unique(dea_imp$age)),
                   sex  = sort(unique(dea_imp$sex))) |>
  left_join(dea_imp, by = c("date", "age", "sex")) |>
  mutate(deaths = replace_na(deaths, 0),
         month = month(date)) |>
  filter(month >= 5, month <= 9) |>
  select(date, age, sex, deaths)




# temperature -------------------------------------------------------------

# import mean daily temperature
tem_imp <- read_delim(path_temp, delim = ";") |>
  mutate(date = as.Date(reference_timestamp, format = "%d.%m.%Y %H:%M")) |>
  rename(temp = ths200d0) |>
  select(date, temp)


# checks
tem_imp |>
  mutate(diff = c(1, diff(date))) |>
  summarize(min_date_diff = min(diff),
            max_date_diff = max(diff),
            temp_NA = sum(is.na(temp)))




# imported data: death and temperature ------------------------------------

# data is limited by death values available
range(dea_imp$date)
range(tem_imp$date)

# imported data
imp <- dea |>
  left_join(tem_imp, by = "date")

# check
nrow(dea)
nrow(imp)

# total (without age and sex)
total <- imp |>
  group_by(date) |>
  summarize(deaths = sum(deaths),
            temp = unique(temp),
            .groups = "drop")




# relative risk (RR) ------------------------------------------------------

# years (only if ten years with past data available)
years <- (year(min(total$date))+9) : year(max(total$date))

# data: relative risk
RR_dat <- map_dfr(years, function(y) {
  RR_fun(
    date = total$date,
    temp = total$temp,
    death = total$deaths,
    year_eval = y
  )
})



# exposure-response relationship ------------------------------------------

ggplot(RR_dat, aes(x = temp, y = RR)) +
  geom_ribbon(aes(ymin = RR_low, ymax = RR_high), alpha = 0.2) +
  geom_line(linewidth = 1) +
  facet_wrap(~ period) +
  geom_hline(yintercept = 1) +
  geom_vline(aes(xintercept = MMT), linetype = "dotted") +
  labs(x = "mean daily temperature (°C)",
       y = "relative risk") +
  theme_minimal()



# Monte Carlo (MC) simulation ----------------------------------------------

# RR stats per period and temperature
RR_stats <- RR_dat |>
  group_by(period, temp) |>
  summarize(RR_mean = mean(RR),
            RR_sd = mean(RR_se),
            RR_low = min(RR_low),
            RR_high = max(RR_high),
            .groups = "drop")

# MC (approx. 3 seconds)
t0 <- Sys.time()
RR_mc <- RR_stats |>
  mutate(sim = purrr::pmap(
    list(RR_mean, RR_sd, RR_low, RR_high),
    \(m, s, lo, hi) rtruncnorm(1000, lo, hi, m, s)
  )) |>
  unnest(sim)
Sys.time() - t0

# confidence interval
RR_ci_mc <- RR_mc |>
  group_by(period, temp) |>
  summarize(
    RR_low_mc  = quantile(sim, 0.025, na.rm = TRUE),
    RR_high_mc = quantile(sim, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# initial RR data with MC confidence interval
RR_dat_mc <- RR_dat |>
  left_join(RR_ci_mc, by = c("period", "temp"))

# check
nrow(RR_dat)
nrow(RR_dat_mc)


# exposure-response relationship with MC confidence interval --------------

ggplot(RR_dat_mc, aes(x = temp, y = RR)) +
  geom_ribbon(aes(ymin = RR_low_mc, ymax = RR_high_mc), alpha = 0.2) +
  geom_line(linewidth = 1) +
  facet_wrap(~ period) +
  geom_hline(yintercept = 1) +
  geom_vline(aes(xintercept = MMT), linetype = "dotted") +
  labs(x = "mean daily temperature (°C)",
       y = "relative risk") +
  theme_minimal()




# attributable fraction and attributable deaths ---------------------------

AF_AD <- RR_dat_mc |>
  filter(year_pred == 1) |>
  mutate(AF = pmax(0, if_else(temp > MMT, (RR - 1) / RR, 0)),
         AF_low = pmax(0, if_else(temp > MMT, (RR_low_mc - 1) / RR_low_mc, 0)),
         AF_high = pmax(0, if_else(temp > MMT, (RR_high_mc - 1) / RR_high_mc, 0)),
         AD = AF * deaths,
         AD_low = AF_low * deaths,
         AD_high = AF_high * deaths)

# plot(ecdf(AF_AD$AF))



# heat-related mortality by year ------------------------------------------

mort_y <- AF_AD |>
  group_by(year) |>
  summarise(across(c(AD, AD_low, AD_high),
                   ~ sum(.x, na.rm = TRUE)),
            .groups = "drop")

ggplot(mort_y) +
  geom_col(aes(x = factor(year), y = AD), fill = qual12[6]) +
  geom_errorbar(aes(x = factor(year), ymin = AD_low, ymax = AD_high),
                width = 0.2, col = "grey40") +
  labs(x = "year", y = "heat-related deaths per year") +
  theme_minimal()



# heat-related mortality by year and heat intensity -----------------------

# temperature classes
temp_classes <- c("moderate", "hot", "very hot")

with_temp_classes <- AF_AD |>
  mutate(temp_class = factor(case_when(
    temp < 25 ~ temp_classes[1],
    (temp >= 25) & (temp < 27) ~ temp_classes[2],
    temp >= 27 ~ temp_classes[3]),
    levels = temp_classes)
  )

# mortality by year and temperature class
mort_yt <- with_temp_classes |>
  group_by(year, temp_class) |>
  summarise(AD = sum(AD), .groups = "drop")



ggplot() +
  geom_col(data = mort_yt,
           aes(x = factor(year), y = AD, fill = fct_rev(temp_class))) +
  geom_errorbar(data = mort_y,
                aes(x = factor(year), ymin = AD_low, ymax = AD_high),
                width = 0.2, col = "grey40") +
  scale_fill_manual(values = qual12[c(5, 11, 12)]) +
  labs(x = "year", y = "heat-related deaths per year", fill = "") +
  theme_minimal()





# proportion of heat related deaths ---------------------------------------

prop_y <- dea |>
  mutate(year = year(date)) |>
  filter(year >= min(mort_y$year)) |>
  group_by(year) |>
  summarize(total = sum(deaths), .groups = "drop") |>
  left_join(mort_y, by = "year") |>
  mutate(percent = pmin(100, AD / total * 100),
         percent_low = pmin(100, AD_low  / total * 100),
         percent_high = pmin(100, AD_high  / total * 100))


ggplot(prop_y) +
  geom_col(aes(x = factor(year), y = percent), fill = qual12[6]) +
  geom_errorbar(aes(x = factor(year), ymin = percent_low, ymax = percent_high),
                width = 0.2, col = "grey40") +
  labs(x = "year", y = "proportion of heat-related deaths in %") +
  theme_minimal()



# selected years (daily values) -------------------------------------------

# Why these years? high values (2015, 2018), and last year (2025)
years_selected <- AF_AD |>
  filter(year %in% c(2015, 2018, 2025)) |>
  mutate(plot_date = as.Date(paste0("2000-", format(date, "%m-%d"))))


p1 <- ggplot(years_selected) +
  geom_line(aes(x = plot_date, y = temp), col = qual12[1]) +
  geom_line(aes(x = plot_date, y = MMT), col = qual12[1], linetype = "dashed") +
  facet_wrap(~ factor(year)) +
  labs(x = "", y = "mean daily temperature (°C)") +
  theme_minimal()

p2 <- ggplot(years_selected) +
  geom_line(aes(x = plot_date, y = AD), col = qual12[5]) +
  geom_ribbon(aes(x = plot_date, ymin = AD_low, ymax = AD_high),
              fill = qual12[5], alpha = 0.2) +
  facet_wrap(~ factor(year)) +
  labs(x = "", y = "daily heat-related deaths") +
  theme_minimal()

p1 / p2



# MMT by year -------------------------------------------------------------

MMT_y <- AF_AD |>
  group_by(year) |>
  summarise(MMT = unique(MMT), .groups = "drop")

ggplot(MMT_y) +
  geom_col(aes(x = factor(year), y = MMT), fill = qual12[1]) +
  labs(x = "year", y = "minimum mortality temperature (°C)") +
  theme_minimal()

