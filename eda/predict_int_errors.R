pacman::p_load(tidyverse, fixest, sf)

if(Sys.info()['user'] == 'tombearpark'){
  dir <- '/Users/tombearpark/Library/CloudStorage/Dropbox/clim_measure/data/'
} else {
  dir <- "X"
}


# load files --------------------------------------------------------------
list.files(dir) %>% str_subset("daymean")

dist     <- read_csv(paste0(dir, "GHCNh_pairwisedistances_daymean_2011-2020_180plus.csv"))
stations <- read_csv(paste0(dir, "GHCNh_stationlist_daymean_2011-2020_180plus.csv"))

# Number of stations within 500, and minimum distance to another station
dist.vals <- tibble(
  dist[,1], 
  dist_to_closest_stn = apply(dist[, -1], 1, function(x) min(x[x>0])), 
  num_within_500 = apply(dist[, -1], 1, function(x) sum(x < 500)), 
  dist_within_500 = apply(dist[, -1], 1, function(x) sum(x*(x < 500) / num_within_500))
)

# start with a single year ------------------------------------------------

years <- 2011:2020

df.vals <- map_dfr(years, function(year)
  read_csv(file.path(dir, paste0("GHCNh_t2m_daymean_", year, ".csv"))) %>% 
    mutate(year = !!year))

df.int  <- map_dfr(years, function(year)
  read_csv(file.path(dir, paste0("GHCNh_t2m_daymean_", year, "_interpolated.csv"))) %>% 
    mutate(year = !!year))

df <- left_join(df.vals, df.int, by = join_by(date, Station_ID, year)) %>% 
  left_join(stations, by = join_by(Station_ID)) %>% 
  mutate(err = t2m - t2m_interpolated)
gc()

ggplot(df) + geom_histogram(aes(x=err))

# funcs -------------------------------------------------------------------

plt <- function(df.mod, var, scales = 'free'){ 
  df.mod %>% 
    pivot_longer(cols = contains('err')) %>% 
    ggplot(aes(x = !!sym(var), y = value)) + geom_point() + geom_smooth() + 
    facet_wrap(~name, scales=scales)
}

# describe interpolation error by cross sectional covariates -------------------

df.mod <- df %>% 
  filter(!is.na(err)) %>% 
  group_by(Station_ID) %>% 
  summarize(mabs_err = mean(abs(err)),
            mean_sqrd_err = mean(err^2),
            clim = mean(t2m),
            n = n()) %>% 
  left_join(stations)  %>% 
  left_join(dist.vals)

plt(df.mod, "Elevation")
plt(df.mod, "Latitude")
plt(df.mod, "Longitude")
plt(df.mod, "n")
plt(df.mod, "dist_to_closest_stn")
plt(df.mod, "num_within_500")
plt(df.mod, "dist_within_500")

feols(
  log(mabs_err) ~ 
    log(Elevation) + 
    log(abs(Latitude)) + 
    log(abs(Longitude)) +
    log(n) + 
    log(dist_to_closest_stn) + 
    log(num_within_500) + 
    log(dist_within_500)+
    log(clim),
  df.mod, 
  vcov = 'hetero'
) %>% 
  etable(fitstat = 'r2')

# include year ------------------------------------------------------------

df.mod <- df %>% 
  filter(!is.na(err)) %>% 
  group_by(Station_ID, year) %>% 
  summarize(mabs_err = mean(abs(err)),
            mean_sqrd_err = mean(err^2),
            clim = mean(t2m),
            n = n()) %>% 
  left_join(stations)  %>% 
  left_join(dist.vals)

feols(
  log(mabs_err) ~ 
    log(Elevation) + 
    log(abs(Latitude)) + 
    log(abs(Longitude)) +
    log(n) + 
    log(dist_to_closest_stn) + 
    log(num_within_500) + 
    log(clim), 
  df.mod, 
  vcov = 'hetero', 
  split = ~year
) %>% 
  coefplot()

feols(
  log(mabs_err) ~ 
    log(Elevation) + 
    log(abs(Latitude)) + 
    log(abs(Longitude)) +
    log(n) + 
    log(dist_to_closest_stn) + 
    log(num_within_500) + 
    log(clim) + 
    i(year), 
  df.mod, 
  vcov = 'hetero', 
) %>% 
  coefplot(drop = 'Const')

# by month of year --------------------------------------------------------

df.mod <- df %>% 
  filter(!is.na(err)) %>% 
  mutate(month = lubridate::month(date)) %>%
  group_by(Station_ID, month) %>% 
  summarize(mabs_err = mean(abs(err)),
            mean_sqrd_err = mean(err^2),
            clim = mean(t2m),
            n = n()) %>% 
  left_join(stations)  %>% 
  left_join(dist.vals)

feols(
  log(mabs_err) ~ 
    log(Elevation) + 
    log(abs(Latitude)) + 
    log(abs(Longitude)) +
    log(n) + 
    log(dist_to_closest_stn) + 
    log(num_within_500) + 
    log(clim), 
  df.mod, 
  vcov = 'hetero', 
  split = ~month
) %>% 
  coefplot()

feols(
  log(mabs_err) ~ 
    log(Elevation) + 
    log(abs(Latitude)) + 
    log(abs(Longitude)) +
    log(n) + 
    log(dist_to_closest_stn) + 
    log(num_within_500) + 
    log(clim) + 
    i(month), 
  df.mod, 
  vcov = 'hetero', 
) %>% 
  coefplot(drop = 'Const')
