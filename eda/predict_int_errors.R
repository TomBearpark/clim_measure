pacman::p_load(tidyverse, fixest, sf)

if(Sys.info()['user'] == 'tombearpark'){
  dir <- '/Users/tombearpark/Library/CloudStorage/Dropbox/clim_measure/data/'
} else {
  dir <- "EDIT THIS PATH TO YOUR DATA FOLDER"
}

# load files --------------------------------------------------------------
# Check file names
list.files(dir) %>% str_subset("daymean")

dist     <- read_csv(paste0(dir, "GHCNh_pairwisedistances_daymean_2011-2020_180plus.csv"))
stations <- read_csv(paste0(dir, "GHCNh_stationlist_daymean_2011-2020_180plus.csv"))

# Some functions of distance that seem intuitive: 
  # Number of stations within 500, 
  # Minimum distance to another station, 
  # Average distance to other stations that are within 500km
dist.vals <- tibble(
  dist[,1], 
  dist_to_closest_stn = apply(dist[, -1], 1, function(x) min(x[x>0])), 
  num_within_500 = apply(dist[, -1], 1, function(x) sum(x < 500)), 
  dist_within_500 = apply(dist[, -1], 1, function(x) sum(x*(x < 500) / num_within_500))
)

stations <- stations %>% left_join(dist.vals)

# load data ---------------------------------------------------------------
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

collapse <- function(df, group_vars, stations){
  df %>% 
    filter(!is.na(err)) %>% 
    group_by(across(all_of(group_vars))) %>% 
    summarize(mabs_err = mean(abs(err)),
              mean_sqrd_err = mean(err^2),
              clim = mean(t2m),
              n = n()) %>% 
    left_join(stations)  
}

# Prediction variables
covs <- c('Elevation', 'Latitude', 'Longitude', 'n', 
          'dist_to_closest_stn', 'num_within_500', 'dist_within_500')
covs.trans <- c("log(Elevation)", "log(abs(Latitude))", "log(abs(Longitude))", "log(n)", 
                "log(dist_to_closest_stn)", "log(num_within_500)", "log(dist_within_500)")

gen.ff <- function(covs){
  paste0("log(mabs_err) ~ ", paste0(covs, collapse = " + ")) %>% 
    as.formula()
}
# describe interpolation error by cross sectional covariates -------------------

df.mod <-collapse(df, "Station_ID", stations)

plt(df.mod, "Elevation")
plt(df.mod, "Latitude")
plt(df.mod, "Longitude")
plt(df.mod, "n")
plt(df.mod, "dist_to_closest_stn")
plt(df.mod, "num_within_500")
plt(df.mod, "dist_within_500")

feols(gen.ff(covs.trans), df.mod, vcov = 'hetero') %>% 
  etable(fitstat = 'r2')

# include year ------------------------------------------------------------

df.mod.y <- collapse(df, c("Station_ID", "year"), stations)

feols(gen.ff(covs.trans), df.mod.y, vcov = 'hetero', split= ~year) %>% 
  coefplot()

feols(gen.ff(c(covs.trans, "i(year)")), df.mod.y, vcov = 'hetero') %>% 
  coefplot()

# by month of year --------------------------------------------------------

df.mod.m <- collapse(df %>% mutate(month = month(date)), 
                     c("Station_ID", "month"), stations)

feols(gen.ff(covs.trans), df.mod.m, vcov = 'hetero', split= ~month) %>% 
  coefplot()

feols(gen.ff(c(covs.trans, "i(month)")), df.mod.m, vcov = 'hetero') %>% 
  etable()

# panel stuff -------------------------------------------------------------
df <- df %>% filter(!is.na(err))
gc()
m.err <- feols(err ~ 1 | date + Station_ID, df, vcov="hetero")
m.t   <- feols(t2m ~ 1 | date + Station_ID, df, vcov="hetero")

df$res.err <- resid(m.err)  
df$res.t   <- resid(m.t)

feols(res.err ~ res.t, df)
