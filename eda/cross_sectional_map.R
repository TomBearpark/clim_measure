pacman::p_load(tidyverse, fixest, sf, rnaturalearth, terra, elevatr)

if(Sys.info()['user'] == 'tombearpark'){
  dir <- '/Users/tombearpark/Library/CloudStorage/Dropbox/clim_measure/data/'
} else {
  dir <- "EDIT THIS PATH TO YOUR DATA FOLDER"
}
theme_set(theme_classic())

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

# Convert stations into a sf object for plotting

world <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")

stations.sf <- stations %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) 
stations.sf %>% 
  ggplot() + 
  geom_sf(data =world, fill=NA)+
  geom_sf(color = 'red')

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

df %>% arrange(-err)
df %>% filter(t2m==-999)
df %>% filter(t2m > -999, !is.na(t2m_interpolated)) %>% arrange(t2m)

ggplot(df) + geom_histogram(aes(x=err))
ggplot(df %>% filter(abs(err)<20)) + geom_histogram(aes(x=err))

df <- df %>% 
  mutate(err = if_else(
    t2m_interpolated < -50 | t2m == -999 | Elevation < 0, NA, err
  ))

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
covs <- c('Elevation', 'Latitude', 'Longitude', 
          #'n', 
          'dist_to_closest_stn', 'num_within_500', 'dist_within_500')
covs.trans <- c("log(Elevation)", "log(abs(Latitude))", "log(abs(Longitude))", 
                # "log(n)", 
                "log(dist_to_closest_stn)", "log(num_within_500)", "log(dist_within_500)")

gen.ff <- function(covs, yvar="log_mabs_err"){
  paste0(yvar, " ~ ", paste0(covs, collapse = " + ")) %>% 
    as.formula()
}

# describe interpolation error by cross sectional covariates -------------------

df.mod <-collapse(df, "Station_ID", stations) %>% 
  filter(Elevation >0) 

ggplot(df.mod) + geom_histogram(aes(x = mabs_err))

plt(df.mod, "Elevation")
plt(df.mod, "Latitude")
plt(df.mod, "Longitude")
# plt(df.mod, "n")
plt(df.mod, "dist_to_closest_stn")
plt(df.mod, "num_within_500")
plt(df.mod, "dist_within_500")

m.reg <- feols(gen.ff(covs.trans, "mabs_err"), df.mod, vcov = 'hetero')

m.reg %>% etable(fitstat = 'r2')

left_join(stations.sf, df.mod) %>% 
  ggplot() + 
  geom_sf(data =world, fill=NA, alpha=.5)+
  geom_sf(aes(color = mabs_err)) +
  scale_color_viridis_c()

# try ML model ------------------------------------------------------------

# ML model
pacman::p_load(tidymodels, vip)

pred.data <- df.mod %>% select(-Station_ID, -Station_name, -mean_sqrd_err, 
                               -n, -clim, -any_of(c("lm", "ml")), 
                               -n_days) %>% 
  mutate(abs_long = abs(Longitude), abs_lat = abs(Latitude))

data_split <- initial_split(pred.data, prop = 0.8)
train_data <- training(data_split)
test_data  <- testing(data_split)

rec <- recipe(mabs_err ~ ., data = train_data) %>%
  step_normalize(all_predictors()) %>%
  step_impute_median(all_predictors())

model_spec <- rand_forest() %>%
  set_engine("ranger", importance="impurity") %>%
  set_mode("regression")

workflow <- workflow() %>%
  add_model(model_spec) %>%
  add_recipe(rec)

model_fit <-  fit(workflow, data = train_data)
final_fit <-  fit(model_fit, data = pred.data)

df.mod$ml <- predict(final_fit, new_data = pred.data)$.pred
df.mod$lm <- predict(m.reg, newdata = pred.data)

final_fit %>% extract_fit_parsnip() %>% vip(num_features = 20)

df.mod %>% 
  summarize(MLrmse = sqrt(mean((mabs_err - ml)^2)), 
            LMrmse = sqrt(mean((mabs_err - lm)^2)))

# Looks like we can already get a 90% R2 with this simple(ish) model
etable(feols(mabs_err ~ lm, df.mod, vcov = 'hetero'), 
       feols(mabs_err ~ ml, df.mod, vcov = 'hetero')
       )

df.mod %>% 
  pivot_longer(cols = c("lm", "ml")) %>%
  ggplot() + geom_point(aes(x = mabs_err, y = value)) + geom_abline() + 
  facet_wrap(~name) + 
  xlab("Error") + ylab("Predicted error")


# create prediction grid --------------------------------------------------

# plot(world['sovereignt'])
# unique(world$sovereignt) %>%sort
uk <- world %>% filter(sovereignt ==  "United Kingdom" , 
                        type == "Country")  %>% 
  select(sovereignt)
ggplot(uk)+geom_sf()

uk_stations <- stations.sf %>% 
  st_intersection(uk)

# Elevation resolution here is 3, which is 13km. This is clearly very insufficient
elev.in <- get_elev_raster(uk, z=3, prj = st_crs(uk), clip='locations')

elev <- as.data.frame(elev.in, xy=T) %>% 
  as_tibble() %>% 
  rename(elev = 3) %>% 
  na.omit()  %>%
  mutate(id = 1:nrow(.))

elev %>% 
  ggplot() + 
  geom_tile(aes(x = x, y = y, fill = elev)) +
  geom_sf(data = uk_stations, color = 'red') +
  scale_fill_viridis_c()

# For each point in the raster, find distance to closest station
elev.sf <- st_as_sf(elev, coords = c("x", "y"), crs = st_crs(uk))

# Compare elevation in the two datasets
st_join(uk_stations, elev.sf, join = st_nearest_feature) %>% 
  ggplot() + geom_point(aes(x = Elevation, y = elev))+
  geom_abline()


dist <- st_distance(elev.sf, uk_stations)
# Convert to miles
dist <- units::set_units(dist, "km")
elev$dist_to_closest_stn <- apply(dist, 1, min)

elev %>% 
  ggplot() + 
  geom_tile(aes(x = x, y = y, fill = log(dist_to_closest_stn)))+
  geom_sf(data = uk_stations, color = 'red') +
  scale_fill_viridis_c()

# For each point in the raster, find number of stations within 500km
elev$num_within_500 <- apply(dist, 1, function(x) sum(x < 500))

# For each point in the raster, find average distance to other stations within 500km
elev$dist_within_500 <- apply(dist, 1, function(x) sum(x*(x < 500))) / elev$num_within_500

elev$Latitude <- st_coordinates(elev.sf)[,2]
elev$Longitude <- st_coordinates(elev.sf)[,1]

elev <- elev %>%  mutate(abs_long = abs(Longitude), abs_lat = abs(Latitude)) %>% 
  rename(Elevation = elev)


elev %>% 
  pivot_longer(cols = c(Elevation, Latitude, Longitude, 
                        dist_to_closest_stn, num_within_500, dist_within_500)) %>% 
  group_by(name) %>% 
  mutate(value = rescale(value))%>% 
  ggplot() + 
  geom_tile(aes(x = x, y = y, fill = value)) + 
  geom_sf(data = uk_stations, color = 'red', size = .1)+
  facet_wrap(~name)+
  scale_fill_viridis_c()

# Predict error for each point in the raster
elev$ml <- predict(final_fit, new_data = elev)$.pred
elev$lm <- predict(m.reg, newdata = elev)

elev %>% 
  ggplot() + 
  geom_abline(color = 'red')+
  geom_point(aes(x = ml, y = lm))

elev %>% 
  pivot_longer(cols = c(ml, lm)) %>% 
  ggplot() + 
  geom_sf(data = uk %>% select(geometry))+
  geom_tile(aes(x = Longitude, y = Latitude, fill = value)) +
  scale_fill_viridis_c() +
  facet_wrap(~name)
