library(here)
library(sf)
library(data.table)
library(stringr)

# bring in complete species list and focal species list (clark rushing 2019-2020 papers)
species_list <- fread(here('Data/BBS/SpeciesList.csv'))
focal_species <- fread(here('Data/BBS/focal_species.csv'))

# read in route data
routes <- fread(here('Data/BBS/routes.csv'))
# read in weather data
weather <- fread(here("Data/BBS/weather.csv"))


# get AOU codes from species list to filter observations
focal_species <- merge(focal_species, species_list[,c('genus_species', "AOU")], 
                       by.x = 'Latin name', by.y = 'genus_species')


# read in all observation data
files <- list.files(path = here('Data/BBS'), pattern = 'fifty', full.names = T)
bbs_observations <- do.call(rbind, lapply(files, fread))
rm(files)



# Keep routes that meet official BBS criteria (listed in RunType.PDF file within BBS folder)
# RouteTypeDetailID (Routes.txt) == 1
# RPID (RunProtocolID) (observations, weather.txt) == 101
# QualityCurrentID (weather.txt) == 1
# RouteTypeID == 1 (roadside == 1, water = 2, off-road = 3)
routes <- routes[RouteTypeDetailID == 1 & RouteTypeID == 1,]
bbs_observations <- bbs_observations[RPID == 101,]
weather <- weather[RPID == 101,]
weather <- weather[QualityCurrentID == 1,]

# join route data to observations to get lat/long, etc
# this will drop observations from routes that are not routetypedetailid == 1
bbs_observations <- merge(bbs_observations, routes, by = c('StateNum', 'Route', 'CountryNum'))

# join weather data to observations for detection information
bbs_observations <- merge(bbs_observations, 
                          weather[, c('RouteDataID', 'Month', 
                                      'Day', 'StartTemp', 'EndTemp',
                                      'StartWind', 'EndWind',
                                      'StartSky', 'EndSky',
                                      'StartTime', 'EndTime',
                                      'Assistant')],
                          by = 'RouteDataID')

# clean up
rm(routes, weather)


# summarize counts across stops
# first get stop columns
stops <- which(str_detect(names(bbs_observations), 'Stop'))

# get total counts along a route
bbs_observations[, route_counts := rowSums(.SD), .SDcols = stops]

# get presence/absence
bbs_observations[, (stops) := lapply(.SD, FUN = function(x) ifelse(x > 0, 1, 0)), .SDcols = stops]
bbs_observations[, route_binom := rowSums(.SD), .SDcols = stops]

# drop the stop columns
#bbs_observations <- bbs_observations[, .SD, .SDcols = !stops]
rm(stops)
