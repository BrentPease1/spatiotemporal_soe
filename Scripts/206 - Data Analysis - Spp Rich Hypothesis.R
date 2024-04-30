library(here)
library(sf)
library(data.table)
library(stringr)
library(exactextractr)
library(unmarked)
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# Bring in BBS observations
source(here('Scripts/01 - Data Prep - Get Species Observations.R'))

# Bring in MODIS Data
source(here('Scripts/02 - Data Prep - Get MODIS Data.R'))


# Bring in IUCN range maps
iucn <- st_read(here('Data/IUCN Range Shapefiles/avian_species/data_0.shp'))
# cerulean warbler not in master set so adding in separately 
cewa <- st_read(here('Data/IUCN Range Shapefiles/Steophaga_cerulea'))
# combine
iucn <- rbind(iucn, cewa)
rm(cewa)
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# fix names
names(landcover) <- paste0('y', 2001:2018)
# set up annual loop here
max_lc_year <- names(landcover) %>% 
  str_extract("[0-9]{4}") %>% 
  as.integer() %>% 
  max()

min_lc_year <- names(landcover) %>% 
  str_extract("[0-9]{4}") %>% 
  as.integer() %>% 
  min()


# create output holder
holder <- data.frame(AOU = 0, year = 0, species_observed = 0, total_indv_counted = 0)
counter = 0

for(fs in 1:nrow(focal_species)){
  

  
  # get another unique ID
  bbs_observations[, state_route_grp := .GRP, by = .(StateNum, Route)]
  
  # get locations where focal species detected
  confirmed_locations <- bbs_observations[AOU %in% focal_species[fs,AOU],unique(state_route_grp)]

  # subset all observations to just confirmed_locations
  bbs_species <- bbs_observations[state_route_grp %in% confirmed_locations,]
  
  # get annual data
  for(yy in min_lc_year:max_lc_year){
    cat('species', fs, 'of', nrow(focal_species), 'year', yy, 'starting.\n\n')
    
    # subset down to single year
    bbs_annual <- bbs_species[Year == yy,]
  
    # get the unique number of observed species across all sites where focal_species was detected
    observed_richness <- bbs_annual[,length(unique(AOU))]
    total_counted <- bbs_annual[, sum(route_counts)]
    
    # store results
    counter = counter + 1
    holder[counter,'AOU'] <- focal_species[fs,AOU]
    holder[counter,'year'] <- yy
    holder[counter,'species_observed'] <- observed_richness
    holder[counter,'total_indv_counted'] <- total_counted
    
    
    cat('species', fs, 'of', nrow(focal_species), 'year', yy, 'completed.\n\n')
    
  }
}
saveRDS(object = holder, file = here('Results/SoE_BBS/Species_Richness/observed_spp_richness.RDS'))

