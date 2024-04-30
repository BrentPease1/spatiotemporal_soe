library(here)
library(sf)
library(data.table)
library(stringr)
library(exactextractr)
library(unmarked)
library(landscapemetrics)
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
# iucn <- dplyr::bind_rows(iucn, cewa)
# rm(cewa)
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



for(fs in 1:nrow(focal_species)){
  
  # create output holder
  holder <- data.frame(AOU = 0, Year = 0, state_route_grp = 0, developed = 0, BCR = 0)
  counter = 0
  
  # get locations where focal species detected
  bbs_species <- bbs_observations[AOU %in% focal_species[fs,AOU],]
  
  # get routes where focal was not detected; keep one row per routeDataID and make count = 0
  no_species <- bbs_observations[!(RouteDataID %in% bbs_species$RouteDataID),]
  no_species <- no_species[!duplicated(RouteDataID)]
  # make all stop data zero
  stops <- which(str_detect(names(no_species), 'Stop'))
  no_species[, (stops) := lapply(.SD, FUN = function(x) ifelse(x > 0, 0, 0)), .SDcols = stops]
  no_species[, route_counts := 0]
  no_species[, route_binom := 0]
  
  no_species[, AOU := focal_species[fs,AOU]]
  
  # bring detect/no-detect back together
  bbs_species <- rbindlist(list(bbs_species, no_species))
  
  rm(no_species)
  
  # get another unique ID
  bbs_species[, state_route_grp := .GRP, by = .(StateNum, Route)]
  
  for(yy in min_lc_year:max_lc_year){
    cat('species', fs, 'of', nrow(focal_species), 'year', yy, 'starting.\n\n')
    
    # subset down to single year
    bbs_annual <- bbs_species[Year == yy,]
    
    
    
    # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    # get focal species IUCN range map------------------------------------------
    # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    if(focal_species[fs, 'Latin name'] == "Setophaga cerulea"){
      focal_iucn <- cewa %>%
        filter(LEGEND == "Extant (breeding)" | LEGEND == "Extant (resident)") %>%
        st_transform(crs = 'ESRI:102008')
    } else{
      focal_iucn <- iucn %>% 
        filter(BINOMIAL == focal_species[fs, `Latin name`]) %>%
        filter(LEGEND == "Extant (breeding)" | LEGEND == "Extant (resident)") %>%
        st_transform(crs = 'ESRI:102008')
    }

    # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    # reduce routes down to IUCN --------------------------------------------
    # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    
    # make bbs_species spatial
    tmp <- st_as_sf(bbs_annual, coords = c('Longitude', 'Latitude'), crs = 4326) %>%
      st_transform(., crs = "ESRI:102008") %>%
      filter(!duplicated(geometry))
    
    tb <- sapply(st_intersects(tmp, focal_iucn), function(z) if (length(z)==0) NA_integer_ else z[1])
    tmp <- tmp[!is.na(tb),]
    rm(tb)
    
    # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    # get environmental information around each BBS route-----------------------
    # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    
    # get bbs_species in meters
    tmp <- st_transform(tmp, crs = "ESRI:102008")
    
    buff <- st_buffer(tmp, dist = 5000)
    
    # align CRS of landcover and buff
    buff <- st_transform(buff, crs = st_crs(landcover))
    
    
    # loop through buffers of each route, crop landcover data, summarize patch size
    for(b in 1:nrow(buff)){
      counter = counter + 1
      this_buff <- buff[b,]
      cropped_lc <- crop(x = landcover[[paste0("y", yy)]],y = this_buff)
      
      
      # make forest raster
      ## from-to-becomes
      # classify the values into three groups 
      # all values >= 0 and <= 0.25 become 1, etc.
      m <- c(0, 12, 0, #0 to 12 -> 0
             13, 13, 1, # class 13 -> 1
             14, 255, 0) # greater than is 0
      rclmat <- matrix(m, ncol=3, byrow=TRUE)
      devel_lc <- reclassify(cropped_lc, rclmat, include.lowest=TRUE)
      
      # calculate proportion forest, store results, and add to holder
      prop_devel <- sum(values(devel_lc)) / length(values(devel_lc))
      prop_devel <- tibble( AOU = this_buff$AOU,Year = this_buff$Year, state_route_grp = this_buff$state_route_grp, 
                            developed = prop_devel,BCR = buff$BCR)
      holder <- bind_rows(holder, prop_devel)
      
    }
    cat('species', fs, 'of', nrow(focal_species), 'year', yy, 'completed.\n\n')
    
  }
  
  save(holder, file = paste0(here('Results/SoE_BBS/prop_devel'), '/SOE_BBS_prop_devel_species_AOU_',focal_species[fs,AOU], '.Rdata'))
  
}





