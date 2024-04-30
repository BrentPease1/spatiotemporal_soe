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
cewa$BINOMIAL <- 'Setophaga cerulea'
# combine

iucn <- rbind(iucn %>% select(BINOMIAL, LEGEND), cewa %>% select(BINOMIAL, LEGEND))
rm(cewa)

# transform Bird Conservation Regions
bcr <- bcr %>%
  st_transform(crs = 'ESRI:102008')
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# fix landcover names
names(landcover) <- paste0('y', 2001:2018)

# set up species loop

for(fs in 5:nrow(focal_species)){
  
  # create output holder
  holder <- data.frame(AOU = 0, BCR = 0, radius = 0, AIC = 0, converged = 0)
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
  
  # get observations restricted to 2001 - 2018
  bbs_species <- bbs_species[Year >= 2001 & Year <= 2018,]
  
  
  # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  # get focal species IUCN range map------------------------------------------
  # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  focal_iucn <- iucn %>% 
    filter(BINOMIAL == focal_species[fs, `Latin name`]) %>%
    filter(LEGEND == "Extant (breeding)" | LEGEND == "Extant (resident)") %>%
    st_transform(crs = 'ESRI:102008')
  
  # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  # reduce routes down to IUCN --------------------------------------------
  # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  
  # make bbs_species spatial
  tmp <- st_as_sf(bbs_species, coords = c('Longitude', 'Latitude'), crs = 4326) %>%
    st_transform(., crs = "ESRI:102008") %>%
    filter(!duplicated(geometry))
  
  tb <- sapply(st_intersects(tmp, focal_iucn), function(z) if (length(z)==0) NA_integer_ else z[1])
  tmp <- tmp[!is.na(tb),]
  rm(tb)
  
  # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  # identify BCRs with routes
  # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
  
  tb <- sapply(st_intersects(bcr, tmp), function(z) if (length(z)==0) NA_integer_ else z[1])
  focal_bcr <- bcr[!is.na(tb),]
  rm(tb)
  
  # loop through BCRs calculate SoE
  these_bcrs <- unique(focal_bcr$BCR)
  
  for(r in these_bcrs){
    cat('species', fs, 'of', nrow(focal_species), 'BCR', which(these_bcrs == r),'of', length(these_bcrs), 'starting.\n\n')
    
    # isolate the shapefile for the focal BCR
    one_bcr <- focal_bcr %>%
      filter(BCR == r)
    
    # get routes within loop bcr
    tb <- sapply(st_intersects(tmp, one_bcr), function(z) if (length(z)==0) NA_integer_ else z[1])
    tmp_bcr <- tmp[!is.na(tb),]
    rm(tb)
    
    # skip if too few of routes
    if(nrow(tmp_bcr) < 10){
      "no good"
    } else{
      # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      # get environmental information around each BBS route-----------------------
      # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
      
      # get bbs_species in meters (just in case)
      tmp_bcr <- st_transform(tmp_bcr, crs = "ESRI:102008")
      
      # generate a range of buffers around bbs routes (meters)
      radii <- seq(from = 100, to = 5000, by = 100)
      
      
      # loop through values, create buffers, extract landcover
      for(radius in radii){
        
        buffs <- st_buffer(tmp_bcr, dist = radius)
        
        # align CRS of landcover and buffs
        buffs <- st_transform(buffs, crs = st_crs(landcover))
        
        
        # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
        # get landcover information in buffers--------------------------------------
        # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
        lc_extract_pred <- landcover[[paste0("y", 2018)]] %>% 
          exact_extract(buffs, progress = T) %>% 
          map(~ count(., landcover = value)) %>% 
          tibble(state_route_grp = buffs$state_route_grp, data = .) %>% 
          unnest(data)
        
        # calculate the percent for each landcover class
        pland_pred <- lc_extract_pred %>% 
          count(state_route_grp, landcover) %>% 
          group_by(state_route_grp) %>% 
          mutate(pland = n / sum(n)) %>% 
          ungroup() %>% 
          select(-n) %>% 
          # remove NAs after tallying so pland is relative to total number of cells
          filter(!is.na(landcover))
        
        
        # convert names to be more descriptive
        lc_names <- tibble(landcover = 0:15,
                           lc_name = c("pland_00_water", 
                                       "pland_01_evergreen_needleleaf", 
                                       "pland_02_evergreen_broadleaf", 
                                       "pland_03_deciduous_needleleaf", 
                                       "pland_04_deciduous_broadleaf", 
                                       "pland_05_mixed_forest",
                                       "pland_06_closed_shrubland", 
                                       "pland_07_open_shrubland", 
                                       "pland_08_woody_savanna", 
                                       "pland_09_savanna", 
                                       "pland_10_grassland", 
                                       "pland_11_wetland", 
                                       "pland_12_cropland", 
                                       "pland_13_urban", 
                                       "pland_14_mosiac", 
                                       "pland_15_barren"))
        
        pland_pred <- pland_pred %>% 
          inner_join(lc_names, by = "landcover") %>% 
          arrange(landcover) %>% 
          select(-landcover)
        
        # tranform to wide format, filling in implicit missing values with 0s
        pland_pred <- pland_pred %>% 
          pivot_wider(names_from = lc_name, 
                      values_from = pland, 
                      values_fill = list(pland = 0)) %>% 
          mutate(year = 2018) %>% 
          select(state_route_grp, year, everything())
        
        # join to buffers
        buffs <- buffs %>%
          left_join(pland_pred %>% select(-year), by = "state_route_grp")
        
        # make sure columns are there
        # happens if specific landcover doesn't exist within buffers
        if(!("pland_12_cropland" %in% names(buffs))){
          buffs$pland_12_cropland <- 0
        }
        if(!("pland_04_deciduous_broadleaf" %in% names(buffs))){
          buffs$pland_04_deciduous_broadleaf <- 0
        }
        if(!("pland_05_mixed_forest" %in% names(buffs))){
          buffs$pland_05_mixed_forest <- 0
        }
        if(!("pland_13_urban" %in% names(buffs))){
          buffs$pland_13_urban <- 0
        }
        # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
        # prep columns for model fitting
        # -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
        
        # get julian date
        j1 <- do.call(paste, list(buffs$Month, buffs$Day, buffs$Year))
        j2 <- as.Date(j1, format=c("%m %d %Y"))
        j3 <- as.numeric(format(j2, "%j"))
        buffs <- buffs %>%
          mutate(julian = j3)
        rm(j1, j2, j3)
        
        # get max wind and duration
        buffs <- buffs %>%
          mutate(maxWind = apply(buffs[,c('StartWind', 'EndWind')] %>% st_drop_geometry, 1, FUN = function(x) max(x))) %>%
          mutate(duration = EndTime - StartTime)
        
        # fix assistant column
        buffs <- buffs %>%
          mutate(Assistant = ifelse(Assistant == 'NULL', 0, Assistant))
        
        # average temperature during observations
        tm <- buffs %>%
          st_drop_geometry() %>%
          mutate(StartTemp = as.numeric(StartTemp), 
                 EndTemp = as.numeric(EndTemp))
        tm <- tm %>%
          mutate(meanTemp = apply(tm[,c('StartTemp', 'EndTemp')], 1, FUN = function(x) mean(x, na.rm = T)))
        
        buffs <- buffs %>%
          mutate(meanTemp = tm$meanTemp)
        
        rm(tm)
        
        # get ready for unmarked insanity (repping dectection covs 50 times)
        julian <- matrix(data = rep(buffs$julian, times = 50), nrow = length(buffs$julian), ncol = 50, byrow = F)
        maxWind <- matrix(data = rep(as.character(buffs$maxWind), times = 50), nrow = length(buffs$maxWind), ncol = 50, byrow = F)
        duration <- matrix(data = rep(buffs$duration, times = 50), nrow = length(buffs$duration), ncol = 50, byrow = F)
        assistant <- matrix(data = rep(as.character(buffs$Assistant), times = 50), nrow = length(buffs$Assistant), ncol = 50, byrow = F)
        meanTemp <- matrix(data = rep(buffs$meanTemp, times = 50), nrow = length(buffs$meanTemp), ncol = 50, byrow = F)
        
        
        # summarize counts across stops
        # first get stop columns
        stops <- which(str_detect(names(buffs), 'Stop'))
        
        # drop the stop columns
        y <- buffs[, stops] %>%
          st_drop_geometry() %>%
          as.matrix()
        
        if(any(names(pland_pred) == 'pland_05_mixed_forest')){
          # package it up
          umf <- unmarkedFrameOccu(
            y = y,                                            
            siteCovs = data.frame(decid = buffs$pland_04_deciduous_broadleaf, 
                                  mixed = buffs$pland_05_mixed_forest,
                                  crop = buffs$pland_12_cropland,
                                  urban = buffs$pland_13_urban),  
            obsCovs = list(julian = julian, 
                           maxWind = maxWind,
                           duration = duration,
                           assistant = assistant,
                           meanTemp = meanTemp))        
          # summary(umf)
          
          # test the water
          #occ1 <- occu(~1~1, data=umf)
          
          # load er up
          occ.mod <- occu(~scale(julian) + maxWind + scale(duration) + assistant + scale(meanTemp)
                          ~ decid + mixed + crop + urban,
                          data = umf)
        } else{
          # package it up
          umf <- unmarkedFrameOccu(
            y = y,                                            
            siteCovs = data.frame(decid = buffs$pland_04_deciduous_broadleaf,
                                  crop = buffs$pland_12_cropland,
                                  urban = buffs$pland_13_urban),  
            obsCovs = list(julian = julian, 
                           maxWind = maxWind,
                           duration = duration,
                           assistant = assistant,
                           meanTemp = meanTemp))        
          # summary(umf)
          
          # test the water
          #occ1 <- occu(~1~1, data=umf)
          
          # load er up
          occ.mod <- occu(~scale(julian) + maxWind + scale(duration) + assistant + scale(meanTemp)
                          ~ decid + crop + urban,
                          data = umf)
        }
        
        if(any(is.na(occ.mod@estimates[1]@covMat))){
          
        }
        
        # store results
        counter = counter + 1
        holder[counter,'AOU'] <- focal_species[fs,AOU]
        holder[counter,'BCR'] <- r
        holder[counter, 'radius'] <- radius
        holder[counter, 'AIC'] <- occ.mod@AIC
        holder[counter, 'converged'] <- ifelse(any(is.na(occ.mod@estimates[1]@covMat)), 0, 1)
        
        cat('species', fs, 'of', nrow(focal_species), 'bcr', which(these_bcrs == r), 'of', length(these_bcrs), 'radius', 
            which(radii == radius), 'of', length(radii), 'completed.\n\n')
        
        save(holder, file = paste0(here('Results/SoE_BBS/BCRs/Jan2024'), '/SOE_BBS_holder_species_AOU_',focal_species[fs,AOU], '.Rdata'))
        
      }
    }
  }
}





