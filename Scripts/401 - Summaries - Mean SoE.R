# Summary statistics of mean scales of effect
library(here)
library(dplyr)
library(data.table)
library(sf)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# TIME #####

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --



# get SoE for all species
# read in results files
soe_results <- list.files(path=here("Results/SoE_BBS/Year"), pattern = '*.Rdata', full.names = T) %>%
  purrr::map_df(~ get(load(file = .x)))

# identify SoE for all species
soe_results <- soe_results %>% group_by(AOU, year) %>%
  filter(AIC == min(AIC))

# get radius in KM
soe_results <- soe_results %>%
  mutate(rad_km = radius / 1000)

# weird duplicates
inds <- soe_results %>%
  group_indices(AOU, year)

soe_results$index <- inds

soe_results <- soe_results[!duplicated(soe_results$index),]
rm(inds)

# bring in complete species list and focal species list 
species_list <- fread(here('Data/BBS/SpeciesList.csv'))
focal_species <- fread(here('Data/BBS/focal_species.csv'))

# get AOU codes from species list to filter observations
focal_species <- merge(focal_species, species_list[,c('genus_species', "AOU")], 
                       by.x = 'Latin name', by.y = 'genus_species')
rm(species_list)

soe_results <- merge(soe_results, focal_species)

soe_results <- rename(soe_results, common_name = `Common name`)

soe_results <- soe_results %>%
  mutate(common_name = gsub("'", "", common_name)) %>%
  mutate(common_name = gsub("<92>", "", common_name))

soe_results$year <- soe_results$year - 2001


mean_temporal <- soe_results %>%
  group_by(common_name) %>%
  summarize(mean_soe = round(mean(rad_km),1),
            sd_soe = round(sd(rad_km), 1),
            min_soe = round(min(rad_km), 1),
            max_soe = round(max(rad_km),1)) %>%
  arrange(mean_soe) %>%
  ungroup() %>%
  rename(`Common Name` = common_name, `Mean SoE` = mean_soe, 
         `SD` = sd_soe, `Min` = min_soe, `Max`= max_soe)
readr::write_csv(mean_temporal, here('Results/soe_over_time.csv'))
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# SPACE #####

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# get SoE for all species
# read in results files
soe_results <- list.files(path=here("Results/SoE_BBS/BCRs"), pattern = '*.Rdata', full.names = T) %>%
  purrr::map_df(~ get(load(file = .x)))

# identify SoE for all species
soe_results <- soe_results %>% group_by(AOU, BCR) %>%
  filter(AIC == min(AIC))

# get radius in KM
soe_results <- soe_results %>%
  mutate(rad_km = radius / 1000)

# weird duplicates
inds <- soe_results %>%
  group_indices(AOU, BCR)

soe_results$index <- inds

soe_results <- soe_results[!duplicated(soe_results$index),]
rm(inds)


# bring in complete species list and focal species list (clark rushing 2019-2020 papers)
species_list <- fread(here('Data/BBS/SpeciesList.csv'))
focal_species <- fread(here('Data/BBS/focal_species.csv'))

# get AOU codes from species list to filter observations
focal_species <- merge(focal_species, species_list[,c('genus_species', "AOU")], 
                       by.x = 'Latin name', by.y = 'genus_species')
rm(species_list)

soe_results <- merge(soe_results, focal_species)

soe_results <- rename(soe_results, common_name = `Common name`)

soe_results <- soe_results %>%
  mutate(common_name = gsub("'", "", common_name)) %>%
  mutate(common_name = gsub("<92>", "", common_name))


# bring in BCRS
bcr <- st_read(here('Data/Bird Conservation Regions/bcr_terrestrial_shape/BCR_Terrestrial_master.shp'))

bcr <- bcr %>%
  filter(BCR %in% unique(soe_results$BCR))


mean_spatial <- soe_results %>%
  group_by(common_name) %>%
  summarize(mean_soe = round(mean(rad_km),1),
            sd_soe = round(sd(rad_km), 1),
            min_soe = round(min(rad_km), 1),
            max_soe = round(max(rad_km),1)) %>%
  arrange(mean_soe) %>%
  ungroup() %>%
  rename(`Common Name` = common_name, `Mean SoE` = mean_soe, 
         `SD` = sd_soe, `Min` = min_soe, `Max`= max_soe)

readr::write_csv(mean_spatial, here('Results/soe_over_space.csv'))
