# compare hypotheses for species with SoE trends
library(here)
library(dplyr)
library(nlme)
library(data.table)
library(sf)
library(tmap)
library(leaflet)
library(mapedit) # map edit has editMap to interactively create a polygon (bounding box)

# EDIT: DEC 2023: Dr. Rob Fletcher recommended Hand-wing index in place of wing length
hwi <- read.csv(here('Data/Global-HWI-v1.1/catherinesheard-Global-HWI-1533081/Dataset HWI 2020-04-10.csv'))


# AVONET morphological traits
traits <- read.csv(here('Data/AVONET/ELEData/ELEData/TraitData/AVONET_Raw_Data.csv'))

# some times a given species is reported more than once, so I am using average values here
traits <- traits %>%
  group_by(Species1_BirdLife) %>%
  summarize(mean.wing = mean(Wing.Length, na.rm = T),
            mean.tarsus = mean(Tarsus.Length, na.rm = T),
            mean.kipps = mean(Kipps.Distance, na.rm = T))

# bring in complete species list and focal species list (clark rushing 2019-2020 papers)
species_list <- fread(here('Data/BBS/SpeciesList.csv'))
focal_species <- fread(here('Data/BBS/focal_species.csv'))

# get AOU codes from species list to filter observations
focal_species <- merge(focal_species, species_list[,c('genus_species', "AOU")], 
                       by.x = 'Latin name', by.y = 'genus_species')
rm(species_list)

# add in info from AviBase, Bird et al 2020
focal_species$body.mass <- c(21.6, 2159, 285, 12.6, 14, 14.16, 50.1, 19.44, 18.99, 69.5, #last is M. carolinus
                             71.6, 48.5, 19.9, 14.69, 40.03, 29.13, 9.99, 5.8, 14.3, 9.04,  # last is S. cerulea
                             10.54, 7.64, 9.69, 27.5, 10.2, 26.18, 12.5,18.96, 68.8, # last is T rufum
                             8.74, 18, 11.4)
focal_species$clutch <- c(6, 2, 4.5, 2.5, 4.5, 4.5, 3.5, 4, 3, 5, #last is M. carolinus
                          5.5, 4, 5, 3.5, 3.5, 4, 5.5, 3, 5, 4, # last is S. cerulea
                          4, 4, 4, 4.5, 5, 4, 4.5, 4.5,4, # last is T rufum
                          4.5,4, 4)

focal_species <- merge(focal_species, y = traits[, c('mean.wing', 'Species1_BirdLife')], by.x = 'Latin name', by.y = 'Species1_BirdLife')
focal_species <- merge(focal_species, y = hwi, by.x = 'Latin name', by.y = 'IUCN.name')


# read in results files
soe_results <- list.files(path=here("Results/SoE_BBS/BCRs/Jan2024"), pattern = '*.Rdata', full.names = T) %>%
  purrr::map_df(~ get(load(file = .x)))

# remove models that have not converged
soe_results <- soe_results %>% 
  filter(converged == 1)

# weird duplicates
inds <- soe_results %>%
  group_indices(AOU, BCR, radius)
soe_results$index <- inds
soe_results <- soe_results[!duplicated(soe_results$index),]
soe_results$index <- NULL
rm(inds)

# get radius in KM
soe_results <- soe_results %>%
  mutate(rad_km = radius / 1000) 

# Create a summary table of group counts
inds <- soe_results %>%
  group_indices(AOU, BCR)
soe_results$group_index <- inds
rm(inds)

group_counts <- soe_results %>% 
  group_by(group_index) %>% 
  summarize(count = n())

# Filter the original data frame based on the count condition
soe_results <- soe_results %>% 
  filter(group_index %in% group_counts$group_index[group_counts$count >= 10])


# identify SoE for all species and calculate weights
soe_results <- soe_results %>% 
  group_by(AOU, BCR) %>%
  mutate(minAIC = min(AIC)) %>%
  mutate(delAIC = AIC - minAIC) %>%
  mutate(relLik = exp(-0.5 * delAIC)) %>%
  mutate(aicweight = relLik / sum(relLik)) %>%
  ungroup()

# Sample and summarize within groups
soe_results <- soe_results %>%
  group_by(AOU, BCR) %>%
  summarise(
    y_med = median(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight))),
    y_sd = sd(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight))),
    y_min = min(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight))),
    y_max = max(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight)))
  )



# Grab species information - add it to Soe_results
soe_results <- merge(soe_results, focal_species)

# deal with column and bird names
soe_results <- rename(soe_results, common_name = `Common name`, latin_name = `Latin name`)
soe_results$common_name[soe_results$common_name == "Swainson\x92s warbler"] <- 'Swainsons warbler'


# bring in BCRS
bcr <- st_read(here('Data/Bird Conservation Regions/bcr_terrestrial_shape/BCR_Terrestrial_master.shp'))

bcr <- bcr %>%
  filter(BCR %in% unique(soe_results$BCR))
#bcr <- bcr[st_is_valid(bcr),]



# plotting
tmap_style('natural')
data('World')
World <- World %>%
  filter(continent == 'North America')

# get bounding box

# Create a basic map
m <- leaflet() %>% addTiles()
# Enable map editing
m <- editMap(m, layerId = "rectangle", options = mapeditOptions(
  featureGroup = "rectangle",
  mode = "edit",
  draw = list(
    rectangle = TRUE
  ),
  edit = list(
    remove = TRUE
  )
))
# the bounding box is stored within the m object - see below
us_bb <- st_bbox(m$drawn$geometry)


out <- list()
counter = 0
for(a in unique(soe_results$AOU)){
  counter = counter + 1
  spp <- soe_results %>%
    filter(AOU == a)
  
  t <- merge(bcr, spp[, c('AOU', 'y_med', 'BCR')], by = 'BCR')
  p <- tm_shape(World, bbox = us_bb) +
    tm_polygons() +
    tm_shape(t) +
    tm_polygons('y_med', palette = 'Reds', title = 'Scale of Effect (km)',
               breaks = c(seq(0,6,by = 1))) +
    tm_layout(panel.labels = spp %>% filter(!duplicated(common_name)) %>% pull(common_name),
               panel.label.size = 0.5) +
    tmap_options(check.and.fix = T) +
    tm_legend(show = F)
  out[[counter]] <- p
}

p <- tm_shape(World, bbox = us_bb) +
  tm_polygons() +
  tm_shape(t) +
  tm_polygons('y_med', palette = 'Reds', title = 'Scale of Effect (km)',
              breaks = c(seq(0,6,by = 1))) +
  tm_layout(panel.labels = spp %>% filter(!duplicated(common_name)) %>% pull(common_name),
            panel.label.size = 0.5, legend.only = T) +
  tmap_options(check.and.fix = T)
out[[counter+1]] <- p
final <- tmap_arrange(ncol = 6, nrow = 6,
             out)
tmap_save(tm = final, filename = here('Results/Figures/scales_over_space_v02_jan2024.png'))

# save summaries for supplementary info

hold <- soe_results %>%
  group_by(AOU) %>%
  mutate(n_bcr = length(unique(BCR))) %>%
  ungroup() %>%
  filter(!duplicated(AOU)) %>%
  arrange(y_med) %>%
  mutate(y_sd = round(y_sd, 1)) %>%
  select(common_name, y_med, y_sd, y_min, y_max, n_bcr)

write.csv(hold, row.names = F, file = here('Results/soe_over_space_jan2024_tables2.csv'))
