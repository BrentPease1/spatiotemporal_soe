# Summary statistics of mean scales of effect
library(here)
library(dplyr)
library(data.table)
library(sf)
library(ggplot2)
library(lme4)
library(colorspace)
library(ggeffects)
library(patchwork)
library(effects)
library(RColorBrewer)
glmerControl(optCtrl = list(maxfun = 100000, tolPwrss = 1e-4, optimizer = "Nelder_Mead"))
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

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

# TIME #####

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

# get SoE for all species
# read in results files
soe_results <- list.files(path=here("Results/SoE_BBS/Year"), pattern = '*.Rdata', full.names = T) %>%
  purrr::map_df(~ get(load(file = .x)))

# weird duplicates
inds <- soe_results %>%
  group_indices(AOU, year, radius)
soe_results$index <- inds
soe_results <- soe_results[!duplicated(soe_results$index),]
soe_results$index <- NULL

# get radius in KM
soe_results <- soe_results %>%
  mutate(rad_km = radius / 1000)

# identify SoE for all species and calculate weights
soe_results <- soe_results %>% 
  group_by(AOU, year) %>%
  mutate(minAIC = min(AIC)) %>%
  mutate(delAIC = AIC - minAIC) %>%
  mutate(relLik = exp(-0.5 * delAIC)) %>%
  mutate(aicweight = relLik / sum(relLik)) %>%
  ungroup()

# Sample and summarize within groups
soe_results <- soe_results %>%
  group_by(AOU, year) %>%
  summarise(
    y_med = median(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight))),
    y_sd = sd(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight)))
  )


# Grab species information - add it to Soe_results
soe_results <- merge(soe_results, focal_species)

# deal with column and bird names
soe_results <- rename(soe_results, common_name = `Common name`, latin_name = `Latin name`)
soe_results$common_name[soe_results$common_name == "Swainson\x92s warbler"] <- 'Swainsons warbler'

# make factors 
soe_results$year_factor <- factor(soe_results$year - 2001)
soe_results <- soe_results %>%
  mutate(class_fac = factor(Class))

soe_results <- soe_results %>% 
  mutate(hwi_scaled = scale(HWI))

# run some models :)
soe_weights <- as.vector(1 / soe_results$y_sd^2)
glme0 <- glmer(y_med ~ 1 + (1 | year_factor), 
               data = soe_results, 
               family = Gamma(link = "log"), 
               weights = soe_weights)
glme1 <- glmer(y_med ~ hwi_scaled + (1 | year_factor), 
               data = soe_results, 
               family = Gamma(link = "log"), 
               weights = soe_weights)
glme2 <- glmer(y_med ~ Body.mass..log. + (1 | year_factor), 
               data = soe_results, 
               family = Gamma(link = "log"), 
               weights = soe_weights)
glme3 <- glmer(y_med ~ -1 + hwi_scaled*class_fac + (1 | year_factor), 
               data = soe_results, 
               family = Gamma(link = "log"), 
               weights = soe_weights)
glme4 <- glmer(y_med ~ Body.mass..log.*class_fac + (1 | year_factor), 
               data = soe_results, 
               family = Gamma(link = "log"), 
               weights = soe_weights)
glme5 <- glmer(y_med ~ class_fac + (1 | year_factor), 
               data = soe_results, 
               family = Gamma(link = "log"), 
               weights = soe_weights)

AICcmodavg::aictab(cand.set = list("~1" = glme0, "hwi_scaled" = glme1, "body.mass" = glme2,
                                   "hwi*class" = glme3, "body*class" = glme4, "class" = glme5))

predictions <- ggpredict(glme3, terms = c('hwi_scaled', 'class_fac'))

p1 <- ggplot(predictions, aes(x, predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, linetype = 2) +
  labs(x = "Scaled HWI", y = "Scale of Effect (km)", colour = "Winter Geography") +
  theme_ggeffects() +
  theme(axis.title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  coord_cartesian(ylim = c(0, NA)) + scale_color_brewer(palette = "Dark2")


# SPACE ####
# get SoE for all species
# read in results files
bcr_results <- list.files(path=here("Results/SoE_BBS/BCRs"), pattern = '*.Rdata', full.names = T) %>%
  purrr::map_df(~ get(load(file = .x)))

# weird duplicates
inds <- bcr_results %>%
  group_indices(AOU, radius, BCR)
bcr_results$index <- inds
bcr_results <- bcr_results[!duplicated(bcr_results$index),]
bcr_results$index <- NULL

# get radius in KM
bcr_results <- bcr_results %>%
  mutate(rad_km = radius / 1000)

# identify SoE for all species and calculate weights
bcr_results <- bcr_results %>% 
  group_by(AOU, BCR) %>%
  mutate(minAIC = min(AIC)) %>%
  mutate(delAIC = AIC - minAIC) %>%
  mutate(relLik = exp(-0.5 * delAIC)) %>%
  mutate(aicweight = relLik / sum(relLik)) %>%
  ungroup()

# Sample and summarize within groups
bcr_results <- bcr_results %>%
  group_by(AOU, BCR) %>%
  summarise(
    y_med = median(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight))),
    y_sd = sd(replicate(1000, sample(rad_km, replace = TRUE, prob = aicweight)))
  )


# Grab species information - add it to bcr_results
bcr_results <- merge(bcr_results, focal_species)

# deal with column and bird names
bcr_results <- rename(bcr_results, common_name = `Common name`, latin_name = `Latin name`)
bcr_results$common_name[bcr_results$common_name == "Swainson\x92s warbler"] <- 'Swainsons warbler'

# make factors 
bcr_results <- bcr_results %>%
  mutate(class_fac = factor(Class))

bcr_results <- bcr_results %>% 
  mutate(hwi_scaled = scale(HWI))

bcr_results <- bcr_results %>% 
  mutate(BCR_factor = factor(BCR))

# run some models :)
soe_weights <- as.vector(1 / bcr_results$y_sd^2)
glme0 <- glmer(y_med ~ 1 + (1 | BCR_factor), 
               data = bcr_results, 
               family = Gamma(link = "log"), 
               weights = soe_weights)
glme1 <- glmer(y_med ~ hwi_scaled + (1 | BCR_factor), 
               data = bcr_results, 
               family = Gamma(link = "log"), 
               weights = soe_weights)
glme2 <- glmer(y_med ~ Body.mass..log. + (1 | BCR_factor), 
               data = bcr_results, 
               family = Gamma(link = "log"), 
               weights = soe_weights)
glme3 <- glmer(y_med ~ -1 + hwi_scaled*class_fac + (1 | BCR_factor), 
               data = bcr_results, 
               family = Gamma(link = "log"), 
               weights = soe_weights)
glme4 <- glmer(y_med ~ -1 + Body.mass..log.*class_fac + (1 | BCR_factor), 
               data = bcr_results, 
               family = Gamma(link = "log"), 
               weights = soe_weights)
glme5 <- glmer(y_med ~ class_fac + (1 | BCR_factor), 
               data = bcr_results, 
               family = Gamma(link = "log"), 
               weights = soe_weights)

AICcmodavg::aictab(cand.set = list("~1" = glme0, "hwi_scaled" = glme1, "body.mass" = glme2,
                                   "hwi*class" = glme3, "body*class" = glme4, "class" = glme5))

predictions <- ggpredict(glme4, terms = c('Body.mass..log.', 'class_fac'))

p2 <- ggplot(predictions, aes(x, predicted, color = group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, linetype = 2) +
  labs(x = "Log Body Mass (g)", y = "Scale of Effect (km)", colour = "Winter Geography") +
  theme_ggeffects() +
  theme(axis.title = element_text(size = 18),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  coord_cartesian(ylim = c(0, NA)) + scale_color_brewer(palette = "Paired")

p1 + p2 + plot_annotation(tag_levels = "A")
ggsave(filename = here('Results/Figures/trait_based_dec2023.jpg'), plot = last_plot(), device = 'jpg', width = 16, height = 8, units = 'in', dpi = 300)
