library(sf)
library(raster)
library(MODIS) #devtools::install_github("MatMatt/MODIS")
library(exactextractr)
library(viridis)
library(tidyverse)
library(here)
library(lubridate)
library(fasterize)

modisCRS <- paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                  "+a=6371007.181 +b=6371007.181 +units=m +no_defs")


bcr <- read_sf(here("Data/Bird Conservation Regions/bcr_terrestrial_shape/BCR_Terrestrial_master.shp"))  %>% 
  st_transform(crs = modisCRS)

if(!file.exists(here('Data/MODIS/modis/modis_mcd12q1_umd_2001.tif'))){



begin_year <- format(min(ebird$observation_date), "%Y.01.01")

tifs <- runGdal(product = "MCD12Q1", collection = "006", SDSstring = "01", 
                extent = bcr, 
                begin = transDate(begin = '2001.01.01', end = '2001.01.01')$beginDOY, 
                end = transDate(begin = '2018.01.01', end = '2018.01.01')$beginDOY, 
                outDirPath = here('Data/MODIS'), job = "modis",
                MODISserverOrder = "LPDAAC") %>% 
  pluck("MCD12Q1.006") %>% 
  unlist()

# rename tifs to have more descriptive names
new_names <- format(as.Date(names(tifs)), "%Y") %>% 
  sprintf("modis_mcd12q1_umd_%s.tif", .) %>% 
  file.path(dirname(tifs), .)
file.rename(tifs, new_names)

# load the landcover data
landcover <- list.files(here("Data/MODIS/modis"), "^modis_mcd12q1_umd", 
                        full.names = TRUE) %>% 
  stack()

# label layers with year
landcover <- names(landcover) %>% 
  str_extract("(?<=modis_mcd12q1_umd_)[0-9]{4}") %>% 
  paste0("y", .) %>% 
  setNames(landcover, .)

} else{
  # load the landcover data
  landcover <- list.files(here("Data/MODIS/modis"), "^modis_mcd12q1_umd", 
                          full.names = TRUE) %>% 
    stack()
  
  # label layers with year
  landcover <- names(landcover) %>% 
    str_extract("(?<=modis_mcd12q1_umd_)[0-9]{4}") %>% 
    paste0("y", .) %>% 
    setNames(landcover, .)
}
