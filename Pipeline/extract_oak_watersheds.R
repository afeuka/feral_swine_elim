library(terra)
library(exactextractr)
library(tidyverse)

oak_mo <- rast("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/gap2011_mo/gap_2011_mo_oak.tif")

if(abund_scale=="watersheds"){
  
  oak_ssg <- exact_extract(oak_mo, study_site_grid, coverage_area = TRUE)
  
  oak_siteid <- data.frame(oak_km =unlist(lapply(1:length(oak_ssg),function(i){
    sum(oak_ssg[[i]]$coverage_area[oak_ssg[[i]]$value==31])/1e6
  })))
  
  oak_siteid$SiteID <- study_site_grid$SiteID
  oak_siteid <- oak_siteid %>% 
    left_join(study_site_grid %>% st_drop_geometry() %>% select(SiteID,area_km))
  
  oak_siteid$prop_oak <- oak_siteid$oak_km/oak_siteid$area_km
  
}  
saveRDS(oak_siteid,"C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/feral_swin_elim/Data/Model Ready Data/oak_siteid.rds")