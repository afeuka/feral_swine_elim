### Title: Extract NLCD data from watersheds
### Author: Abbey Feuka
### Date: 19MAR25
### Notes:
library(sf)
library(terra)
library(tidyverse)

occ_scale <- "watersheds"

study_site_grid <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/huc10_siteIDs_cond_EA.shp")
study_site_grid <- study_site_grid %>% rename(elim_area_idx=elm_r_d,
                                              Area_Name=Area_Nm)

nlcd <- rast("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/NLCD_MO/NLCD_2019_Land_Cover_L48_20210604_jj2pnCpncP4xHk8rj4R5.tiff")
study_site <- st_union(study_site_grid)

study_site_v <- st_transform(study_site,crs(nlcd))
study_site_v <- vect(study_site_v)
nlcd.crop <- crop(nlcd,study_site_v)

nlcd.crop <- project(nlcd.crop,"epsg:3601")
nlcd.crop[is.nan(nlcd.crop)] <- 0

nlcd.water <- nlcd.crop
values(nlcd.water)[values(nlcd.water)==11] <- 1
values(nlcd.water)[values(nlcd.water)!=1] <- 0

study_site_grid <- st_transform(study_site_grid,"epsg:3601")
# plot(study_site_grid)
water <- terra::extract(nlcd.water,study_site_grid,fun=mean,na.rm=T) %>% 
  rename(SiteID=ID,water=Layer_1) %>% 
  right_join(study_site_grid,by="SiteID") %>% 
  st_as_sf()# %>% select(-c(neighbors)) 

st_write(water,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/so_mo_water_nlcd2019_",occ_scale,".shp"),append=F)

nlcd.developed <- nlcd.crop
values(nlcd.developed)[values(nlcd.developed)%in%c(21:24)] <- 1
values(nlcd.developed)[values(nlcd.developed)!=1] <- 0

developed <- terra::extract(nlcd.developed,
                            study_site_grid,fun=mean) %>% 
  rename(SiteID=ID,developed=Layer_1) %>% 
  right_join(study_site_grid,by="SiteID") %>% 
  st_as_sf()# %>% select(-c(neighbors)) 
# select(-c(neighbors,FID,SITE,geometry))
st_write(developed,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/so_mo_developed_nlcd2019_",occ_scale,".shp"),append=F)

nlcd.forest <- nlcd.crop
values(nlcd.forest)[values(nlcd.forest)%in%c(41:43)] <- 1
values(nlcd.forest)[values(nlcd.forest)!=1] <- 0

forest <- terra::extract(nlcd.forest,study_site_grid,fun=mean) %>%
  rename(SiteID=ID,forest=Layer_1) %>%
  right_join(study_site_grid,by="SiteID") %>%
  st_as_sf() #%>% select(-c(neighbors)) 
st_write(forest,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/so_mo_forest_nlcd2019_",occ_scale,".shp"),append=F)

nlcd.wetland <- nlcd.crop
values(nlcd.wetland)[values(nlcd.wetland)%in%c(90,95)] <- 1
values(nlcd.wetland)[values(nlcd.wetland)!=1] <- 0

wetland <- terra::extract(nlcd.wetland,study_site_grid,fun=mean) %>% 
  rename(SiteID=ID,wetland=Layer_1) %>% 
  right_join(study_site_grid,by="SiteID") %>% 
  # select(-c(neighbors,FID,SITE,geometry))
  st_as_sf() #%>% select(-c(neighbors)) 
st_write(wetland,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/so_mo_wetland_nlcd2019_",occ_scale,".shp"),append=F)

nlcd.agri <- nlcd.crop
values(nlcd.agri)[values(nlcd.agri)%in%c(81,82)] <- 1
values(nlcd.agri)[values(nlcd.agri)!=1] <- 0

agri <- terra::extract(nlcd.agri,study_site_grid,fun=mean) %>% 
  rename(SiteID=ID,agri=Layer_1) %>% 
  right_join(study_site_grid,by="SiteID") %>% 
  st_as_sf() #%>% select(-c(neighbors)) 
st_write(agri,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/so_mo_agri_nlcd2019_",occ_scale,".shp"),append=F)

nlcd_siteid <- full_join(water %>% st_drop_geometry(),
                         developed%>% st_drop_geometry()) %>% 
  full_join(agri%>% st_drop_geometry()) %>% 
  full_join(wetland%>% st_drop_geometry()) %>% 
  full_join(forest%>% st_drop_geometry()) %>%
  full_join(study_site_grid) %>% 
  st_drop_geometry()#%>% dplyr::select(-neighbors)

# save(nlcd_siteid,file=paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Data/Landscape Covariates/nlcd_",occ_scale,".RData"))

saveRDS(nlcd_siteid,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/feral_swine_elim/Data/nlcd_",occ_scale,".rds"))
