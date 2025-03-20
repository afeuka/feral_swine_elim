### Title: Quantifying aerial ops effort area via track data
### Author: Abbey Feuka
### Date: 21JUN24
### Notes:

library(sf)
library(tidyverse)

setwd("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/Aerial Ops Tracks")

#yearly aerial ops track file names -----------------------
ao_files <- c("./Raw Data/FY2021 AERIAL OPS TRACK DATA.csv",
              "./Raw Data/FY2022 MO AERIAL OPS TRACK DATA.csv",
              "./Raw Data/FY2023 MO AERIAL OPS TRACK DATA.csv",
              "./Raw Data/FY24_MOFlightTracks.gdb")

#function to id last character ------------------------
substrRight <- function(x, n=1){
  substr(x, nchar(x)-n+1, nchar(x))
}

#align files --------------------
ao <- list()

for(i in 1:length(ao_files)){
  if(grepl("csv",ao_files[i])){
    ao[[i]] <- read.csv(ao_files[i])
    colnames(ao[[i]]) <- gsub("\\.","_",colnames(ao[[i]]))
    ao[[i]] <- st_as_sf(ao[[i]],coords=c("Longitude_decimal_","Latitude_decimal_"),crs="epsg:4326")
  } else {
    ao[[i]] <- st_read(ao_files[i])
  }
  
  # if(sum(grepl("Longitud",colnames(ao[[i]])))>0){
    ao[[i]] <- ao[[i]] %>% 
      dplyr::select(-c(grep("Longitud",colnames(ao[[i]])),grep("Latitud",colnames(ao[[i]]))))
  # }
  
  if("DateTime_L"%in%colnames(ao[[i]])){
    ao[[i]] <- ao[[i]] %>% rename(DateTime_Local=DateTime_L,
                                  DateTime_UTC=DateTime_U)
  }

  colnames(ao[[i]])[substrRight(colnames(ao[[i]]))=="_"]<-
  substr(colnames(ao[[i]])[substrRight(colnames(ao[[i]]))=="_"],1,
         nchar(colnames(ao[[i]])[substrRight(colnames(ao[[i]]))=="_"])-1)

  if("Time"%in%colnames(ao[[i]])){
    ao[[i]] <- ao[[i]] %>% dplyr::select(-Time)
  }
  if("Shape"%in%colnames(ao[[i]])){
    ao[[i]] <- ao[[i]] %>% mutate(geometry=Shape) %>% st_drop_geometry() %>% st_as_sf()
  }
  if("Descriptio"%in%colnames(ao[[i]])){
    ao[[i]] <- ao[[i]] %>% rename(Description=Descriptio)
  }
  if("Altitude_f"%in%colnames(ao[[i]])){
    ao[[i]] <- ao[[i]] %>% rename(Altitude_ft=Altitude_f)
  }
  if("Speed_knot"%in%colnames(ao[[i]])){
    ao[[i]] <- ao[[i]] %>% rename(Speed_knots=Speed_knot)
  }
  if("OBJECTID"%in%colnames(ao[[i]])){
    ao[[i]] <- ao[[i]] %>% dplyr::select(-OBJECTID)
  }
  if("DateTime_UTC"%in%colnames(ao[[i]])){
    ao[[i]] <- ao[[i]] %>% dplyr::select(-DateTime_UTC)
    # ao[[i]]$DateTime_UTC <- as.POSIXct(ao[[i]]$DateTime_UTC,format="%Y-%m-%d %H:%M:%S",tz="UTC")
  }
  ao[[i]]$DateTime_Local <- as.POSIXct(ao[[i]]$DateTime_Local,
                                       # format="%Y-%m-%d %H:%M:%S",
                                       tryFormats=c("%Y-%m-%d %H:%M:%S",
                                                    "%m/%d/%Y %H:%M:%S"),
                                       tz="CST6CDT")
  if(!("Altitude_ft"%in%colnames(ao[[i]]))){
    ao[[i]]$Altitude_ft <- NA
  }
  ao[[i]] <- ao[[i]] %>% dplyr::select("Aircraft","Altitude_ft","DateTime_Local","Speed_knots",
                                "Bearing","geometry")
  }

ao <- do.call("rbind",ao)
ao <- ao[!duplicated(ao),]

#filter out ferrying locations (out of state and speeds >25 knots) -----------------
library(tigris)
options(tigris_use_cache = TRUE)
states_sf <- states()
states_nb <- states_sf %>% filter(STUSPS%in%c("MO","KY","TN","IL","IN","AR","OK"))

ao <- ao %>% st_intersection(states_nb %>% st_transform(st_crs(ao)))
ao$year <- lubridate::year(ao$DateTime_Local)

ao <- ao %>% filter(STUSPS%in%c("MO","AR")) %>% 
  filter(Speed_knots<25) %>% 
  dplyr::select(Aircraft,Altitude_ft,DateTime_Local,Speed_knots,Bearing,
                STUSPS,year,geometry)

# summary(ao$Altitude_ft)
# hist(ao$Altitude_ft) #mean altitude is 1000 ft

#take mean of multiple points taken at same time
ao_coords <- st_coordinates(ao)
ao_coords <- cbind.data.frame(ao %>% st_drop_geometry(),
                              ao_coords)
ao_coords_clean <- ao_coords %>% 
  group_by(Aircraft,DateTime_Local) %>% 
  summarise(n=n(),
            Speed_knots=ifelse(n>1,mean(Speed_knots),Speed_knots),
            Bearing=ifelse(n>1,mean(Bearing),Bearing),
            x = ifelse(n>1,mean(X),X),
            y = ifelse(n>1,mean(Y),Y))

ao_clean_sf <- st_as_sf(ao_coords_clean,coords=c("x","y"),crs=st_crs(4326))

ggplot()+
  geom_sf(data=ao,col="black",alpha=0.5)+
  geom_sf(data=ao_clean_sf,col="red",alpha=0.5)

library(amt)
study_site_grid <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/huc10_siteIDs_cond_EA.shp")

ao_coords_clean <- tibble(ao_coords_clean)
ao_coords_clean$day <- as.Date(ao_coords_clean$DateTime_Local)
ao_coords_clean$aircraft_day <- paste(ao_coords_clean$Aircraft,ao_coords_clean$day)

# make amt tracks by aircraft day
ao_amt <- ao_coords_clean |> make_track(.x=x,.y=y,.t=DateTime_Local,
                     id=aircraft_day,crs=st_crs(4326))
ao_amt <- transform_coords(ao_amt,(st_crs(study_site_grid)))

#nest by aircraft day
ao_amt_nest <- ao_amt |> nest(data=-"id")

#calculate steps and step length by aircraft day
ao_amt_nest <- ao_amt_nest |>
  mutate(steps = map(data, function(x) x |>
                       track_resample(rate = minutes(1), tolerance = seconds(15)) |>
                       steps_by_burst()),
         sl =map(data, function(x) x |>
                         step_lengths()))

#filter out large steps - if desired 
# sl_df <- ao_amt_nest %>%
#   select(id,steps,sl) %>%
#   unnest(cols=steps)

# ao_amt_nest <- ao_amt_nest %>%
#   select(id,steps,sl) %>%
#   unnest(cols=steps) %>%
#   filter(is.na(x2_))
  # filter(sl_<quantile(sl_df$sl_,prob=0.95)) %>%
  # nest(data=-"id")

ao_amt_steps <- ao_amt_nest %>%  
  unnest(cols=steps)

#buffer distance from trajectory = tan(angle of shooting) * altitude
buff_dist <- tan(60*pi/180) * 304.8 
units(buff_dist) <- "meters"
mat_list <- mls_list <- buff_list <- list()
for(i in 1:length(ao_amt_nest$steps)){
  if(nrow(ao_amt_nest$steps[[i]])>0){
    mat_list[[i]] <- list()
    for(j in 1:nrow(ao_amt_nest$steps[[i]])){
      mat_list[[i]][[j]] <- as.data.frame(matrix(c(ao_amt_nest$steps[[i]]$x1_[j],
                                                   ao_amt_nest$steps[[i]]$x2_[j],
                                                   ao_amt_nest$steps[[i]]$y1_[j],
                                                   ao_amt_nest$steps[[i]]$y2_[j]),
                                                 2,2))
      colnames(mat_list[[i]][[j]]) <- c("x","y")
      mat_list[[i]][[j]]$burst <- j
      mat_list[[i]][[j]]$dt <- ao_amt_nest$steps[[i]]$dt_[j]
      
      mat_list[[i]][[j]] <- mat_list[[i]][[j]] %>% 
        st_as_sf(coords=c('x','y'),crs=st_crs(study_site_grid)) %>% 
        group_by(burst) %>% 
        summarise(do_union=FALSE,
                  dt=unique(dt)) %>% 
        st_cast("LINESTRING")
    }
  mat_list[[i]] <- do.call("rbind",mat_list[[i]]) %>% 
    st_cast("MULTILINESTRING")
  mat_list[[i]]$aircraft_day <- ao_amt_nest$id[[i]]
  
  mat_list[[i]]$Date <- as.Date(sapply(1:nrow(mat_list[[i]]), function(j) {
    str_split(mat_list[[i]]$aircraft_day[j],"ER ")[[1]][2]
  }))
  mat_list[[i]]$Aircraft <- sapply(1:nrow(mat_list[[i]]), function(j) {
    str_flatten(str_split(mat_list[[i]]$aircraft_day[j]," ")[[1]][1:3])
  })
  
  buff_list[[i]] <- st_buffer(mat_list[[i]],dist=buff_dist) #1000 ft distance
  }
}

ao_buff <- do.call("rbind",buff_list)

#convert seconds to minutes
ao_buff$dt[ao_buff$dt>2] <- ao_buff$dt[ao_buff$dt>2]/60

#combine all bursts per day 
ao_buff_day <- ao_buff %>% 
  group_by(Aircraft,Date) %>% 
  summarise(geometry=st_union(geometry),
            dt_hr=sum(dt)/60,
            year=unique(lubridate::year(Date)))

st_write(ao_buff_day,
         paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/Aerial Ops Tracks/Cleaned/FY21_24_MOFlightTracks_clean_",as.numeric(round(buff_dist)),"m_buff.shp"),append=F)


#load saved data---------
ao_buff_day <- st_read(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/Aerial Ops Tracks/Cleaned/FY21_24_MOFlightTracks_clean_",as.numeric(round(buff_dist)),"m_buff.shp"))

ggplot()+
  # geom_sf(data=states_nb %>% filter(STUSPS=="MO"))+
  geom_sf(data=study_site_grid,fill="white")+
  geom_sf(data=study_site_grid %>% 
            group_by(Area_Name) %>% 
            summarise(geometry=st_union(geometry)),
          aes(col=Area_Name),fill="transparent",lwd=2) +
  geom_sf(data=ao_buff_day,col="blue") +
  guides(fill="none")+
  facet_wrap(.~year)

ggplot()+
  # geom_sf(data=states_nb %>% filter(STUSPS=="MO"))+
  geom_sf(data=study_site_grid,fill="white")+
  # geom_sf(data=study_site_grid %>% 
  #           group_by(Area_Name) %>% 
            # summarise(geometry=st_union(geometry)),
          # aes(col=Area_Name),fill="transparent",lwd=2) +
  geom_sf(data=ao_buff_day %>% filter(year==2023),col="blue") +
  guides(fill="none")

#merge with watershed data ----------------------------
study_site_grid <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/huc10_siteIDs_cond_EA.shp")
study_site_grid <- study_site_grid %>% rename(elim_area_idx=elm_r_d,
                                              Area_Name=Area_Nm,
                                              site_area_km=area_km)
# ao_buff_day_ssg <- ao_buff_day %>% st_intersection(study_site_grid %>% 
#                                         st_transform(st_crs(ao_buff_day)))
# ao_buff_day_ssg <- ao_buff_day_ssg %>% 
#   dplyr::select(Aircraft,Date,SiteID,Area_Name,
#                 elim_area_idx,site_area_km,geometry)

ao_buff_day_ssg <- st_intersection(study_site_grid %>% 
                                     st_transform(st_crs(ao_buff_day)), 
                                  ao_buff_day)

intersect_ssg <- ao_buff_day_ssg %>% 
  mutate(intersect_area = st_area(.),
         intersect_area_km=intersect_area/1e6) %>%   
  dplyr::select(SiteID,Date,intersect_area_km) %>%   
  st_drop_geometry() %>% 
  left_join(study_site_grid %>% st_drop_geometry()) %>% 
  mutate(prop_site_flown = as.numeric(intersect_area_km/site_area_km))

ao_ssg_day <- ao_buff_day_ssg %>% st_drop_geometry() %>% 
  left_join(intersect_ssg) %>% as.data.frame()
ao_ssg_day %>% filter(dt_hr>100)
ao_buff %>% filter(Aircraft=="N332AB-ER" & Date==as.Date("2021-12-08"))

hist(ao_ssg_day$dt_hr)
hist(ao_ssg_day$prop_site_flown)
plot(ao_ssg_day$dt_hr,ao_ssg_day$prop_site_flown)

ao_ssg_day<- ao_ssg_day %>% rename(flight_time_hr=dt_hr)
write.csv(ao_ssg_day,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Model Ready Data/fy21_24_flight_time_ssg_clean_",as.numeric(round(buff_dist)),"m_buff.csv"))

## intersect with elim areas ----------------
ea_merge <- study_site_grid %>% 
  group_by(Area_Name) %>% 
  summarise(geometry=st_union(geometry),
            ea_area_km=st_area(geometry)/1e6)

# ao_ssg <- ao_ssg %>% left_join(ea_merge %>% st_drop_geometry())

#calculate proportion of EA covered 
ao_buff_day_ea <- st_intersection(ea_merge %>% st_transform(st_crs(ao_buff)), 
                ao_buff_day)

intersect_ea <- ao_buff_day_ea %>% 
  mutate(intersect_area = st_area(.),
         intersect_area_km = intersect_area/1e6) %>%   
  dplyr::select(Area_Name,Date,intersect_area_km) %>%   
  st_drop_geometry() %>% 
  left_join(ea_merge %>% st_drop_geometry()) %>% 
  mutate(prop_ea_flown = intersect_area_km/ea_area_km)
# hist(intersect_pct$prop_ea_flown)

#add flight time
ao_ea_day <- ao_buff_day_ea %>% st_drop_geometry() %>% 
  left_join(intersect_ea) %>% as.data.frame()

hist(ao_ea_day$dt_hr)
hist(ao_ea_day$prop_ea_flown)
plot(ao_ea_day$dt_hr,ao_ea_day$prop_ea_flown)

ao_ea_day<- ao_ea_day %>% rename(flight_time_hr=dt_hr)

write.csv(ao_ea_day,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Model Ready Data/fy21_24_flight_time_ea_clean_",as.numeric(round(buff_dist)),"m_buff.csv"))

#merge with farm bill ----------------------------
fb <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/Farm Bill/mo_fb_counties_zones.shp")

ao_buff_fb <- st_intersection(fb %>% st_transform(st_crs(ao_buff_day)), 
                                           ao_buff_day)

st_write(ao_buff_fb,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/Aerial Ops Tracks/Cleaned/FY21_24_MOFlightTracks_clean_fb_",as.numeric(round(buff_dist)),"m_buff.shp",append=F))

## intersect with farm bill counties
fb$county_area_km <- st_area(fb)/1e6

ao_buff_day_count <- st_intersection(fb %>% st_transform(st_crs(ao_buff)), 
                                  ao_buff_day)

#calculate proportion of EA covered 
intersect_count <- ao_buff_day_count %>% 
  mutate(intersect_area = st_area(.),
         intersect_area_km=intersect_area/1e6) %>%   
  dplyr::select(zone,CNTY_NAME,Date,intersect_area_km) %>%   
  st_drop_geometry() %>% 
  left_join(fb %>% st_drop_geometry()) %>% 
  mutate(prop_county_flown = as.numeric(intersect_area_km/county_area_km))

#add flight time
ao_county_day <- ao_buff_day_count %>% st_drop_geometry() %>% 
  left_join(intersect_count) %>% as.data.frame()

hist(ao_county_day$dt_hr)
hist(ao_county_day$prop_county_flown)
plot(ao_county_day$dt_hr,ao_county_day$prop_county_flown)

ao_county_day<- ao_county_day %>% rename(flight_time_hr=dt_hr)

write.csv(ao_county_day,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Farm Bill/MIS Data/fy21_24_flight_time_fb_county_zone_clean_",as.numeric(round(buff_dist)),"m_buff.csv"))


## intersect mcps with farm bill zones
fb_merge <- fb %>% group_by(zone) %>% 
  summarise(geometry=st_union(geometry),
            fb_area_km=st_area(geometry)/1e6)

ao_buff_day_fb <- st_intersection(fb_merge %>% st_transform(st_crs(ao_buff)), 
                                  ao_buff_day)

#calculate proportion of EA covered 
intersect_pct_fb <- ao_buff_day_fb %>% 
  mutate(intersect_area = st_area(.),
         intersect_area_km=intersect_area/1e6) %>%   
  dplyr::select(zone,Date,intersect_area_km) %>%   
  st_drop_geometry() %>% 
  left_join(fb_merge %>% st_drop_geometry()) %>% 
  mutate(prop_fb_flown = as.numeric(intersect_area_km/fb_area_km))

#add flight time
ao_fb_day <- ao_buff_day_fb %>% st_drop_geometry() %>% 
  left_join(intersect_pct_fb) %>% as.data.frame()

hist(ao_fb_day$dt_hr)
hist(ao_fb_day$prop_fb_flown)
plot(ao_fb_day$dt_hr,ao_fb_day$prop_fb_flown)

ao_fb_day<- ao_fb_day %>% rename(flight_time_hr=dt_hr)

write.csv(ao_fb_day,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Farm Bill/MIS Data/fy21_24_flight_time_fb_zone_clean_",as.numeric(round(buff_dist)),"m_buff.csv"))



#merge with county data ----------------------------
study_site_grid <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/so_mo_counties.shp")
study_site_grid <- study_site_grid %>% rename(county_area_km=area_km)
# ao_buff_day_ssg <- ao_buff_day %>% st_intersection(study_site_grid %>% 
#                                         st_transform(st_crs(ao_buff_day)))
# ao_buff_day_ssg <- ao_buff_day_ssg %>% 
#   dplyr::select(Aircraft,Date,SiteID,Area_Name,
#                 elim_area_idx,site_area_km,geometry)

ao_buff_day_counties <- st_intersection(study_site_grid %>% 
                                          st_transform(st_crs(ao_buff_day)), 
                                        ao_buff_day)

intersect_counties <- ao_buff_day_counties %>% 
  mutate(intersect_area = st_area(.),
         intersect_area_km=intersect_area/1e6) %>%   
  dplyr::select(COUNTY,Date,intersect_area_km) %>%   
  st_drop_geometry() %>% 
  left_join(study_site_grid %>% st_drop_geometry()) %>% 
  mutate(prop_county_flown = as.numeric(intersect_area_km/county_area_km))

ao_counties_day <- ao_buff_day_counties %>% st_drop_geometry() %>% 
  left_join(intersect_counties) %>% as.data.frame()
ao_counties_day %>% filter(dt_hr>100)
ao_buff %>% filter(Aircraft=="N332AB-ER" & Date==as.Date("2021-12-08"))

hist(ao_counties_day$dt_hr)
hist(ao_counties_day$prop_county_flown)
plot(ao_counties_day$dt_hr,ao_counties_day$prop_county_flown)

ao_counties_day<- ao_counties_day %>% rename(flight_time_hr=dt_hr)
write.csv(ao_counties_day,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Model Ready Data/fy21_24_flight_time_counties_clean_",as.numeric(round(buff_dist)),"m_buff.csv"))


#merge with missouri property data ----------------------------
prp_nat <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/National Density/MIS_Pull_01JUL2024/fs_national_property_01JUL2024_clean.shp")

prp_nat <- prp_nat %>% rename(PropID=AGRP_PR) %>% 
  filter(ST_NAME=="MISSOURI")

prp_nat <- prp_nat[!duplicated(prp_nat),]

prp_nat_coords <- st_coordinates(prp_nat)
prp_nat$long <- prp_nat_coords[,1]
prp_nat$lat <- prp_nat_coords[,2]

#multiple parcels in same property
prp_nat %>% st_drop_geometry() %>% 
  group_by(PropID) %>% 
  summarise(n=n(),
            lat=length(unique(lat)),
            long=length(unique(long))) %>% 
  filter(n>1)

prp_nat <- prp_nat %>%
  st_drop_geometry() %>% 
  group_by(PropID) %>% 
  summarise(prop_area_acres = sum(unique(PRPS_QT)),
            lat_prop=unique(lat),
            long_prop=unique(long))

prp_nat_sf <- prp_nat %>% 
  st_as_sf(coords=c("long_prop","lat_prop"),crs="epsg:4326") %>% 
  st_transform(crs='epsg:26915') %>% 
  mutate(prop_area_km=prop_area_acres * 0.00404686,
         geometry=st_buffer(geometry,dist=sqrt(prop_area_km*1e6/pi))) %>% 
  st_as_sf()

#check
# prp_nat_sf %>% mutate(prop_area_buff=st_area(geometry)/1e6)

ao_buff_day_prop <- st_intersection(prp_nat_sf %>% 
                                    st_transform(st_crs(ao_buff_day)), 
                                        ao_buff_day)

intersect_prop <- ao_buff_day_prop %>% 
  mutate(intersect_area = st_area(.),
         intersect_area_km=intersect_area/1e6) %>%   
  dplyr::select(PropID,Date,intersect_area_km) %>%   
  st_drop_geometry() %>% 
  left_join(prp_nat_sf %>% st_drop_geometry()) %>%
  mutate(prop_property_flown = as.numeric(intersect_area_km/prop_area_km))


ao_prop_day <- ao_buff_day_prop %>% st_drop_geometry() %>% 
  left_join(intersect_prop) %>% as.data.frame()
ao_prop_day %>% filter(dt_hr>100)

hist(ao_prop_day$dt_hr)
hist(ao_prop_day$prop_property_flown)
plot(ao_prop_day$dt_hr,ao_prop_day$prop_property_flown)

ao_prop_day<- ao_prop_day %>% rename(flight_time_hr=dt_hr)
write.csv(ao_prop_day,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Model Ready Data/fy21_24_flight_time_properites_clean_",as.numeric(round(buff_dist)),"m_buff.csv"))
