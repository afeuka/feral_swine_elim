### Title: Data prep functions for MO dynamic occupancy model
### Author: Abbey Feuka
### Date: 05DEC2023
### Notes: numbering weeks strictly within month

#grid systematic baiting by watershed
grid_sysbait_take <- function(study_site_grid, #occupancy grid/sites
                              period) {#"month"or "season"
  require(tidyverse)
  require(sf)
  
  #section ------------------------------------
  sec <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/MO_Public_Land_Survey_System/MO_Public_Land_Survey_System.shp",quiet=T)
  sec <- st_transform(sec,st_crs(study_site_grid))
  sec <- st_intersection(sec,study_site_grid)
  
  #clean systematic baiting -------------------------------
  sys <- read.csv("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Systematic Baiting/Ops_Sept2020_Dec2022.csv")
  sys1 <- read.csv("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Systematic Baiting/Ops_Jan2023_Jun2023.csv")
  sys1<- sys1 %>% select(-c(GlobalID..,EditDate,Editor))
  sys <- full_join(sys,sys1) 
  sys2 <- read.csv("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Systematic Baiting/Ops_Jul2023_Nov2023.csv")
  sys2 <- sys2 %>% select(-c(GlobalID..))
  sys <- full_join(sys,sys2 %>% mutate(Lat=as.character(Lat),Long=as.character(Long))) 
  
  ##remove NA geometries -------------------------------
  sys$Lat <- as.numeric(sys$Lat)
  sys$Long <- as.numeric(sys$Long)
  sys <- subset(sys, !is.na(Lat))
  
  ##make date objects ---------------------------------
  sys$Baiting.Start.Date <- as.Date(sys$Baiting.Start.Date,format="%m/%d/%Y %H:%M")
  sys$Date.of.Hot.Bait <- as.Date(sys$Date.of.Hot.Bait,format="%m/%d/%Y %H:%M")
  sys$Complete.Date <- as.Date(sys$Complete.Date,format="%m/%d/%Y %H:%M")
  sys$Trap.Start.Date <- as.Date(sys$Trap.Start.Date,format="%m/%d/%Y %H:%M")
  sys$Site.Visit.Date <- as.Date(sys$Site.Visit.Date,format="%m/%d/%Y %H:%M")
  
  sys$Site.Visit.Reason[sys$Site.Visit.Reason=="<Null>"] <- NA
  sys <- subset(sys,!(Site.Visit.Reason%in%c("Domestic Observed","Misinformation")))
  
  ##use southern mo only ---------------------------------
  sys_sf <- st_as_sf(sys,coords=c("Long","Lat"),crs="epsg:4326")
  sys_sf <- sys_sf[!(sys_sf$OBJECTID..%in%sys$OBJECTID..[which(sys$Lat>39)]),]
  sys_sf <- st_transform(sys_sf,st_crs(study_site_grid))
  sys_sf <- sys_sf %>% st_intersection(sec) 
  sys_sf <- sys_sf %>% select(Site.Status,Baiting.Start.Date,Date.of.Hot.Bait,
                              Trap.Start.Date,Complete.Date,LABEL_,SiteID,geometry) %>% 
    rename(section = LABEL_)

  #filter to study period
  sys_sf <- sys_sf %>% filter((Trap.Start.Date>=as.Date("2020-09-01",format="%Y-%m-%d") |
                                is.na(Trap.Start.Date)))
  
  sys_sf <- sys_sf %>% filter((Baiting.Start.Date>=as.Date("2020-09-01",format="%Y-%m-%d") |
                                 is.na(Baiting.Start.Date)))

  #set up whole site data frame ---------------------
  fy_mo <- c(10,1,4,7)
  max_yr <- max(year(sys_sf$Trap.Start.Date),year(sys_sf$Baiting.Start.Date),year(sys_sf$Date.of.Hot.Bait),na.rm=T)
  det_yr <- 2020:max_yr
  per_start <- as.Date(sapply(1:length(fy_mo),function(i){
      sapply(1:length(det_yr),function(j){
        paste0(det_yr[j],"-",fy_mo[i],"-01")})}))
  per_start <- per_start[order(per_start)]
  per_start <- per_start[-which(per_start<=as.Date("2020-09-01"))]
  subper <- c(sapply(1:length(per_start),function(i){
      if(i<length(per_start)){
        seq(per_start[i],per_start[i+1],by=10)
      } else {
        seq(per_start[i],per_start[i] %m+% months(3),by=10)
      }}))
  subper <- as.Date(subper,origin="1970-01-01")

  subper<- data.frame(subper_start=subper)
  subper$subper_end <- subper$subper_start+days(9)

  subper$year <- year(subper$subper_start)
  subper$month <- month(subper$subper_start)
  subper$fy <- subper$year
  subper$fy[subper$month%in%c(10,11,12)] <- subper$year[subper$month%in%c(10,11,12)]+1

  seasons <- data.frame(month=1:12,
                          season=c("Winter","Winter","Winter",
                                   "Spring","Spring","Spring",
                                   "Summer","Summer","Summer",
                                   "Fall","Fall","Fall"))
  seasons$scode <- as.numeric(as.factor(seasons$season))
  seasons$scode[seasons$season=="Fall"] <- 2
  seasons$scode[seasons$season=="Winter"] <- 3
  seasons$scode[seasons$season=="Spring"] <- 4
  seasons$scode[seasons$season=="Summer"] <- 1 #first time period is in september

  subper$season <- seasons$scode[month(subper$subper_start)]
  subper$period <- paste(subper$season,subper$year)
  subper$period <- as.numeric(as.factor(subper$period))
  subper$period <- as.numeric(factor(subper$period,levels=unique(subper$period)))

  subper$subper_end_new <-NA
  #cut overlapping periods
  for(i in 1:nrow(subper)){
      if(i<nrow(subper)){
        if(subper$subper_end[i]>subper$subper_start[i+1]){
          subper$subper_end_new[i] <- subper$subper_start[i+1]-1
        } else {
          subper$subper_end_new[i] <- subper$subper_end[i]
        }
      } else {
        subper$subper_end_new[i] <- subper$subper_end[i]
      }
  }
  subper <- subper %>%
    mutate(subper_end_new = as.Date(subper_end_new)) %>%
    dplyr::select(-subper_end) %>%
    rename(subper_end=subper_end_new)
  
  subper <- cbind.data.frame(subper,matrix(0,nrow=nrow(subper),ncol=length(unique(study_site_grid$SiteID))))
  subper <- as.data.frame(subper)
  
  sys_sf <- sys_sf %>% 
    filter(Baiting.Start.Date<Trap.Start.Date | 
             is.na(Baiting.Start.Date) | 
             is.na(Trap.Start.Date)) 
  
  sys_df <- sys_sf %>% st_drop_geometry()

  # hog detections occurring three days prior to trapping
  sys_df$Hog.Detection.Sys <- rep(NA,nrow(sys_df))
  for(i in 1:nrow(sys_df)){
    if(!is.na(sys_df$Trap.Start.Date[i]) & is.na(sys_df$Date.of.Hot.Bait[i])){
      sys_df$Hog.Detection.Sys[i] <- sys_df$Trap.Start.Date[i] - 3
      sys_df$Date.of.Hot.Bait[i] <- as.Date(sys_df$Hog.Detection.Sys[i],origin="1970-01-01") 
    }
    if(!is.na(sys_df$Date.of.Hot.Bait[i])){
      sys_df$Hog.Detection.Sys[i] <- sys_df$Date.of.Hot.Bait[i] 
    }
  }
  sys_df$Hog.Detection.Sys <- as.Date(sys_df$Hog.Detection.Sys,origin="1970-01-01")
  
  ##clean dates  ---------------------------------------------
  for(i in 1:nrow(sys_df)){
    #no complete date but baiting start date
    if((!is.na(sys_df$Baiting.Start.Date[i]) & is.na(sys_df$Complete.Date[i]))){
      sys_df$Complete.Date[i] <- sys_df$Baiting.Start.Date[i] + days(90)
    }
    #no baiting start date but complete date
    if(!is.na(sys_df$Complete.Date[i]) & is.na(sys_df$Baiting.Start.Date[i])){
      sys_df$Baiting.Start.Date[i] <- sys_df$Complete.Date[i] - days(90)
    }
    # #no complete date but trap start date
    # if(!is.na(sys_df$Trap.Start.Date[i]) & is.na(sys_df$Complete.Date[i])){
    #   sys_df$Complete.Date[i] <- sys_df$Trap.Start.Date[i] + days(90)
    # }
    if(!is.na(sys_df$Baiting.Start.Date[i]) & !is.na(sys_df$Hog.Detection.Sys[i])){
      if(sys_df$Hog.Detection.Sys[i]<sys_df$Baiting.Start.Date[i]){
        sys_df$Baiting.Start.Date[i] <- sys_df$Hog.Detection.Sys[i]
      }
    }else if(!is.na(sys_df$Trap.Start.Date[i]) & !is.na(sys_df$Hog.Detection.Sys[i])){
      if(sys_df$Hog.Detection.Sys[i]<sys_df$Trap.Start.Date[i]){
        sys_df$Trap.Start.Date[i] <- sys_df$Hog.Detection.Sys[i]
      }
    }
    # #complete date is earlier than baiting start
    # if(!(is.na(sys_df$Baiting.Start.Date[i]) | is.na(sys_df$Complete.Date[i]))){
    #   if(sys_df$Baiting.Start.Date[i] > sys_df$Complete.Date[i]){
    #     sys_df$Complete.Date[i] <- sys_df$Baiting.Start.Date[i] + days(90)
    #   }
    # }
  }

  sys_df <- sys_df %>% filter(Baiting.Start.Date<Complete.Date) 
  
  #correct complete date for end of study
  sys_df$Complete.Date[sys_df$Complete.Date>as.Date("11-13-2023",format="%m-%d-%Y") & 
                         !is.na(sys_df$Complete.Date)] <- 
    as.Date("11-13-2023",format="%m-%d-%Y")

  #remove days with no start or end date
  sys_df <- subset(sys_df,!(is.na(Baiting.Start.Date) & is.na(Complete.Date)))

  #remove duplicates
  sys_df <- sys_df[!duplicated(sys_df),]

  #add hot baits to capture history ------------------
  for(i in 1:nrow(sys_df)){
    subper[which(subper$subper_start<=sys_df$Hog.Detection.Sys[i] & 
                   subper$subper_end>=sys_df$Hog.Detection.Sys[i]),
              as.character(sys_df$SiteID[i])] <- 1
  }
  
  #removal data ----------------------------------
  rem <- read.csv("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Removal/Take/MO_ops_take_cleaned.csv")
  rem$Date <- as.Date(rem$Date)
  rem_sf <- st_as_sf(rem,coords=c("Long","Lat"),crs="epsg:4326")

  ## grid to study_site ------------------------------------------
  rem_sf <- st_transform(rem_sf, st_crs(study_site_grid))
  rem_sf_trap <- st_intersection(rem_sf,study_site_grid) %>%
    st_intersection(sec) %>%
    select(Date,Total,Males,Females,Adults,Sub_adults,Method,NonBreed_M,
           NonBreed_F,Adult_M,Adult_F,CY,SiteID,LABEL_,geometry) %>%
    rename(section = LABEL_) %>% filter(Method=="Trap") %>%
    filter(Date>=as.Date("01-01-2014",format="%m-%d-%Y"))

  #remove duplicates
  #might have different lat/longs but in same section = same trap
  #removes differences in individual hog counts (m/f, adult/juv), uses totals only
  rem_sf_trap <- rem_sf_trap[!duplicated(rem_sf_trap),]
  rem_df_trap <- rem_sf_trap %>% st_drop_geometry()
  
  # align systematic baiting and take data ------------------------------------
  sys_df$detection_trap <- rep(NA,nrow(sys_df))
  sys_df_orig <- sys_df
  for(i in 1:nrow(sys_df_orig)){
    rem_same <- which(rem_df_trap$Date>= sys_df_orig$Baiting.Start.Date[i] &
                        rem_df_trap$Date<= sys_df_orig$Complete.Date[i] &
                        rem_df_trap$section==sys_df_orig$section[i])

    # x[i] <- length(rem_same)
    if(sum(rem_same)>0){
      for(j in 1:length(rem_same)){
        x <- sys_df_orig[i,]
        x$detection_trap <- rem_df_trap$Date[rem_same][j]
        sys_df <- rbind(sys_df,x)
      }
    }
  }
  sys_df$detection_trap <- as.Date(sys_df$detection_trap, origin="1970-01-01")

  sys_df <- sys_df %>% rename(detection_sys=Hog.Detection.Sys) %>%
    pivot_longer(cols=c("detection_sys","detection_trap"),names_to="det_typ",
                 values_to="det_date")
  sys_df <- as.data.frame(sys_df)

  #add trapping to capture history ------------------
  for(i in 1:nrow(sys_df)){
    if(sys_df$det_typ[i]=="detection_trap" & !is.na(sys_df$det_date[i])){
      subper[which(subper$subper_start<=sys_df$det_date[i] &
                     subper$subper_end>=sys_df$det_date[i]),
             as.character(sys_df$SiteID[i])] <- 1
    }
  }
  
  # effort in trap days------------------------------------
  # j<-7
  # i <- which(subper_eff$subper_start>=samp$Baiting.Start.Date[1]&
  #              subper_eff$subper_end<=samp$Complete.Date[1])[1]
  subper_eff <- subper
  subper_eff[,as.character(unique(study_site_grid$SiteID))] <- 0
  # sys_df_trap <- sys_df %>% filter(det_typ=="detection_trap")
  # sys_df_trap <- sys_df_trap[!duplicated(sys_df_trap),]
  for(j in 1:length(unique(study_site_grid$SiteID))){
    samp <- sys_df[which(sys_df$SiteID==j),]
    if(nrow(samp)>0){
      for(i in 1:nrow(subper_eff)){
        idx <- which(samp$Baiting.Start.Date<=subper_eff$subper_start[i] &
                samp$Complete.Date>=subper_eff$subper_end[i])
        if(length(idx)>0){
          temp <- samp[idx,]
          temp$Baiting.Start.Date[temp$Baiting.Start.Date<=subper_eff$subper_start[i]] <- subper_eff$subper_start[i]
          temp$Complete.Date[temp$Complete.Date>=subper_eff$subper_end[i]] <- subper_eff$subper_end[i]
          temp$days <- temp$Complete.Date-temp$Baiting.Start.Date +1 #include first day
          #1 trap per site * number of days/trap = trap days
          subper_eff[i,as.character(j)] <- sum(as.numeric(temp$days[!duplicated(temp$section)]))
        }
      }
    }
  }
  
  subper <- subper %>% pivot_longer(cols=as.character(study_site_grid$SiteID),
                          names_to="site_idx",values_to="detection")
  subper_eff <- subper_eff %>% pivot_longer(cols=as.character(study_site_grid$SiteID),
                                    names_to="site_idx",values_to="trap_nights")
  
  subper <- subper %>% left_join(subper_eff)
  
  #change detections with no trap nights to 0
  subper$detection[subper$detection>subper$trap_nights] <- 0

  list(subper=subper)
}
# subper <- sys$subper

#grid removals by elimination area
grid_removals <- function(study_site_grid,
                          subper, #output from grid_sys_effort, by subperiod (10 days)
                          period#season" or "month"
                          ){
  require(tidyverse)
  require(lubridate)
  require(sf)
  require(terra)
  
  #section ------------------------------------
  sec <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/MO_Public_Land_Survey_System/MO_Public_Land_Survey_System.shp",quiet=T)
  sec <- st_transform(sec,st_crs(study_site_grid))
  sec <- st_intersection(sec,study_site_grid)
  
  #removal data ----------------------------------
  ## collector data
  rem <- read.csv("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Removal/Take/MO_ops_take_cleaned.csv")
  
  ## mis data 
  # rem <- read.csv("./take_mo_all_dates.csv")
  
  rem$Date <- as.Date(rem$Date)
  rem_sf <- st_as_sf(rem,coords=c("Long","Lat"),crs="epsg:4326")
  
  rem_sf <- st_transform(rem_sf, st_crs(study_site_grid))
  rem_site <- st_intersection(rem_sf,study_site_grid) %>% 
    select(Date,SiteID,Total,Males,Females,Adults,Sub_adults,Method,NonBreed_M,
           NonBreed_F,Adult_M,Adult_F,EA,County) %>% 
    filter(year(Date)>=2014)
  
  rem_sf <- rem_sf %>% select(Date,Total,Method,geometry) %>% 
    st_intersection(study_site_grid)
  
  #intersect with section
  rem_site <- rem_site %>% 
    st_intersection(sec %>% 
                      rename(section=LABEL_) %>% 
                      select(section,geometry))
  #remove duplicates
  #might have different lat/longs but in same section = same trap
  #removes differences in individual hog counts (m/f, adult/juv), uses totals only
  rem_site <- rem_site[!duplicated(rem_site),]
  rem_site <- rem_site %>% select(-EA) %>% 
    st_drop_geometry() %>% 
    left_join(study_site_grid %>% 
                st_drop_geometry() %>% 
                select(SiteID,Area_Name,area_km) %>% 
                rename(site_area_km=area_km))
  
  rem_site <- rem_site %>% 
    left_join(study_site_grid %>% 
                group_by(Area_Name) %>% 
                summarise(geometry=st_union(geometry),
                          ea_area_km=as.numeric(st_area(geometry))/1e6) %>% 
                st_drop_geometry()) %>% 
    select(Date,SiteID,Total,Method,section,Area_Name,site_area_km,ea_area_km) %>% 
    filter(Method!="Night Shoot")
  
  period_dates <- subper %>% group_by(period) %>% 
    summarise(per_start=min(subper_start),
              per_end=max(subper_end))
  
  rem_site <- rem_site %>% filter(Date>=min(period_dates$per_start) &
                        Date<=max(period_dates$per_end))
  rem_site$period <- NA
  for(i in 1:nrow(rem_site)){
    rem_site$period[i] <- period_dates$period[period_dates$per_start<=rem_site$Date[i] &
                                                 period_dates$per_end>=rem_site$Date[i]]
  }
  
  list(rem_site=rem_site)
}

# subper <- sys$subper
# rem_site <- rem$rem_site

grid_effort <- function(subper,#output from grid_sysbaittake,by subperiod
                        rem_site, #output from grid_removals
                        study_site_grid,#individual site boundaries
                        grid_typ#"counties" or "watersheds"
                        ){ 
  
  #load effort data ------------------------------
  #data with point locations
  if(grid_typ=="watersheds"){
    load("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/MIS/effort_2014_2023_clean_locs.RData")
  } else {
    eff <- read.csv("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/MIS/MIS Pulls/mo_fs_effort_Jan2014_Dec2023_q5.csv")
    eff <- eff %>% rename(lat =PRPU_N_LAT,lon=PRPU_E_LONG)
    
    eff$lon <- str_trim(eff$lon)
    eff$lat <- str_trim(eff$lat)
    
    eff$lon[grep("--",eff$lon)] <-
      gsub("--","-",eff$lon[grep("--",eff$lon)])
    eff$lon[grep(" ",eff$lon)] <-
      gsub(" ","",eff$lon[grep(" ",eff$lon)])
    eff$lon[!grepl("-",eff$lon)] <-
      paste0("-",eff$lon[!grepl("-",eff$lon)])
    
    eff$lat[grep(" ",eff$lat)] <-
      gsub(" ","",eff$lat[grep(" ",eff$lat)])
    
    
    samp <- eff$lon[which(is.na(unlist(lapply(eff$lon,as.numeric))))]
    samp.idx <- which(is.na(unlist(lapply(eff$lon,as.numeric))))
    eff$lon[samp.idx] <-
      sapply(1:length(samp),function(i){
        paste0(str_split_fixed(samp[i],pattern='\\.',3)[1],
               ".",str_split_fixed(samp[i],pattern='\\.',3)[2])})
    
    samp <- eff$lat[which(is.na(unlist(lapply(eff$lat,as.numeric))))]
    samp.idx <- which(is.na(unlist(lapply(eff$lat,as.numeric))))
    eff$lat[samp.idx] <-
      sapply(1:length(samp),function(i){
        paste0(str_split_fixed(samp[i],pattern='\\.',3)[1],
               ".",str_split_fixed(samp[i],pattern='\\.',3)[2])})
    
    eff$lon <- as.numeric(eff$lon)
    eff$lat <- as.numeric(eff$lat)
    
    eff_sf <- st_as_sf(eff %>% filter(!is.na(lat)),coords=c("lon","lat"),crs="epsg:4326")
    
    #southern missouri watersheds
    smws <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/huc10_siteIDs_cond_EA.shp")
    
    eff_sf <- eff_sf %>% st_transform(st_crs(smws)) %>% st_intersection(smws)
    save(eff_sf,file = "C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/MIS/effort_2014_2023_clean_locs.RData")
  }
  
  eff <- eff_sf
  eff <- eff %>% rename(date=WT_WORK_DATE,method=CMP_NAME)
  
  ## only use MIS effort for traps, ground shooting, systematic baiting
  eff$method[eff$method %in% c("TRAPS, LIVE, FERAL HOGS",
                               "TRAPS, BODY GRIP","TRAPS, CORRAL",
                               "TRAPS, CAGE")] <- "Trap"
  eff$method[eff$method %in% c("FIREARMS")] <- "Ground"
  eff$method[eff$method %in% c("HELICOPTER")] <- "Aerial"
  # eff$method[eff$method %in% c("MONITORING CAMERA")] <- "Camera"
  eff$method[eff$method %in% c("PREBAIT","BAIT STATION")] <- "SysBait"
  eff <- eff %>% filter(method %in% c("Trap","Ground","Aerial",
                                      "Camera",
                                      "SysBait"))
  eff$date <- as.Date(eff$date,format="%d-%b-%y")
  
  eff <- eff %>% filter(UOM_NAME%in%c("MINUTES","HOURS","HOBBS METER"))
  
  eff$time_hr <-  eff$WTM_QTY
  eff$time_hr[eff$UOM_NAME=="MINUTES"] <- eff$WTM_QTY[eff$UOM_NAME=="MINUTES"]/60

  #only hobbs meter readings for aerial ops
  eff <- eff %>% filter(!(method=="Aerial" & UOM_NAME!="HOBBS METER"))
  
  eff$UOM_NAME[eff$UOM_NAME=="HOBBS METER"] <- "HOURS"
  
  eff<- eff %>% filter(USET_NAME%in%c("APPLIED/USED","CHECKED","RESET","SET"))
  
  eff <- eff[!duplicated(eff),]

  #systematic baiting effort hours ---------------
  subper$eff_hrs <- rep(0,nrow(subper)) 
  for(i in 1:nrow(subper)){
    x <- eff[eff$SiteID==subper$site_idx[i] &
               eff$date>=subper$subper_start[i] &
               eff$date<=subper$subper_end[i],]
    
    if(nrow(x)>0){
      if(any(x$method%in%c("SysBait","Trap"))){ #systematic baiting effort 
        subper$eff_hrs[i] <- sum(subset(x,method=="SysBait" | method=="Trap")$time_hr,na.rm=T)
      }
    }
  }

  ##add NAs for missing hours ------------------------
  subper$eff_hrs[subper$eff_hrs==0 & subper$trap_nights>0] <- NA
  
  ##standardize effort by area --------------------
  subper <- subper %>% 
    left_join(study_site_grid %>% mutate(site_idx=as.character(SiteID)))
  subper$eff_hrs_km <- subper$eff_hrs/subper$area_km
  
  ##add proportion of area covered by feral swine distribution ----------------------
  load("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Model Ready Data/NFSP Watershed Overlap/ssg_nfsp_all.RData")
  
  subper <- subper %>% 
    left_join(ssg_nfsp %>% select(SiteID,prp_nfs,year) %>% 
                rename(fy=year) %>% 
                mutate(site_idx=as.character(as.numeric(SiteID))) %>% 
                dplyr::select(-SiteID))
  
  subper$trap_nights_km <- subper$trap_nights/subper$area_km
  
  # removal effort -------------------------------------------------
  ##elimination area - areas ---------------
  ea <- study_site_grid %>% 
    rename(elim_area=Area_Name) %>% 
    group_by(elim_area) %>% 
    summarise(geometry=st_union(geometry),
              area_km=as.numeric(st_area(geometry))/1e6)
  ea$elim_area[ea$elim_area=="0"] <- NA

  rem_day_ea <-
    rem_site %>%  
    group_by(Date=floor_date(Date,"day"),period, Method,Area_Name) %>% 
    summarise(tot_rem=sum(Total),
              ea_area_km=unique(ea_area_km))
  
  ##clip removals to systematic baiting -----------------------
  rem_day_ea <- rem_day_ea %>% 
    filter(Date>=min(subper$subper_start) &
             Date<=max(subper$subper_end))
  
  ##clip systematic baiting to removals --------------------
  subper <- subper %>% 
    filter(subper_start>=min(rem_day_ea$Date) &
             subper_end<=max(rem_day_ea$Date))%>% 
    st_drop_geometry() %>% 
    select(-c(season,month,n)) %>%
    rename(site_area_km=area_km)

  ##effective area and time ------------------------
  # rem_site_ea <- rem_site %>%
  #   left_join(ao %>%
  #               mutate(Date=as.Date(day,format="%Y-%m-%d"),
  #                      Area_Name=as.character(Area_Name)) %>%
  #               select(Area_Name,Date,flight_time,prop_ea_flown,
  #                      intersect_area,Method)) %>%
  #   rename(eff_hrs=flight_time,
  #          prop_ea_impact= prop_ea_flown) %>%
  #   mutate(effect_area_km=intersect_area/1e6)
  
  eff_day <- eff %>% 
    st_drop_geometry() %>% 
    group_by(Date=floor_date(date,"day"),
             WT_AGRPROP_ID,method) %>% 
    summarise(tot_hrs=sum(time_hr),#/length(unique(WT_EMP_ID)),
              num_events=n(),
              Area_Name=unique(Area_Name),
              SiteID=unique(SiteID)) %>% 
    left_join(study_site_grid %>% st_drop_geometry() %>% 
                rename(site_area_km=area_km)
              %>% select(-n)) %>% 
    left_join(ea %>% st_drop_geometry() %>% 
                          rename(Area_Name=elim_area,
                                 ea_area_km=area_km))
  

  period_dates <- subper %>% group_by(period) %>% 
    summarise(per_start=min(subper_start),
              per_end=max(subper_end))
  
  eff_day <- eff_day %>% filter(Date>=min(period_dates$per_start) &
                                   Date<=max(period_dates$per_end)) %>% 
    filter(!is.na(Date))
  
  eff_day$period <- NA
  for(i in 1:nrow(eff_day)){
    eff_day$period[i] <- period_dates$period[period_dates$per_start<=eff_day$Date[i] &
                                                period_dates$per_end>=eff_day$Date[i]]
  }
  
  ### join aerial gps data -----------------------------
  ao <- read.csv("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Model Ready Data/fy21_24_flight_time_ea_clean.csv")
  ao$method<-"Aerial"
  
  eff_day <- eff_day %>%
    left_join(ao %>% rename(Date=day,
                            prop_ea_impact=prop_ea_flown) %>% 
                mutate(Date=as.Date(Date),
                       Area_Name=as.character(Area_Name),
                       effect_area_km=as.numeric(intersect_area)/1e6) %>% 
                select(Area_Name,Date,method,flight_time,effect_area_km,prop_ea_impact))

  #remove MIS flight times and replace with GPS tracks--------------------
  eff_day$tot_hrs[eff_day$method=="Aerial"] <- NA
  eff_day$tot_hrs[eff_day$method=="Aerial"] <- eff_day$flight_time[eff_day$method=="Aerial"]
  
  # OR assume helicopter covered entire watershed -------------------
  # eff_day <- eff_day %>% mutate(effect_area_km=site_area_km)
  
  ## ground = 5 km2 ------------
  eff_day$effect_area_km[eff_day$method=="Ground"] <- 5
  
  ## trap = 6.7 km2 (McCrae et al 2020)----------------------------
  eff_day$effect_area_km[eff_day$method=="Trap"] <- 6.7

  # remove misisng property areas
  # eff_day <- eff_day %>% filter(!is.na(effect_area_km))
  
  #remove systematic baiting ---------------
  eff_day <- eff_day %>% filter(method!="SysBait")
  
  #multiply by number of events
  eff_day$effect_area_km[eff_day$method=="Ground"] <- eff_day$num_events[eff_day$method=="Ground"]*eff_day$effect_area_km[eff_day$method=="Ground"]
  eff_day$effect_area_km[eff_day$method=="Trap"] <- eff_day$num_events[eff_day$method=="Trap"]*eff_day$effect_area_km[eff_day$method=="Trap"]
  
  # summarise by elim area-------------------
  eff_day_ea <- eff_day %>% 
    group_by(Date,period,method,Area_Name) %>% 
    summarise(tot_hrs=sum(tot_hrs),
              num_events=sum(num_events),
              effect_area_km=max(effect_area_km),
              ea_area_km=unique(ea_area_km)) 
    
  # calculate proportional area of impact on elimination area
  eff_day_ea$prop_ea_impact <- eff_day_ea$effect_area_km/eff_day_ea$ea_area_km

  eff_day_ea$eff_area_hrs <- eff_day_ea$tot_hrs/eff_day_ea$effect_area_km
  eff_day_ea$eff_area_events <- eff_day_ea$num_events/eff_day_ea$effect_area_km

  rem_day_ea <- rem_day_ea %>% rename(method=Method)
  rem_day_ea$method[rem_day_ea$method=="Ground Shoot"] <- "Ground"

  rem_eff_ea <- eff_day_ea %>% 
    filter(Date>=min(rem_day_ea$Date) & 
             Date<=max(rem_day_ea$Date)) %>% 
    left_join(rem_day_ea) 
  rem_eff_ea$tot_rem[is.na(rem_eff_ea$tot_rem)] <- 0
  
  rem_eff_ea <- rem_eff_ea %>% group_by(period,method,Area_Name) %>% 
    mutate(pass_idx=1:n()) 

   #grid landscape covariates to watershed -----------------------------
  if(file.exists(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Landscape Covariates/nlcd_",grid_typ,".RData"))){
    load(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Landscape Covariates/nlcd_",grid_typ,".RData"))
  } else {
    library(terra)
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
    
    st_write(water,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/so_mo_water_nlcd2019_",grid_typ,".shp"),append=F)
    
    nlcd.developed <- nlcd.crop
    values(nlcd.developed)[values(nlcd.developed)%in%c(21:24)] <- 1
    values(nlcd.developed)[values(nlcd.developed)!=1] <- 0
    
    developed <- terra::extract(nlcd.developed,
                                study_site_grid,fun=mean) %>% 
      rename(SiteID=ID,developed=Layer_1) %>% 
      right_join(study_site_grid,by="SiteID") %>% 
      st_as_sf()# %>% select(-c(neighbors)) 
    # select(-c(neighbors,FID,SITE,geometry))
    st_write(developed,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/so_mo_developed_nlcd2019_",grid_typ,".shp"),append=F)
    
    nlcd.forest <- nlcd.crop
    values(nlcd.forest)[values(nlcd.forest)%in%c(41:43)] <- 1
    values(nlcd.forest)[values(nlcd.forest)!=1] <- 0
    
    forest <- terra::extract(nlcd.forest,study_site_grid,fun=mean) %>%
      rename(SiteID=ID,forest=Layer_1) %>%
      right_join(study_site_grid,by="SiteID") %>%
      st_as_sf() #%>% select(-c(neighbors)) 
    st_write(forest,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/so_mo_forest_nlcd2019_",grid_typ,".shp"),append=F)
    
    nlcd.wetland <- nlcd.crop
    values(nlcd.wetland)[values(nlcd.wetland)%in%c(90,95)] <- 1
    values(nlcd.wetland)[values(nlcd.wetland)!=1] <- 0
    
    wetland <- terra::extract(nlcd.wetland,study_site_grid,fun=mean) %>% 
      rename(SiteID=ID,wetland=Layer_1) %>% 
      right_join(study_site_grid,by="SiteID") %>% 
      # select(-c(neighbors,FID,SITE,geometry))
      st_as_sf() #%>% select(-c(neighbors)) 
    st_write(wetland,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/so_mo_wetland_nlcd2019_",grid_typ,".shp"),append=F)
    
    nlcd.agri <- nlcd.crop
    values(nlcd.agri)[values(nlcd.agri)%in%c(81,82)] <- 1
    values(nlcd.agri)[values(nlcd.agri)!=1] <- 0
    
    agri <- terra::extract(nlcd.agri,study_site_grid,fun=mean) %>% 
      rename(SiteID=ID,agri=Layer_1) %>% 
      right_join(study_site_grid,by="SiteID") %>% 
      st_as_sf() #%>% select(-c(neighbors)) 
    st_write(agri,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/so_mo_agri_nlcd2019_",grid_typ,".shp"),append=F)
    
    nlcd_siteid <- full_join(water %>% st_drop_geometry(),
                             developed%>% st_drop_geometry()) %>% 
      full_join(agri%>% st_drop_geometry()) %>% 
      full_join(wetland%>% st_drop_geometry()) %>% 
      full_join(forest%>% st_drop_geometry()) %>%
      full_join(study_site_grid) %>% 
      st_drop_geometry()#%>% dplyr::select(-neighbors)

    save(nlcd_siteid,file=paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Landscape Covariates/nlcd_",grid_typ,".RData"))
  }
  
  list(nlcd_siteid=nlcd_siteid,
       sysbait_det_eff=subper,
       rem_eff_ea=rem_eff_ea)
}