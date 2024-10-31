### Title: Data prep functions for MO dynamic occupancy model
### Author: Abbey Feuka
### Date: 05DEC2023
### Notes: numbering weeks strictly within month

#grid systematic baiting by watershed
grid_sysbait_take <- function(study_site_grid, #occupancy grid/sites
                              start_date="2020-09-01", #start date for data in y-m-d string
                              end_date="2023-09-30", #end date for data in y-m-d string
                              period) {#"month"or "season"
  require(tidyverse)
  require(sf)
  require(lubridate)
  
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
  sys <- sys %>% select(Site.Status, Baiting.Start.Date, Date.of.Hot.Bait, Site.Visit.Date, Site.Visit.Reason, Trap.Start.Date,
                        Complete.Date, Lat, Long, CreationDate, Creator, Sounder.Size)
  sys3 <- read.csv("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Systematic Baiting/Ops_July23_Sept24.csv")
  sys3 <- sys3 %>% rename(Site.Status=SiteStatus,
                          Baiting.Start.Date=BaitingStartDate,
                          Date.of.Hot.Bait=DateHotBait,
                          Site.Visit.Date=SiteVisitDate,
                          Site.Visit.Reason=SiteVisitReas,
                          Trap.Start.Date=TrapStartDate,
                          Complete.Date=CompleteDate,
                          Responsible.Trapper=ResposTrapper,
                          CreationDate=CreationDate,
                          Lat=lat,
                          Long=long,
                          Sounder.Size=SounderSize) %>% 
    mutate(Lat=as.character(Lat),
           Long=as.character(Long),
           Sounder.Size=as.character(Sounder.Size)) %>% 
    select(Site.Status, Baiting.Start.Date, Date.of.Hot.Bait, Site.Visit.Date, Site.Visit.Reason, Trap.Start.Date,
           Complete.Date, Lat, Long, CreationDate, Creator, Sounder.Size)
  
  ##make date objects ---------------------------------
  sys$Baiting.Start.Date <- as.Date(sys$Baiting.Start.Date,format="%m/%d/%Y %H:%M")
  sys$Date.of.Hot.Bait <- as.Date(sys$Date.of.Hot.Bait,format="%m/%d/%Y %H:%M")
  sys$Complete.Date <- as.Date(sys$Complete.Date,format="%m/%d/%Y %H:%M")
  sys$Trap.Start.Date <- as.Date(sys$Trap.Start.Date,format="%m/%d/%Y %H:%M")
  sys$Site.Visit.Date <- as.Date(sys$Site.Visit.Date,format="%m/%d/%Y %H:%M")
  sys$CreationDate <- as.Date(sys$CreationDate,format="%m/%d/%Y %H:%M:%S")
  
  sys3$Baiting.Start.Date <-as.Date(sys3$Baiting.Start.Date,format="%Y-%m-%d %H:%M:%S")
  sys3$Date.of.Hot.Bait <- as.Date(sys3$Date.of.Hot.Bait,format="%Y-%m-%d %H:%M:%S")
  sys3$Complete.Date <- as.Date(sys3$Complete.Date,format="%Y-%m-%d %H:%M:%S")
  sys3$Trap.Start.Date <- as.Date(sys3$Trap.Start.Date,format="%Y-%m-%d %H:%M:%S")
  sys3$Site.Visit.Date <- as.Date(sys3$Site.Visit.Date,format="%Y-%m-%d %H:%M:%S")
  sys3$CreationDate <- as.Date(sys3$CreationDate,format="%Y-%m-%d %H:%M:%S")
  
  # sys3 %>% filter(Baiting.Start.Date>=as.Date("2023-07-01")) %>% arrange(Baiting.Start.Date)
  # sys %>% filter(Baiting.Start.Date>=as.Date("2023-07-01")) %>% arrange(Baiting.Start.Date)
  
  sys <- full_join(sys,sys3)
  
  ##remove NA geometries -------------------------------
  sys$Lat <- as.numeric(sys$Lat)
  sys$Long <- as.numeric(sys$Long)
  sys <- subset(sys, !is.na(Lat))
  
  sys$Site.Visit.Reason[sys$Site.Visit.Reason=="<Null>"] <- NA
  sys <- subset(sys,!(Site.Visit.Reason%in%c("Domestic Observed","Misinformation")))
  
  ##use southern mo only ---------------------------------
  sys <- sys[which(sys$Lat<39),]
  sys <- sys[!duplicated(sys),]
  
  sys_sf <- st_as_sf(sys,coords=c("Long","Lat"),crs="epsg:4326")
  sys_sf <- st_transform(sys_sf,st_crs(study_site_grid)) %>% 
    st_intersection(sec) %>% 
    select(Site.Status,Baiting.Start.Date,Date.of.Hot.Bait,
           Trap.Start.Date,Complete.Date,LABEL_,SiteID,geometry) %>% 
    rename(section = LABEL_)
  
  #filter to study period
  start_date <- as.Date(start_date,format="%Y-%m-%d")
  sys_sf <- sys_sf %>% filter((Trap.Start.Date>=start_date |is.na(Trap.Start.Date)))
  sys_sf <- sys_sf %>% filter((Baiting.Start.Date>=start_date|is.na(Baiting.Start.Date)))
  
  #set up whole site data frame ---------------------
  if(month(start_date)%in%c(10,11,12)){
    fy_mo <- c(10,1,4,7)
  } else if(month(start_date)%in%c(1,2,3)){
    fy_mo <- c(1,4,7,10)
  }else if(month(start_date)%in%c(4,5,6)){
    fy_mo <- c(4,7,10,1)
  }else if(month(start_date)%in%c(7,8,9)){
    fy_mo <- c(7,10,1,4)
  }
  
  max_yr <- max(year(sys_sf$Trap.Start.Date),year(sys_sf$Baiting.Start.Date),year(sys_sf$Date.of.Hot.Bait),na.rm=T)
  det_yr <- year(start_date):max_yr
  per_start <- as.Date(sapply(1:length(fy_mo),function(i){
    sapply(1:length(det_yr),function(j){
      paste0(det_yr[j],"-",fy_mo[i],"-01")})}))
  per_start <- per_start[order(per_start)]
  per_start <- per_start[-which(per_start<=as.Date(start_date))]
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
  
  if(month(start_date)%in%c(10,11,12)){
    seasons$scode[seasons$season=="Fall"] <- 1
    seasons$scode[seasons$season=="Winter"] <- 2
    seasons$scode[seasons$season=="Spring"] <- 3
    seasons$scode[seasons$season=="Summer"] <- 4
  } else if(month(start_date)%in%c(1,2,3)){
    seasons$scode[seasons$season=="Fall"] <- 4
    seasons$scode[seasons$season=="Winter"] <- 1
    seasons$scode[seasons$season=="Spring"] <- 2
    seasons$scode[seasons$season=="Summer"] <- 3
  }else if(month(start_date)%in%c(4,5,6)){
    seasons$scode[seasons$season=="Fall"] <- 3
    seasons$scode[seasons$season=="Winter"] <- 4
    seasons$scode[seasons$season=="Spring"] <- 1
    seasons$scode[seasons$season=="Summer"] <- 2
  }else if(month(start_date)%in%c(7,8,9)){
    seasons$scode[seasons$season=="Fall"] <- 2
    seasons$scode[seasons$season=="Winter"] <- 3
    seasons$scode[seasons$season=="Spring"] <- 4
    seasons$scode[seasons$season=="Summer"] <- 1
  }
  
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
  end_date <- as.Date(end_date,format="%Y-%m-%d")
  sys_df$Complete.Date[sys_df$Complete.Date>end_date & 
                         !is.na(sys_df$Complete.Date)] <- end_date
  
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
  rem$Date <- as.POSIXct(rem$Date,format="%Y-%m-%d %H:%M:%S",tz="CST6CDT")
  rem_sf <- st_as_sf(rem,coords=c("Long","Lat"),crs="epsg:4326")
  
  ## grid to study_site ------------------------------------------
  rem_sf <- st_transform(rem_sf, st_crs(study_site_grid))
  rem_sf_trap <- st_intersection(rem_sf,study_site_grid) %>%
    st_intersection(sec) %>%
    select(Date,Total,Males,Females,Adults,Sub_adults,Method,NonBreed_M,
           NonBreed_F,Adult_M,Adult_F,SiteID,LABEL_,geometry) %>%
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
  sys_df$detection_trap <- as.POSIXct(sys_df$detection_trap, tz="CST6CDT", origin="1970-01-01")
  
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
  
  #filter to study period 
  subper <- subper %>% filter(subper_start<=as.Date(end_date))
  
  list(sysbait_det_eff=subper)
  
}
# sysbait_det_eff <- sys$sysbait_det_eff

#grid removals by elimination area
grid_removals <- function(study_site_grid,
                          sysbait_det_eff, #output from grid_sys_effort, by subperiod (10 days)
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
  
  rem$Date <- as.POSIXct(rem$Date,tz="CST6CDT",format="%Y-%m-%d %H:%M:%S")
  rem_sf <- st_as_sf(rem,coords=c("Long","Lat"),crs="epsg:4326")
  
  rem_sf <- st_transform(rem_sf, st_crs(study_site_grid))
  rem_site_sf <- st_intersection(rem_sf,study_site_grid) %>% 
    select(Date,SiteID,Total,Males,Females,Adults,Sub_adults,Method,NonBreed_M,
           NonBreed_F,Adult_M,Adult_F,County) %>% 
    filter(year(Date)>=2020)
  
  rem_sf <- rem_sf %>% select(Date,Total,Method,geometry) %>% 
    st_intersection(study_site_grid)
  
  #intersect with section
  rem_site_sf <- rem_site_sf %>% 
    st_intersection(sec %>% 
                      rename(section=LABEL_) %>% 
                      select(section,geometry))
  #remove duplicates
  #might have different lat/longs but in same section = same trap
  #removes differences in individual hog counts (m/f, adult/juv), uses totals only
  rem_site_sf <- rem_site_sf[!duplicated(rem_site_sf),]
  rem_site_sf <- rem_site_sf %>%
    # st_drop_geometry() %>%
    left_join(study_site_grid %>% 
                st_drop_geometry() %>% 
                select(SiteID,Area_Name,area_km) %>% 
                rename(site_area_km=area_km))
  
  rem_site_sf <- rem_site_sf %>% 
    left_join(study_site_grid %>% 
                group_by(Area_Name) %>% 
                summarise(geometry=st_union(geometry),
                          ea_area_km=as.numeric(st_area(geometry))/1e6) %>% 
                st_drop_geometry()) %>% 
    select(Date,SiteID,Total,Method,section,Area_Name,site_area_km,ea_area_km) %>% 
    filter(Method!="Night Shoot")
  
  period_dates <- sysbait_det_eff %>% group_by(period) %>% 
    summarise(per_start=min(subper_start),
              per_end=max(subper_end))
  
  rem_site_sf <- rem_site_sf %>% filter(Date>=min(period_dates$per_start) &
                                          Date<=max(period_dates$per_end))
  rem_site_sf$period <- NA
  for(i in 1:nrow(rem_site_sf)){
    rem_site_sf$period[i] <- period_dates$period[period_dates$per_start<=as.Date(rem_site_sf$Date[i]) &
                                                   period_dates$per_end>=as.Date(rem_site_sf$Date[i])]
  }
  
  list(rem_site_sf=rem_site_sf)
}

# sysbait_det_eff <- sys$sysbait_det_eff
# rem_site_sf <- rem$rem_site_sf

grid_effort <- function(sysbait_det_eff,#output from grid_sysbaittake,by subperiod
                        rem_site_sf, #output from grid_removals
                        study_site_grid,#individual site boundaries
                        occ_scale,#"counties" or "watersheds"
                        abund_scale#Watersheds" or "elim_areas"
){ 
  
  #load effort data ------------------------------
  #data with point locations
  if(file.exists("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/MIS/effort_2020_2024_clean_locs.RData")){
    load("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/MIS/effort_2020_2024_clean_locs.RData")
  } else {
    eff <- read.csv("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/MIS/MIS Pulls/mo_fs_effort_Sep2020_Sep2024_q6.csv")
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
    save(eff_sf,file = "C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/MIS/effort_2020_2024_clean_locs.RData")
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
  eff <- eff %>% filter(method %in% c("Trap","Ground","Aerial","Camera","SysBait"))
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
  sysbait_det_eff$eff_hrs <- rep(0,nrow(sysbait_det_eff)) 
  for(i in 1:nrow(sysbait_det_eff)){
    x <- eff[eff$SiteID==sysbait_det_eff$site_idx[i] &
               eff$date>=sysbait_det_eff$subper_start[i] &
               eff$date<=sysbait_det_eff$subper_end[i],]
    
    if(nrow(x)>0){
      if(any(x$method%in%c("SysBait","Trap"))){ #systematic baiting effort 
        sysbait_det_eff$eff_hrs[i] <- sum(subset(x,method=="SysBait" | method=="Trap")$time_hr,na.rm=T)
      }
    }
  }
  
  ##add NAs for missing hours ------------------------
  sysbait_det_eff$eff_hrs[sysbait_det_eff$eff_hrs==0 & sysbait_det_eff$trap_nights>0] <- NA
  
  ##standardize effort by area --------------------
  sysbait_det_eff <- sysbait_det_eff %>% 
    left_join(study_site_grid %>% mutate(site_idx=as.character(SiteID)))
  sysbait_det_eff$eff_hrs_km <- sysbait_det_eff$eff_hrs/sysbait_det_eff$area_km
  
  ##add proportion of area covered by feral swine distribution ----------------------
  ## LAGGED a year
  load("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Model Ready Data/NFSP Watershed Overlap/ssg_nfsp_all.RData")
  
  ssg_nfsp <- rbind.data.frame(ssg_nfsp,
                   data.frame(SiteID=1:length(unique(ssg_nfsp$SiteID)),
                              full_ar=NA,intrsc_=NA,prp_nfs=NA,year=2024))
  
  ssg_nfsp$prp_nfs_lag <- NA
  for(i in 1:nrow(ssg_nfsp)){
    if(ssg_nfsp$year[i]>2015){
      ssg_nfsp$prp_nfs_lag[i] <- ssg_nfsp$prp_nfs[ssg_nfsp$SiteID==ssg_nfsp$SiteID[i] &
                                          ssg_nfsp$year==(ssg_nfsp$year[i]-1)]
    }
  }

  sysbait_det_eff <- sysbait_det_eff %>% 
    left_join(ssg_nfsp %>% select(SiteID,prp_nfs,prp_nfs_lag,year) %>% 
                rename(fy=year) %>% 
                mutate(site_idx=as.character(as.numeric(SiteID))) %>% 
                dplyr::select(-SiteID))
  
  sysbait_det_eff$trap_nights_km <- sysbait_det_eff$trap_nights/sysbait_det_eff$area_km
  
  # removal effort -------------------------------------------------
  ##elimination area - areas ---------------
  ea <- study_site_grid %>% 
    group_by(Area_Name) %>% 
    summarise(geometry=st_union(geometry),
              area_km=as.numeric(st_area(geometry))/1e6)
  
  # rem_day_site <-
  #   rem_site_sf %>%  
  #   group_by(Date=floor_date(Date,"day"),period,Method,Area_Name) %>% 
  #   summarise(tot_rem=sum(Total),
  #             n_events=n(),
  #             ea_area_km=unique(ea_area_km))
  
  ##clip removals to systematic baiting -----------------------
  # rem_day_site <- rem_day_site %>% 
  #   filter(Date>=min(sysbait_det_eff$subper_start) &
  #            Date<=max(sysbait_det_eff$subper_end))
  
  ##effective area and time ------------------------
  # rem_site_ea <- rem_site_sf %>%
  #   left_join(ao %>%
  #               mutate(Date=as.Date(day,format="%Y-%m-%d"),
  #                      Area_Name=as.character(Area_Name)) %>%
  #               select(Area_Name,Date,flight_time_hr,prop_ea_flown,
  #                      intersect_area,Method)) %>%
  #   rename(eff_hrs=flight_time_hr,
  #          prop_ea_impact= prop_ea_flown) %>%
  #   mutate(effect_area_km=intersect_area/1e6)
  
  eff_day <- eff %>% 
    st_drop_geometry() %>% 
    group_by(Date=floor_date(date,"day"),
             WT_AGRPROP_ID,method) %>% 
    summarise(tot_hrs=sum(time_hr),#/length(unique(WT_EMP_ID)),
              num_events=n(),
              Area_Name=unique(Area_Nm),
              SiteID=unique(SiteID)) %>% 
    left_join(study_site_grid %>% st_drop_geometry() %>% 
                rename(site_area_km=area_km)%>% 
                select(-n)) %>% 
    left_join(ea %>% st_drop_geometry() %>% 
                rename(ea_area_km=area_km))
  
  period_dates <- sysbait_det_eff %>% group_by(period) %>% 
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
  if(abund_scale=="watersheds"){
    ao <- read.csv("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Model Ready Data/fy21_24_flight_time_ssg_clean_1000ft_buff.csv")
    ao$method<-"Aerial"
    
    ao <- ao %>% 
      rename(prop_site_impact=prop_site_flown) %>% 
      mutate(Date=as.Date(Date),
             Area_Name=as.character(Area_Name),
             effect_area_km=as.numeric(intersect_area_km)) %>% 
      select(SiteID,Area_Name,Date,method,flight_time_hr,effect_area_km,prop_site_impact) %>% 
      filter(Date<=as.Date(max(sysbait_det_eff$subper_end)))
  } else {
    ao <- read.csv("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Model Ready Data/fy21_24_flight_time_ea_clean_1000ft_buff.csv")
    ao$method<-"Aerial"
    
    ao <- ao %>% 
      rename(prop_site_impact=prop_ea_flown) %>% 
      mutate(Date=as.Date(Date),
             Area_Name=as.character(Area_Name),
             effect_area_km=as.numeric(intersect_area_km)) %>% 
      select(Area_Name,Date,method,flight_time_hr,effect_area_km,prop_site_impact) %>% 
      filter(Date<=as.Date(max(sysbait_det_eff$subper_end)))
  }
  hist(ao$flight_time_hr)
  ao$period <- NA
  for(i in 1:nrow(ao)){
    ao$period[i] <- period_dates$period[period_dates$per_start<=ao$Date[i] &
                                          period_dates$per_end>=ao$Date[i]]
  }
  
  # eff_day <- eff_day %>%
  #   left_join(ao %>% rename(Date=Date,
  #                           prop_ea_impact=prop_ea_flown) %>% 
  #               mutate(Date=as.Date(Date),
  #                      Area_Name=as.character(Area_Name),
  #                      effect_area_km=as.numeric(intersect_area)/1e6) %>% 
  #               select(Area_Name,Date,method,flight_time_hr,effect_area_km,prop_ea_impact))
  
  #remove MIS flight times and replace with GPS tracks--------------------
  # eff_day$tot_hrs[eff_day$method=="Aerial"] <- NA
  # eff_day$tot_hrs[eff_day$method=="Aerial"] <- eff_day$flight_time_hr[eff_day$method=="Aerial"]
  
  #calculate area buffers for ground and trap --------------
  rem_site_trap <- rem_site_sf %>% filter(Method=="Trap") %>% 
    mutate(geometry=st_buffer(geometry,dist=sqrt(6.7e6/3.14)))
  rem_site_ground <- rem_site_sf %>% filter(Method=="Ground Shoot") %>% 
    mutate(geometry=st_buffer(geometry,dist=sqrt(5e6/3.14)))
  rem_site_aerial <- rem_site_sf %>% filter(Method=="Aerial")
  rem_site_buff <- rbind(rem_site_aerial,rem_site_trap,rem_site_ground)
  
  # ggplot()+
  #   geom_sf(data=study_site_grid,fill="transparent")+
  #   geom_sf(data=ea,aes(col=Area_Name),lwd=2)+
  #   geom_sf(data=rem_site_buff %>% 
  #             group_by(Method,Date=floor_date(as.Date(Date),"day")) %>% 
  #             summarise(geometry=st_union(geometry)) %>% 
  #             filter(Date>=as.Date("2024-05-01") & 
  #                      Date<=as.Date("2024-05-30")),aes(fill=Method),alpha=0.5)
  # 
  # ggplot()+
  #   geom_sf(data=study_site_grid,fill="transparent")+
  #   geom_sf(data=ea,aes(col=Area_Name),lwd=2)+
  #   geom_sf(data=sys_sf %>% filter(Trap.Start.Date>=as.Date("2024-05-01") & 
  #                                    Trap.Start.Date<=as.Date("2024-05-30") &
  #                                    !is.na(Trap.Start.Date)))
  
  #merge all trap/ground footprints to calculate effective area --------------
  if(abund_scale=="watersheds"){
    rem_day_site <- rem_site_buff %>% 
      group_by(Date=floor_date(as.Date(Date),"day"),period,SiteID,Method) %>% 
      summarise(tot_rem=sum(Total),
                # n_events=n(),
                Area_Name=unique(Area_Name),
                site_area_km=unique(site_area_km),
                ea_area_km=unique(ea_area_km),
                effect_area_km=as.numeric(st_area(st_union(geometry))/1e6)) %>% 
      rename(method=Method) %>% 
      filter(Date>=min(sysbait_det_eff$subper_start) &
               Date<=max(sysbait_det_eff$subper_end)) %>% 
      st_drop_geometry()
  } else {
    rem_day_site <- rem_site_buff %>% 
      group_by(Date=floor_date(as.Date(Date),"day"),period,Area_Name,Method) %>% 
      summarise(tot_rem=sum(Total),
                # n_events=n(),
                ea_area_km=unique(ea_area_km),
                effect_area_km=as.numeric(st_area(st_union(geometry))/1e6)) %>% 
      rename(method=Method) %>% 
      filter(Date>=min(sysbait_det_eff$subper_start) &
               Date<=max(sysbait_det_eff$subper_end)) %>% 
      st_drop_geometry()
  }
  
  rem_day_site$effect_area_km[rem_day_site$method=="Aerial"] <- NA
  
  ##clip systematic baiting to removals --------------------
  sysbait_det_eff <-sysbait_det_eff %>% 
    filter(subper_start>=min(rem_day_site$Date) &
             subper_end<=max(rem_day_site$Date))%>% 
    st_drop_geometry() %>% 
    select(-c(season,month,n)) %>%
    rename(site_area_km=area_km)
  
  #merge with aerial ops track data ----------------
  if(abund_scale=="watersheds"){
    rem_day_site <- rem_day_site %>% 
      full_join(ao,by=c("Date","period","method","SiteID","Area_Name")) %>% 
      rename(tot_hrs=flight_time_hr)
  } else {
    rem_day_site <- rem_day_site %>% 
      full_join(ao,by=c("Date","period","method","Area_Name")) %>% 
      rename(tot_hrs=flight_time_hr)
  }
  
  for(i in 1:nrow(rem_day_site)){
    if(is.na(rem_day_site$effect_area_km.x[i])){
      rem_day_site$effect_area_km.x[i] <-  rem_day_site$effect_area_km.y[i]
    }
  }

  
  rem_day_site <- rem_day_site %>% 
    filter(!is.na(effect_area_km.x)) %>% 
    rename(effect_area_km=effect_area_km.x) %>% 
    select(-effect_area_km.y) %>% 
    st_drop_geometry()
  
  rem_day_site$method[rem_day_site$method=="Ground Shoot"] <- "Ground"
  
  #remove systematic baiting ---------------
  eff_day <- eff_day %>% filter(method!="SysBait")
  
  # summarise by site from property  -------------------
  if(abund_scale=="watersheds"){
    eff_day_site <- eff_day %>% 
      group_by(Date,period,method,SiteID) %>% 
      summarise(tot_hrs=sum(tot_hrs),
                num_events=sum(num_events)) 
  } else {
    eff_day_site <- eff_day %>% 
      group_by(Date,period,method,Area_Name) %>% 
      summarise(tot_hrs=sum(tot_hrs),
                num_events=sum(num_events)) 
  }

  #merge removal and effort ------------
  if(abund_scale=="watersheds"){
    rem_eff_site <- eff_day_site %>%
      filter(Date>=min(rem_day_site$Date) & 
               Date<=max(rem_day_site$Date)) %>% 
      full_join(rem_day_site,by=c("Date","period","method","SiteID"))
  } else {
    rem_eff_site <- eff_day_site %>%
      filter(Date>=min(rem_day_site$Date) & 
               Date<=max(rem_day_site$Date)) %>% 
      full_join(rem_day_site,by=c("Date","period","method","Area_Name"))
  }
  
  for(i in 1:nrow(rem_eff_site)){
    if(is.na(rem_eff_site$tot_hrs.x[i])){
      rem_eff_site$tot_hrs.x[i] <-  rem_eff_site$tot_hrs.y[i]
    }
  }
  
  #remove removal observations with no effort associated with them
  rem_eff_site$tot_rem[is.na(rem_eff_site$tot_rem)] <- 0
  
  rem_eff_site <- 
    rem_eff_site %>% 
    filter(!is.na(tot_hrs.x)) %>% 
    rename(tot_hrs=tot_hrs.x) %>% 
    select(-tot_hrs.y) %>%
    filter(!is.na(Area_Name))

  ea_df <- ea %>% st_drop_geometry()
  
  #add ea areas and proportion of ea impacted by removal method
  if(abund_scale=="watersheds"){
    for(i in 1:nrow(rem_eff_site)){
      if(is.na(rem_eff_site$site_area_km[i])){
        rem_eff_site$site_area_km[i] <- study_site_grid$area_km[study_site_grid$SiteID==rem_eff_site$SiteID[i]]
      }

      
      if(is.na(rem_eff_site$prop_site_impact[i])){
        rem_eff_site$prop_site_impact[i] <- rem_eff_site$effect_area_km[i]/rem_eff_site$site_area_km[i]
      }
    }
  } else {
    for(i in 1:nrow(rem_eff_site)){
      if(is.na(rem_eff_site$ea_area_km[i])){
        rem_eff_site$ea_area_km[i] <- ea_df$area_km[ea_df$Area_Name==rem_eff_site$Area_Name[i]]
      }
      
      if(is.na(rem_eff_site$prop_site_impact[i])){
        rem_eff_site$prop_site_impact[i] <- rem_eff_site$effect_area_km[i]/rem_eff_site$ea_area_km[i]
      }
    }
  }

  #remove observations with no area estimate
  rem_eff_site <- rem_eff_site %>% filter(!is.na(effect_area_km))
  rem_eff_site %>% filter(tot_hrs>24) %>% as.data.frame()
  #remove effort observations over 24 hours per day
  if(abund_scale=="watersheds"){
    rem_eff_site <- rem_eff_site %>% 
      group_by(period,method,SiteID) %>% 
      mutate(pass_idx=1:n())# %>% 
      # filter(tot_hrs>24)
  } else {
    rem_eff_site <- rem_eff_site %>% 
      group_by(period,method,Area_Name) %>% 
      mutate(pass_idx=1:n()) #%>% 
      # filter(tot_hrs<24)
  }
  
  #calculate effective hours spent removing per area
  rem_eff_site$effect_area_hrs <- rem_eff_site$tot_hrs/rem_eff_site$effect_area_km
  
  #grid landscape covariates to watershed -----------------------------
  if(file.exists(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Landscape Covariates/nlcd_",occ_scale,".RData"))){
    load(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Landscape Covariates/nlcd_",occ_scale,".RData"))
  } else {
    nlcd_siteid
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
    
    save(nlcd_siteid,file=paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Landscape Covariates/nlcd_",occ_scale,".RData"))
  }
  
  list(nlcd_siteid=nlcd_siteid,
       sysbait_det_eff=sysbait_det_eff,
       rem_eff_site=rem_eff_site)
}

