### Title: Data prep functions for MO dynamic occupancy model
### Author: Abbey Feuka
### Date: 05DEC2023
### Notes: numbering weeks strictly within month

# start_date="2020-09-01"
# end_date=paste0("20",end_fy,"-12-31")

#grid systematic baiting by watershed
grid_sysbait_take <- function(study_site_grid, #occupancy grid/sites
                              start_date="2020-09-01", #start date for data in y-m-d string
                              end_date="2023-09-30", #end date for data in y-m-d string
                              period,
                              data_dir) {#"month"or "season"
  require(tidyverse)
  require(sf)
  require(lubridate)
  
  #section ------------------------------------
  sec <- st_read(file.path(data_dir,"MO_Public_Land_Survey_System","MO_Public_Land_Survey_System.shp"),quiet=T)
  sec <- st_transform(sec,st_crs(study_site_grid))
  sec <- st_intersection(sec,study_site_grid)
  
  #collated MO ops/sys bait data
  sys <- readRDS(file.path(data_dir,"MO_sysbait_clean.rds"))
  
  ##remove NA geometries -------------------------------
  sys$Lat <- as.numeric(sys$Lat)
  sys$Long <- as.numeric(sys$Long)
  sys <- subset(sys, !is.na(Lat))
  
  sys$SiteVisitReason[sys$SiteVisitReason=="<Null>"] <- NA
  sys <- subset(sys,!(SiteVisitReason%in%c("Domestic Observed","Misinformation")))
  
  ##use southern mo only ---------------------------------
  sys <- sys[which(sys$Lat<39),]
  sys <- sys[!duplicated(sys),]
  
  sys_sf <- st_as_sf(sys,coords=c("Long","Lat"),crs="epsg:4326")
  sys_sf <- st_transform(sys_sf,st_crs(study_site_grid)) %>% 
    st_intersection(sec) %>% 
    select(SiteStatus,BaitingStartDate,DateHotBait,
           TrapStartDate,CompleteDate,LABEL_,SiteID,geometry) %>% 
    rename(section = LABEL_)
  
  #filter to study period
  start_date <- as.Date(start_date,format="%Y-%m-%d")
  sys_sf <- sys_sf %>% filter((TrapStartDate>=start_date |is.na(TrapStartDate)))
  sys_sf <- sys_sf %>% filter((BaitingStartDate>=start_date|is.na(BaitingStartDate)))
  
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
  
  max_yr <- max(year(sys_sf$TrapStartDate),year(sys_sf$BaitingStartDate),
                year(sys_sf$DateHotBait),na.rm=T)
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
    filter(BaitingStartDate<TrapStartDate | 
             is.na(BaitingStartDate) | 
             is.na(TrapStartDate)) 
  
  sys_df <- sys_sf %>% st_drop_geometry()
  
  # hog detections occurring three days prior to trapping
  sys_df$HogDetectionSys <- rep(NA,nrow(sys_df))
  for(i in 1:nrow(sys_df)){
    if(!is.na(sys_df$TrapStartDate[i]) & is.na(sys_df$DateHotBait[i])){
      sys_df$HogDetectionSys[i] <- sys_df$TrapStartDate[i] - 3
      sys_df$DateHotBait[i] <- as.Date(sys_df$HogDetectionSys[i],origin="1970-01-01") 
    }
    if(!is.na(sys_df$DateHotBait[i])){
      sys_df$HogDetectionSys[i] <- sys_df$DateHotBait[i] 
    }
  }
  sys_df$HogDetectionSys <- as.Date(sys_df$HogDetectionSys,origin="1970-01-01")
  
  ##clean dates  ---------------------------------------------
  for(i in 1:nrow(sys_df)){
    #no complete date but baiting start date
    if((!is.na(sys_df$BaitingStartDate[i]) & is.na(sys_df$CompleteDate[i]))){
      sys_df$CompleteDate[i] <- sys_df$BaitingStartDate[i] + days(90)
    }
    #no baiting start date but complete date
    if(!is.na(sys_df$CompleteDate[i]) & is.na(sys_df$BaitingStartDate[i])){
      sys_df$BaitingStartDate[i] <- sys_df$CompleteDate[i] - days(90)
    }
    # #no complete date but trap start date
    # if(!is.na(sys_df$TrapStartDate[i]) & is.na(sys_df$CompleteDate[i])){
    #   sys_df$CompleteDate[i] <- sys_df$TrapStartDate[i] + days(90)
    # }
    if(!is.na(sys_df$BaitingStartDate[i]) & !is.na(sys_df$HogDetectionSys[i])){
      if(sys_df$HogDetectionSys[i]<sys_df$BaitingStartDate[i]){
        sys_df$BaitingStartDate[i] <- sys_df$HogDetectionSys[i]
      }
    }else if(!is.na(sys_df$TrapStartDate[i]) & !is.na(sys_df$HogDetectionSys[i])){
      if(sys_df$HogDetectionSys[i]<sys_df$TrapStartDate[i]){
        sys_df$TrapStartDate[i] <- sys_df$HogDetectionSys[i]
      }
    }
    # #complete date is earlier than baiting start
    # if(!(is.na(sys_df$BaitingStartDate[i]) | is.na(sys_df$CompleteDate[i]))){
    #   if(sys_df$BaitingStartDate[i] > sys_df$CompleteDate[i]){
    #     sys_df$CompleteDate[i] <- sys_df$BaitingStartDate[i] + days(90)
    #   }
    # }
  }
  
  sys_df <- sys_df %>% filter(BaitingStartDate<CompleteDate) 
  
  #correct complete date for end of study
  end_date <- as.Date(end_date,format="%Y-%m-%d")
  sys_df$CompleteDate[sys_df$CompleteDate>end_date & 
                         !is.na(sys_df$CompleteDate)] <- end_date
  
  #remove days with no start or end date
  sys_df <- subset(sys_df,!(is.na(BaitingStartDate) & is.na(CompleteDate)))
  
  #remove duplicates
  sys_df <- sys_df[!duplicated(sys_df),]
  
  #add hot baits to capture history ------------------
  for(i in 1:nrow(sys_df)){
    subper[which(subper$subper_start<=sys_df$HogDetectionSys[i] & 
                   subper$subper_end>=sys_df$HogDetectionSys[i]),
           as.character(sys_df$SiteID[i])] <- 1
  }
  
  #removal data ----------------------------------
  rem <- readRDS(file.path(data_dir,"MO_take_clean.rds"))
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
    rem_same <- which(rem_df_trap$Date>= sys_df_orig$BaitingStartDate[i] &
                        rem_df_trap$Date<= sys_df_orig$CompleteDate[i] &
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
  
  sys_df <- sys_df %>% rename(detection_sys=HogDetectionSys) %>%
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
  # i <- which(subper_eff$subper_start>=samp$BaitingStartDate[1]&
  #              subper_eff$subper_end<=samp$CompleteDate[1])[1]
  subper_eff <- subper
  subper_eff[,as.character(unique(study_site_grid$SiteID))] <- 0
  # sys_df_trap <- sys_df %>% filter(det_typ=="detection_trap")
  # sys_df_trap <- sys_df_trap[!duplicated(sys_df_trap),]
  
  for(j in 1:length(unique(study_site_grid$SiteID))){
    samp <- sys_df[which(sys_df$SiteID==j),]
    if(nrow(samp)>0){
      for(i in 1:nrow(subper_eff)){
        idx <- which(samp$BaitingStartDate<=subper_eff$subper_start[i] &
                       samp$CompleteDate>=subper_eff$subper_end[i])
        if(length(idx)>0){
          temp <- samp[idx,]
          temp$BaitingStartDate[temp$BaitingStartDate<=subper_eff$subper_start[i]] <- subper_eff$subper_start[i]
          temp$CompleteDate[temp$CompleteDate>=subper_eff$subper_end[i]] <- subper_eff$subper_end[i]
          temp$days <- temp$CompleteDate-temp$BaitingStartDate +1 #include first day
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
  
  sysbait_det_eff <- subper
  
  list(sysbait_det_eff=sysbait_det_eff,
       sec=sec)
  
}
# sysbait_det_eff <- sys$sysbait_det_eff

#grid removals by elimination area
grid_removals <- function(sec,
                          study_site_grid,
                          sysbait_det_eff, #output from grid_sys_effort, by subperiod (10 days)
                          period,#season" or "month"
                          data_dir
){
  require(tidyverse)
  require(lubridate)
  require(sf)
  # require(terra)
  
  #removal data ----------------------------------
  ## collector data
  rem <- readRDS(file.path(data_dir,"MO_take_clean.rds"))
  
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
                        abund_scale,#Watersheds" or "elim_areas"
                        pass_scale, #"day" or "sys"
                        data_dir
){ 
  require(ROracle)
  require(keyring)
  
  #load effort data ------------------------------
  #data with point locations
  # if(file.exists("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Data/MIS/effort_2020_2024_clean_locs.RData")){
  #   # load("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Data/MIS/effort_2020_2024_clean_locs.RData")
  # } else {
  
  #effort query
  source(file.path(data_dir,"hours_query_mo.R"))
  
  #connect to MIS
  host <- "aap00va0ora407.usda.net "
  port <- 1521
  service <- "wsmis.usda.net" 
  drv <- dbDriver("Oracle")
  
  connect.string <- paste(
    
    "(DESCRIPTION=",
    
    "(ADDRESS=(PROTOCOL=tcp)(HOST=", host, ")(PORT=", port, "))",
    
    "(CONNECT_DATA=(SERVICE_NAME=", service, ")))", sep = "")
  
  
  con <- dbConnect(drv, 
                   username = "wsdata", 
                   password = key_get("mis","wsdata"), 
                   dbname = connect.string)
  
  eff=dbSendQuery(con,hours_query)
  eff=fetch(eff)
  eff=as.data.frame(eff)
  dbDisconnect(con)

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
  smws <- st_read(file.path(data_dir,"HUC10_MO","huc10_siteIDs_cond_EA.shp"))
  
  eff_sf <- eff_sf %>% st_transform(st_crs(smws)) %>% st_intersection(smws)
  
  eff <- eff_sf
  eff <- eff %>% rename(Date=WT_WORK_DATE,Method=CMP_NAME)
  
  ## only use MIS effort for traps, ground shooting, systematic baiting
  eff$Method[eff$Method %in% c("TRAPS, LIVE, FERAL HOGS",
                               "TRAPS, BODY GRIP","TRAPS, CORRAL",
                               "TRAPS, CAGE")] <- "Trap"
  eff$Method[eff$Method %in% c("FIREARMS")] <- "Ground"
  eff$Method[eff$Method %in% c("HELICOPTER")] <- "Aerial"
  eff$Method[eff$Method %in% c("PREBAIT","BAIT STATION")] <- "SysBait"
  eff <- eff %>% filter(Method %in% c("Trap","Ground","Aerial","SysBait"))
  eff$Date <- as.Date(eff$Date,format="%d-%b-%y")
  
  eff <- eff %>% filter(UOM_NAME%in%c("MINUTES","HOURS","HOBBS METER"))
  
  eff$time_hr <-  eff$WTM_QTY
  eff$time_hr[eff$UOM_NAME=="MINUTES"] <- eff$WTM_QTY[eff$UOM_NAME=="MINUTES"]/60
  
  #only hobbs meter readings for aerial ops
  eff <- eff %>% filter(!(Method=="Aerial" & UOM_NAME!="HOBBS METER"))
  
  eff$UOM_NAME[eff$UOM_NAME=="HOBBS METER"] <- "HOURS"
  
  eff<- eff %>% filter(USET_NAME%in%c("APPLIED/USED","CHECKED","RESET","SET"))
  
  eff <- eff[!duplicated(eff),]
  
  #systematic baiting effort hours ---------------
  sysbait_det_eff$eff_hrs <- rep(0,nrow(sysbait_det_eff)) 
  for(i in 1:nrow(sysbait_det_eff)){
    x <- eff[eff$SiteID==sysbait_det_eff$site_idx[i] &
               eff$Date>=sysbait_det_eff$subper_start[i] &
               eff$Date<=sysbait_det_eff$subper_end[i],]
    
    if(nrow(x)>0){
      if(any(x$Method%in%c("SysBait","Trap"))){ #systematic baiting effort 
        sysbait_det_eff$eff_hrs[i] <- sum(subset(x,Method=="SysBait" | Method=="Trap")$time_hr,na.rm=T)
      }
    }
  }
  
  ##add NAs for missing hours ------------------------
  sysbait_det_eff$eff_hrs[sysbait_det_eff$eff_hrs==0 & sysbait_det_eff$trap_nights>0] <- NA
  
  ##standardize effort by area --------------------
  sysbait_det_eff <- sysbait_det_eff %>% 
    left_join(study_site_grid %>% mutate(site_idx=as.character(SiteID)))
  sysbait_det_eff$eff_hrs_km <- sysbait_det_eff$eff_hrs/sysbait_det_eff$area_km
  
  sysbait_det_eff$trap_nights_km <- sysbait_det_eff$trap_nights/sysbait_det_eff$area_km
  
  # removal effort -------------------------------------------------
  ##elimination area - areas ---------------
  ea <- study_site_grid %>% 
    group_by(Area_Name) %>% 
    summarise(geometry=st_union(geometry),
              area_km=as.numeric(st_area(geometry))/1e6)
  
  eff_day <- eff %>% 
    st_drop_geometry() %>% 
    group_by(Date=floor_date(Date,"day"),
             WT_AGRPROP_ID,Method) %>% 
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
  
  #calculate area buffers for ground and trap --------------
  rem_site_trap <- rem_site_sf %>% filter(Method=="Trap") %>% 
    mutate(geometry=st_buffer(geometry,dist=sqrt(6.7e6/3.14)))
  rem_site_ground <- rem_site_sf %>% filter(Method=="Ground Shoot") %>% 
    mutate(geometry=st_buffer(geometry,dist=sqrt(5e6/3.14)))
  rem_site_aerial <- rem_site_sf %>% filter(Method=="Aerial")
  rem_site_buff <- rbind(rem_site_aerial,rem_site_trap,rem_site_ground)

  #merge all trap/ground footprints to calculate effective area --------------
  if(abund_scale=="watersheds"){
    if(pass_scale=="day"){
      rem_pass_site <- rem_site_buff %>% 
        group_by(Date=floor_date(as.Date(Date),"day"),period,SiteID,Method) %>% 
        summarise(tot_rem=sum(Total),
                  # n_events=n(),
                  Area_Name=unique(Area_Name),
                  site_area_km=unique(site_area_km),
                  ea_area_km=unique(ea_area_km),
                  effect_area_km=as.numeric(st_area(st_union(geometry))/1e6)) %>% 
        filter(Date>=min(sysbait_det_eff$subper_start) &
                 Date<=max(sysbait_det_eff$subper_end)) %>% 
        st_drop_geometry()
    } else if(pass_scale=="sys"){
      
      period_dates <- sysbait_det_eff %>% group_by(period) %>% 
        summarise(per_start=min(subper_start),
                  per_end=max(subper_end))
      
      subper_dates <- list()
      for(t in 1:nrow(period_dates)){
        n<-max(unique(sysbait_det_eff$subper_end-sysbait_det_eff$subper_start)) +1
        from <- as.Date(period_dates$per_start[t])
        to <- as.Date(period_dates$per_end[t])
        by = paste(n, "days")
        subper_starts <- seq(from, to - n + 1, by = by)
        subper_ends <- seq(from + n - 1, to, by = by)
        last <- length(subper_starts)
        remaining <- to - subper_ends[last]
  
        if(remaining != 0){
          subper_starts <- c(subper_starts, subper_ends[last] + 1)
          subper_ends <- c(subper_ends, to)
        }
        subper_dates[[t]] <- data.frame(per_start=period_dates$per_start[t],
                   per_end=period_dates$per_end[t],
                   subper_starts, subper_ends,
                   subper_id = 1:length(subper_starts))
      }
      subper_dates <- do.call("rbind.data.frame",subper_dates)
      period_dates <- period_dates %>% left_join(subper_dates)
      
      #remove singleton subperiods, combine into larger
      period_dates$ndays <- period_dates$subper_ends-period_dates$subper_starts
      singles <- period_dates %>% filter(ndays<=1)

      for(i in 1:nrow(period_dates)){
        for(j in 1:nrow(singles)){
          if(period_dates$period[i]==singles$period[j] &
             period_dates$per_end[i]>=singles$subper_starts[j]){
            period_dates$subper_ends[period_dates$period==singles$period[j] &
                               period_dates$subper_id==(n-1)] <- 
              period_dates$subper_ends[period_dates$period==singles$period[j] &
                                 period_dates$subper_id==n]
          }
        }
      }
      period_dates$ndays <- period_dates$subper_ends-period_dates$subper_starts
      period_dates <- period_dates %>% filter(ndays>1)
      
      rem_site_buff$subper_id <-rem_site_buff$subper_start <- rem_site_buff$subper_end <- NA
      for(i in 1:nrow(rem_site_buff)){
        samp <- as.Date(rem_site_buff$Date[i])
        per_samp <- period_dates[period_dates$subper_starts<= samp &
                       period_dates$subper_ends>= samp,]
        rem_site_buff$subper_id[i] <- per_samp$subper_id
        rem_site_buff$subper_start[i] <- per_samp$subper_starts
        rem_site_buff$subper_end[i] <- per_samp$subper_ends
      }
      rem_site_buff$subper_end <- as.Date(rem_site_buff$subper_end)
      rem_site_buff$subper_start <- as.Date(rem_site_buff$subper_start)
      
      rem_pass_site <- rem_site_buff %>% 
        group_by(subper_id,period,SiteID,Method,subper_start,subper_end) %>% 
        summarise(tot_rem=sum(Total),
                  Area_Name=unique(Area_Name),
                  site_area_km=unique(site_area_km),
                  ea_area_km=unique(ea_area_km),
                  effect_area_km=as.numeric(st_area(st_union(geometry))/1e6)) %>% 
        filter(subper_start>=min(sysbait_det_eff$subper_start) &
                 subper_end<=max(sysbait_det_eff$subper_end)) %>% 
        st_drop_geometry()
    }
    
  } else {
    rem_pass_site <- rem_site_buff %>% 
      group_by(Date=floor_date(as.Date(Date),"day"),period,Area_Name,Method) %>% 
      summarise(tot_rem=sum(Total),
                # n_events=n(),
                ea_area_km=unique(ea_area_km),
                effect_area_km=as.numeric(st_area(st_union(geometry))/1e6)) %>% 
      filter(Date>=min(sysbait_det_eff$subper_start) &
               Date<=max(sysbait_det_eff$subper_end)) %>% 
      st_drop_geometry()
  }

  rem_pass_site$effect_area_km[rem_pass_site$Method=="Aerial"] <- NA
  
  ##clip systematic baiting to removals --------------------
  sysbait_det_eff <-sysbait_det_eff %>% 
    filter(subper_start>=min(rem_pass_site$subper_start) &
             subper_end<=max(rem_pass_site$subper_end)) %>% 
    st_drop_geometry() %>% 
    select(-c(season,month,n)) %>%
    rename(site_area_km=area_km)
  
  #merge with aerial ops track data ----------------
  if(abund_scale=="watersheds"){
    ao <- read.csv(file.path(data_dir,paste0("fy21_",end_fy,"_flight_time_ssg_clean_528m_buff.csv")))
    ao$Method<-"Aerial"
    
    ao <- ao %>% 
      rename(prop_site_impact=prop_site_flown) %>% 
      mutate(Date=as.Date(Date),
             Area_Name=as.character(Area_Name),
             effect_area_km=as.numeric(intersect_area_km)) %>% 
      select(SiteID,Area_Name,Date,Method,flight_time_hr,effect_area_km,prop_site_impact) %>% 
      filter(Date<=as.Date(max(sysbait_det_eff$subper_end)))
    
    ao$period <- NA
    for(i in 1:nrow(ao)){
      ao$period[i] <- unique(period_dates$period[period_dates$per_start<=ao$Date[i] &
                                            period_dates$per_end>=ao$Date[i]])
    }
    
    if(pass_scale=="sys"){

      ao$subper_id <-ao$subper_start <- ao$subper_end <- NA
      for(i in 1:nrow(ao)){
        samp <- as.Date(ao$Date[i])
        per_samp <- period_dates[period_dates$subper_starts<= samp &
                                   period_dates$subper_ends>= samp,]
        ao$subper_id[i] <- per_samp$subper_id
        ao$subper_start[i] <- per_samp$subper_starts
        ao$subper_end[i] <- per_samp$subper_ends
      }
      
      ao <- ao %>% group_by(SiteID,Area_Name,Method,period,subper_id,subper_start,subper_end) %>% 
        summarise(flight_time_hr=sum(flight_time_hr,na.rm=T),
                  effect_area_km=sum(effect_area_km,na.rm=T)) %>% 
        mutate(subper_start=as.Date(subper_start),
               subper_end=as.Date(subper_end)) %>% 
        left_join(study_site_grid %>% st_drop_geometry()) %>% 
        mutate(prop_site_impact=effect_area_km/area_km)
    }

  } else {
    ao <-read.csv(file.path(data_dir,paste0("fy21_",end_fy,"_flight_time_ssg_ea_clean_528m_buff.csv")))
    ao$Method<-"Aerial"
    
    ao <- ao %>% 
      rename(prop_site_impact=prop_ea_flown) %>% 
      mutate(Date=as.Date(Date),
             Area_Name=as.character(Area_Name),
             effect_area_km=as.numeric(intersect_area_km)) %>% 
      select(Area_Name,Date,Method,flight_time_hr,effect_area_km,prop_site_impact) %>% 
      filter(Date<=as.Date(max(sysbait_det_eff$subper_end)))
  }

  
  if(abund_scale=="watersheds"){
    if(pass_scale=="day"){
      rem_pass_site <- rem_pass_site %>% 
        full_join(ao,by=c("Date","period","Method","SiteID","Area_Name")) %>% 
        rename(tot_hrs=flight_time_hr)
      
    } else if(pass_scale=="sys"){
      rem_pass_site <- rem_pass_site %>% 
        full_join(ao,by=c("subper_id","subper_start","subper_end","period","Method","SiteID","Area_Name")) %>% 
        rename(tot_hrs=flight_time_hr)
    }

  } else {
    rem_pass_site <- rem_pass_site %>% 
      full_join(ao,by=c("Date","period","Method","Area_Name")) %>% 
      rename(tot_hrs=flight_time_hr)
  }
  
  for(i in 1:nrow(rem_pass_site)){
    if(is.na(rem_pass_site$effect_area_km.x[i])){
      rem_pass_site$effect_area_km.x[i] <-  rem_pass_site$effect_area_km.y[i]
    }
  }
  
  rem_pass_site <- rem_pass_site %>% 
    filter(!is.na(effect_area_km.x)) %>% 
    rename(effect_area_km=effect_area_km.x) %>% 
    select(-effect_area_km.y) %>% 
    st_drop_geometry()
  
  rem_pass_site$Method[rem_pass_site$Method=="Ground Shoot"] <- "Ground"
  
  #remove systematic baiting ---------------
  eff_day <- eff_day %>% filter(Method!="SysBait")
  
  # summarise by site from property  -------------------
  if(abund_scale=="watersheds"){
    if(pass_scale=="day"){
      eff_pass_site <- eff_day %>% 
        group_by(Date,period,Method,SiteID) %>% 
        summarise(tot_hrs=sum(tot_hrs),
                  num_events=sum(num_events)) 
    } else if(pass_scale=="sys"){
      
      eff_day$subper_id <-eff_day$subper_start <- eff_day$subper_end <- NA
      for(i in 1:nrow(eff_day)){
        samp <- as.Date(eff_day$Date[i])
        per_samp <- period_dates[period_dates$subper_starts<= samp &
                                   period_dates$subper_ends>= samp,]
        eff_day$subper_id[i] <- per_samp$subper_id
        eff_day$subper_start[i] <- per_samp$subper_starts
        eff_day$subper_end[i] <- per_samp$subper_ends
      }
      
      eff_pass_site <- eff_day %>% 
        group_by(subper_id,subper_start,subper_end,period,Method,SiteID) %>% 
        summarise(tot_hrs=sum(tot_hrs),
                  num_events=sum(num_events)) %>% 
        mutate(subper_start=as.Date(subper_start),
               subper_end=as.Date(subper_end))
    }

  } else {
    eff_pass_site <- eff_day %>% 
      group_by(Date,period,Method,Area_Name) %>% 
      summarise(tot_hrs=sum(tot_hrs),
                num_events=sum(num_events)) 
  }

  #merge removal and effort ------------
  if(abund_scale=="watersheds"){
    if(pass_scale=="day"){
      rem_eff_site <- eff_pass_site %>%
        filter(Date>=min(rem_pass_site$Date) & 
                 Date<=max(rem_pass_site$Date)) %>% 
        full_join(rem_pass_site,by=c("Date","period","Method","SiteID"))
      
    } else if(pass_scale=="sys"){
      rem_eff_site <- eff_pass_site %>%
        filter(subper_start>=min(rem_pass_site$subper_start) & 
                 subper_end<=max(rem_pass_site$subper_end)) %>% 
        full_join(rem_pass_site,by=c("subper_id","subper_start","subper_end","period","Method","SiteID"))
      
    }

  } else {
    rem_eff_site <- eff_pass_site %>%
      filter(Date>=min(rem_pass_site$Date) & 
               Date<=max(rem_pass_site$Date)) %>% 
      full_join(rem_pass_site,by=c("Date","period","Method","Area_Name"))
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
  

  if(pass_scale=="day"){
    #remove effort observations over 24 hours per day
    rem_eff_site %>% filter(tot_hrs>24) %>% as.data.frame()
  }
  
  #add pass idx
  if(abund_scale=="watersheds"){
    rem_eff_site <- rem_eff_site %>% 
      group_by(period,Method,SiteID) %>% 
      mutate(pass_idx=1:n())# %>% 
    # filter(tot_hrs>24)
  } else {
    rem_eff_site <- rem_eff_site %>% 
      group_by(period,Method,Area_Name) %>% 
      mutate(pass_idx=1:n()) #%>% 
    # filter(tot_hrs<24)
  }

  #calculate effective hours spent removing per area
  rem_eff_site$effect_area_hrs <- rem_eff_site$tot_hrs/rem_eff_site$effect_area_km
  
  list(sysbait_det_eff=sysbait_det_eff,
       rem_eff_site=rem_eff_site)
  
}

grid_covs <- function(occ_scale,
                      abund_scale,
                      sysbait_det_eff,
                      study_site_grid,
                      data_dir
                      ){
  
  #grid landscape covariates to watershed -----------------------------
  if(file.exists(file.path(data_dir,paste0("nlcd_",occ_scale,".rds")))){
    nlcd_siteid <- readRDS(file.path(data_dir,paste0("nlcd_",occ_scale,".rds")))
  } else {
    #TODO source extract_nlcd_watersheds
  }

  #oak cover -----------
  if(file.exists(file.path(data_dir,"oak_siteid.rds"))){
    oak_siteid <- readRDS(file.path(data_dir,"oak_siteid.rds"))
  } else {
    #TODO: source extract_oak_watersheds
  }
   
  ##add proportion of area covered by feral swine distribution ----------------------
  ## LAGGED a year
  if(file.exists(file.path(data_dir,paste0("ssg_nfsp_2015_",
                        max(year(sysbait_det_eff$subper_end)),".rds")))){
    nfsp_siteid <- readRDS(file.path(data_dir,paste0("ssg_nfsp_2015_",
                                                  max(year(sysbait_det_eff$subper_end)),".rds")))
  } else {
    #TODO source watershed_nfsp_overlap.R to update years
  }
  
  list(nlcd_siteid=nlcd_siteid,
       oak_siteid=oak_siteid,
       nfsp_siteid=nfsp_siteid
  )
}
