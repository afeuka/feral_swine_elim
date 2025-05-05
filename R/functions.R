### Title: Data prep functions for MO dynamic occupancy model
### Author: Abbey Feuka
### Date: 05DEC2023
### Notes: numbering weeks strictly within month

format_ssg <- function(ssg){
  
  st_read(ssg) %>% 
    rename(elim_area_idx=elm_r_d,
           Area_Name=Area_Nm)
}

make_elim_areas <- function(ssg_formatted){
  ssg_formatted %>%
    group_by(Area_Name,elim_area_idx) %>% 
    summarise(geometry=st_union(geometry))
}

format_sec <- function(sec,
                       ssg_formatted){
  
  st_read(sec) %>% 
    st_transform(st_crs(ssg_formatted)) %>% 
    st_intersection(ssg_formatted)
}

#grid systematic baiting by watershed
grid_sysbait <- function(sysbait_clean,
                         take_clean,
                         sec_formatted,
                         ssg_formatted,
                         start_date="2020-10-01", #start date for data in y-m-d string
                         end_date="2024-09-30" #end date for data in y-m-d string
                         ){
  
  ##remove NA geometries -------------------------------
  sysbait_clean$Lat <- as.numeric(sysbait_clean$Lat)
  sysbait_clean$Long <- as.numeric(sysbait_clean$Long)
  sysbait_clean <- subset(sysbait_clean, !is.na(Lat))

  sysbait_clean$SiteVisitReason[sysbait_clean$SiteVisitReason=="<Null>"] <- NA
  sysbait_clean <- subset(sysbait_clean,!(SiteVisitReason%in%c("Domestic Observed","Misinformation")))

  ##use southern mo only ---------------------------------
  sysbait_clean <- sysbait_clean[which(sysbait_clean$Lat<39),]
  sysbait_clean <- sysbait_clean[!duplicated(sysbait_clean),]

  sys_sf <- st_as_sf(sysbait_clean,coords=c("Long","Lat"),crs="epsg:4326")
  sys_sf <- st_transform(sys_sf,st_crs(ssg_formatted)) %>%
    st_intersection(sec_formatted) %>%
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
  per_start <- per_start[-which(per_start<as.Date(start_date))]
  sysbait_grid <- c(sapply(1:length(per_start),function(i){
    if(i<length(per_start)){
      seq(per_start[i],per_start[i+1],by=10)
    } else {
      seq(per_start[i],per_start[i] %m+% months(3),by=10)
    }}))
  sysbait_grid <- as.Date(sysbait_grid,origin="1970-01-01")

  sysbait_grid<- data.frame(subper_start=sysbait_grid)
  sysbait_grid$subper_end <- sysbait_grid$subper_start+days(9)

  sysbait_grid$year <- year(sysbait_grid$subper_start)
  sysbait_grid$month <- month(sysbait_grid$subper_start)
  sysbait_grid$fy <- sysbait_grid$year
  sysbait_grid$fy[sysbait_grid$month%in%c(10,11,12)] <- sysbait_grid$year[sysbait_grid$month%in%c(10,11,12)]+1

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

  sysbait_grid$season <- seasons$scode[month(sysbait_grid$subper_start)]
  sysbait_grid$period <- paste(sysbait_grid$season,sysbait_grid$year)
  sysbait_grid$period <- as.numeric(as.factor(sysbait_grid$period))
  sysbait_grid$period <- as.numeric(factor(sysbait_grid$period,levels=unique(sysbait_grid$period)))

  sysbait_grid$subper_end_new <-NA
  #cut overlapping periods
  for(i in 1:nrow(sysbait_grid)){
    if(i<nrow(sysbait_grid)){
      if(sysbait_grid$subper_end[i]>sysbait_grid$subper_start[i+1]){
        sysbait_grid$subper_end_new[i] <- sysbait_grid$subper_start[i+1]-1
      } else {
        sysbait_grid$subper_end_new[i] <- sysbait_grid$subper_end[i]
      }
    } else {
      sysbait_grid$subper_end_new[i] <- sysbait_grid$subper_end[i]
    }
  }
  sysbait_grid <- sysbait_grid %>%
    mutate(subper_end_new = as.Date(subper_end_new)) %>%
    dplyr::select(-subper_end) %>%
    rename(subper_end=subper_end_new)

  sysbait_grid <- cbind.data.frame(sysbait_grid,matrix(0,nrow=nrow(sysbait_grid),ncol=length(unique(ssg_formatted$SiteID))))
  sysbait_grid <- as.data.frame(sysbait_grid)

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
    sysbait_grid[which(sysbait_grid$subper_start<=sys_df$HogDetectionSys[i] &
                   sysbait_grid$subper_end>=sys_df$HogDetectionSys[i]),
           as.character(sys_df$SiteID[i])] <- 1
  }

  #removal data ----------------------------------
  take_clean$Date <- as.POSIXct(take_clean$Date,format="%Y-%m-%d %H:%M:%S",tz="CST6CDT")
  rem_sf <- st_as_sf(take_clean,coords=c("Long","Lat"),crs="epsg:4326")

  ## grid to study_site ------------------------------------------
  rem_sf <- st_transform(rem_sf, st_crs(ssg_formatted))
  rem_sf_trap <- st_intersection(rem_sf,ssg_formatted) %>%
    st_intersection(sec_formatted) %>%
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
      sysbait_grid[which(sysbait_grid$subper_start<=sys_df$det_date[i] &
                     sysbait_grid$subper_end>=sys_df$det_date[i]),
             as.character(sys_df$SiteID[i])] <- 1
    }
  }

  # effort in trap days------------------------------------
  # j<-7
  # i <- which(subper_eff$subper_start>=samp$BaitingStartDate[1]&
  #              subper_eff$subper_end<=samp$CompleteDate[1])[1]
  subper_eff <- sysbait_grid
  subper_eff[,as.character(unique(ssg_formatted$SiteID))] <- 0
  # sys_df_trap <- sys_df %>% filter(det_typ=="detection_trap")
  # sys_df_trap <- sys_df_trap[!duplicated(sys_df_trap),]

  for(j in 1:length(unique(ssg_formatted$SiteID))){
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

  sysbait_grid <- sysbait_grid %>% pivot_longer(cols=as.character(ssg_formatted$SiteID),
                                    names_to="site_idx",values_to="detection")
  subper_eff <- subper_eff %>% pivot_longer(cols=as.character(ssg_formatted$SiteID),
                                            names_to="site_idx",values_to="trap_nights")

  sysbait_grid <- sysbait_grid %>% left_join(subper_eff)

  #change detections with no trap nights to 0
  sysbait_grid$detection[sysbait_grid$detection>sysbait_grid$trap_nights] <- 0

  #filter to study period
  sysbait_grid <- sysbait_grid %>% filter(subper_start<=as.Date(end_date))
  
}

#grid removals by elimination area
grid_removals <- function(take_clean,
                          sysbait_grid,
                          sec_formatted,
                          ssg_formatted
){
  
  take_clean$Date <- as.POSIXct(take_clean$Date,tz="CST6CDT",format="%Y-%m-%d %H:%M:%S")
  rem_sf <- st_as_sf(take_clean,coords=c("Long","Lat"),crs="epsg:4326")

  rem_sf <- st_transform(rem_sf, st_crs(ssg_formatted))
  rem_site_sf <- st_intersection(rem_sf,ssg_formatted) %>%
    select(Date,SiteID,Total,Males,Females,Adults,Sub_adults,Method,NonBreed_M,
           NonBreed_F,Adult_M,Adult_F,County) %>%
    filter(year(Date)>=2020)

  rem_sf <- rem_sf %>% select(Date,Total,Method,geometry) %>%
    st_intersection(ssg_formatted)

  #intersect with section
  rem_site_sf <- rem_site_sf %>%
    st_intersection(sec_formatted %>%
                      rename(section=LABEL_) %>%
                      select(section,geometry))
  #remove duplicates
  #might have different lat/longs but in same section = same trap
  #removes differences in individual hog counts (m/f, adult/juv), uses totals only
  rem_site_sf <- rem_site_sf[!duplicated(rem_site_sf),]
  rem_site_sf <- rem_site_sf %>%
    # st_drop_geometry() %>%
    left_join(ssg_formatted %>%
                st_drop_geometry() %>%
                select(SiteID,Area_Name,area_km) %>%
                rename(site_area_km=area_km))

  rem_site_sf <- rem_site_sf %>%
    left_join(ssg_formatted %>%
                group_by(Area_Name) %>%
                summarise(geometry=st_union(geometry),
                          ea_area_km=as.numeric(st_area(geometry))/1e6) %>%
                st_drop_geometry()) %>%
    select(Date,SiteID,Total,Method,section,Area_Name,site_area_km,ea_area_km) %>%
    filter(Method!="Night Shoot")

  period_dates <- sysbait_grid %>% group_by(period) %>%
    summarise(per_start=min(subper_start),
              per_end=max(subper_end))

  rem_site_sf <- rem_site_sf %>% filter(Date>=min(period_dates$per_start) &
                                          Date<=max(period_dates$per_end))
  rem_site_sf$period <- NA
  for(i in 1:nrow(rem_site_sf)){
    rem_site_sf$period[i] <- period_dates$period[period_dates$per_start<=as.Date(rem_site_sf$Date[i]) &
                                                   period_dates$per_end>=as.Date(rem_site_sf$Date[i])]
  }
  rem_site_sf
}

mis_effort_query <- function(take_grid){
  
  end_date <- 
  
  paste("SELECT M2_AGREEMENT_PROPERTY.AGRP_PRP_ID,",
        "M2_WORK_TASK.WT_AGRPROP_ID,",
        "M2_WORK_TASK_UOM.WTM_WT_ID,",
        "M2_WORK_TASK.WT_EMP_ID,",
        "M2_WORK_TASK.WT_WORK_DATE,",
        "M2_WORK_TASK_UOM.WTM_QTY,",
        "M2_UNIT_OF_MEASURE.UOM_NAME,",
        "M2_DAMAGE_AGENT.ID,",
        "M2_DAMAGE_AGENT.DA_NAME,",
        "M2_COMPONENT.CMP_NAME,",
        "M2_COMPONENT.CMP_TYPE,",
        "M2_WORK_TASK.WT_STATE,",
        "M2_COUNTY.CNTY_NAME,",
        "M2_PROPERTY.PRP_NAME,",
        "M2_PROPERTY_UOL.PRPU_N_LAT,",
        "M2_PROPERTY_UOL.PRPU_E_LONG,",
        "M2_USE_TYPE.USET_NAME",
        "FROM M2_WORK_TASK",
        "INNER JOIN M2_CONFLICT_WORK",
        "ON M2_WORK_TASK.ID = M2_CONFLICT_WORK.CNFW_WT_ID",
        "INNER JOIN M2_DAMAGE_AGENT",
        "ON M2_DAMAGE_AGENT.ID = M2_CONFLICT_WORK.CNFW_DA_ID",
        "INNER JOIN M2_WORK_TASK_UOM",
        "ON M2_WORK_TASK.ID = M2_WORK_TASK_UOM.WTM_WT_ID",
        "INNER JOIN M2_UNIT_OF_MEASURE",
        "ON M2_UNIT_OF_MEASURE.ID = M2_WORK_TASK_UOM.WTM_UOM_ID",
        "LEFT JOIN M2_AGREEMENT_PROPERTY",
        "ON M2_AGREEMENT_PROPERTY.ID = M2_WORK_TASK.WT_AGRPROP_ID",
        "LEFT JOIN M2_WORK_TASK_COMPONENT",
        "ON M2_WORK_TASK_COMPONENT.WTC_WT_ID = M2_WORK_TASK.ID",
        "LEFT JOIN M2_COMPONENT",
        "ON M2_COMPONENT.ID = M2_WORK_TASK_COMPONENT.WTC_CMP_ID",
        "LEFT JOIN M2_PROPERTY",
        "ON M2_PROPERTY.ID = M2_AGREEMENT_PROPERTY.AGRP_PRP_ID",
        "LEFT JOIN M2_COUNTY",
        "ON M2_PROPERTY.PRP_CNTY_ID = M2_COUNTY.ID",
        "LEFT JOIN M2_PROPERTY_UOL",
        "ON M2_PROPERTY_UOL.PRPU_PRP_ID = M2_PROPERTY.ID",
        "LEFT JOIN M2_USE_TYPE",
        "ON M2_USE_TYPE.ID = M2_WORK_TASK_COMPONENT.WTC_USET_ID",
        "WHERE M2_DAMAGE_AGENT.ID=8 ",
        "AND WT_WORK_DATE>='01-SEP-20'",
        paste(
          " AND WT_WORK_DATE <= to_date('",
          as.Date(max(take_grid$Date,na.rm=T)),
          "', 'YYYY-MM-DD') ",
          sep=""
        ),
        " AND WT_STATE='MO'",
        "AND CMP_NAME!='SNARES, NECK'")
}

pull_effort_mis <- function(effort_query,
                            ssg_formatted){
  
  #connect to MIS
  host <- "xxxxxxxxx"
  port <- 1521
  service <- "xxxxxxx"
  drv <- dbDriver("Oracle")
  
  connect.string <- paste(
    
    "(DESCRIPTION=",
    
    "(ADDRESS=(PROTOCOL=tcp)(HOST=", host, ")(PORT=", port, "))",
    
    "(CONNECT_DATA=(SERVICE_NAME=", service, ")))", sep = "")
  
  
  con <- dbConnect(drv,
                   username = "wsdata",
                   password = key_get("mis","wsdata"),
                   dbname = connect.string)
  
  eff=dbSendQuery(con,effort_query)
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
  
  eff_sf <- eff_sf %>% st_transform(st_crs(ssg_formatted)) %>% 
    st_intersection(ssg_formatted)
  
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
  
  eff
}

grid_sysbait_effort <- function(effort_raw, 
                                sysbait_grid,
                                ssg_formatted){

  sysbait_grid$eff_hrs <- rep(0,nrow(sysbait_grid))
  for(i in 1:nrow(sysbait_grid)){
    x <- effort_raw[effort_raw$SiteID==sysbait_grid$site_idx[i] &
                      effort_raw$Date>=sysbait_grid$subper_start[i] &
                      effort_raw$Date<=sysbait_grid$subper_end[i],]
    
    if(nrow(x)>0){
      if(any(x$Method%in%c("SysBait","Trap"))){ #systematic baiting effort
        sysbait_grid$eff_hrs[i] <- sum(subset(x,Method=="SysBait" | Method=="Trap")$time_hr,na.rm=T)
      }
    }
  }
  
  ##add NAs for missing hours ------------------------
  sysbait_grid$eff_hrs[sysbait_grid$eff_hrs==0 & sysbait_grid$trap_nights>0] <- NA
  
  ##standardize effort by area --------------------
  sysbait_grid <- sysbait_grid %>%
    left_join(ssg_formatted %>% mutate(site_idx=as.character(SiteID)))
  sysbait_grid$eff_hrs_km <- sysbait_grid$eff_hrs/sysbait_grid$area_km
  
  sysbait_grid$trap_nights_km <- sysbait_grid$trap_nights/sysbait_grid$area_km
  
  sysbait_grid
}

grid_take_effort <- function(sysbait_effort,
                             take_grid,
                             effort_raw,
                             aerial_ops_tracks,
                             ssg_formatted){
  
  abund_scale <- "watersheds"
  pass_scale <- "sys"
  
  # removal effort -------------------------------------------------
  ##elimination area - areas ---------------
  ea <- ssg_formatted %>%
    group_by(Area_Name) %>%
    summarise(geometry=st_union(geometry),
              area_km=as.numeric(st_area(geometry))/1e6)

  eff_day <- effort_raw %>%
    st_drop_geometry() %>%
    group_by(Date=floor_date(Date,"day"),
             WT_AGRPROP_ID,Method) %>%
    summarise(tot_hrs=sum(time_hr),#/length(unique(WT_EMP_ID)),
              num_events=n(),
              Area_Name=unique(Area_Name),
              SiteID=unique(SiteID)) %>%
    left_join(ssg_formatted %>% st_drop_geometry() %>%
                rename(site_area_km=area_km)%>%
                select(-n)) %>%
    left_join(ea %>% st_drop_geometry() %>%
                rename(ea_area_km=area_km))

  period_dates <- sysbait_effort %>% group_by(period) %>%
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
  rem_site_trap <- take_grid %>% filter(Method=="Trap") %>%
    mutate(geometry=st_buffer(geometry,dist=sqrt(6.7e6/3.14)))
  rem_site_ground <- take_grid %>% filter(Method=="Ground Shoot") %>%
    mutate(geometry=st_buffer(geometry,dist=sqrt(5e6/3.14)))
  rem_site_aerial <- take_grid %>% filter(Method=="Aerial")
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
        filter(Date>=min(sysbait_effort$subper_start) &
                 Date<=max(sysbait_effort$subper_end)) %>%
        st_drop_geometry()
    } else if(pass_scale=="sys"){

      period_dates <- sysbait_effort %>% group_by(period) %>%
        summarise(per_start=min(subper_start),
                  per_end=max(subper_end))

      subper_dates <- list()
      for(t in 1:nrow(period_dates)){
        n<-max(unique(sysbait_effort$subper_end-sysbait_effort$subper_start)) +1
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
        filter(subper_start>=min(sysbait_effort$subper_start) &
                 subper_end<=max(sysbait_effort$subper_end)) %>%
        st_drop_geometry()
    }

  } else {
    rem_pass_site <- rem_site_buff %>%
      group_by(Date=floor_date(as.Date(Date),"day"),period,Area_Name,Method) %>%
      summarise(tot_rem=sum(Total),
                # n_events=n(),
                ea_area_km=unique(ea_area_km),
                effect_area_km=as.numeric(st_area(st_union(geometry))/1e6)) %>%
      filter(Date>=min(sysbait_effort$subper_start) &
               Date<=max(sysbait_effort$subper_end)) %>%
      st_drop_geometry()
  }

  rem_pass_site$effect_area_km[rem_pass_site$Method=="Aerial"] <- NA

  #merge with aerial ops track data ----------------
  if(abund_scale=="watersheds"){
    ao <- aerial_ops_tracks 
    
    ao$Method<-"Aerial"

    ao <- ao %>%
      rename(prop_site_impact=prop_site_flown) %>%
      mutate(Date=as.Date(Date),
             Area_Name=as.character(Area_Name),
             effect_area_km=as.numeric(intersect_area_km)) %>%
      select(SiteID,Area_Name,Date,Method,flight_time_hr,effect_area_km,prop_site_impact) %>%
      filter(Date<=as.Date(max(sysbait_effort$subper_end)))

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
        left_join(ssg_formatted %>% st_drop_geometry()) %>%
        mutate(prop_site_impact=effect_area_km/area_km)
    }

  } else {
    # ao <- aerial_ops_tracks
    # ao$Method<-"Aerial"
    #
    # ao <- ao %>%
    #   rename(prop_site_impact=prop_ea_flown) %>%
    #   mutate(Date=as.Date(Date),
    #          Area_Name=as.character(Area_Name),
    #          effect_area_km=as.numeric(intersect_area_km)) %>%
    #   select(Area_Name,Date,Method,flight_time_hr,effect_area_km,prop_site_impact) %>%
    #   filter(Date<=as.Date(max(sysbait_effort$subper_end)))
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
    # rem_pass_site <- rem_pass_site %>%
    #   full_join(ao,by=c("Date","period","Method","Area_Name")) %>%
    #   rename(tot_hrs=flight_time_hr)
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
    # eff_pass_site <- eff_day %>%
    #   group_by(Date,period,Method,Area_Name) %>%
    #   summarise(tot_hrs=sum(tot_hrs),
    #             num_events=sum(num_events))
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
    # rem_eff_site <- eff_pass_site %>%
    #   filter(Date>=min(rem_pass_site$Date) &
    #            Date<=max(rem_pass_site$Date)) %>%
    #   full_join(rem_pass_site,by=c("Date","period","Method","Area_Name"))
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

  # ea_df <- ea %>% st_drop_geometry()

  #add ea areas and proportion of ea impacted by removal method
  if(abund_scale=="watersheds"){
    for(i in 1:nrow(rem_eff_site)){
      if(is.na(rem_eff_site$site_area_km[i])){
        rem_eff_site$site_area_km[i] <- ssg_formatted$area_km[ssg_formatted$SiteID==rem_eff_site$SiteID[i]]
      }


      if(is.na(rem_eff_site$prop_site_impact[i])){
        rem_eff_site$prop_site_impact[i] <- rem_eff_site$effect_area_km[i]/rem_eff_site$site_area_km[i]
      }
    }
  } else {
    # for(i in 1:nrow(rem_eff_site)){
    #   if(is.na(rem_eff_site$ea_area_km[i])){
    #     rem_eff_site$ea_area_km[i] <- ea_df$area_km[ea_df$Area_Name==rem_eff_site$Area_Name[i]]
    #   }
    #
    #   if(is.na(rem_eff_site$prop_site_impact[i])){
    #     rem_eff_site$prop_site_impact[i] <- rem_eff_site$effect_area_km[i]/rem_eff_site$ea_area_km[i]
    #   }
    # }
  }

  #remove observations with no area estimate
  rem_eff_site <- rem_eff_site %>% filter(!is.na(effect_area_km))

  if(pass_scale=="day"){
    #remove effort observations over 24 hours per day
    rem_eff_site<- rem_eff_site %>% filter(tot_hrs<24) %>% as.data.frame()
  }

  #add pass idx
  if(abund_scale=="watersheds"){
    rem_eff_site <- rem_eff_site %>%
      group_by(period,Method,SiteID) %>%
      mutate(pass_idx=1:n())# %>%
    # filter(tot_hrs>24)
  } else {
    # rem_eff_site <- rem_eff_site %>%
    #   group_by(period,Method,Area_Name) %>%
    #   mutate(pass_idx=1:n()) #%>%
    # # filter(tot_hrs<24)
  }

  #calculate effective hours spent removing per area
  rem_eff_site$effect_area_hrs <- rem_eff_site$tot_hrs/rem_eff_site$effect_area_km

  rem_eff_site
}

clean_aerial_tracks <- function(raw_ao_folder,
                                ssg_formatted){
  
  #read files in folder
  ao_files <- list.files(raw_ao_folder)
  # ao_files <- list.files('C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/feral_swine_elim/data/raw/aerial_ops_tracks')
  
  #function to id last character ------------------------
  substrRight <- function(x, n=1){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  #align files --------------------
  ao <- list()

  for(i in 1:length(ao_files)){
    if(grepl("csv",ao_files[i])){
      ao[[i]] <- read.csv(file.path(raw_ao_folder,ao_files[i]))
      colnames(ao[[i]]) <- gsub("\\.","_",colnames(ao[[i]]))
      ao[[i]] <- st_as_sf(ao[[i]],coords=c("Longitude_decimal_","Latitude_decimal_"),crs="epsg:4326")
    } else {
      ao[[i]] <- st_read(file.path(raw_ao_folder,ao_files[i]))
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
    ao[[i]] <- ao[[i]] %>% dplyr::select("Aircraft","Altitude_ft",
                                         "DateTime_Local","Speed_knots",
                                         "Bearing","geometry")
  }
  
  ao <- do.call("rbind",ao)
  ao <- ao[!duplicated(ao),]
  
  #filter out ferrying locations (out of state and speeds >25 knots) -----------------
  options(tigris_use_cache = TRUE)
  states_sf <- tigris::states()
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
  
  ao_coords_clean <- tibble(ao_coords_clean)
  ao_coords_clean$day <- as.Date(ao_coords_clean$DateTime_Local)
  ao_coords_clean$aircraft_day <- paste(ao_coords_clean$Aircraft,ao_coords_clean$day)
  
  # make amt tracks by aircraft day
  ao_amt <- ao_coords_clean |> make_track(.x=x,.y=y,.t=DateTime_Local,
                                          id=aircraft_day,crs=st_crs(4326))
  ao_amt <- transform_coords(ao_amt,(st_crs(ssg_formatted)))
  
  #nest by aircraft day
  ao_amt_nest <- ao_amt |> nest(data=-"id")
  
  #calculate steps and step length by aircraft day
  ao_amt_nest <- ao_amt_nest |>
    mutate(steps = map(data, function(x) x |>
                         track_resample(rate = minutes(1), tolerance = seconds(15)) |>
                         steps_by_burst()),
           sl =map(data, function(x) x |>
                     step_lengths()))
  
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
          st_as_sf(coords=c('x','y'),crs=st_crs(ssg_formatted)) %>% 
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
  
  ao_buff_day_ssg <- st_intersection(ssg_formatted %>% 
                                       st_transform(st_crs(ao_buff_day)), 
                                     ao_buff_day)
  
  intersect_ssg <- ao_buff_day_ssg %>% 
    mutate(intersect_area = st_area(.),
           intersect_area_km=intersect_area/1e6) %>%   
    dplyr::select(SiteID,Date,intersect_area_km) %>%   
    st_drop_geometry() %>% 
    left_join(ssg_formatted %>% st_drop_geometry()) %>% 
    mutate(prop_site_flown = as.numeric(intersect_area_km/area_km))
  
  ao_ssg_day <- ao_buff_day_ssg %>% st_drop_geometry() %>% 
    left_join(intersect_ssg) %>% as.data.frame()
  ao_ssg_day %>% filter(dt_hr>100)
  ao_buff %>% filter(Aircraft=="N332AB-ER" & Date==as.Date("2021-12-08"))
  
  ao_ssg_day<- ao_ssg_day %>% rename(flight_time_hr=dt_hr)
  
  ao_ssg_day

}

clean_sysbait <- function(#raw_sys_sep20_dec22,
                          #raw_sys_jan23_jun23,
                          #raw_sys_jul23_nov23,
                          #raw_sys_jul23_sep24,
                          #raw_sys_oct24_dec24,
                          #raw_systake_jan25_mar25_gdb
                          raw_sysbait_folder){
  
  sys_files <- sort(list.files(raw_sysbait_folder))
  
  sys_list <- list()

  #clean systematic baiting -------------------------------
  for(i in 1:length(sys_files)){
    if(grepl("csv",sys_files[i])){
      if(grepl("Sept2020",sys_files[i])){ #sys
        sys_list[[i]] <- read.csv(file.path(raw_sysbait_folder,sys_files[i]))
        
      } else if(grepl("Jan2023",sys_files[i])){ #sys1
        sys_list[[i]] <- read.csv(file.path(raw_sysbait_folder,sys_files[i]))
        sys_list[[i]] <- sys_list[[i]] %>% select(-c(GlobalID..,EditDate,Editor))
        
        sys <- full_join(sys_list[[i-1]],sys_list[[i]])
        
      } else if(grepl("Nov2023",sys_files[i])){ #sys 2
        sys_list[[i]] <- read.csv(file.path(raw_sysbait_folder,sys_files[i]))
        sys_list[[i]] <- sys_list[[i]] %>% 
          select(-c(GlobalID..))%>% 
          mutate(Lat=as.character(Lat),Long=as.character(Long))
        
        sys <- full_join(sys,sys_list[[i]]) 
        
        sys <- sys %>% select(Site.Status, Baiting.Start.Date, Date.of.Hot.Bait, Site.Visit.Date, Site.Visit.Reason, Trap.Start.Date,
                              Complete.Date, Lat, Long, CreationDate, Creator, Sounder.Size) %>% 
          rename(SiteStatus=Site.Status,
                 BaitingStartDate=Baiting.Start.Date,
                 DateHotBait=Date.of.Hot.Bait,
                 SiteVisitDate=Site.Visit.Date,
                 SiteVisitReason=Site.Visit.Reason,
                 TrapStartDate=Trap.Start.Date,
                 CompleteDate=Complete.Date,
                 SounderSize=Sounder.Size) %>% 
          mutate(SounderSize=as.integer(SounderSize))
        
        sys$BaitingStartDate <- as.Date(sys$BaitingStartDate,format="%m/%d/%Y %H:%M")
        sys$DateHotBait <- as.Date(sys$DateHotBait,format="%m/%d/%Y %H:%M")
        sys$CompleteDate <- as.Date(sys$CompleteDate,format="%m/%d/%Y %H:%M")
        sys$TrapStartDate <- as.Date(sys$TrapStartDate,format="%m/%d/%Y %H:%M")
        sys$SiteVisitDate <- as.Date(sys$SiteVisitDate,format="%m/%d/%Y %H:%M")
        sys$CreationDate <- as.Date(sys$CreationDate,format="%m/%d/%Y %H:%M")
        
      } else if(grepl("Sept24",sys_files[i])){ #sys 3
        sys_list[[i]] <- read.csv(file.path(raw_sysbait_folder,sys_files[i]))
        sys_list[[i]] <- sys_list[[i]] %>% rename(Lat=lat,
                                                  Long=long,
                                                  SiteVisitReason=SiteVisitReas) %>% 
          mutate(Lat=as.character(Lat),
                 Long=as.character(Long),
                 SounderSize) %>% 
          select(SiteStatus, BaitingStartDate, DateHotBait, SiteVisitDate, SiteVisitReason, TrapStartDate,
                 CompleteDate, Lat, Long, CreationDate, Creator, SounderSize)
        
        sys_list[[i]]$BaitingStartDate <-as.Date(sys_list[[i]]$BaitingStartDate,format="%Y-%m-%d %H:%M:%S")
        sys_list[[i]]$DateHotBait <- as.Date(sys_list[[i]]$DateHotBait,format="%Y-%m-%d %H:%M:%S")
        sys_list[[i]]$CompleteDate <- as.Date(sys_list[[i]]$CompleteDate,format="%Y-%m-%d %H:%M:%S")
        sys_list[[i]]$TrapStartDate <- as.Date(sys_list[[i]]$TrapStartDate,format="%Y-%m-%d %H:%M:%S")
        sys_list[[i]]$SiteVisitDate <- as.Date(sys_list[[i]]$SiteVisitDate,format="%Y-%m-%d %H:%M:%S")
        sys_list[[i]]$CreationDate <- as.Date(sys_list[[i]]$CreationDate,format="%Y-%m-%d %H:%M:%S")
        
        sys <- full_join(sys,sys_list[[i]]) 
        
      } else if(grepl("Oct24",sys_files[i])){ #sys 4
        sys_list[[i]] <- read.csv(file.path(raw_sysbait_folder,sys_files[i]))
        
        sys_list[[i]] <- sys_list[[i]] %>% select(SiteStatus, BaitingStartDate, DateHotBait, SiteVisitDate, SiteVisitReason, TrapStartDate,
                                                  CompleteDate, Lat, Long, CreationDate, Creator, SounderSize) %>% 
          mutate(Lat=as.character(Lat),
                 Long=as.character(Long))
        sys_list[[i]]$BaitingStartDate <-as.Date(sys_list[[i]]$BaitingStartDate,format="%Y-%m-%d")
        sys_list[[i]]$DateHotBait <- as.Date(sys_list[[i]]$DateHotBait,format="%Y-%m-%d")
        sys_list[[i]]$CompleteDate <- as.Date(sys_list[[i]]$CompleteDate,format="%Y-%m-%d")
        sys_list[[i]]$TrapStartDate <- as.Date(sys_list[[i]]$TrapStartDate,format="%Y-%m-%d")
        sys_list[[i]]$SiteVisitDate <- as.Date(sys_list[[i]]$SiteVisitDate,format="%Y-%m-%d")
        sys_list[[i]]$CreationDate <- as.Date(sys_list[[i]]$CreationDate,format="%Y-%m-%d")
        
        sys <- full_join(sys,sys_list[[i]]) 
        
      }
    } else if(grepl("gdb",sys_files[i])){
      if(grepl("Q1CY2025",sys_files[i])){ #sys5
        sys_list[[i]] <- st_read(file.path(raw_sysbait_folder,sys_files[i]),
                                 layer="Q12025_FeralHogOps")
        
        sys_list[[i]] <- sys_list[[i]] %>% st_transform("epsg:4326")
        sys_coords <- st_coordinates(sys_list[[i]])
        sys_list[[i]]$Lat <- sys_coords[,2]
        sys_list[[i]]$Long <- sys_coords[,1]
        sys_list[[i]] <- sys_list[[i]] %>% st_drop_geometry()
        
        sys_list[[i]] <- sys_list[[i]] %>% 
          rename(SiteVisitReason=SiteVisitReas) %>% 
          select(SiteStatus, BaitingStartDate, DateHotBait, 
                 SiteVisitDate, SiteVisitReason, TrapStartDate,
                 CompleteDate, Lat, Long, CreationDate, Creator, SounderSize) %>% 
          mutate(Lat=as.character(Lat),
                 Long=as.character(Long))
        
        sys_list[[i]]$BaitingStartDate <-as.Date(sys_list[[i]]$BaitingStartDate,format="%Y-%m-%d")
        sys_list[[i]]$DateHotBait <- as.Date(sys_list[[i]]$DateHotBait,format="%Y-%m-%d")
        sys_list[[i]]$CompleteDate <- as.Date(sys_list[[i]]$CompleteDate,format="%Y-%m-%d")
        sys_list[[i]]$TrapStartDate <- as.Date(sys_list[[i]]$TrapStartDate,format="%Y-%m-%d")
        sys_list[[i]]$SiteVisitDate <- as.Date(sys_list[[i]]$SiteVisitDate,format="%Y-%m-%d")
        sys_list[[i]]$CreationDate <- as.Date(sys_list[[i]]$CreationDate,format="%Y-%m-%d")
        
        sys <- full_join(sys,sys_list[[i]])
      } else { #new gdbs
        layer_name <- st_layers(file.path(raw_sysbait_folder,sys_files[i]))
        layer_name <- layer_name[which(grepl("ops",layer_name) |
                                         grepl("Ops",layer_name) |
                                         grepl("OPS",layer_name))]
        sys_list[[i]] <- st_read(file.path(raw_sysbait_folder,sys_files[i]),
                                 layer=layer_name)
        
        sys_list[[i]] <- sys_list[[i]] %>% st_transform("epsg:4326")
        sys_coords <- st_coordinates(sys_list[[i]])
        sys_list[[i]]$Lat <- sys_coords[,2]
        sys_list[[i]]$Long <- sys_coords[,1]
        sys_list[[i]] <- sys_list[[i]] %>% st_drop_geometry()
        
        sys_list[[i]] <- sys_list[[i]] %>% 
          rename(SiteVisitReason=SiteVisitReas) %>% 
          select(SiteStatus, BaitingStartDate, DateHotBait, 
                 SiteVisitDate, SiteVisitReason, TrapStartDate,
                 CompleteDate, Lat, Long, CreationDate, Creator, SounderSize) %>% 
          mutate(Lat=as.character(Lat),
                 Long=as.character(Long))
        
        sys_list[[i]]$BaitingStartDate <-as.Date(sys_list[[i]]$BaitingStartDate,format="%Y-%m-%d")
        sys_list[[i]]$DateHotBait <- as.Date(sys_list[[i]]$DateHotBait,format="%Y-%m-%d")
        sys_list[[i]]$CompleteDate <- as.Date(sys_list[[i]]$CompleteDate,format="%Y-%m-%d")
        sys_list[[i]]$TrapStartDate <- as.Date(sys_list[[i]]$TrapStartDate,format="%Y-%m-%d")
        sys_list[[i]]$SiteVisitDate <- as.Date(sys_list[[i]]$SiteVisitDate,format="%Y-%m-%d")
        sys_list[[i]]$CreationDate <- as.Date(sys_list[[i]]$CreationDate,format="%Y-%m-%d")
        
        sys <- full_join(sys,sys_list[[i]])
      }
    } else {
      stop("New data not csv or gdb, edit clean_sysbait() target to read new files.")
    }
  }#end file loop

  sys
}

clean_take <- function(#raw_take_sep20_dec22,
                       #raw_take_jan23_jun23,
                       #raw_take_jul23_nov23,
                       #raw_take_jul23_sep24,
                       #raw_take_oct24_dec24,
                       #raw_systake_jan25_mar25_gdb,
                       raw_take_folder,
                       ssg_formatted){
  
  #character substituion function
  charsub <- function(x){gsub("\\.","_",x)}

  take_files <- sort(list.files(raw_take_folder))

  take_list <- list()

  #clean take-------------------------------
  for(i in 1:length(take_files)){
    if(grepl("csv",take_files[i])){
      if(grepl("Sept2020",take_files[i])){
        take_list[[i]] <- read.csv(file.path(raw_take_folder,take_files[i]))
        take_list[[i]] <- as.data.frame(take_list[[i]] %>% select(-c(OBJECTID..,Shape..,GlobalID..)))
        take_list[[i]] <- take_list[[i]][,!grepl("UAV",colnames(take_list[[i]]))]
        take_list[[i]] <- take_list[[i]][,!grepl("Minutes",colnames(take_list[[i]]))]
        take_list[[i]] <- as.data.frame(take_list[[i]] %>%
                                          rename(Date=Take.Date,
                                                 Total=Total.Take,
                                                 Adult_M=Adult.Male.Take,
                                                 Adult_F=Adult.Female.Take,
                                                 Method=Method.Type,
                                                 NonBreed_F=Non.Breeding.Female.Take,
                                                 NonBreed_M=Non.Breeding.Male.Take))
        take_list[[i]]$Males <- take_list[[i]]$Adult_M+take_list[[i]]$NonBreed_M
        take_list[[i]]$Females <- take_list[[i]]$Adult_F+take_list[[i]]$NonBreed_F
        take_list[[i]]$Adults <- take_list[[i]]$Adult_M+take_list[[i]]$Adult_F
        take_list[[i]]$Sub_adults <- take_list[[i]]$NonBreed_F+take_list[[i]]$NonBreed_M #sub-adults = sub-adults and juveniles
        take_list[[i]] <- take_list[[i]] %>%
          mutate(Trapper.Name=toupper(Trapper.Name),
                 created_user=toupper(created_user),
                 last_edited_user=toupper(last_edited_user),
                 Comments=toupper(Comments))
        take_list[[i]]$Date <- as.POSIXct(take_list[[i]]$Date,format="%m/%d/%Y %H:%M",tz="CST6CDT")
        
      } else if(grepl("Jan2023",take_files[i]) ){
        take_list[[i]] <- read.csv(file.path(raw_take_folder,take_files[i]))
        take_list[[i]] <- as.data.frame(take_list[[i]] %>% select(-c(OBJECTID..,Shape..,GlobalID..)))
        take_list[[i]] <- take_list[[i]] %>% select(-grep("UAV",colnames(take_list[[i]])))
        take_list[[i]] <- as.data.frame(take_list[[i]] %>% rename(Method=Method.Type,
                                                                  Total=Total.Take,
                                                                  Date=Take.Date,
                                                                  Adult_M=Adult.Male.Take,
                                                                  Adult_F=Adult.Female.Take,
                                                                  NonBreed_M=Non.Breeding.Male.Take,
                                                                  NonBreed_F=Non.Breeding.Female.Take))
        take_list[[i]]$Total <- as.numeric(take_list[[i]]$Total)
        take_list[[i]]$Adult_M <- as.numeric(take_list[[i]]$Adult_M)
        take_list[[i]]$Adult_F <- as.numeric(take_list[[i]]$Adult_F)
        take_list[[i]]$NonBreed_M <- as.numeric(take_list[[i]]$NonBreed_M)
        take_list[[i]]$NonBreed_F <- as.numeric(take_list[[i]]$NonBreed_F)
        
        take_list[[i]]$Males <- take_list[[i]]$Adult_M+take_list[[i]]$NonBreed_M
        take_list[[i]]$Females <- take_list[[i]]$Adult_F+take_list[[i]]$NonBreed_F
        take_list[[i]]$Adults <- take_list[[i]]$Adult_M+take_list[[i]]$Adult_F
        take_list[[i]]$Sub_adults <- take_list[[i]]$NonBreed_F+take_list[[i]]$NonBreed_M #sub-adults = sub-adults and juveniles
        take_list[[i]] <- take_list[[i]] %>%
          select(-Minutes.of.Effort) %>%
          mutate(Trapper.Name=toupper(Trapper.Name),
                 created_user=toupper(created_user),
                 last_edited_user=toupper(last_edited_user),
                 Comments=toupper(Comments))
        take_list[[i]]$Date <- as.POSIXct(take_list[[i]]$Date,format="%m/%d/%Y %H:%M",tz="CST6CDT")
        
        take <- full_join(take_list[[i-1]],take_list[[i]])
        colnames(take) <- unlist(lapply(colnames(take),charsub))
        
      } else if(grepl("Nov2023",take_files[i])){
        take_list[[i]] <- read.csv(file.path(raw_take_folder,take_files[i]))
        colnames(take_list[[i]]) <- unlist(lapply(colnames(take_list[[i]]),charsub))
        charsub_ <- function(x){gsub("__","_",x)}
        colnames(take_list[[i]]) <- unlist(lapply(colnames(take_list[[i]]),charsub_))
        
        take_list[[i]] <- take_list[[i]] %>%
          rename(Date=Take_Date,
                 Total=Total_Take,
                 Adult_M=Adult_Male_Take,
                 Adult_F=Adult_Female_Take,
                 Method=Method_Type,
                 NonBreed_F=Non_Breeding_Female_Take,
                 NonBreed_M=Non_Breeding_Male_Take) %>%
          select(-c(GlobalID_,UAV_Used,Minutes_of_UAV_Flight,
                    Primary_Shooter_UAV_,Secondary_Shooter_UAV_,Pilot_Name_UAV_,Minutes_of_Effort,
                    OBJECTID_,Shape_)) %>%
          mutate(Trapper_Name=toupper(Trapper_Name),
                 Comments=toupper(Comments))
        
        take_list[[i]]$Total <- as.numeric(take_list[[i]]$Total)
        take_list[[i]]$Adult_M <- as.numeric(take_list[[i]]$Adult_M)
        take_list[[i]]$Adult_F <- as.numeric(take_list[[i]]$Adult_F)
        take_list[[i]]$NonBreed_M <- as.numeric(take_list[[i]]$NonBreed_M)
        take_list[[i]]$NonBreed_F <- as.numeric(take_list[[i]]$NonBreed_F)
        
        take_list[[i]]$Males <- take_list[[i]]$Adult_M+take_list[[i]]$NonBreed_M
        take_list[[i]]$Females <- take_list[[i]]$Adult_F+take_list[[i]]$NonBreed_F
        take_list[[i]]$Adults <- take_list[[i]]$Adult_M+take_list[[i]]$Adult_F
        take_list[[i]]$Sub_adults <- take_list[[i]]$NonBreed_F+take_list[[i]]$NonBreed_M #sub-adults = sub-adults and juveniles
        
        take_list[[i]]$Date <- as.POSIXct(take_list[[i]]$Date,format="%m/%d/%Y %H:%M",tz="CST6CDT")
        # colnames(take_list[[i]])[!colnames(take_list[[i]])%in%colnames(take)]
        
        take <- full_join(take,take_list[[i]])
        
      } else if(grepl("Sept24",take_files[i])){
        take_list[[i]] <- read.csv(file.path(raw_take_folder,take_files[i]))
        take_list[[i]] <- take_list[[i]] %>% rename(Date=Take_Date,
                                                    Total=Total_Take,
                                                    Adult_M=Adult_Male_Take,
                                                    Adult_F=Adult_Female_Take,
                                                    NonBreed_F=NonBreeding_Female_Take,
                                                    NonBreed_M=NonBreeding_Male_Take,
                                                    Long=long,
                                                    Lat=lat,
                                                    Feral_Domestic_Pot_belly=FeralDomestic_Potbelly,
                                                    Method=Method_Type) %>%
          select(-c(X,GlobalID,UAV_Used,Minutes_UAV_Flight,
                    Primary_Shooter_UAV,Secondary_Shooter_UAV,Pilot_Name_UAV,Minutes_of_Effort)) %>%
          mutate(Trapper_Name=toupper(Trapper_Name),
                 created_user=toupper(created_user),
                 last_edited_user=toupper(last_edited_user),
                 Comments=toupper(Comments))
        take_list[[i]]$Date <- as.POSIXct(take_list[[i]]$Date,format="%Y-%m-%d %H:%M:%S",tz="CST6CDT")
        
        take <- full_join(take,take_list[[i]])
        
      } else if(grepl("Oct24",take_files[i])){
        take_list[[i]] <- st_read(file.path(raw_take_folder,take_files[i]),
                                  layer="FeralHogTake_Oct24_Dec24")
        take_list[[i]] <- take_list[[i]] %>%
          rename(Trapper_Name=Trapper_Name,
                 Method=Method_Type,
                 Date=Take_Date,
                 Total=Total_Take,
                 Adult_M=Adult_Male_Take,
                 Adult_F=Adult_Female_Take,
                 NonBreed_M=NonBreeding_Male_Take ,
                 NonBreed_F=NonBreeding_Female_Take ,
                 Feral_Domestic_Pot_belly=FeralDomestic_Potbelly) %>%
          st_transform(crs="epsg:4326") %>%
          mutate(created_date=as.character(created_date),
                 last_edited_date=as.character(last_edited_date))
        
        take_coords <- st_coordinates(take_list[[i]])
        take_list[[i]]$Lat <- take_coords[,2]
        take_list[[i]]$Long <- take_coords[,1]
        take_list[[i]]$Date <- as.POSIXct(take_list[[i]]$Date,format="%Y-%m-%d %H:%M:%S",tz="CST6CDT")
        
        take_list[[i]] <- take_list[[i]] %>%
          select(colnames(take_list[[i]])[colnames(take_list[[i]])%in%colnames(take_list[[i-1]])])
        
        take <- full_join(take,take_list[[i]])
        
      }
    } else if(grepl("gdb",take_files[i])){
      if(grepl("Q1CY2025",take_files[i])){
        take_list[[i]] <- st_read(file.path(raw_take_folder,take_files[i]),
                                  layer="Q12025_FeralHogTake")
        take_list[[i]] <- take_list[[i]] %>%
          rename(Trapper_Name=Trapper_Name,
                 Method=Method_Type,
                 Date=Take_Date,
                 Total=Total_Take,
                 Adult_M=Adult_Male_Take,
                 Adult_F=Adult_Female_Take,
                 NonBreed_M=NonBreeding_Male_Take ,
                 NonBreed_F=NonBreeding_Female_Take ,
                 Feral_Domestic_Pot_belly=FeralDomestic_Potbelly) %>%
          st_transform(crs="epsg:4326") %>%
          mutate(created_date=as.character(created_date),
                 last_edited_date=as.character(last_edited_date))
        
        take_coords <- st_coordinates(take_list[[i]])
        take_list[[i]]$Lat <- take_coords[,2]
        take_list[[i]]$Long <- take_coords[,1]
        take_list[[i]]$Date <- as.POSIXct(take_list[[i]]$Date,format="%Y-%m-%d %H:%M:%S",tz="CST6CDT")
        
        take_list[[i]] <- take_list[[i]] %>%
          select(colnames(take_list[[i]])[colnames(take_list[[i]])%in%colnames(take_list[[i-1]])])
        
        take <- full_join(take,take_list[[i]])
        
      } else { #other gdb
        layer_name <- st_layers(file.path(raw_take_folder,take_files[i]))
        layer_name <- layer_name[which(grepl("take",layer_name) |
                                         grepl("Take",layer_name) |
                                         grepl("TAKE",layer_name))]
        
        take_list[[i]] <- st_read(file.path(raw_take_folder,take_files[i]),
                                  layer=layer_name)
        take_list[[i]] <- take_list[[i]] %>%
          rename(Trapper_Name=Trapper_Name,
                 Method=Method_Type,
                 Date=Take_Date,
                 Total=Total_Take,
                 Adult_M=Adult_Male_Take,
                 Adult_F=Adult_Female_Take,
                 NonBreed_M=NonBreeding_Male_Take ,
                 NonBreed_F=NonBreeding_Female_Take ,
                 Feral_Domestic_Pot_belly=FeralDomestic_Potbelly) %>%
          st_transform(crs="epsg:4326") %>%
          mutate(created_date=as.character(created_date),
                 last_edited_date=as.character(last_edited_date))
        
        take_coords <- st_coordinates(take_list[[i]])
        take_list[[i]]$Lat <- take_coords[,2]
        take_list[[i]]$Long <- take_coords[,1]
        take_list[[i]]$Date <- as.POSIXct(take_list[[i]]$Date,format="%Y-%m-%d %H:%M:%S",tz="CST6CDT")
        
        take_list[[i]] <- take_list[[i]] %>%
          select(colnames(take_list[[i]])[colnames(take_list[[i]])%in%colnames(take_list[[i-1]])])
        
        take <- full_join(take,take_list[[i]])
      }
    }else {
      stop("New data not csv or gdb, edit clean_take() target to read new files.")
    }
  }


  take$Method[grep("Aerial",take$Method)] <- "Aerial"
  take$Method[grep("Shoot",take$Method)] <- "Ground Shoot"
  take <- subset(take,Method!="Detection Dogs" & Method!="Snares")
  take <- take[!duplicated(take),]

  #remove wrong timestamps
  take$date_notime <- as.Date(take$Date)
  take <- take %>%
    filter(year(date_notime)>=2020 & year(date_notime)<=year(Sys.Date())) %>%
    filter(Method!="<Null>")

  #remove points outside of study area
  ssg_formatted <- ssg_formatted %>% st_union()

  take$id <- 1:nrow(take)
  take_sf <- subset(take,!is.na(Lat))
  take_sf <- st_as_sf(take_sf,coords=c("Long","Lat"),crs="epsg:4326")
  take_sf <- st_intersection(take_sf,st_transform(ssg_formatted,st_crs(take_sf)))
  take <- take[which(take$id %in% take_sf$id),]
  take <- take %>% select(-Shape)
  take
  
}

extract_nlcd <- function(ssg_formatted,
                         nlcd_raster){
  
  nlcd <- rast(nlcd_raster)
  
  occ_scale <- "watersheds"
  study_site <- st_union(ssg_formatted)
  
  study_site_v <- st_transform(study_site,crs(nlcd))
  study_site_v <- vect(study_site_v)
  nlcd.crop <- crop(nlcd,study_site_v)
  
  nlcd.crop <- project(nlcd.crop,"epsg:3601")
  nlcd.crop[is.nan(nlcd.crop)] <- 0
  
  nlcd.water <- nlcd.crop
  values(nlcd.water)[values(nlcd.water)==11] <- 1
  values(nlcd.water)[values(nlcd.water)!=1] <- 0
  
  ssg_formatted <- st_transform(ssg_formatted,"epsg:3601")

  water <- terra::extract(nlcd.water,ssg_formatted,fun=mean,na.rm=T) %>% 
    rename(SiteID=ID,water=Layer_1) %>% 
    right_join(ssg_formatted,by="SiteID") %>% 
    st_as_sf()# %>% select(-c(neighbors)) 
  
  nlcd.developed <- nlcd.crop
  values(nlcd.developed)[values(nlcd.developed)%in%c(21:24)] <- 1
  values(nlcd.developed)[values(nlcd.developed)!=1] <- 0
  
  developed <- terra::extract(nlcd.developed,
                              ssg_formatted,fun=mean) %>% 
    rename(SiteID=ID,developed=Layer_1) %>% 
    right_join(ssg_formatted,by="SiteID") %>% 
    st_as_sf()# %>% select(-c(neighbors)) 
  # select(-c(neighbors,FID,SITE,geometry))
 
  nlcd.forest <- nlcd.crop
  values(nlcd.forest)[values(nlcd.forest)%in%c(41:43)] <- 1
  values(nlcd.forest)[values(nlcd.forest)!=1] <- 0
  
  forest <- terra::extract(nlcd.forest,ssg_formatted,fun=mean) %>%
    rename(SiteID=ID,forest=Layer_1) %>%
    right_join(ssg_formatted,by="SiteID") %>%
    st_as_sf() #%>% select(-c(neighbors)) 
  
  nlcd.wetland <- nlcd.crop
  values(nlcd.wetland)[values(nlcd.wetland)%in%c(90,95)] <- 1
  values(nlcd.wetland)[values(nlcd.wetland)!=1] <- 0
  
  wetland <- terra::extract(nlcd.wetland,ssg_formatted,fun=mean) %>% 
    rename(SiteID=ID,wetland=Layer_1) %>% 
    right_join(ssg_formatted,by="SiteID") %>% 
    # select(-c(neighbors,FID,SITE,geometry))
    st_as_sf() #%>% select(-c(neighbors)) 
 
  nlcd.agri <- nlcd.crop
  values(nlcd.agri)[values(nlcd.agri)%in%c(81,82)] <- 1
  values(nlcd.agri)[values(nlcd.agri)!=1] <- 0
  
  agri <- terra::extract(nlcd.agri,ssg_formatted,fun=mean) %>% 
    rename(SiteID=ID,agri=Layer_1) %>% 
    right_join(ssg_formatted,by="SiteID") %>% 
    st_as_sf() #%>% select(-c(neighbors)) 
  
  nlcd_siteid <- full_join(water %>% st_drop_geometry(),
                           developed%>% st_drop_geometry()) %>% 
    full_join(agri%>% st_drop_geometry()) %>% 
    full_join(wetland%>% st_drop_geometry()) %>% 
    full_join(forest%>% st_drop_geometry()) %>%
    full_join(ssg_formatted) %>% 
    st_drop_geometry()#%>% dplyr::select(-neighbors)
  
  nlcd_siteid$site_idx <- as.numeric(as.factor(nlcd_siteid$SiteID))
  
  nlcd_siteid
}

extract_oak <- function(ssg_formatted,
                        oak_raster){
  
  oak_mo <- rast(oak_raster)
  abund_scale <- "watersheds"
  
  if(abund_scale=="watersheds"){
    
    oak_ssg <- exact_extract(oak_mo, ssg_formatted, coverage_area = TRUE)
    
    oak_siteid <- data.frame(oak_km =unlist(lapply(1:length(oak_ssg),function(i){
      sum(oak_ssg[[i]]$coverage_area[oak_ssg[[i]]$value==31])/1e6
    })))
    
    oak_siteid$SiteID <- ssg_formatted$SiteID
    oak_siteid <- oak_siteid %>% 
      left_join(ssg_formatted %>% st_drop_geometry() %>% select(SiteID,area_km))
    
    oak_siteid$prop_oak <- oak_siteid$oak_km/oak_siteid$area_km
  }  
  
  oak_siteid$prop_oak_sc <- scale(oak_siteid$prop_oak)
  
  oak_siteid
}

occ_data_prep <- function(sysbait_effort,
                          ssg_formatted){
 
  sysbait_effort$site_idx_occ <- as.numeric(sysbait_effort$site_idx)
  
  sysbait_effort <- sysbait_effort %>% 
    rename(period_idx=period) %>% 
    left_join(ssg_formatted %>% 
                st_drop_geometry()) %>% 
    group_by(period_idx,site_idx_occ) %>% 
    summarise(detections=sum(detection),
              nweeks = length(unique(subper_start)),
              per_start=min(subper_start),
              trap_nights_km = sum(trap_nights_km),
              elim_area = max(as.numeric(Area_Name))) %>% 
    arrange(period_idx,site_idx_occ)
  
  sysbait_effort$trap_nights_km_sc <- scale(sysbait_effort$trap_nights_km)
  
  sysbait_effort 
}

# take_effort <- tar_read(take_effort)
# dat_occ <- tar_read(dat_occ)
rem_data_prep <- function(take_effort,
                          ssg_formatted,
                          dat_occ){

  ##filter effort outliers ------------
  eff_cutoff <- quantile(take_effort$effect_area_hrs,prob=0.95)
  dat_rem <- take_effort %>% 
    filter(effect_area_hrs<eff_cutoff) %>% 
    rename(period_idx = period,
           method=Method) 

  dat_rem$site_idx_rem <- as.numeric(as.factor(dat_rem$SiteID))
  dat_rem$site_idx_occ <- as.numeric(dat_rem$SiteID)
  
  ##add 0 aerial data for later periods (only occurs sesaonally) -----
  last_per_meth <- dat_rem %>% 
    filter(period_idx==max(dat_rem$period_idx)) %>% 
    ungroup() %>% 
    reframe(method=unique(method))
  
  if(!"Aerial"%in%last_per_meth$method){
    a_pers <- dat_rem %>% ungroup() %>%
      filter(method=="Aerial") %>%
      select(period_idx) %>% 
      distinct
    a_pers_n <- which(!(1:max(dat_rem$period_idx)%in%a_pers$period_idx)) #no aerial ops during these periods

    sids <- sort(unique(dat_rem$site_idx_rem))
    
    site_idx_lookup <- dat_rem %>% 
      ungroup() %>% 
      select(SiteID,site_idx_occ,site_idx_rem) %>% 
      distinct() %>% 
      arrange(site_idx_rem) %>% 
      drop_na
    
    # samp <- dat_rem[1:(length(a_pers_n)*length(sids)),] #use existing data frame structure and clear fields
    samp <- data.frame(
      subper_id =1,
      subper_start =NA,
      subper_end = NA,
      period_idx = rep(a_pers_n,length(sids)),
      method= "Aerial",
      tot_hrs =0, 
      num_events=0,
      tot_rem =0,
      ea_area_km=NA,
      effect_area_km = 0,
      prop_site_impact=0,
      pass_idx = 1,
      effect_area_hrs =0,
      site_idx_rem = sort(rep(sids,length(a_pers_n)))
    )
    samp <- samp %>% 
      left_join(site_idx_lookup)%>% 
      left_join(ssg_formatted %>% 
                  st_drop_geometry())
    
    # samp$site_idx_rem <- sort(rep(unique(sids),length(a_pers_n)))
    dat_rem <- dat_rem %>% 
      full_join(samp)
  }
  
  #fill in missing SiteIDs from added aerial 0's
  # dat_rem$SiteID[is.na(dat_rem$SiteID)] <- 
  #   sapply(1:length(dat_rem$SiteID[is.na(dat_rem$SiteID)]),function(x){
  #     site_idx_lookup$SiteID[site_idx_lookup$site_idx_rem==
  #                              dat_rem$site_idx_rem[is.na(dat_rem$SiteID)][x]]
  #   })
  
  dat_rem
}

align_sites <- function(dat_rem){
  dat_rem %>% 
    ungroup() %>% 
    select(SiteID,site_idx_occ,site_idx_rem) %>% 
    distinct() %>% 
    arrange(site_idx_rem) %>% 
    drop_na
}

site_lat_longs <- function(ssg_formatted){
  
  lat_long_siteid <- ssg_formatted
  lat_long_siteid$lat <- st_coordinates(st_centroid(lat_long_siteid))[,2]
  lat_long_siteid$long <- st_coordinates(st_centroid(lat_long_siteid))[,1]
  lat_long_siteid$lat_sc <- scale(lat_long_siteid$lat)
  lat_long_siteid$long_sc <- scale(lat_long_siteid$long)
  
  lat_long_siteid
}

nimble_fun_avail <- function(){
  avail_fun <- nimbleFunction(
    run = function(pass = double(0), #removal pass number
                   # type = double(0), #removal type index
                   gamma = double(0), #area of removal impact/area of site (time t)
                   gamma_past = double(1), #area of removal impact/area of site (time t-1)
                   theta = double(0), #capture probability (time t)
                   theta_past = double(1)) { #capture probability (time)
      returnType(double(0))
      if(pass!=1) {
        return(gamma*theta*prod((1-gamma_past)+gamma_past*(1-theta_past)))
      } else {
        return(gamma*theta)
      }
    })
  avail_fun
}

# model_inputs_list <- function(dat_occ,
#                               dat_rem, 
#                               avail_fun,
#                               ssg_formatted, 
#                               oak_siteid,
#                               nlcd_siteid,
#                               lat_long_siteid,
#                               site_idx_lookup,
#                               chain_idx){
#   
#   nChains <- max(chain_idx)
#   data_list <- list(dat_occ=dat_occ,
#                     dat_rem=dat_rem, 
#                     avail_fun=avail_fun,
#                     ssg_formatted=ssg_formatted, 
#                     oak_siteid=oak_siteid,
#                     nlcd_siteid=nlcd_siteid,
#                     lat_long_siteid=lat_long_siteid,
#                     site_idx_lookup=site_idx_lookup)
#   
#   full_data_list <- rep(list(data_list),nChains)
#   names(full_data_list)<-rep("model_inputs",nChains)
#   full_data_list
# }
make_chain_idx <- function(nChains){
  seq(1:nChains)
}

# model_data_prep <- function(nlcd_siteid,
#                             lat_long_siteid,
#                             oak_siteid,
#                             site_idx_lookup,
#                             dat_occ,
#                             dat_rem,
#                             eff_weeks){
#   
#   #covariates
#   develop <- data.frame(develop=nlcd_siteid$developed)
#   develop$develop_sc <- scale(develop$develop)[,1]
#   
#   agri <- data.frame(agri=nlcd_siteid$agri)
#   agri$agri_sc <- scale(agri$agri)[,1]
#   
#   lat_long_siteid <- lat_long_siteid %>% filter(SiteID%in%site_idx_lookup$SiteID)
#   oak_siteid <- oak_siteid %>% filter(SiteID%in%site_idx_lookup$SiteID)
#   rem_covs <- lat_long_siteid %>% left_join(oak_siteid)
#   
#   ##indices for nimble ---------------------------------
#   nsites_occ <- length(unique(dat_occ$site_idx_occ))
#   nsites_rem <- length(unique(dat_rem$site_idx_rem))
#   nperiods <- length(unique(dat_occ$period_idx))
#   
#   #split by removal type 
#   dat_aerial <- dat_rem %>% filter(method=="Aerial")
#   dat_trap <- dat_rem %>% filter(method=="Trap")
#   
#   gamma_a <- array(0,dim=c(max(dat_aerial$pass_idx),nsites_rem,nperiods))
#   gamma_t <- array(0,dim=c(max(dat_trap$pass_idx),nsites_rem,nperiods))
#   for(i in 1:nsites_rem){
#     for(t in 1:nperiods){
#       for(j in 1:max(dat_aerial$pass_idx,dat_trap$pass_idx)){
#         x <- dat_aerial$prop_site_impact[dat_aerial$pass_idx==j &
#                                            dat_aerial$site_idx_rem==i &
#                                            dat_aerial$period_idx==t]
#         z <- dat_trap$prop_site_impact[dat_trap$pass_idx==j &
#                                          dat_trap$site_idx_rem==i &
#                                          dat_trap$period_idx==t]
#         if(length(x)>0){
#           gamma_a[j,i,t] <- max(x) 
#         }
#         if(length(z)>0){
#           gamma_t[j,i,t] <- max(z)
#         }
#       }
#     }
#   }
#   
#   dat_aerial$effect_area_hrs_sc <- scale(dat_aerial$effect_area_hrs)
#   dat_trap$effect_area_hrs_sc <- scale(dat_trap$effect_area_hrs)
#   
#   ##removal matrix ------------------
#   remtot <- dat_rem %>% group_by(site_idx_rem,period_idx) %>%
#     reframe(removal=sum(tot_rem))
#   
#   yrem <- matrix(0,nsites_rem,nperiods)
#   site_idx <- sort(unique(dat_rem$site_idx_rem))
#   for(j in 1:nsites_rem){
#     for(t in 1:nperiods){
#       x <- remtot$removal[remtot$site_idx_rem==site_idx[j] & remtot$period_idx==t]
#       if(length(x)>0){
#         yrem[j,t] <- x
#       }}}
#   
#   #mean trapping effort  -------------------------------
#   mn_te <- (((0.5/2.59)*45)-attr(dat_occ$trap_nights_km_sc,"scaled:center"))/
#     attr(dat_occ$trap_nights_km_sc,"scaled:scale") 
#   
#   mod_covs <- 
#     list(rem_covs=rem_covs,
#         nsites_rem=nsites_rem,
#         nsites_occ=nsites_occ,
#         nperiods=nperiods,
#         dat_aerial=dat_aerial,
#         dat_trap=dat_trap,
#         gamma_a=gamma_a,
#         gamma_t=gamma_t,
#         yrem=yrem,
#         mn_te=mn_te,
#         eff_weeks=eff_weeks)
#   mod_covs
# }
# library(nimble)
fit_zibinom_mod <- function(dat_occ,
                            dat_rem,
                            avail_fun,
                            ssg_formatted,
                            site_idx_lookup,
                            nlcd_siteid,
                            lat_long_siteid,
                            oak_siteid,
                            # eff_weeks,
                            niter,
                            thin,
                            burnProp,
                            chain_idx){
  
  #covariates
  develop <- data.frame(develop=nlcd_siteid$developed)
  develop$develop_sc <- scale(develop$develop)[,1]
  
  agri <- data.frame(agri=nlcd_siteid$agri)
  agri$agri_sc <- scale(agri$agri)[,1]
  
  lat_long_siteid <- lat_long_siteid %>% filter(SiteID%in%site_idx_lookup$SiteID)
  oak_siteid <- oak_siteid %>% filter(SiteID%in%site_idx_lookup$SiteID)
  rem_covs <- lat_long_siteid %>% left_join(oak_siteid)
  
  ##indices for nimble ---------------------------------
  nsites_occ <- length(unique(dat_occ$site_idx_occ))
  nsites_rem <- length(unique(dat_rem$site_idx_rem))
  nperiods <- length(unique(dat_occ$period_idx))
  
  #split by removal type 
  dat_aerial <- dat_rem %>% filter(method=="Aerial")
  dat_trap <- dat_rem %>% filter(method=="Trap")
  
  gamma_a <- array(0,dim=c(max(dat_aerial$pass_idx),nsites_rem,nperiods))
  gamma_t <- array(0,dim=c(max(dat_trap$pass_idx),nsites_rem,nperiods))
  for(i in 1:nsites_rem){
    for(t in 1:nperiods){
      for(j in 1:max(dat_aerial$pass_idx,dat_trap$pass_idx)){
        x <- dat_aerial$prop_site_impact[dat_aerial$pass_idx==j &
                                           dat_aerial$site_idx_rem==i &
                                           dat_aerial$period_idx==t]
        z <- dat_trap$prop_site_impact[dat_trap$pass_idx==j &
                                         dat_trap$site_idx_rem==i &
                                         dat_trap$period_idx==t]
        if(length(x)>0){
          gamma_a[j,i,t] <- max(x) 
        }
        if(length(z)>0){
          gamma_t[j,i,t] <- max(z)
        }
      }
    }
  }
  
  dat_aerial$effect_area_hrs_sc <- scale(dat_aerial$effect_area_hrs)
  dat_trap$effect_area_hrs_sc <- scale(dat_trap$effect_area_hrs)
  
  ##removal matrix ------------------
  remtot <- dat_rem %>% group_by(site_idx_rem,period_idx) %>%
    reframe(removal=sum(tot_rem))
  
  yrem <- matrix(0,nsites_rem,nperiods)
  site_idx <- sort(unique(dat_rem$site_idx_rem))
  for(j in 1:nsites_rem){
    for(t in 1:nperiods){
      x <- remtot$removal[remtot$site_idx_rem==site_idx[j] & remtot$period_idx==t]
      if(length(x)>0){
        yrem[j,t] <- x
      }}}
  
  #mean trapping effort  -------------------------------
  # mn_te <- (((0.5/2.59)*45)-attr(dat_occ$trap_nights_km_sc,"scaled:center"))/
  #   attr(dat_occ$trap_nights_km_sc,"scaled:scale") 

  ## model specification ------------------------
  ZIbinomcode <- nimbleCode({
    #occupancy regression coefficients
    for(i in 1:nbeta){
      beta[i] ~ dlogis(0,1)
    }
    
    #lambda regression coefficients
    for(i in 1:nbeta_lam){
      beta_lam[i] ~ dnorm(0,5)
    }
    
    for(t in 1:nperiods){
      beta0[t] ~ dnorm(beta[1],sd=sd_beta0)
    }
    sd_beta0 ~ dgamma(10,10)
    
    for(i in 1:2){
      alpha[i] ~ dlogis(0,1)
      delta_t[i] ~ dlogis(0,1)
      delta_a[i] ~ dlogis(0,1)
    }
    
    sd_pdet ~ dgamma(10,10)
    sd_theta_t ~ dgamma(10,10)
    sd_theta_a ~ dgamma(10,10)
    sd_lam ~ dgamma(1,20)
    # lambda ~ dnorm(1,0.1)
    
    p_n1 ~ dbeta(0.2,12)
    r_n1 ~ dgamma(6,1) 
    
    
    #mean effort detection probability
    # mu_pdet_mn <- alpha[1] + alpha[2]*mn_te
    # logit(p_sys) ~ dnorm(mu_pdet_mn,sd=sd_pdet)
    
    
    for(i in 1:nsites_rem){
      N[i,1] ~ dnegbin(size=r_n1,prob=p_n1)
      
      for(t in 1:nperiods){
        N_latent[i,t] <- N[i,t] * z_site[site_idx_occ_rem[i],t]
      }
      
      #abundace/removal process
      mu_lam[i]<- beta_lam[1] + beta_lam[2]*oak[i] + 
        beta_lam[3]*y[i] + beta_lam[4]*x[i]
      log(lambda[i]) ~ dnorm(mu_lam[i],sd=sd_lam)
      
      for(t in 2:nperiods){
        
        mu_intens[i,t] <- lambda[i] * (N[i,t-1] - yrem_mat[i,t-1])
        N[i,t] ~ dpois(mu_intens[i,t])
      }
    }
    
    for(k in 1:nsites_occ){ #watershed loop
      for(t in 1:nperiods){
        logit(psi[k,t]) <- beta0[t] + beta[2]*develop[k] + beta[3]*agri[k]
        
        #occupancy/detection
        pabs[k,t] <- 1-psi[k,t]
        z_site[k,t] ~ dbern(psi[k,t])
        
        #effective detection given effort
        # pstar_site[k,t] <- 1-(1-p_sys)^eff_weeks
        
        #p(elimination) given effort
        # pelim[k,t] <- pabs[k,t]/(pabs[k,t] + psi[k,t]*(1-pstar_site[k,t]))
      }
    }
    
    #abundance
    ##aerial
    for(i in 1:nsamp_aerial){
      mu_theta_a[i] <- delta_a[1] + delta_a[2]*eff_area_a[i]
      logit(theta_a[pass_idx_a[i],site_idx_a[i],period_idx_a[i]]) ~ 
        dnorm(mu_theta_a[i],sd=sd_theta_a)
      
      pip_a[i] <-
        avail_fun(pass=pass_idx_a[i],
                  gamma=gamma_a[pass_idx_a[i],site_idx_a[i],period_idx_a[i]],
                  gamma_past=gamma_a[1:(pass_idx_a[i]-1),site_idx_a[i],period_idx_a[i]],
                  theta=theta_a[pass_idx_a[i],site_idx_a[i],period_idx_a[i]],
                  theta_past=theta_a[1:(pass_idx_a[i]-1),site_idx_a[i],period_idx_a[i]])
      
      yrem_aerial[i] ~ dbinom(size=N_latent[site_idx_a[i],period_idx_a[i]], prob=pip_a[i])
      yrem_pred_aerial[i] ~ dbinom(size=N_latent[site_idx_a[i],period_idx_a[i]], prob=pip_a[i])
    }
    
    ##traps
    for(i in 1:nsamp_trap){
      #removal probability
      mu_theta_t[i] <- delta_t[1] + delta_t[2]*eff_area_t[i]
      logit(theta_t[pass_idx_t[i],site_idx_t[i],period_idx_t[i]]) ~
        dnorm(mu_theta_t[i],sd=sd_theta_t)
      
      pip_t[i] <-
        avail_fun(pass=pass_idx_t[i],
                  gamma=gamma_t[pass_idx_t[i],site_idx_t[i],period_idx_t[i]],
                  gamma_past=gamma_t[1:(pass_idx_t[i]-1),site_idx_t[i],period_idx_t[i]],
                  theta=theta_t[pass_idx_t[i],site_idx_t[i],period_idx_t[i]],
                  theta_past=theta_t[1:(pass_idx_t[i]-1),site_idx_t[i],period_idx_t[i]])
      
      #abundance likelihood
      yrem_trap[i] ~ dbinom(size=N_latent[site_idx_t[i],period_idx_t[i]], prob=pip_t[i])
      yrem_pred_trap[i] ~ dbinom(size=N_latent[site_idx_t[i],period_idx_t[i]], prob=pip_t[i])
    }
    
    #occupancy
    for(i in 1:nsamp_occ){
      #detection probability
      mu_pdet[i] <- alpha[1] + alpha[2]*trap_nights_km[i]
      logit(pdet[i]) ~ dnorm(mu_pdet[i],sd=sd_pdet)
      
      #occupancy likelihood
      pocc[i] <- pdet[i]*z_site[site_idx_occ[i],period_idx_occ[i]]
      yocc[i] ~ dbinom(prob=pocc[i],size=nweeks[i])
      yocc_pred[i] ~ dbinom(prob=pocc[i],size=nweeks[i])
    }
  })#nimblecode end
  
  nbeta <- 3
  nbeta_lam <- 4
  # nsites_occ <- mod_covs$nsites_occ
  # nsites_rem <- mod_covs$nsites_rem
  # nperiods <- mod_covs$nperiods
  # 
  # dat_trap <- mod_covs$dat_trap
  # dat_aerial <- mod_covs$dat_aerial
  # 
  # gamma_t <- mod_covs$gamma_t
  # gamma_a <- mod_covs$gamma_a
  # 
  # yrem <- mod_covs$yrem
  # rem_covs <- mod_covs$rem_covs
  # develop <- mod_covs$develop
  # agri <- mod_covs$agri
  # 
  # eff_weeks <- mod_covs$eff_weeks
  # mn_te <- mod_covs$mn_te
  
  ## nimble lists  -------------------------------------------------
  modDat <- list(yrem_trap=dat_trap$tot_rem,
                 eff_area_t=dat_trap$effect_area_hrs_sc[,1],
                 gamma_t=gamma_t,
                 yrem_aerial=dat_aerial$tot_rem,
                 eff_area_a=dat_aerial$effect_area_hrs_sc[,1],
                 gamma_a=gamma_a,
                 yocc=dat_occ$detections,
                 nweeks=dat_occ$nweeks,
                 trap_nights_km=dat_occ$trap_nights_km_sc[,1],
                 yrem_mat=yrem,
                 oak=rem_covs$prop_oak_sc[,1],
                 y=rem_covs$lat_sc[,1],
                 x=rem_covs$long_sc[,1]
  )
  ### constants -----------------------------
  const <- list(develop=develop$develop_sc,
                agri=agri$agri_sc,
                nbeta=nbeta,
                nbeta_lam=nbeta_lam,
                nperiods=nperiods,
                nsamp_occ=nrow(dat_occ),
                nsites_occ=nsites_occ,
                site_idx_occ=dat_occ$site_idx_occ,
                site_idx_occ_rem=site_idx_lookup$site_idx_occ, 
                period_idx_occ=dat_occ$period_idx,
                nsites_rem=nsites_rem,
                # mn_te=mn_te,
                # eff_weeks=eff_weeks,
                site_idx_t=dat_trap$site_idx_rem,
                site_idx_a=dat_aerial$site_idx_rem,
                period_idx_t=dat_trap$period_idx,
                period_idx_a=dat_aerial$period_idx,
                pass_idx_t=dat_trap$pass_idx,
                pass_idx_a=dat_aerial$pass_idx,
                nsamp_trap=nrow(dat_trap),
                nsamp_aerial=nrow(dat_aerial)
  )

  ### initial values -----------------------------
  inits <- list(
    beta=rnorm(nbeta,c(0,0.5),0.1),
    beta_lam=rnorm(nbeta_lam,c(-0.05,0.01,-0.02,0),0.005),
    sd_lam=rnorm(1,0.005,0.0001),
    N = matrix(NA,nsites_rem,nperiods),
    z_site = matrix(1,nsites_occ,nperiods),
    p_n1 = rnorm(1,0.005,0.0001),
    r_n1 = rnorm(1,6,0.1),
    theta_a = array(0,dim=dim(gamma_a)),
    theta_t = array(0,dim=dim(gamma_t))
  )
  
  inits$N[,1] <- apply(yrem,1,max)*20
  
  ## nimble configuration -------------------------------------------------
  mod <- nimbleModel(code = ZIbinomcode,
                     data = modDat,
                     constants = const,
                     inits = inits)
  
  mcmc_conf <- configureMCMC(mod,enableWAIC=F)
  
  monitors <- c("lambda",
                "beta",
                "beta0",
                "beta_lam",
                "alpha",
                "sd_pdet",
                "delta_t",
                "delta_a",
                "N",
                "pabs",
                "pocc",
                "p_n1",
                "r_n1",
                "N_latent",
                "pip_a",
                "pip_t",
                "yrem_pred_aerial",
                "yrem_pred_trap",
                "yocc_pred"
  )

  mcmc_conf$setMonitors(monitors)
  
  Cmcmc <- buildMCMC(mcmc_conf,resetFuncitions=T) #uncompiled MCMC
  Cmod <- compileNimble(mod,Cmcmc) #compiled model

  samples <- runMCMC(Cmod$Cmcmc,
                     niter = niter,
                     thin= thin,
                     nburn = niter*burnProp,
                     setSeed=chain_idx,
                     nchains = 1)
  
  samples <- as.data.frame(samples)
  samples$chain_idx <- chain_idx
  samples
  
}

# fit_zibinom_mod <- function(Cmod,
#                             niter,
#                             thin,
#                             burnProp,
#                             chain_idx){
#   
# 
# }

model_validation <- function(dat_occ,
                             dat_rem,
                             samples,
                             out_dir){
  
  if(class(samples)=="list"){
    nChains <- length(samples)
    nmcmc <- nrow(samples[[1]])
  } else {
    nChains <- 1
    nmcmc <- nrow(samples)
  }
  
  #split by removal type 
  dat_aerial <- dat_rem %>% filter(method=="Aerial")
  dat_trap <- dat_rem %>% filter(method=="Trap")
  
  #occupancy ----------------------
  if(nChains>1){
    psi <- 1- do.call("rbind",lapply(1:nChains,function(i){
      samples[[i]][,grep("pabs",colnames(samples[[i]]))]}))
  } else {
    psi <- 1- samples[,grep("pabs",colnames(samples))]
  }
  psi_mn <- colMeans(psi)
  
  psi_md <- sapply(1:ncol(psi),function(i)quantile(psi[,i],prob=0.5))
  
  psi_sum <- data.frame(idx=names(psi_mn),mn=psi_mn)
  psi_sum$site_idx <- sapply(1:nrow(psi_sum),function(x){
    as.numeric(gsub("\\[","",gsub("pabs","",unlist(str_split(psi_sum[x,"idx"],","))[1])))})
  
  psi_sum$period_idx <- sapply(1:nrow(psi_sum),function(x){
    as.numeric(gsub("\\]","",unlist(str_split((psi_sum[x,"idx"]),","))[2]))})
  
  ## AUC for prevalance process ----------------
  if(nChains>1){
    pocc <- do.call("rbind",lapply(1:nChains,function(i){
      samples[[i]][grep("pocc",colnames(samples[[i]]))]}))
  } else {
    pocc <- samples[,grep("pocc",colnames(samples))]
  }
  pocc_mn <- colMeans(pocc)
  pocc_md <- sapply(1:ncol(pocc),function(i)quantile(pocc[,i],prob=0.5))
  
  pocc_sum <- data.frame(idx=names(pocc_mn),mn=pocc_mn,md=pocc_md)
  
  ypred_vec <- rbinom(nrow(dat_occ),dat_occ$nweeks,pocc_sum$mn)
  # ypred_vec_bin <- ifelse(ypred_vec>0,1,0)
  ybin <- ifelse(dat_occ$detections>0,1,0)
  # auc_prev_bin <- as.numeric(pROC::roc(ybin,ypred_vec_bin)$auc)
  auc_prev_pocc <- as.numeric(pROC::roc(ybin,pocc_sum$mn)$auc)
  
  ## AUC for occupancy process ---------------------
  psi_mn_vec <- numeric(nrow(dat_occ))
  for(i in 1:nrow(dat_occ)){
    psi_mn_vec[i] <- psi_sum$mn[psi_sum$site_idx==dat_occ$site_idx_occ[i] & 
                                  psi_sum$period_idx==dat_occ$period_idx[i]]
  }
  auc_occ <- suppressMessages(as.numeric(pROC::roc(ybin,psi_mn_vec)$auc))
  # auc_prev_bin
  
  #bayesian p-value ----------------------
  pocc <- as.matrix(pocc)
  occ_ll <- dev_ypred_occ <- matrix(NA,nrow(dat_occ),nmcmc*nChains)
  
  if(nChains>1){
    ypred_occ <- as.matrix(do.call("rbind",lapply(1:nChains,function(i){samples[[i]][grep("yocc_pred",colnames(samples[[i]]))]})))
  } else {
    ypred_occ <- samples[,grep("yocc_pred",colnames(samples))]
  }
  
  for(i in 1:nrow(dat_occ)){
    occ_ll[i,] <- dbinom(dat_occ$detections[i],dat_occ$nweeks[i],pocc[,i],log=T)
    dev_ypred_occ[i,] <- dbinom(ypred_occ[,i],dat_occ$nweeks[i],pocc[,i],log=T)
  }
  
  if(nChains>1){
    N_latent <- as.matrix(do.call("rbind",lapply(1:nChains,function(i){
      samples[[i]][,grepl("N_latent",colnames(samples[[i]]))]})))
    pip_a <- as.matrix(do.call("rbind",lapply(1:nChains,function(i){
      samples[[i]][,grepl("pip_a",colnames(samples[[i]]))]})))
    ypred_rem_a <- as.matrix(do.call("rbind",lapply(1:nChains,function(i){samples[[i]][,grep("yrem_pred_aerial",colnames(samples[[i]]))]})))
    pip_t <- as.matrix(do.call("rbind",lapply(1:nChains,function(i){
      samples[[i]][,grepl("pip_t",colnames(samples[[i]]))]})))
    ypred_rem_t <- as.matrix(do.call("rbind",lapply(1:nChains,function(i){samples[[i]][,grep("yrem_pred_trap",colnames(samples[[i]]))]})))
    pip_g <- as.matrix(do.call("rbind",lapply(1:nChains,function(i){
      samples[[i]][,grepl("pip_g",colnames(samples[[i]]))]})))
    ypred_rem_g <- as.matrix(do.call("rbind",lapply(1:nChains,function(i){samples[[i]][,grep("yrem_pred_ground",colnames(samples[[i]]))]})))
  } else {
    N_latent <- samples[,grepl("N_latent",colnames(samples))]
    pip_a <- samples[,grepl("pip_a",colnames(samples))]
    ypred_rem_a <- samples[,grepl("yrem_pred_aerial",colnames(samples))]
    pip_t <- samples[,grepl("pip_t",colnames(samples))]
    ypred_rem_t <- samples[,grepl("yrem_pred_trap",colnames(samples))]
    pip_g <- samples[,grepl("pip_g",colnames(samples))]
    ypred_rem_g <- samples[,grepl("yrem_pred_ground",colnames(samples))]
  }
  
  samp_N_a <- rem_ll_a <-dev_ypred_rem_a <- matrix(NA,nrow(dat_aerial),nmcmc*nChains)
  samp_N_t <- rem_ll_t <-dev_ypred_rem_t <- matrix(NA,nrow(dat_trap),nmcmc*nChains)
  for(i in 1:nrow(dat_aerial)){
    samp_N_a[i,] <- N_latent[,paste0("N_latent[",dat_aerial$site_idx_rem[i],", ",
                                     dat_aerial$period_idx[i],"]")]
    rem_ll_a[i,] <- dbinom(x=dat_aerial$tot_rem[i],size=samp_N_a[i,],
                           prob=pip_a[,i],log=T)
    dev_ypred_rem_a[i,] <- dbinom(x=ypred_rem_a[,i],size=samp_N_a[i,],
                                  prob=pip_a[,i],log=T)
  }
  for(i in 1:nrow(dat_trap)){
    samp_N_t[i,] <- N_latent[,paste0("N_latent[",dat_trap$site_idx_rem[i],", ",
                                     dat_trap$period_idx[i],"]")]
    rem_ll_t[i,] <- dbinom(x=dat_trap$tot_rem[i],size=samp_N_t[i,],
                           prob=pip_t[,i],log=T)
    dev_ypred_rem_t[i,] <- dbinom(x=ypred_rem_t[,i],size=samp_N_t[i,],
                                  prob=pip_t[,i],log=T)
  }
  
  rem_ll <- rbind(rem_ll_a,rem_ll_t)
  dev_ypred_rem <- rbind(dev_ypred_rem_a,dev_ypred_rem_t)
  
  dev_y <- colSums(occ_ll) + colSums(rem_ll,na.rm=T)
  dev_ypred <- colSums(dev_ypred_occ) + colSums(dev_ypred_rem,na.rm = T)
  
  pVal_ll <- sum(dev_ypred>dev_y)/(nmcmc*nChains)
  
  pVal_occ <- sum(colSums(dev_ypred_occ)>colSums(occ_ll))/(nmcmc*nChains)
  pVal_rem <- sum(colSums(dev_ypred_rem,na.rm=T)>colSums(rem_ll,na.rm=T))/(nmcmc*nChains)
  
  ggplot()+
    geom_histogram(aes(x=colSums(dev_ypred_occ),fill="predictions"),alpha=0.5)+
    geom_histogram(aes(x=colSums(occ_ll),fill="data"),alpha=0.5)+
    xlab('loglikelihood - occupancy')+
    ggtitle(paste("pVal =",round(pVal_occ,2)))
  ggsave(filename = file.path(out_dir,"validation","pval_occ_ll.jpeg"),
         width=7,height=5,units="in",device="jpeg")

  ggplot()+
    geom_histogram(aes(x=colSums(dev_ypred_rem,na.rm=T),fill="predictions"),alpha=0.5)+
    geom_histogram(aes(x=colSums(rem_ll,na.rm=T),fill="data"),alpha=0.5)+
    xlab('loglikelihood - abundance')+
    ggtitle(paste("pVal =",round(pVal_rem,2)))
  ggsave(filename = file.path(out_dir,"validation","pval_rem_ll.jpeg"),
         width=7,height=5,units="in",device="jpeg")
  
  #mean squared error -----------------------------
  ypred_occ <- as.matrix(ypred_occ)
  mse_occ <- mse_pred_occ <- rep(NA,nmcmc*nChains)
  for(i in 1:(nmcmc*nChains)){
    mse_occ[i] <- mean((dat_occ$detections - dat_occ$nweeks*pocc[i,])^2)
    mse_pred_occ[i] <- mean((ypred_occ[i,] - dat_occ$nweeks*pocc[i,])^2)
  }
  
  ypred_rem_a <- as.matrix(ypred_rem_a)
  ypred_rem_t <- as.matrix(ypred_rem_t)
  N_latent <- as.matrix(N_latent)
  mse_rem_a <- mse_pred_rem_a <- matrix(NA,nrow(dat_aerial),nmcmc*nChains)
  mse_rem_t <- mse_pred_rem_t <- matrix(NA,nrow(dat_trap),nmcmc*nChains)
  for(i in 1:nrow(dat_aerial)){
    exp_N <- samp_N_a[i,]*pip_a[,i]
    mse_rem_a[i,] <- (dat_aerial$tot_rem[i] - exp_N)^2
    mse_pred_rem_a[i,] <- (ypred_rem_a[,i] - exp_N)^2
  }
  for(i in 1:nrow(dat_trap)){
    exp_N <- samp_N_t[i,]*pip_t[,i]
    mse_rem_t[i,] <- (dat_trap$tot_rem[i] - exp_N)^2
    mse_pred_rem_t[i,] <- (ypred_rem_t[,i] - exp_N)^2
  }
  
  mse_rem <- colMeans(rbind(mse_rem_a,mse_rem_t),na.rm=T)
  mse_pred_rem <- colMeans(rbind(mse_pred_rem_a,mse_pred_rem_t),na.rm=T)
  
  pVal_occ_mse <- sum(mse_pred_occ>mse_occ)/(nmcmc*nChains)
  pVal_rem_mse <- sum(mse_pred_rem>mse_rem)/(nmcmc*nChains)
  
  ggplot()+
    geom_histogram(aes(x=mse_pred_occ,fill="predictions"),alpha=0.5)+
    geom_histogram(aes(x=mse_occ,fill="data"),alpha=0.5)+
    xlab('mse - occupancy')+
    ggtitle(paste("pVal =",round(pVal_occ_mse,2)))
  ggsave(filename = file.path(out_dir,"validation","pval_occ_mse.jpeg"),
         width=7,height=5,units="in",device="jpeg")

  ggplot()+
    geom_histogram(aes(x=mse_pred_rem,fill="predictions"),alpha=0.5)+
    geom_histogram(aes(x=mse_rem,fill="data"),alpha=0.5)+
    xlab('mse - abundance')+
    ggtitle(paste("pVal =",round(pVal_rem_mse,2)))
  ggsave(filename = file.path(out_dir,"validation","pval_rem_mse.jpeg"),
         width=7,height=5,units="in",device="jpeg")
  
  
  #p-value mean and var stats -------------------
  ypred_rem_mn <- rowMeans(cbind(ypred_rem_a,ypred_rem_t))
  dat_rem_mn <- mean(c(dat_aerial$tot_rem,dat_trap$tot_rem))
  pVal_rem_mn <- sum(ypred_rem_mn>dat_rem_mn)/(nmcmc*nChains)
  ggplot()+geom_histogram(aes(x=ypred_rem_mn))+
    geom_vline(xintercept=dat_rem_mn,col="red",lwd=1.5)+
    xlab("mean - removal")+
    ggtitle(paste("pVal =",round(pVal_rem_mn,2)))
  ggsave(filename = file.path(out_dir,"validation","pval_rem_mn.jpeg"),
         width=7,height=5,units="in",device="jpeg")

  ypred_occ_mn <- rowMeans(ypred_occ)
  dat_occ_mn <- mean(dat_occ$detections)
  pVal_occ_mn <- sum(ypred_occ_mn>dat_occ_mn)/(nmcmc*nChains)
  ggplot()+geom_histogram(aes(x=ypred_occ_mn))+
    geom_vline(xintercept=dat_occ_mn,col="red",lwd=1.5)+
    xlab("mean - occupancy")+
    ggtitle(paste("pVal =",round(pVal_occ_mn,2)))
  ggsave(filename = file.path(out_dir,"validation","pval_occ_mn.jpeg"),
         width=7,height=5,units="in",device="jpeg")

  ## standard deviation ------------------
  ypred_rem_var <- sqrt(apply(cbind(ypred_rem_a,ypred_rem_t),1,var))
  dat_rem_var <- sqrt(var(c(dat_aerial$tot_rem,dat_trap$tot_rem)))
  pVal_rem_var <- sum(ypred_rem_var>dat_rem_var)/(nmcmc*nChains)
  ggplot()+geom_histogram(aes(x=ypred_rem_var))+
    geom_vline(xintercept=dat_rem_var,col="red",lwd=1.5)+
    xlab("sd - removal")+
    ggtitle(paste("pVal =",round(pVal_rem_var,2)))
  ggsave(filename = file.path(out_dir,"validation","pval_rem_sd.jpeg"),
         width=7,height=5,units="in",device="jpeg")

  ypred_occ_var <- sqrt(apply(ypred_occ,1,var))
  dat_occ_var <- sqrt(var(dat_occ$detections))
  pVal_occ_var <- sum(ypred_occ_var>dat_occ_var)/(nmcmc*nChains)
  ggplot()+geom_histogram(aes(x=ypred_occ_var))+
    geom_vline(xintercept=dat_occ_var,col="red",lwd=1.5)+
    xlab("sd - occupancy")+
    ggtitle(paste("pVal =",round(pVal_occ_var,2)))
  ggsave(filename = file.path(out_dir,"validation","pval_occ_sd.jpeg"),
         width=7,height=5,units="in",device="jpeg")
  
  val_tab <- data.frame("model"=c("joint",
                       rep(c("occ","rem"),4),
                       "occ_det",
                       "occ_occ"),
             "validation_meth"=c(rep("pVal",9),
                                 rep("AUC",2)),
             "stat" =c(rep("log_liklihood",3),
                       rep("mean",2),
                       rep("var",2),
                       rep("mse",2),
                       rep(NA,2)),
             "value"=c(  pVal_ll,
                         pVal_occ,
                         pVal_rem,
                         
                         pVal_occ_mn,
                         pVal_rem_mn,
                         
                         pVal_occ_var,
                         pVal_rem_var,
                        
                         pVal_occ_mse,
                         pVal_rem_mse,
                         
                         auc_prev_pocc,
                         auc_occ))
  
  write.csv(val_tab,file.path(out_dir,"validation","val_table.csv"))
  
  val_tab
  
}

# extract_nfsp <- function(ssg_formatted){
#   
#   mo <- tigris::states()
#   mo <- mo %>% filter(NAME=="Missouri")
#   yrs<-15:23
#   
#   #crop to study site
#   ssg <- ssg_formatted %>% mutate(full_area=st_area(.)) %>% st_drop_geometry()
#   i<-length(yrs)
#   ssg_nfsp_l <- list()
#   for(i in 1:length(yrs)){
#     x <- st_read(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/NFSP Mapping/Missouri/fy",yrs[i],"_mo.shp"))
#     x <- st_transform(x,st_crs(ssg_formatted))
#     samp <- st_intersection(x,ssg_formatted) %>% mutate(intersect_area=st_area(.)) %>% 
#       group_by(SiteID) %>% summarise(intersect_area=sum(intersect_area))
#     
#     y <- ssg %>% left_join(samp) %>% mutate(prop_nfsp=intersect_area/full_area) 
#     z <- y %>% select(-geometry) %>% left_join(ssg_formatted)
#     z$prop_nfsp[is.na(z$prop_nfsp)] <- 0
#     z$prop_nfsp <- as.numeric(z$prop_nfsp)
#     ssg_nfsp_l[[i]] <- z
#   }
#   
#   # q15 <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri Data/Model Ready Data/NFSP Watershed Overlap/watersheds_nfsp_overlap_fy15.shp")
#   # q16 <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri Data/Model Ready Data/NFSP Watershed Overlap/watersheds_nfsp_overlap_fy16.shp")
#   # q17 <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri Data/Model Ready Data/NFSP Watershed Overlap/watersheds_nfsp_overlap_fy17.shp")
#   # q18 <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri Data/Model Ready Data/NFSP Watershed Overlap/watersheds_nfsp_overlap_fy18.shp")
#   # q19 <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri Data/Model Ready Data/NFSP Watershed Overlap/watersheds_nfsp_overlap_fy19.shp")
#   # q20 <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri Data/Model Ready Data/NFSP Watershed Overlap/watersheds_nfsp_overlap_fy20.shp")
#   # q21 <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri Data/Model Ready Data/NFSP Watershed Overlap/watersheds_nfsp_overlap_fy21.shp")
#   # q22 <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri Data/Model Ready Data/NFSP Watershed Overlap/watersheds_nfsp_overlap_fy22.shp")
#   # q23 <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri Data/Model Ready Data/NFSP Watershed Overlap/watersheds_nfsp_overlap_fy23.shp")
#   
#   ssg_nfsp <- ssg_nfsp_l[[1]]
#   for(i in 2:length(ssg_nfsp_l)){
#     ssg_nfsp <- ssg_nfsp %>% 
#       st_drop_geometry() %>% 
#       full_join(ssg_nfsp_l[[i]] %>% 
#                   st_drop_geometry()) 
#   }
#   ssg_nfsp$prp_nfs[is.na(ssg_nfsp$prp_nfs)]<- 0
#   
#   ssg_nfsp
# }
