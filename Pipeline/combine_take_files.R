###Title: Combining multiple take data sheets into one
###Author: Abbey Feuka
###Date: 19042023
###Notes: creates MO_ops_take_cleaned.csv and MO_ops_trapping_clean.csv
#######################
library(tidyverse)
library(lubridate)
library(sf)

setwd("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Data/Take")

##############
# #ops
# ops1 <- read.csv("Ops_Oct2015_May2019.csv")
# colnames(ops1) <- tolower(colnames(ops1))
# ops1 <- ops1[,!grepl("object",colnames(ops1))]
# 
# #ops1
# colnames(ops1)[colnames(ops1)=="shape.."] <- "shape"
# colnames(ops1)[colnames(ops1)=="date"] <- "field_date"
# colnames(ops1)[colnames(ops1)=="creationda"] <- "date_entered"
# colnames(ops1)[colnames(ops1)=="editdate"] <- "date_edited"
# colnames(ops1)[colnames(ops1)=="location_n"] <- "location_info"
# colnames(ops1)[colnames(ops1)=="methods"] <- "trap_status"
# # ops1$coords <- paste(ops1$lat,ops1$long)
# ops1 <- as.data.frame(ops1 %>% select(-c(globalid,shape)))
# 
# #ops2
# ops2 <- read.csv("Ops_May2019_Sept2020.csv")
# colnames(ops2) <- tolower(colnames(ops2))
# ops2 <- ops2[,!grepl("object",colnames(ops2))]
# 
# colnames(ops2)[colnames(ops2)=="shape.."] <- "shape"
# colnames(ops2)[colnames(ops2)=="created_da"] <- "date_entered"
# colnames(ops2)[colnames(ops2)=="last_edite"] <- "date_edited"
# colnames(ops2)[colnames(ops2)=="date_remov"] <- "date_removed"
# colnames(ops2)[colnames(ops2)=="date_preba"] <- "date_prebait"
# colnames(ops2)[colnames(ops2)=="date_trap_"] <- "date_trapset"
# colnames(ops2)[colnames(ops2)=="trap_statu"] <- "trap_status"
# colnames(ops2)[colnames(ops2)=="location_i"] <- "location_info"
# ops2 <- ops2[,!grepl("1",colnames(ops2))]
# ops2 <- as.data.frame(ops2 %>% select(-c(globalid,shape,created_us)))
# 
# ops <- full_join(ops1,ops2,by=colnames(ops2)[which(colnames(ops2)%in%colnames(ops1))])
# 
# #ops3
# ops3 <- read.csv("Ops_Sept2020_Dec2022.csv")
# colnames(ops3) <- tolower(colnames(ops3))
# ops3 <- ops3[,!grepl("object",colnames(ops3))]
# 
# colnames(ops3)[colnames(ops3)=="shape.."] <- "shape"
# colnames(ops3)[colnames(ops3)=="site.status"] <- "site_status"
# colnames(ops3)[colnames(ops3)=="baiting.start.date"] <- "date_prebait"
# colnames(ops3)[colnames(ops3)=="date.of.hot.bait"] <- "date_hotbait" 
# colnames(ops3)[colnames(ops3)=="site.visit.date"] <- "date_sitevisit"
# colnames(ops3)[colnames(ops3)=="site.visit.reason"] <- "reason_sitevisit"
# colnames(ops3)[colnames(ops3)=="trap.start.date"] <- "date_trapset"
# colnames(ops3)[colnames(ops3)=="complete.date"] <- "date_completed"
# colnames(ops3)[colnames(ops3)=="responsible.trapper"] <- "employee"
# colnames(ops3)[colnames(ops3)=="comments"] <- "location_info"
# ops3$lat <- as.numeric(ops3$lat)
# ops3$long <- as.numeric(ops3$long)
# ops3$site_status <- NA #all marked as "complete"
# ops3 <- as.data.frame(ops3 %>% select(-c(date_sitevisit,reason_sitevisit)) %>% 
#   pivot_longer(cols=c("date_prebait","date_hotbait","date_trapset","date_completed"),
#                names_to="trap_status",values_to="field_date"))
# ops3$trap_status[ops3$trap_status=="date_prebait"] <- "Pre_Baiting"
# ops3$trap_status[ops3$trap_status=="date_trapset"] <- "Trap_Set"
# ops3$trap_status[ops3$trap_status=="date_hotbait"] <- "Hot_Baiting"
# ops3$trap_status[ops3$trap_status=="date_completed"] <- "Trap_Inactive"
# 
# ops <- full_join(ops,ops3,by=colnames(ops3)[which(colnames(ops3)%in%colnames(ops))])
# 
# ops <- ops %>% select(-c(site_status,shape))
# 
# #remove points outside of missouri
# ops$id <- seq(1:nrow(ops))
# ops.sf <- subset(ops,!is.na(lat))
# ops.sf <- st_as_sf(ops.sf,coords=c("long","lat"),crs="epsg:4326")
# mo <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri Data/missouri.shp")
# mo <- st_transform(mo,st_crs(take_sf))
# ops.sf <- st_intersection(ops.sf,mo)
# ops <- ops[which(ops.sf$id %in% ops$id),]
# 
# #remove wrong timestamps
# ops$field_date <- as.Date(ops$field_date,format="%m/%d/%Y")
# range(ops$field_date,na.rm=T)
# # write.csv(ops,"./MO_ops_trapping_clean.csv")

##############
#take
# take1 <- read.csv("Take_Oct2015_Sept2021.csv")
# take1 <- as.data.frame(take1 %>% select(-c(OBJECTID..,Shape..,POINT_X,POINT_Y)))
# take1 <- as.data.frame(take1 %>% 
#                          rename(Method=Methods,
#                                 Total=Take))
# take1$Sub_adults <- take1$Sub_adults +take1$Juveniles  #sub-adults = sub-adults and juveniles

take2 <- read.csv("Take_Sept2020_Dec2022.csv")
take2 <- as.data.frame(take2 %>% select(-c(OBJECTID..,Shape..,GlobalID..)))
take2 <- take2[,!grepl("UAV",colnames(take2))]
take2 <- take2[,!grepl("Minutes",colnames(take2))]
take2 <- as.data.frame(take2 %>% 
                         rename(Date=Take.Date,
                                Total=Total.Take,
                                Adult_M=Adult.Male.Take,
                                Adult_F=Adult.Female.Take,
                                Method=Method.Type,
                                NonBreed_F=Non.Breeding.Female.Take,
                                NonBreed_M=Non.Breeding.Male.Take))
take2$Males <- take2$Adult_M+take2$NonBreed_M
take2$Females <- take2$Adult_F+take2$NonBreed_F
take2$Adults <- take2$Adult_M+take2$Adult_F
take2$Sub_adults <- take2$NonBreed_F+take2$NonBreed_M #sub-adults = sub-adults and juveniles
take2 <- take2 %>% 
  mutate(Trapper.Name=toupper(Trapper.Name),
         created_user=toupper(created_user),
         last_edited_user=toupper(last_edited_user),
         Comments=toupper(Comments))
take2$Date <- as.POSIXct(take2$Date,format="%m/%d/%Y %H:%M",tz="CST6CDT")

take3 <- read.csv("Take_Jan2023_Jun2023.csv")
take3 <- as.data.frame(take3 %>% select(-c(OBJECTID..,Shape..,GlobalID..)))
take3 <- take3 %>% select(-grep("UAV",colnames(take3)))
take3 <- as.data.frame(take3 %>% rename(Method=Method.Type,
                                Total=Total.Take,
                                Date=Take.Date,
                                Adult_M=Adult.Male.Take,
                                Adult_F=Adult.Female.Take,
                                NonBreed_M=Non.Breeding.Male.Take,
                                NonBreed_F=Non.Breeding.Female.Take))
take3$Total <- as.numeric(take3$Total)
take3$Adult_M <- as.numeric(take3$Adult_M)
take3$Adult_F <- as.numeric(take3$Adult_F)
take3$NonBreed_M <- as.numeric(take3$NonBreed_M)
take3$NonBreed_F <- as.numeric(take3$NonBreed_F)

take3$Males <- take3$Adult_M+take3$NonBreed_M
take3$Females <- take3$Adult_F+take3$NonBreed_F
take3$Adults <- take3$Adult_M+take3$Adult_F
take3$Sub_adults <- take3$NonBreed_F+take3$NonBreed_M #sub-adults = sub-adults and juveniles
take3 <- take3 %>% 
  select(-Minutes.of.Effort) %>% 
  mutate(Trapper.Name=toupper(Trapper.Name),
         created_user=toupper(created_user),
         last_edited_user=toupper(last_edited_user),
         Comments=toupper(Comments))
take3$Date <- as.POSIXct(take3$Date,format="%m/%d/%Y %H:%M",tz="CST6CDT")

take <- full_join(take2,take3)

charsub <- function(x){gsub("\\.","_",x)}
colnames(take) <- unlist(lapply(colnames(take),charsub))

take4 <- read.csv("./Take_Jul2023_Nov2023.csv")
colnames(take4) <- unlist(lapply(colnames(take4),charsub))
charsub_ <- function(x){gsub("__","_",x)}
colnames(take4) <- unlist(lapply(colnames(take4),charsub_))

take4 <- take4 %>% 
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

take4$Total <- as.numeric(take4$Total)
take4$Adult_M <- as.numeric(take4$Adult_M)
take4$Adult_F <- as.numeric(take4$Adult_F)
take4$NonBreed_M <- as.numeric(take4$NonBreed_M)
take4$NonBreed_F <- as.numeric(take4$NonBreed_F)

take4$Males <- take4$Adult_M+take4$NonBreed_M
take4$Females <- take4$Adult_F+take4$NonBreed_F
take4$Adults <- take4$Adult_M+take4$Adult_F
take4$Sub_adults <- take4$NonBreed_F+take4$NonBreed_M #sub-adults = sub-adults and juveniles

take4$Date <- as.POSIXct(take4$Date,format="%m/%d/%Y %H:%M",tz="CST6CDT")
# colnames(take4)[!colnames(take4)%in%colnames(take)]

take <- full_join(take,take4)

take5 <- read.csv("Take_July23_Sept24.csv")
take5 <- take5 %>% rename(Date=Take_Date,
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
take5$Date <- as.POSIXct(take5$Date,format="%Y-%m-%d %H:%M:%S",tz="CST6CDT")

take <- full_join(take,take5)

st_layers('C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Data/Original MO Data/MOFeralHogData_Oct24_Dec24.gdb')
take6 <- st_read('C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Data/Original MO Data/MOFeralHogData_Oct24_Dec24.gdb',
                 layer="FeralHogTake_Oct24_Dec24")
take6 <- take6 %>% 
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

take6_coords <- st_coordinates(take6)
take6$Lat <- take6_coords[,2]
take6$Long <- take6_coords[,1]
take6$Date <- as.POSIXct(take6$Date,format="%Y-%m-%d %H:%M:%S",tz="CST6CDT")

take6 <- take6 %>% 
  select(colnames(take6)[colnames(take6)%in%colnames(take5)])

take <- full_join(take,take6)


take$Method[grep("Aerial",take$Method)] <- "Aerial"
take$Method[grep("Shoot",take$Method)] <- "Ground Shoot"
take <- subset(take,Method!="Detection Dogs" & Method!="Snares")
take <- take[!duplicated(take),]

#remove wrong timestamps
take$date_notime <- as.Date(take$Date)
range(take$Date,na.rm=T)
take <- take %>% 
  filter(year(date_notime)>=2020 & year(date_notime)<=year(Sys.Date())) %>% 
  filter(Method!="<Null>")

#remove points outside of missouri
mo <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/missouri.shp")

take$id <- 1:nrow(take)
take_sf <- subset(take,!is.na(Lat))
take_sf <- st_as_sf(take_sf,coords=c("Long","Lat"),crs="epsg:4326")
take_sf <- st_intersection(take_sf,mo)
take <- take[which(take$id %in% take_sf$id),]
take <- take %>% select(-Shape)

ggplot()+geom_sf(data=mo)+
  geom_sf(data=take_sf,aes(col=Method,size=Total),alpha=0.5)


saveRDS(take,"C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/feral_swine_elim/Data/MO_take_clean.rds")

