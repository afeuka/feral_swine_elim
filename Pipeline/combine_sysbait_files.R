###Title: Organizing MO systematic feral swine baiting/trapping data
###Author: Abbey Feuka
###Date: 19MAR25

library(tidyverse)
library(sf)
library(lubridate)

setwd('C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Data/Systematic Baiting')

#clean systematic baiting -------------------------------
sys <- read.csv("./Ops_Sept2020_Dec2022.csv")
sys1 <- read.csv("./Ops_Jan2023_Jun2023.csv")
sys1<- sys1 %>% select(-c(GlobalID..,EditDate,Editor))
sys <- full_join(sys,sys1) 
sys2 <- read.csv("./Ops_Jul2023_Nov2023.csv")
sys2 <- sys2 %>% select(-c(GlobalID..))
sys <- full_join(sys,sys2 %>% mutate(Lat=as.character(Lat),Long=as.character(Long))) 
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
sys$CreationDate <- as.Date(sys$CreationDate,format="%m/%d/%Y %H:%M:%S")

sys3 <- read.csv("./Ops_July23_Sept24.csv")
sys3 <- sys3 %>% rename(Lat=lat,
                        Long=long,
                        SiteVisitReason=SiteVisitReas) %>% 
  mutate(Lat=as.character(Lat),
         Long=as.character(Long),
         SounderSize) %>% 
  select(SiteStatus, BaitingStartDate, DateHotBait, SiteVisitDate, SiteVisitReason, TrapStartDate,
         CompleteDate, Lat, Long, CreationDate, Creator, SounderSize)

sys3$BaitingStartDate <-as.Date(sys3$BaitingStartDate,format="%Y-%m-%d %H:%M:%S")
sys3$DateHotBait <- as.Date(sys3$DateHotBait,format="%Y-%m-%d %H:%M:%S")
sys3$CompleteDate <- as.Date(sys3$CompleteDate,format="%Y-%m-%d %H:%M:%S")
sys3$TrapStartDate <- as.Date(sys3$TrapStartDate,format="%Y-%m-%d %H:%M:%S")
sys3$SiteVisitDate <- as.Date(sys3$SiteVisitDate,format="%Y-%m-%d %H:%M:%S")
sys3$CreationDate <- as.Date(sys3$CreationDate,format="%Y-%m-%d %H:%M:%S")

sys <- full_join(sys,sys3)

sys4 <- read.csv("./Ops_Oct24_Dec24.csv")

sys4 <- sys4 %>% select(SiteStatus, BaitingStartDate, DateHotBait, SiteVisitDate, SiteVisitReason, TrapStartDate,
                        CompleteDate, Lat, Long, CreationDate, Creator, SounderSize) %>% 
  mutate(Lat=as.character(Lat),
         Long=as.character(Long))
sys4$BaitingStartDate <-as.Date(sys4$BaitingStartDate,format="%Y-%m-%d")
sys4$DateHotBait <- as.Date(sys4$DateHotBait,format="%Y-%m-%d")
sys4$CompleteDate <- as.Date(sys4$CompleteDate,format="%Y-%m-%d")
sys4$TrapStartDate <- as.Date(sys4$TrapStartDate,format="%Y-%m-%d")
sys4$SiteVisitDate <- as.Date(sys4$SiteVisitDate,format="%Y-%m-%d")
sys4$CreationDate <- as.Date(sys4$CreationDate,format="%Y-%m-%d")

sys <- full_join(sys,sys4)

saveRDS(sys,'C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/feral_swine_elim/Data/MO_sysbait_clean.rds')