### Title: Hurdle model posterior results (in text)
### Author: Abbey Feuka
### Date: 22AUG24
### Notes: uses output from fit_zi_rem_mod.R

library(tidyverse)
library(sf)

setwd("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble")

##study site--------------------
study_site_grid <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/huc10_siteIDs_cond_EA.shp")
study_site_grid <- study_site_grid %>% rename(elim_area_idx=elm_r_d,
                                              Area_Name=Area_Nm)
## elimination areas -----------------
elim_areas <- study_site_grid %>%
  group_by(Area_Name,elim_area_idx) %>% 
  summarise(geometry=st_union(geometry))

# load samples ---------------------
load("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/ziBinMod_area_04SEP24_logit_det.Rdata")
load("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/Plots/NFSP/posterior_summaries_04SEP24.Rdata")

if(nbeta==3){
  subfolder<-"No NFSP"
} else {
  subfolder<-"NFSP"
}

## removal data -----------
dat_rem <- rbind.data.frame(dat_aerial,dat_trap)
dat_rem_sum <- dat_rem %>% 
  group_by(Area_Name,method,month=floor_date(Date,"month")) %>% 
  summarise(Season=min(Date),
            removal=sum(tot_rem),
            tot_hrs=sum(tot_hrs),
            tot_events=sum(num_events)) %>% 
  filter(month<=max(dat_occ$per_start))

nmcmc <- nrow(samples[[1]])
nChains <- length(samples)

# occupancy regression --------------
beta<- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
  cbind.data.frame(samples[[i]][,grep("beta",colnames(samples[[i]]))],
                   chain=i)})),
  samp=rep(1:nmcmc,nChains))

if(nbeta==3){
  beta_names <- c("Intercept","Development","Agriculture")
} else {
  beta_names <- c("Intercept","Feral swine range","Development","Agriculture")
}

colnames(beta)[1:nbeta] <- beta_names

beta_long <- beta %>% 
  pivot_longer(cols=all_of(1:nbeta),names_to="beta",values_to="value") %>% 
  mutate(chain=factor(chain))

beta_long %>% 
  group_by(beta) %>% 
  summarise(mn=mean(value),
            md=median(value),
            lci=quantile(value,0.025),
            uci=quantile(value,0.975))

# probability of absence ------------------
## watershed values ------------------
# all sites and periods mean and quantile
mean(rowMeans(pabs))
quantile(rowMeans(pabs),c(0.025,0.975))

#maximum and minimum p(absence)
pabs_sum[which.max(pabs_sum$mn),]
pabs_sum[which.min(pabs_sum$mn),]

## elimination area medians -------------
pabs_ea_sum[which.max(pabs_ea_sum$ea_md),]
pabs_ea_sum[which.min(pabs_ea_sum$ea_md),]

# probability of elimination (conditional) ----------------
## watershed values ---------------
# all sites and periods mean and quantile
mean(rowMeans(pelim))
quantile(rowMeans(pelim),c(0.025,0.975))

#maximum and minimum p(absence)
pelim_sum[which.max(pelim_sum$mn),]
pelim_sum[which.min(pelim_sum$mn),]

## elimination area medians ---------------
pelim_ea[which.max(pelim_ea$md),]
pelim_ea[which.min(pelim_ea$md),]

## spatial extent of elimination -----------
pabs_thresh_yr

# effort to 95% elim -------------
#range of sample values
eff_sum <- st_read(paste0("./Model outputs/Plots/",subfolder,"/eff_sum_sf.shp"))
eff_sum <- eff_sum %>% st_drop_geometry() %>% 
  rename(period_idx=perd_dx,
         site_idx=site_dx,
         eff_mn=mn,
         eff_md=md,
         eff_lci=lci,
         eff_uci=uci) %>% 
  select(eff_mn,eff_md,eff_lci,eff_uci,site_idx,period_idx)
eff_sum <- pabs_sum %>% left_join(eff_sum)
eff_sum %>% filter(mn>=0.4 & mn<=0.6)
# eff_sum[which.min(eff_sum$mn),]
# eff_sum[which.max(eff_sum$mn),]

# abundance -------------------
# density estimate min/max for EAs 
N_sum_sf<- st_read(paste0("./Model outputs/Plots/",subfolder,"/N_sum_sf.shp"))
N_sum_4 <- N_sum_sf %>% st_drop_geometry() %>% filter(Area_Nm==4)
N_sum_4[which.min(N_sum_4$mn_dens),]
N_sum_4[which.max(N_sum_4$mn_dens),]

N_sum_6 <- N_sum_sf %>% st_drop_geometry() %>% filter(Area_Nm==6)
N_sum_6[which.min(N_sum_6$mn_dens),]
N_sum_6[which.max(N_sum_6$mn_dens),]

## fy trends ---------------
N_yr_sum

# lambda --------------------
if(!exists("lambda")){
  lambda <- do.call("rbind",lapply(1:nChains,function(i){
    samples[[i]][grepl("lambda",colnames(samples[[i]]))]}))
}
lambda <- lambda[,c("lambda[5]","lambda[7]")]
data.frame(mn=colMeans(lambda),
           lci=sapply(1:2,function(i)quantile(lambda[,i],probs=c(0.025))),
           uci=sapply(1:2,function(i)quantile(lambda[,i],probs=c(0.975))))

# removal probability ----------------
# standardize effort and compare
# a_mn <- max(dat_aerial$prop_ea_impact)*(1-(1-p_rem[,"Aerial"])^max(dat_aerial$eff_area_hrs[dat_aerial$eff_area_hrs!=0]))
# t_mn <- max(dat_trap$prop_ea_impact)*(1-(1-p_rem[,"Trap"])^max(dat_trap$eff_area_hrs[dat_trap$eff_area_hrs!=0]))
# 
# rem_mn <- cbind.data.frame(aerial=a_mn,trap=t_mn,ground=g_mn)
# rem_df <- rem_mn %>% pivot_longer(cols=1:3,names_to="rem_typ",
#                                   values_to="value") %>% 
#   group_by(rem_typ) %>% 
#   summarise(mn=mean(value),
#             lci=quantile(value,0.025),
#             uci=quantile(value,0.975))
# rem_df$mn[1]/rem_df$mn[2]
# rem_df$mn[1]/rem_df$mn[3]
# range(elim_areas %>% filter(Area_Name!="0") %>% st_area()/1e6)

