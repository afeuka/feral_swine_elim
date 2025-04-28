### Title: MO Feral Swine Manuscript Stats
### Author: Abbey Feuka
### Date: 16APR25
### Notes: Use after running targets pipeline

library(targets)
library(tidyverse)
library(sf)

#removal data summary
dat_rem <- tar_read(dat_rem)
dat_rem %>% 
  filter(method!="Ground") %>% 
  group_by(method) %>% 
  summarise(tot_rem=sum(tot_rem),
            tot_hrs=sum(tot_hrs))

#occ data summary
study_site_grid <- tar_read(ssg_formatted)

dat_occ <- tar_read(dat_occ)
dat_occ %>% 
  left_join(study_site_grid %>% st_drop_geometry() %>% rename(site_idx_occ=SiteID)) %>% 
  ungroup() %>% 
  summarise(ndet=sum(detections),
            ntrapnights=sum(trap_nights_km*area_km))

#model evaluation
modval <- read.csv(file.path('output','validation','val_table.csv'))
modval %>% 
  filter(validation_meth=='pVal'&stat=='mean')
modval %>% 
  filter(validation_meth=='pVal'&stat=='var')
modval %>% 
  filter(validation_meth=='AUC')

#pabs
#range in estimates
load(file.path('output','samples','posterior_summaries_10APR25.RData'))
pabs_sum[which.min(pabs_sum$md),]
pabs_sum[which.max(pabs_sum$md),]

#ea level medians
pabs_sum_ea %>% 
  filter(Area_Name!="Outside EAs") %>% 
  group_by(Area_Name) %>% 
  summarise(min=min(md),
            max=max(md)) %>% 
  arrange(min)

#highest pabs
pabs_sum_ea %>% 
  filter(Area_Name=="Elimination Area 1") %>% 
  arrange(md) %>% 
  select(-mn)

#lowest pabs
pabs_sum_ea %>% 
  filter(Area_Name=="Elimination Area 4") %>% 
  arrange(md) %>% 
  select(-mn)

#area above 95
area95 <- read.csv(file.path('output','tables','area_above_0.95_table.csv'))
area95 %>% select(fy,md,lci,uci) %>% round(0)
area95 %>% select(fy,md_prop,lci_prop,uci_prop) %>% round(3)
area95 %>% select(fy,md_prop,lci_prop,uci_prop) %>% round(3) %>% arrange(md_prop)

#beta summary
trace_out <- tar_read(trace_plots)
trace_out$beta_sum %>% select(-mn)

trace_out$beta_lam_sum %>% select(-mn) %>% 
  mutate(md=round(md,2),
         lci=round(lci,2),
         uci=round(uci,2))

#density trends 
ntab <- read.csv(file.path('output','tables','N_change_table.csv'))
ntab <- ntab %>% 
  select(SiteID,per_start,md,lci,uci) %>% 
  mutate(md=round(md*100,2),
         lci=round(lci*100,2),
         uci=round(uci*100,2)) %>% 
  arrange(md)

# removal needed 
abund <- tar_read(abundance_plots)
rem_need <- abund$N_sum %>% 
  select(SiteID,per_start,md_lambda_prop_dens,lci_lambda_prop_dens,uci_lambda_prop_dens) %>% 
  mutate(md_lambda_prop_dens=round(md_lambda_prop_dens,2),
         lci_lambda_prop_dens=round(lci_lambda_prop_dens,2),
         uci_lambda_prop_dens=round(uci_lambda_prop_dens,2))
rem_need[which.min(rem_need$md_lambda_prop_dens),]
rem_need[which.max(rem_need$md_lambda_prop_dens),]

#removal prob
rem_prob <- tar_read(removal_prob_plots)
theta_gamma <- rem_prob$rem_df %>%  
  group_by(rem_typ,effect_area_hrs,prop_site_impact) %>% 
  summarise(md_rem_prob=md*prop_site_impact,
            lci_rem_prob=lci*prop_site_impact,
            uci_rem_prob=uci*prop_site_impact)
theta_gamma

days_attempted <- 45
theta_gamma %>% 
  summarise(md=md_rem_prob*45,
            lci=lci_rem_prob*45,
            uci=uci_rem_prob*45)




N_sum <- abund$N_sum
N_sum %>% 
  filter(SiteID==43) %>% 
  select(md_N,md_N_no_rem,lci_N,lci_N_no_rem)

samples <- tar_read(samples)
site_idx_lookup <- tar_read(site_idx_lookup)
idx<- site_idx_lookup[site_idx_lookup$SiteID==43,]

colnames(samples[,grepl("lambda",colnames(samples))])
plot(lam,typ="l")
