### Title: Gridding MO Camera, systematic baiting, and removal data
# and fitting to dynamic occupancy model
### Author: Abbey Feuka
### Date: 24072023
### Notes:

# setup ------------------------------------
library(sf)
library(tidyverse)
library(nimble)
library(snow)

# setwd("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/feral_swine_elim")

#data cleaning/firmatting functions
source("./Functions/data_functions_occ_rem_dynamic.R")

#where output files are sent
input_dir <- file.path('./Data','input')
output_dir <- file.path('./Data','model_ready')

end_fy <- "24"
occ_scale <- "watersheds"
abund_scale <- "watersheds"
pass_scale <- "sys" 
period <- "season"
study_site_grid <- st_read(file.path(input_dir,
                                     "HUC10_MO","huc10_siteIDs_cond_EA.shp")) %>% 
  rename(elim_area_idx=elm_r_d,
         Area_Name=Area_Nm)

if(file.exists(file.path(output_dir,"mo_sysbait_rem_2020-10-01_2024-12-31.RData"))){
  load(file.path(output_dir,"mo_sysbait_rem_2020-10-01_2024-12-31.RData"))
} else {
  
  ## systematic baiting  ---------------------------
  sys <- grid_sysbait_take(study_site_grid = study_site_grid,
                           start_date="2020-09-01", 
                           end_date=paste0("20",end_fy,"-12-31"),
                           period=period,
                           data_dir=input_dir)
  ## removal  ---------------------------
  rem <- grid_removals(sec=sys$sec,
                       study_site_grid = study_site_grid,
                       sysbait_det_eff = sys$sysbait_det_eff,
                       period=period,
                       data_dir=input_dir)
  
  eff <- grid_effort(sysbait_det_eff=sys$sysbait_det_eff,
                     rem_site_sf = rem$rem_site_sf,
                     study_site_grid=study_site_grid,
                     occ_scale=occ_scale,
                     abund_scale=abund_scale,
                     pass_scale=pass_scale,
                     data_dir=input_dir)
  
  cov <- grid_covs(occ_scale=occ_scale,
                   abund_scale=abund_scale,
                   sysbait_det_eff=sys$sysbait_det_eff,
                   study_site_grid=study_site_grid,
                   data_dir=input_dir)
  
  sysbait_det_eff <- eff$sysbait_det_eff
  rem_eff_site <- eff$rem_eff_site
  
  nlcd_siteid <- cov$nlcd_siteid
  oak_siteid <- cov$oak_siteid
  nfsp_siteid <- cov$ssg_nfsp

  save(nlcd_siteid,oak_siteid,nfsp_siteid,sysbait_det_eff,rem_eff_site,
       file=file.path(output_dir,
             paste0("mo_sysbait_rem_",
                    min(sysbait_det_eff$subper_start),"_",
                    max(sysbait_det_eff$subper_end),
                    ".RData")
                   ))
}

# max(sysbait_det_eff$subper_start)
# max(rem_eff_site$Date)

# hist(rem_eff_site$effect_area_hrs)
# hist(rem_eff_site$prop_site_impact)

#fit model -------------------
source("./Functions/avail_fun.R")
source("./Functions/fit_zi_rem_occ_mod_binom_dynamic.R")

##single chain -----------------
mcmc.out <- fit_zi_rem_occ(sysbait_det_eff = sysbait_det_eff,
                           rem_eff_site = rem_eff_site,
                           study_site_grid=study_site_grid,
                           oak_siteid=oak_siteid,
                           abund_scale=abund_scale,
                           eff_weeks =10,
                           monitors = NA,
                           niter=100000,
                           thin=5,
                           burnProp=0.8,
                           nChains=1)

samples<- mcmc.out$samples
dat_occ <- mcmc.out$dat_occ
dat_aerial <- mcmc.out$dat_aerial
dat_trap <- mcmc.out$dat_trap
dat_ground <- mcmc.out$dat_ground
nsites <- mcmc.out$nsites
nea <- mcmc.out$nea
nperiods <- mcmc.out$nperiods
nbeta <- mcmc.out$nbeta
site_idx_lookup <- mcmc.out$site_idx_lookup
oak_siteid <- mcmc.out$oak_siteid

# save(samples,dat_occ,dat_rem,nsites,nperiods,nbeta,
#      file = "./Model outputs/ziBinMod_determ_psi.Rdata")

##parallel --------------------
nChains <- 8
niter <- 500000
burnProp <- 0.80
thin <- 10

mcmcPar <- function(j){
  samp <- fit_zi_rem_occ(dat_occ = sysbait_det_eff,
                         dat_rem = rem_eff_site,
                         study_site_grid=study_site_grid,
                         oak_siteid=oak_siteid,
                         nlcd_siteid=nlcd_siteid,
                         nfsp_siteid=nfsp_siteid,
                         eff_weeks=4,
                         abund_scale=abund_scale,
                         monitors=NA, 
                         niter=niter, 
                         thin=thin,
                         burnProp=burnProp,
                         nChains=1)
}

cl <- makeCluster(nChains, "SOCK")
clusterExport(cl, list("sysbait_det_eff","rem_eff_site","nlcd_siteid",
                       "abund_scale","oak_siteid",
                       "study_site_grid","niter","burnProp","thin",
                       "fit_zi_rem_occ","avail_fun"))

system.time(
  parSamples<- clusterApply(cl, 1:nChains, mcmcPar)
)
stopCluster(cl)

samples <- list(parSamples[[1]]$samples,
                parSamples[[2]]$samples,
                parSamples[[3]]$samples,
                parSamples[[4]]$samples,
                parSamples[[5]]$samples,
                parSamples[[6]]$samples,
                parSamples[[7]]$samples,
                parSamples[[8]]$samples
                )
samples <- lapply(1:nChains,function(i)parSamples[[i]]$samples)

dat_occ <- parSamples[[1]]$dat_occ
dat_trap <- parSamples[[1]]$dat_trap
dat_aerial <- parSamples[[1]]$dat_aerial
nsites_rem <- parSamples[[1]]$nsites_rem
nsites_occ <- parSamples[[1]]$nsites_occ
# nea <- parSamples[[1]]$nea
nperiods <- parSamples[[1]]$nperiods
nbeta <- parSamples[[1]]$nbeta
agri <- parSamples[[1]]$agri
develop <- parSamples[[1]]$develop
mn_te <- parSamples[[1]]$mn_te
site_idx_lookup <- parSamples[[1]]$site_idx_lookup
rem_covs <- parSamples[[1]]$rem_covs

save(samples,dat_occ,nsites_occ,nsites_rem,nperiods,nbeta,agri,develop,
     dat_trap,dat_aerial,mn_te,site_idx_lookup,rem_covs,
     file=paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/mo_zibinom",
                 min(sysbait_det_eff$subper_start),"_",
                 max(sysbait_det_eff$subper_end),"_v1.Rdata"))
