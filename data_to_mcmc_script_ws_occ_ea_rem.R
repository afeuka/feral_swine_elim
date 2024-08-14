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

setwd("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/feral_swine_elim")
source("./Functions/data_functions_ws_occ_ea_rem.R")
grid_typ <- "watersheds"
period <- "season"
study_site_grid <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/huc10_siteIDs_cond_EA.shp")
study_site_grid <- study_site_grid %>% rename(elim_area_idx=elm_r_d,
                                              Area_Name=Area_Nm)
if(file.exists(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Model Ready Data/sysbait_10day_season_nlcd_neweff.RData"))){
  load(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Model Ready Data/sysbait_10day_season_nlcd_neweff.RData"))
} else {
  
  ## systematic baiting  ---------------------------
  sys <- grid_sysbait_take(study_site_grid = study_site_grid,
                           period=period)
  ## removal  ---------------------------
  rem <- grid_removals(study_site_grid = study_site_grid,
                       sysbait_det_eff = sys$sysbait_det_eff,
                       period=period)
  
  eff <- grid_effort(sysbait_det_eff=sys$sysbait_det_eff,
                     rem_site = rem$rem_site,
                     study_site_grid=study_site_grid,
                     grid_typ=grid_typ)
  
  sysbait_det_eff <- eff$sysbait_det_eff
  nlcd_siteid <- eff$nlcd_siteid
  rem_eff_ea <- eff$rem_eff_ea

  save(nlcd_siteid,sysbait_det_eff,rem_eff_ea,
       file="C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Model Ready Data/sysbait_10day_season_nlcd_neweff.RData")
}

#fit model -------------------
source("./Functions/avail_fun.R")
source("./Functions/fit_zi_rem_occ_mod_binom_ea_area.R")

##single chain -----------------
mcmc.out <- fit_zi_rem_occ(sysbait_det_eff = sysbait_det_eff,
                           rem_eff_ea = rem_eff_ea,
                           study_site_grid=study_site_grid,
                           elim_prob=0.95,
                           eff_weeks =10,
                           monitors = NA,
                           subset_data=F,
                           niter=10000,
                           thin=1,
                           burnProp=0.5,
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

# save(samples,dat_occ,dat_rem,nsites,nperiods,nbeta,
#      file = "./Model outputs/ziBinMod_determ_psi.Rdata")

##parallel --------------------
nChains <- 3
niter <- 300000
burnProp <- 0.75
thin <- 5

mcmcPar <- function(j){
  samp <- fit_zi_rem_occ(sysbait_det_eff = sysbait_det_eff,
                         rem_eff_ea = rem_eff_ea,
                         study_site_grid=study_site_grid,
                         elim_prob=0.95,
                         eff_weeks=10,
                         monitors=NA, 
                         subset_data=F,
                         niter=niter, 
                         thin=thin,
                         burnProp=burnProp,
                         nChains=1)
}

cl <- makeCluster(nChains, "SOCK")
clusterExport(cl, list("sysbait_det_eff","rem_eff_ea","nlcd_siteid",
                       "study_site_grid","niter","burnProp","thin",
                       "fit_zi_rem_occ","avail_fun"))

system.time(
  parSamples<- clusterApply(cl, 1:nChains, mcmcPar)
)
stopCluster(cl)

samples <- list(parSamples[[1]]$samples,
                parSamples[[2]]$samples,
                parSamples[[3]]$samples)
                # parSamples[[4]]$samples,
                # parSamples[[5]]$samples)
dat_occ <- parSamples[[1]]$dat_occ
dat_trap <- parSamples[[1]]$dat_trap
dat_aerial <- parSamples[[1]]$dat_aerial
dat_ground <- parSamples[[1]]$dat_ground
nsites <- parSamples[[1]]$nsites
nea <- parSamples[[1]]$nea
nperiods <- parSamples[[1]]$nperiods
nbeta <- parSamples[[1]]$nbeta
agri <- parSamples[[1]]$agri
develop <- parSamples[[1]]$develop

save(samples,dat_occ,nsites,nea,nperiods,nbeta,agri,develop,
     dat_trap,dat_ground,dat_aerial,
     file="C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/ziBinMod_area_08AUG24_validation.Rdata")
