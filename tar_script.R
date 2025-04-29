library(targets)

# tar_manifest() 

# tar_visnetwork()
my_visNetwork <- function(physics = FALSE, scale = 1, levelSeperation = 150) {
  
  net <- tar_visnetwork(targets_only = T, exclude = dplyr::ends_with("_file"))
  
  net <- visNetwork::visNodes(net, physics = physics)
  
  net$x$nodes$level <- net$x$nodes$level * scale
  
  net <- visNetwork::visHierarchicalLayout(net, levelSeparation = levelSeperation, direction = "LR")
  
  net 
}

my_visNetwork(scale=0.1, levelSeperation = 1000)

tar_make()

samp <- tar_read(samples)
saveRDS(samp,file.path('output','samples',"mo_fs_elim_samp.rds"))

# save output
trace <- tar_read(trace_plots)
beta_sum <- trace$beta_sum
beta_lam_sum <- trace$beta_lam_sum
lam_sum <- trace$lam_sum

pabs <- tar_read(pabs_plots)
pabs_sum <- pabs$pabs_sum
pabs_sum_fy <- pabs$pabs_sum_fy
pabs_sum_ea <- pabs$pabs_sum_ea
pabs_thresh <- pabs$pabs_thresh
pabs_thresh_yr <- pabs$pabs_thresh_yr
pabs_thresh_ext <- pabs$pabs_thresh_ext

eff_elim_sum <- tar_read(elim_samp_plot)

det_sum <- tar_read(det_plot)

abund <- tar_read(abundance_plots)
N_sum <- abund$N_sum
N_yr_sum <- abund$N_yr_sum
dat_rem_sum <- abund$dat_rem_sum

rem_prob <- tar_read(removal_prob_plots)
rem_df <- rem_prob$rem_df
det_a_sum <- rem_prob$det_a_sum
det_t_sum <- rem_prob$det_t_sum

dat_occ <- tar_read(dat_occ)
dat_rem <- tar_read(dat_rem)

date_text <- "10APR25"

save(beta_sum,beta_lam_sum,lam_sum,
     pabs_sum,pabs_sum_fy,pabs_sum_ea,
     pabs_thresh,pabs_thresh_yr,pabs_thresh_ext,
     eff_elim_sum,
     N_sum,N_yr_sum,dat_rem_sum,
     rem_df,det_sum,det_a_sum,det_t_sum,
     dat_occ,dat_rem,
     file=file.path('output','samples',
                    'posterior_summaries_',
                    date_text,'.RData'))

#save ouptut for shiny
library(sf)
library(tidyverse)
ssg_formatted <- tar_read(ssg_formatted)
pabs_sum_sf<- pabs_sum %>% 
  left_join(ssg_formatted %>% rename(site_idx=SiteID)) %>% 
  st_as_sf() %>% 
  st_simplify()

date_text <- "10APR25"
st_write(pabs_sum_sf,
         'C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Shiny/mo_feral_swine_elim_app/resources/afeuka/mo-fs-elim/model_outputs/pabs_sum_sf_simple_',
         date_text,'.shp')


site_idx_lookup <- tar_read(site_idx_lookup)

samp <- tar_read(samples)
alpha <- samp[,grepl("alpha",colnames(samp)) | 
                grepl("chain",colnames(samp))]
sd_pdet <- samp[,grepl("sd_pdet",colnames(samp)) | 
                  grepl("chain",colnames(samp))]

save(beta_sum,beta_lam_sum,lam_sum,
     pabs_sum,pabs_sum_fy,pabs_sum_ea,
     pabs_thresh,pabs_thresh_yr,pabs_thresh_ext,
     eff_elim_sum,
     N_sum,N_yr_sum,dat_rem_sum,
     rem_df,det_sum,det_a_sum,det_t_sum,
     dat_occ,dat_rem,
     alpha,sd_pdet,site_idx_lookup,
     file='C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Shiny/mo_feral_swine_elim_app/resources/afeuka/mo-fs-elim/model_outputs/posterior_summaries_',
     date_text,'.RData')
