# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) 
library(crew)

input_dir <- file.path('./data','input')
raw_dir <- file.path('./data','raw')
out_dir <- file.path('./output')

mcmc_controller = crew_controller_local(
  name = "mcmc_controller",
  workers = 3,
  seconds_idle = 60
)

# local_controller = crew_controller_local(
#   name = "local_contoller",
#   workers = 1,
#   seconds_idle = 10
# )
# Set target options:
tar_option_set(
  packages = c("sf", "tidyverse", "nimble","snow",
               "lubridate","ROracle","keyring",
               "tigris","amt","terra","exactextractr",
               "pROC","gridExtra"), # packages that your targets use
  controller = mcmc_controller
  # resources = tar_resources(
  #   crew = tar_resources_crew(controller = "local_controller")
  # )
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
tar_source("plot_functions.R")

list(
  
  ## cleaning/formatting raw data --------------
  tar_target(
    name = ssg,
    command = file.path(input_dir,"HUC10_MO","huc10_siteIDs_cond_EA.shp"),
    format="file",
    description = "HUC10 study site grid file" 
  ),
  tar_target(
    name = sec,
    command = file.path(input_dir,"MO_Public_Land_Survey_System",
                        "MO_Public_Land_Survey_System.shp"),
    format="file",
    description = "MO Public Land Survey sections file" 
  ),
  tar_target(
    name = ssg_formatted,
    command = format_ssg(ssg),
    description = "Format HUC10 watershed file for MO" 
  ),
  tar_target(
    name=elim_areas,
    command=make_elim_areas(ssg_formatted),
    description = "Make elimination areas from study site grid"
  ),
  tar_target(
    name = sec_formatted,
    command = format_sec(sec,ssg_formatted),
    description = "Format PLS sections file for MO" 
  ),
  tar_target(
    name = sysbait_grid,
    command = grid_sysbait(sysbait_clean,
                           take_clean,
                           sec_formatted,
                           ssg_formatted,
                           start_date="2020-10-01", #start date for data in y-m-d string
                           end_date="2025-03-30"),
    description = "Grid systematic baiting data" 
  ),
  tar_target(
    name = take_grid,
    command = grid_removals(take_clean,sysbait_grid,sec_formatted,ssg_formatted),
    description = "Grid take/removal data" 
  ),
  tar_target(
    name=effort_query,
    command=mis_effort_query(take_grid),
    description = "Read effort query"
  ),
  tar_target(
    name=effort_raw,
    command=pull_effort_mis(effort_query,ssg_formatted),
    description = "Pull effort data from MIS"
  ),
  tar_target(
    name=sysbait_effort,
    command = grid_sysbait_effort(effort_raw,sysbait_grid,ssg_formatted),
    description="Join effort and sysbait data"
  ),
  tar_target(
    name=take_effort,
    command = grid_take_effort(sysbait_effort,
                               take_grid,
                               effort_raw,
                               aerial_ops_tracks,
                               ssg_formatted),
    description = "Join effort, aerial ops tracks, and take data"
  ),
  # tar_target(
  #   name=nfsp_siteid,
  #   command=file.path(input_dir,"ssg_nfsp.rds"),
  #   format="file",
  #   description="NFSP feral swine range data"
  # ),
  # tar_target(
  #   name=raw_ao_fy21,
  #   command=file.path(raw_dir,"aerial_ops_tracks","FY2021 AERIAL OPS TRACK DATA.csv"),
  #   format="file",
  #   description="FY21 raw aerial ops track data"
  # ),
  # tar_target(
  #   name=raw_ao_fy22,
  #   command=file.path(raw_dir,"aerial_ops_tracks","FY2022 MO AERIAL OPS TRACK DATA.csv"),
  #   format="file",
  #   description="FY22 raw aerial ops track data"
  # ),
  # tar_target(
  #   name=raw_ao_fy23,
  #   command=file.path(raw_dir,"aerial_ops_tracks","FY2023 MO AERIAL OPS TRACK DATA.csv"),
  #   format="file",
  #   description="FY23 raw aerial ops track data"
  # ),
  # tar_target(
  #   name=raw_ao_fy24,
  #   command=file.path(raw_dir,"aerial_ops_tracks","FY24_MOFlightTracks.gdb"),
  #   format="file",
  #   description="FY24 raw aerial ops track data"
  # ),
  # tar_target(
  #   name=raw_ao_dec24_mar25,
  #   command=file.path(raw_dir,"aerial_ops_tracks","FlightTracks_Dec_MarchFY25.gdb"),
  #   format="file",
  #   description="Dec 2024 to Mar 2025 raw aerial ops track data"
  # ),
  tar_target(
    name=raw_ao_folder,
    command=file.path(raw_dir,"aerial_ops_tracks"),
    format="file",
    description="Location of raw aerial ops track data"
  ),
  tar_target(
    name=aerial_ops_tracks,
    command=clean_aerial_tracks(#raw_ao_fy21,
                                #raw_ao_fy22,
                                #raw_ao_fy23,
                                #raw_ao_fy24,
                                #raw_ao_dec24_mar25,
                                raw_ao_folder,
                                ssg_formatted),
    description="Read, clean, and collate raw aerial track data"
  ),
  # tar_target(
  #   name=raw_sys_sep20_dec22,
  #   command=file.path(raw_dir,"sysbait","Ops_Sept2020_Dec2022.csv"),
  #   format="file",
  #   description="Sept20-Dec22 raw sysbait data"
  # ),
  # tar_target(
  #   name=raw_sys_jan23_jun23,
  #   command=file.path(raw_dir,"sysbait","Ops_Jan2023_Jun2023.csv"),
  #   format="file",
  #   description="Jan23-Jun22 raw sysbait data"
  # ),
  # tar_target(
  #   name=raw_sys_jul23_nov23,
  #   command=file.path(raw_dir,"sysbait","Ops_Jul2023_Nov2023.csv"),
  #   format="file",
  #   description="Jul23-Nov23 raw sysbait data"
  # ),
  # tar_target(
  #   name=raw_sys_jul23_sep24,
  #   command=file.path(raw_dir,"sysbait","Ops_July23_Sept24.csv"),
  #   format="file",
  #   description="Jul23-Sept24 raw sysbait data"
  # ),
  # tar_target(
  #   name=raw_sys_oct24_dec24,
  #   command=file.path(raw_dir,"sysbait","Ops_Oct24_Dec24.csv"),
  #   format="file",
  #   description="Oct24-Dec24 raw sysbait data"
  # ),
  tar_target(
    name=raw_sysbait_folder,
    command=file.path(raw_dir,"sysbait"),
    format="file",
    description="Location of raw systematic baiting/ops files"
  ),
  ## data cleaning functions ----------------
  tar_target(
    name=sysbait_clean,
    command=clean_sysbait(#raw_sys_sep20_dec22,
                          #raw_sys_jan23_jun23,
                          #raw_sys_jul23_nov23,
                          #raw_sys_jul23_sep24,
                          #raw_sys_oct24_dec24,
                          #raw_systake_jan25_mar25_gdb,
                          raw_sysbait_folder),
    description="Read, clean and collate raw sysbait data"
  ),
  # tar_target(
  #   name=raw_take_sep20_dec22,
  #   command=file.path(raw_dir,"take","Take_Sept2020_Dec2022.csv"),
  #   format="file",
  #   description="Sept20-Dec22 raw take data"
  # ),
  # tar_target(
  #   name=raw_take_jan23_jun23,
  #   command=file.path(raw_dir,"take","Take_Jan2023_Jun2023.csv"),
  #   format="file",
  #   description="Jan23-Jun22 raw take data"
  # ),
  # tar_target(
  #   name=raw_take_jul23_nov23,
  #   command=file.path(raw_dir,"take","Take_Jul2023_Nov2023.csv"),
  #   format="file",
  #   description="Jul23-Nov23 raw take data"
  # ),
  # tar_target(
  #   name=raw_take_jul23_sep24,
  #   command=file.path(raw_dir,"take","Take_July23_Sept24.csv"),
  #   format="file",
  #   description="Jul23-Sept24 raw take data"
  # ),
  # tar_target(
  #   name=raw_take_oct24_dec24,
  #   command=file.path(raw_dir,"take","MOFeralHogData_Oct24_Dec24.gdb"),
  #   format="file",
  #   description="Oct24-Dec24 raw take data"
  # ),
  # tar_target(
  #   name=raw_systake_jan25_mar25_gdb,
  #   command=file.path(raw_dir,"take","Q1CY2025_FeralHogData.gdb"),
  #   format="file",
  #   description="Oct24-Dec24 raw sys and take data"
  # ),
  tar_target(
    name=raw_take_folder,
    command=file.path(raw_dir,"take"),
    format="file",
    description="Read, cocation of raw take data"
  ),
  tar_target(
    name=take_clean,
    command=clean_take(#raw_take_sep20_dec22,
                       #raw_take_jan23_jun23,
                       #raw_take_jul23_nov23,
                       #raw_take_jul23_sep24,
                       #raw_take_oct24_dec24,
                       #raw_systake_jan25_mar25_gdb,
                       raw_take_folder,
                       ssg_formatted),
    description="Clean and collate raw take data"
  ),
  
  ## data prep for model --------------------
  tar_target(
    name=nlcd_raster,
    command=  file.path(raw_dir,"nlcd","NLCD_2019_Land_Cover_L48_20210604_jj2pnCpncP4xHk8rj4R5.tiff"),
    description="Raw MO NLCD raster file",
    format="file"
  ),
  tar_target(
    name=nlcd_siteid,
    command=extract_nlcd(ssg_formatted,nlcd_raster),
    description = "Extract NLCD data for study site grid"
  ),
  tar_target(
    name=oak_raster,
    command=  file.path(raw_dir,"gap2011_mo","gap_2011_mo_oak.tif"),
    description="Raw MO oak raster file",
    format="file"
  ),
  tar_target(
    name=oak_siteid,
    command=extract_oak(ssg_formatted,oak_raster),
    description = "Extract GAP oak data for study site grid"
  ),
  # tar_target(
  #   name=nfsp_siteid,
  #   command=extract_nfsp(ssg_formatted,...),
  #   description = "Extract feral swine range data for study site grid"
  # )
  tar_target(
    name=dat_occ,
    command=occ_data_prep(sysbait_effort,ssg_formatted),
    description="Format sysbait data for occupancy model"
  ),
  tar_target(
    name=dat_rem,
    command=rem_data_prep(take_effort,dat_occ),
    description="Format take data for removal model"
  ),
  tar_target(
    name=lat_long_siteid,
    command=site_lat_longs(ssg_formatted),
    description = "Extract lat longs from study site grid as covariate"
  ),
  
  ## model fitting functions -----------------
  tar_target(
    name=site_idx_lookup,
    command=align_sites(dat_rem),
    description="Create index key for occupancy and removal sites"
  ),
  tar_target(
    name=avail_fun,
    command=nimble_fun_avail(),
    description = "Make nimbleFunction for availability"
  ),
  tar_target(
    name=chain_idx,
    command=make_chain_idx(nChains=3),
    description="Make chain for parallelization"
  ),
  tar_target(
    name=samples,
    command=fit_zibinom_mod(dat_occ,
                            dat_rem,
                            avail_fun,
                            ssg_formatted,
                            site_idx_lookup,
                            nlcd_siteid,
                            lat_long_siteid,
                            oak_siteid,
                            # eff_weeks=4,
                            chain_idx=chain_idx,
                            niter=100000,
                            thin=5,
                            burnProp=0.8),
    description="Compile and fit occupancy/removal model",
    pattern=map(chain_idx)
    # resources=tar_resources(
    #   crew = tar_resources_crew(controller = "mcmc_controller")
    # )
  ),
  tar_target(
    name=wi_mod_val,
    command=model_validation(dat_occ,dat_rem,samples,out_dir),
    description="Within sample model validation and goodness of fit"
  ),
  
  ## plotting functions --------------------
  tar_target(
    name=trace_plots,
    command=plot_trace(samples,
                       dat_rem,
                       out_dir,
                       chain_idx,
                       site_idx_lookup,
                       lat_long_siteid,
                       oak_siteid),
    description="Make trace plots"
  ),
  tar_target(
    name=pabs_plots,
    command=plot_pabs(samples,dat_occ,ssg_formatted,elim_areas,out_dir),
    description="Make p(absence) plots",
    deployment="main"
  ),
  tar_target(
    name=det_plot,
    command=plot_det(samples,dat_occ,out_dir),
    description="Make occupancy detection curve",
    deployment="main"
  ),
  tar_target(
    name=elim_samp_plot,
    command=plot_elim_samp(samples,dat_occ,out_dir),
    description = "Make number of samples to elimination plot",
    deployment="main"
  ),
  tar_target(
    name=abundance_plots,
    command=plot_abundance(samples,
                           dat_rem,
                           dat_occ,
                           ssg_formatted,
                           elim_areas,
                           site_idx_lookup,
                           out_dir),
    description="Make abundance plots",
    deployment="main"
  ),
  tar_target(
    name=removal_prob_plots,
    command=plot_removal_prob(samples,dat_rem,out_dir),
    description = "Make removal probability plots and dataframe",
    deployment="main"
  ),
  tar_target(
    name=data_plots,
    command=plot_data(abundance_plots,out_dir),
    description="Make raw data plots",
    deployment="main"
  )
)
