### Title: Hurdle model (removal) for southern MO feral swine occupancy/abundance
### Author: Abbey Feuka
### Date: 27NOV2023 edited 22FEB2024
### Notes: 
# 
# source("./Functions/avail_fun.R")
# monitors = NA
# niter=10000
# burnProp=0.1
# nChains=1
# thin=1
# elim_prob=0.95
# eff_weeks=10
# nfsp_reg=T
# abund_scale="properties"

fit_zi_rem_occ <- function(sysbait_det_eff, #output of data_functions_ws_occ_ea_rem.R
                           rem_eff_prop, #output of data_functions_ws_occ_ea_rem.R
                           study_site_grid, #output from mapping huc10 to eas.R
                           # elim_prob=0.95,#threshold for probability of elimination
                           eff_weeks=10,#number of trap nights to determine eliminiation probability
                           nfsp_reg=TRUE,#use nfsp range in occupancy estimates
                           abund_scale,#watersheds or ea
                           monitors, #parameters to save mcmc samples for
                           niter, #number of mcmc iterations
                           thin=1,#thinning interval
                           burnProp,
                           nChains
){
  
  require(nimble)
  require(tidyverse)
  require(sf)
  
  # data setup -------------------------------------------------
  dat_rem <- rem_eff_prop
  range(rem_eff_prop$county_area_km,na.rm=T)
  
  #filter no effect area (for now)
  dat_rem <- dat_rem %>% filter(!is.na(effect_area_km) & !is.na(effect_area_hrs))
  
  ##filter effort outliers ------------
  eff_cutoff <- quantile(dat_rem$effect_area_hrs,prob=0.95,na.rm=T)
  dat_rem <- dat_rem %>% filter(effect_area_hrs<eff_cutoff)
  dat_rem <- dat_rem %>% 
    group_by(period,Method,PropID) %>% 
    mutate(pass_idx=1:n())
  
  ## only include sites with consistent removal ----------
  site_rem <- dat_rem %>%
    group_by(Date,PropID) %>%
    summarise(tot_rem=sum(tot_rem)) %>%
    filter(tot_rem>0) %>%
    group_by(PropID) %>%
    summarise(n=n()) %>%
    arrange(desc(n))
  
  rem_cutoff <- quantile(site_rem$n,prob=0.75)
  site_rem_cutoff <- site_rem %>% filter(n>=rem_cutoff)
  dat_rem <- dat_rem %>% filter(PropID%in%site_rem_cutoff$PropID)
  
  dat_rem <- dat_rem %>% rename(period_idx = period) 
  dat_rem$prop_rem_idx <- as.numeric(as.factor(dat_rem$PropID))
  
  # dat_rem %>% group_by(Method,period_idx,Method,pass_idx) %>% 
  #   summarise(n=n()) %>% 
  #   filter(n>1) 
  
  ##add 0 data for later periods -aerial (ends in march) -------------
  a_pers <- dat_rem %>% ungroup() %>%
    filter(Method=="Aerial") %>%
    select(period_idx) %>% distinct
  a_pers_n <- which(!(1:max(dat_rem$period_idx)%in%a_pers$period_idx)) #no aerial ops during these periods
  
  # dat_rem$site_occ_idx <- as.numeric(dat_rem$SiteID)
  
  sids <- sort(unique(dat_rem$prop_rem_idx))
  
  samp <- dat_rem[1:(length(a_pers_n)*length(sids)),] #use existing data frame structure and clear fields
  samp$Method <- "Aerial"
  samp$period_idx <- rep(a_pers_n,length(sids))
  samp$prop_rem_idx <- sort(rep(unique(sids),length(a_pers_n)))
  samp <- samp %>% 
    ungroup() %>% 
    select(-c('SiteID','PropID','prop_area_km','county_area_km'))%>% 
    left_join(dat_rem %>% 
                ungroup() %>% 
                select(prop_rem_idx,SiteID,PropID,
                       prop_area_km,county_area_km) %>% 
                distinct())
  
  samp$Date <- NA
  samp$tot_rem <- samp$tot_hrs <- samp$effect_area_hrs <-
    samp$prop_property_impact <- samp$effect_area_km <- 0
  samp$pass_idx <- 1
  
  dat_rem <- rbind(dat_rem,samp)
  
  # dat_rem %>% group_by(prop_rem_idx) %>% 
  #   summarise(n=length(unique(prop_area_km))) %>% 
  #   filter(n>1)
  
  ##occupancy data-----------------------
  dat_occ <- sysbait_det_eff %>% 
    # mutate(site_occ_idx = as.numeric(SiteID)) %>% 
    group_by(period,SiteID) %>% 
    summarise(detections=sum(detection),
              nweeks = length(unique(subper_start)),
              per_start=min(subper_start),
              trap_nights_km = sum(trap_nights_km)
    ) %>% 
    rename(period_idx=period) %>% 
    arrange(period_idx,SiteID)
  
  dat_occ$trap_nights_km_sc <- scale(dat_occ$trap_nights_km)
  
  site_idx_lookup <- dat_rem %>% 
    ungroup() %>% 
    select(prop_rem_idx,SiteID) %>% 
    distinct() %>% 
    arrange(prop_rem_idx) %>% 
    drop_na
  
  #fill in missing SiteIDs from added aerial 0's
  dat_rem$SiteID[is.na(dat_rem$SiteID)] <- 
    unlist(sapply(1:length(dat_rem$SiteID[is.na(dat_rem$SiteID)]), function(x){
      site_idx_lookup$SiteID[site_idx_lookup$prop_rem_idx==
                               dat_rem$prop_rem_idx[is.na(dat_rem$SiteID)][x]]
    }))
  
  # #check no duplicated values from adding dummy data
  # dat_rem %>% group_by(prop_rem_idx,period_idx,Method,pass_idx) %>%
  #   summarise(n=n()) %>%
  #   filter(n>1)
  
  #county indices
  dat_rem$county_idx <- as.numeric(as.factor(dat_rem$SiteID))
  
  county_idx <- dat_rem %>% 
    ungroup() %>% 
    select(prop_rem_idx,county_idx,SiteID,prop_area_km,county_area_km) %>% 
    distinct() %>% 
    arrange(prop_rem_idx) %>% 
    mutate(property_prop_county=prop_area_km/county_area_km)
  
  ## map occupancy counties to removal counties
  rem_county_to_occ_site <- dat_occ %>% 
    ungroup() %>% 
    select(SiteID) %>% 
    distinct() %>% 
    left_join(county_idx %>% 
                select(county_idx,SiteID) %>% 
                distinct()) %>% 
    drop_na()
  
  ##covariates -----------------------------
  nbeta <- 3
  nbeta_lam <- 2
  
  nlcd_siteid_orig <- nlcd_siteid
  # nlcd_siteid$site_idx <- as.numeric(as.factor(nlcd_siteid$SiteID))
  
  develop <- data.frame(develop=nlcd_siteid$developed)
  develop$develop_sc <- scale(develop$develop)[,1]
  
  agri <- data.frame(agri=nlcd_siteid$agri)
  agri$agri_sc <- scale(agri$agri)[,1]
  
  ##indices for nimble ---------------------------------
  nsites_occ <- length(unique(dat_occ$SiteID))
  nproperties <- length(unique(dat_rem$prop_rem_idx))
  ncounties <- length(unique(dat_rem$SiteID))
  nperiods <- length(unique(dat_occ$period_idx))
  
  ## removal data-------------------
  #for instances of no aerial ops, set effort area to small number to avoid NaN in model division
  dat_rem$effect_area_km[dat_rem$effect_area_km==0 & dat_rem$tot_hrs==0] <- 0.001
  
  ##split by removal type -------------------
  dat_aerial <- dat_rem %>% filter(Method=="Aerial")
  dat_trap <- dat_rem %>% filter(Method=="Trap")
  
  gamma_a <- array(0,dim=c(max(dat_aerial$pass_idx),nproperties,nperiods))
  gamma_t <- array(0,dim=c(max(dat_trap$pass_idx),nproperties,nperiods))
  for(i in 1:nproperties){
    for(t in 1:nperiods){
      for(j in 1:max(dat_aerial$pass_idx,dat_trap$pass_idx)){
        x <- dat_aerial$prop_property_impact[dat_aerial$pass_idx==j &
                                               dat_aerial$prop_rem_idx==i &
                                               dat_aerial$period_idx==t]
        z <- dat_trap$prop_property_impact[dat_trap$pass_idx==j &
                                             dat_trap$prop_rem_idx==i &
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
  
  ## nfsp covariate ------------------
  if(!exists("nfsp")){
    if(file.exists(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Model Ready Data/NFSP County Overlap/nfsp_lag_2020_",max(year(sysbait_det_eff$subper_start)),".RData"))){
      load(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Model Ready Data/NFSP County Overlap/nfsp_lag_2020_",max(year(sysbait_det_eff$subper_start)),".RData"))
    } else {
      nfsp <- matrix(NA,nsites_occ,nperiods)
      for(i in 1:nsites_occ){
        for(t in 1:nperiods){
          nfsp[i,t] <- unique(sysbait_det_eff$prp_nfs_lag[as.numeric(sysbait_det_eff$site_idx)==i & sysbait_det_eff$period==t])
        }
      }
      nfsp_sc <- scale(nfsp)
      save(nfsp,nfsp_sc,file=paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/Model Ready Data/NFSP County Overlap/nfsp_lag_2020_",max(year(sysbait_det_eff$subper_start)),".RData"))
    }
  }

  ##removal matrix ------------------
  remtot <- dat_rem %>% group_by(prop_rem_idx,period_idx) %>%
    reframe(removal=sum(tot_rem))
  
  yrem <- matrix(0,nproperties,nperiods)
  prop_idx <- sort(unique(dat_rem$prop_rem_idx))
  for(j in 1:nproperties){
    for(t in 1:nperiods){
      x <- remtot$removal[remtot$prop_rem_idx==prop_idx[j] & remtot$period_idx==t]
      if(length(x)>0){
        yrem[j,t] <- x
      }}}
  
  remtot_counties <- dat_rem %>% 
    group_by(SiteID,period_idx) %>%
    reframe(removal=sum(tot_rem))
  count_idx <- sort(unique(dat_rem$SiteID))
  yrem_counties<- matrix(0,ncounties,nperiods)
  for(j in 1:ncounties){
    for(t in 1:nperiods){
      x <- remtot_counties$removal[remtot_counties$SiteID==count_idx[j] & remtot_counties$period_idx==t]
      if(length(x)>0){
        yrem_counties[j,t] <- x
      }}}
  
  #mean trapping effort  -------------------------------
  mn_te <- (((0.5/2.59)*45)-attr(dat_occ$trap_nights_km_sc,"scaled:center"))/
    attr(dat_occ$trap_nights_km_sc,"scaled:scale") 
  
  # fit model ---------------------------
  # range(rowSums(yrem))
  # hist(rnbinom(10000,size=1000,prob=0.1))
  # hist(rgamma(10000,1000,0.1))
  # hist(rbeta(10000,1,2))

  ## model specification -------------------------------------------------
  ZIbinomcode <- nimbleCode({
    
    #priors
    
    #occupancy regression coefficients
    for(i in 1:nbeta){
      beta[i] ~ dlogis(0,1)
    }
    
    #lambda regression coefficients
    for(i in 1:nbeta_lam){
      beta_lam[i] ~ dnorm(0,5)
    }
    
    #seasonal random effects for occupancy
    for(t in 1:nperiods){
      beta0[t] ~ dnorm(beta[1],sd=sd_beta0)
    }
    sd_beta0 ~ dgamma(10,10)
    
    #detection coefficients
    for(i in 1:2){
      alpha[i] ~ dlogis(0,1)
      delta_t[i] ~ dlogis(0,1)
      delta_a[i] ~ dlogis(0,1)
    }
    sd_pdet ~ dgamma(10,10)
    sd_theta_t ~ dgamma(10,10)
    sd_theta_a ~ dgamma(10,10)
    
    #initial abundance
    p_n1 ~ dbeta(1,2)
    r_n1 ~ dgamma(800,0.1) 
    
    #lambda variance 
    sd_lam ~ dgamma(1,20)
    
    #mean effort detection probability
    mu_pdet_mn <- alpha[1] + alpha[2]*mn_te
    logit(p_sys) ~ dnorm(mu_pdet_mn,sd=sd_pdet)
    
    #county process
    for(w in 1:ncounties){
      for(t in 1:nperiods){
        #county-level lambda
        mu_lam[w,t] <- beta_lam[1] + beta_lam[2]*nfsp[w,t]
        log(lambda[w,t]) ~ dnorm(mu_lam[w,t],sd=sd_lam)
      }
    
    }

    #property process
    for(i in 1:nproperties){
      N_prop[i,1] ~ dnegbin(size=r_n1,prob=p_n1)
      
      for(t in 2:nperiods){
        mu_intens[i,t] <- lambda[county_idx_rem[i],t-1] * 
                              (N_prop[i,t-1] - yrem_mat[i,t-1])
        N_prop[i,t] ~ dpois(mu_intens[i,t])
      }
      
      for(t in 1:nperiods){
        N_prop_latent[i,t] <- N_prop[i,t] * z_site[site_idx_occ_rem[i],t]
      }
    }

    #county-level occupancy process
    for(k in 1:nsites_occ){ 
      for(t in 1:nperiods){
      # logit(psi[k,t]) <- beta[1] + beta[2]*develop[k] + beta[3]*agri[k]
        logit(psi[k,t]) <- beta0[t] + beta[2]*develop[k] + beta[3]*agri[k]
        
      #occupancy/detection
        pabs[k,t] <- 1-psi[k,t]
        z_site[k,t] ~ dbern(psi[k,t])

        #effective detection given effort
        pstar_site[k,t] <- 1-(1-p_sys)^eff_weeks
        
        #p(elimination) given effort
        pelim[k,t] <- pabs[k,t]/(pabs[k,t] + psi[k,t]*(1-pstar_site[k,t]))
      }
    }

    #abundance
    ##aerial
    for(i in 1:nsamp_aerial){
      mu_theta_a[i] <- delta_a[1] + delta_a[2]*eff_area_a[i]
      logit(theta_a[pass_idx_a[i],prop_idx_a[i],period_idx_a[i]]) ~ 
        dnorm(mu_theta_a[i],sd=sd_theta_a)
      
      pip_a[i] <-
        avail_fun(pass=pass_idx_a[i],
                  gamma=gamma_a[pass_idx_a[i],prop_idx_a[i],period_idx_a[i]],
                  gamma_past=gamma_a[1:(pass_idx_a[i]-1),prop_idx_a[i],period_idx_a[i]],
                  theta=theta_a[pass_idx_a[i],prop_idx_a[i],period_idx_a[i]],
                  theta_past=theta_a[1:(pass_idx_a[i]-1),prop_idx_a[i],period_idx_a[i]])
      
      yrem_aerial[i] ~ dbinom(size=N_prop_latent[prop_idx_a[i],period_idx_a[i]], prob=pip_a[i])
      yrem_pred_aerial[i] ~ dbinom(size=N_prop_latent[prop_idx_a[i],period_idx_a[i]], prob=pip_a[i])
    }
    
    ##traps
    for(i in 1:nsamp_trap){
      mu_theta_t[i] <- delta_t[1] + delta_t[2]*eff_area_t[i]
      logit(theta_t[pass_idx_t[i],prop_idx_t[i],period_idx_t[i]]) ~
        dnorm(mu_theta_t[i],sd=sd_theta_t)
      
      pip_t[i] <-
        avail_fun(pass=pass_idx_t[i],
                  gamma=gamma_t[pass_idx_t[i],prop_idx_t[i],period_idx_t[i]],
                  gamma_past=gamma_t[1:(pass_idx_t[i]-1),prop_idx_t[i],period_idx_t[i]],
                  theta=theta_t[pass_idx_t[i],prop_idx_t[i],period_idx_t[i]],
                  theta_past=theta_t[1:(pass_idx_t[i]-1),prop_idx_t[i],period_idx_t[i]])
      
      #abundance likelihood
      yrem_trap[i] ~ dbinom(size=N_prop_latent[prop_idx_t[i],period_idx_t[i]], prob=pip_t[i])
      yrem_pred_trap[i] ~ dbinom(size=N_prop_latent[prop_idx_t[i],period_idx_t[i]], prob=pip_t[i])
    }
    
    #occupancy
    for(i in 1:nsamp_occ){
      #detection probability
      mu_pdet[i] <- alpha[1] + alpha[2]*trap_nights_km[i]
      logit(pdet[i]) ~ dnorm(mu_pdet[i],sd=sd_pdet)
      
      #occupancy likelihood
      pocc[i] <- pdet[i]*z_site[site_occ_idx[i],period_idx_occ[i]]
      yocc[i] ~ dbinom(prob=pocc[i],size=nweeks[i])
      yocc_pred[i] ~ dbinom(prob=pocc[i],size=nweeks[i])
    }
  })#nimblecode end


  ## nimble lists  -------------------------------------------------
  modDat <- list(yrem_trap=dat_trap$tot_rem,
                 # prop_property_county=county_idx$property_prop_county,
                 eff_area_t=dat_trap$effect_area_hrs_sc[,1],
                 gamma_t=gamma_t,
                 yrem_aerial=dat_aerial$tot_rem,
                 eff_area_a=dat_aerial$effect_area_hrs_sc[,1],
                 gamma_a=gamma_a,
                 yocc=dat_occ$detections,
                 nweeks=dat_occ$nweeks,
                 trap_nights_km=dat_occ$trap_nights_km_sc[,1],
                 yrem_mat=yrem
                 # yrem_mat_counties=yrem_counties
                 )
  
  
  ### constants -----------------------------
  const <- list(develop=develop$develop_sc,
                agri=agri$agri_sc,
                nfsp=nfsp_sc,
                nbeta=nbeta,
                nbeta_lam=nbeta_lam,
                nperiods=nperiods,
                nsamp_occ=nrow(dat_occ),
                nsites_occ=nsites_occ,
                site_occ_idx=dat_occ$SiteID,
                # prop_idx_occ_rem=site_idx_lookup$prop_rem_idx, #indexes to property site idx over occ idx
                site_idx_occ_rem=site_idx_lookup$SiteID, #indexes to occupancy site idx over rem idx
                # rem_county_to_occ_site=rem_county_to_occ_site$SiteID,
                county_idx_rem=county_idx$county_idx,
                period_idx_occ=dat_occ$period_idx,
                nproperties=nproperties,
                ncounties=ncounties,
                mn_te=mn_te,
                eff_weeks=eff_weeks,
                prop_idx_t=dat_trap$prop_rem_idx,
                # site_idx_t=dat_trap$county_idx,
                prop_idx_a=dat_aerial$prop_rem_idx,
                # site_idx_a=dat_aerial$county_idx,
                period_idx_t=dat_trap$period_idx,
                period_idx_a=dat_aerial$period_idx,
                pass_idx_t=dat_trap$pass_idx,
                pass_idx_a=dat_aerial$pass_idx,
                nsamp_trap=nrow(dat_trap),
                nsamp_aerial=nrow(dat_aerial)
  )
  const <- const[!sapply(const,is.null)]
  
  ### initial values -----------------------------
  inits <- list(beta=rnorm(nbeta,c(0,0.5),0.1),
                beta_lam=rnorm(2,c(0,0.001),0.001),
                lambda=matrix(rnorm(ncounties*nperiods,1,0.05),
                              ncounties,nperiods),
                sd_lam=rnorm(1,0.3,0.01),
                z_site = matrix(1,nsites_occ,nperiods),
                N_prop = matrix(NA,nproperties,nperiods),
                p_n1 = rnorm(1,0.1,0.01),
                r_n1 = rnorm(1,800,50)
  )
  inits <- inits[!sapply(inits,is.null)]
  
  inits$N_prop[,1] <- apply(yrem,1,max)*50
  # inits$N_site[,1] <- apply(yrem_counties,1,max)*1000
  
  ## nimble configuration -------------------------------------------------
  mod <- nimbleModel(code = ZIbinomcode,
                     data = modDat,
                     constants = const,
                     inits = inits)
  # mod$initializeInfo()
  
  mcmc.conf <- configureMCMC(mod,enableWAIC=F)
  
  # if(is.na(monitors)){
  monitors <- c("lambda",
                "beta",
                "beta_lam",
                "beta0",
                "alpha",
                "sd_lam",
                "sd_pdet",
                "delta_t",
                "delta_a",
                "pabs",
                "pelim",
                "pocc",
                "p_n1",
                "r_n1",
                "N_prop"
                # "pip_a",
                # "pip_t",
                # "yrem_pred_aerial",
                # "yrem_pred_trap",
                # "yocc_pred"
  )
  # }
  mcmc.conf$setMonitors(monitors)
  
  Cmcmc <- buildMCMC(mcmc.conf,resetFuncitions=T) #uncompiled MCMC
  Cmod <- compileNimble(mod,Cmcmc) #compiled model
  
  ## run mcmc ----------------------------------- 
  mcmc.out <- runMCMC(Cmod$Cmcmc,
                      niter = 0,
                      thin= thin,
                      nburn = 0,
                      # setSeed = 1:nChains,
                      nchains = nChains)
  
  range(exp(rnorm(10000,c(0,0.01),Cmod$mod$sd_lam)))
  Cmod$mod$N_prop
  Cmod$mod$N_prop[which(is.nan(rowSums(Cmod$mod$N_prop))),]
  hist(Cmod$mod$lambda)
  Cmod$mod$beta_lam
  # yrem[which(is.nan(rowSums(Cmod$mod$N))),]
  # Cmod$mod$p_n1 <- inits$p_n1
  # Cmod$mod$r_n1 <- inits$r_n1

  niter<-50000
  burnProp<-0.5
  thin <-5
  
  samples <- runMCMC(Cmod$Cmcmc,
                     niter = niter,
                     thin= thin,
                     nburn = niter*burnProp,
                     # setSeed = 1:nChains,
                     setSeed=1,
                     nchains = nChains)
  
  samples <- as.data.frame(samples)
  
  #multichain
  if(nChains>1){
    samples1 <- mcmc.out$chain1
    samples2 <- mcmc.out$chain2
    samples3 <- mcmc.out$chain3
    samples <- rbind(samples1,samples2,samples3)
  }
  
  list(samples=samples,
       if(nChains>1){
         samples1=samples1
       },
       if(nChains>1){
         samples2=samples2
       },
       if(nChains>1){
         samples3=samples3
       },
       nimbleMod=ZIbinomcode,
       nbeta=nbeta,
       nperiods=nperiods,
       nsites_occ=nsites_occ,
       nproperties=nproperties,
       mn_te=mn_te,
       eff_weeks=eff_weeks,
       dat_occ=dat_occ,
       dat_trap=dat_trap,
       dat_aerial=dat_aerial,
       develop=develop,
       agri=agri,
       site_idx_lookup=site_idx_lookup
  )
}

