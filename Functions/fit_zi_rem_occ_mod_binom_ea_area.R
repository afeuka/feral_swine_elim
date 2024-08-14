### Title: Hurdle model (removal) for southern MO feral swine occupancy/abundance
### Author: Abbey Feuka
### Date: 27NOV2023 edited 22FEB2024
### Notes: 

# monitors = NA
# niter=10000
# burnProp=0.1
# nChains=1
# thin=1
# elim_prob=0.95
# eff_weeks=10
# subset_data=F

fit_zi_rem_occ <- function(sysbait_det_eff, #output of data_functions_ws_occ_ea_rem.R
                           rem_eff_ea, #output of data_functions_ws_occ_ea_rem.R
                           study_site_grid, #output from mapping huc10 to eas.R
                           elim_prob=0.95,#threshold for probability of elimination
                           eff_weeks=10,#number of trap nights to determine eliminiation probability
                           monitors, #parameters to save mcmc samples for
                           niter, #number of mcmc iterations
                           thin=1,#thinning interval
                           burnProp,
                           nChains,
                           subset_data=F#subset data for testing
){
  
  require(nimble)
  require(tidyverse)
  require(sf)
  
  # data setup -------------------------------------------------
  elim_areas <- study_site_grid %>% 
    group_by(Area_Name) %>% 
    summarise(geometry=st_union(geometry)) %>% 
    mutate(area_m=st_area(.),
           area_km=area_m/1e6,
           elim_area_idx= as.numeric(Area_Name)+1)
  
  sysbait_det_eff$Area_Name[is.na(sysbait_det_eff$Area_Name)] <- 0
  sysbait_det_eff$elim_area_idx <- as.numeric(sysbait_det_eff$Area_Name)+1

  rem_eff_ea$Area_Name[is.na(rem_eff_ea$Area_Name)]<- 0
  rem_eff_ea$elim_area_idx <- as.numeric(rem_eff_ea$Area_Name)+1
  
  rem_eff_ea <- rem_eff_ea %>% filter(!is.na(prop_ea_impact)) %>% 
    filter(!(method=="Ground" & is.na(eff_area_events))) %>% 
    filter(!(method%in%c("Aerial","Trap") & is.na(eff_area_hrs))) %>% 
    filter(!(method=="Aerial" & eff_area_hrs>0.15)) %>% 
    rename(period_idx = period) %>% 
    group_by(period_idx,method,Area_Name) %>% 
    mutate(pass_idx=1:n()) 
  
  ##add 0 data for later periods -aerial -------------
  a_pers <- rem_eff_ea %>% ungroup() %>%  
    filter(method=="Aerial") %>% 
    select(period_idx) %>% distinct
  a_pers_n <- which(!(1:max(rem_eff_ea$period_idx)%in% a_pers$period_idx))
  eids <- unique(elim_areas$elim_area_idx)
  
  samp <- rem_eff_ea[1:(length(a_pers_n)*length(eids)),]
  samp$method <- "Aerial"
  samp$period_idx <- rep(a_pers_n,length(eids))
  samp$elim_area_idx <- sort(rep(1:length(eids),length(a_pers_n)))
  samp$Area_Name <-samp$Date <- NA
  samp$tot_rem <- samp$tot_hrs <- samp$num_events <- samp$eff_area_events <- samp$eff_area_hrs <-
    samp$prop_ea_impact <- samp$effect_area_km <- 0
  samp$pass_idx <- 1
  rem_eff_ea <- rbind(rem_eff_ea,samp)

  ##occupancy data-----------------------
  dat_occ <- sysbait_det_eff %>% 
    left_join(study_site_grid %>% 
                mutate(site_idx=as.character(SiteID)) %>% 
                st_drop_geometry() %>% 
                dplyr::select(-c(SiteID))) %>% 
    group_by(period,site_idx) %>% 
    mutate(site_idx = as.numeric(site_idx)) %>% 
    summarise(detections=sum(detection),
              nweeks = length(unique(subper_start)),
              per_start=min(subper_start),
              trap_nights_km = sum(trap_nights_km),
              elim_area = max(as.numeric(Area_Name)),
              elim_area_idx = max(elim_area_idx)
              ) %>% 
    rename(period_idx=period) %>% 
    arrange(period_idx,site_idx)

  elim_site_idx <- dat_occ %>% group_by(site_idx) %>% 
     summarise(elim_area_idx=as.numeric(unique(elim_area_idx))) %>% 
     arrange(site_idx)

  ##covariates -----------------------------
  nbeta <- 4
  
  nlcd_siteid_orig <- nlcd_siteid
  nlcd_siteid$site_idx <- as.numeric(as.factor(nlcd_siteid$SiteID))
  
  develop <- data.frame(develop=nlcd_siteid$developed)
  develop$develop_sc <- scale(develop$develop)[,1]
  
  agri <- data.frame(agri=nlcd_siteid$agri)
  agri$agri_sc <- scale(agri$agri)[,1]
  
  ##indices for nimble ---------------------------------
  nsites <- length(unique(dat_occ$site_idx))
  nea <- length(unique(dat_occ$elim_area))
  nperiods <- length(unique(dat_occ$period_idx))

  ## removal data-------------------
  ##split by removal type -------------------
  dat_aerial <- rem_eff_ea %>% filter(method=="Aerial")
  dat_trap <- rem_eff_ea %>% filter(method=="Trap")

  gamma_a <- array(0,dim=c(max(dat_aerial$pass_idx),nea,nperiods))
  gamma_t <- array(0,dim=c(max(dat_trap$pass_idx),nea,nperiods))
  for(i in 1:nea){
    for(t in 1:nperiods){
      for(j in 1:max(dat_aerial$pass_idx,dat_trap$pass_idx)){
        x <- dat_aerial$prop_ea_impact[dat_aerial$pass_idx==j &
                                           dat_aerial$elim_area_idx==i &
                                           dat_aerial$period_idx==t]
        z <- dat_trap$prop_ea_impact[dat_trap$pass_idx==j &
                                         dat_trap$elim_area_idx==i &
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

  remtot <- rem_eff_ea %>% group_by(elim_area_idx,period_idx) %>%
    reframe(removal=sum(tot_rem))

  nfsp <- matrix(NA,nsites,nperiods)
  for(i in 1:nsites){
    for(t in 1:nperiods){
      nfsp[i,t] <- unique(sysbait_det_eff$prp_nfs[as.numeric(sysbait_det_eff$site_idx)==i & sysbait_det_eff$period==t])
    }
  }
  nfsp_sc <- scale(nfsp)
  
  ##covariate correlation-------------------
  # samp <- sysbait_det_eff %>% left_join(nlcd_siteid %>% mutate(site_idx=as.character(site_idx)) %>% select(-geometry))
  # samp <- samp %>% select(agri,developed,prp_nfs) %>% distinct()
  # cor(samp$agri,samp$prp_nfs)
  # cor(samp$developed,samp$prp_nfs)
  
  ##removal matrix ------------------
  yrem <- matrix(0,nea,nperiods)
  ea_idx <- sort(unique(rem_eff_ea$elim_area_idx))
  for(j in 1:nea){
    for(t in 1:nperiods){
      x <- remtot$removal[remtot$elim_area_idx==ea_idx[j] & remtot$period_idx==t]
      if(length(x)>0){
        yrem[j,t] <- x
      }
    }
   }

  #one trap per mi2 -------------------------------
  mn_te <- 1/2.59

  if(subset_data){
    samp_idx <- sample(unique(dat_occ$site_idx),20,replace=F)

    nsites <- length(samp_idx)
    
    dat_occ <- dat_occ %>% filter(site_idx%in%samp_idx) %>% 
      mutate(site_idx = as.numeric(as.factor(site_idx))) 
    dat_rem <- dat_rem %>% filter(site_idx%in%samp_idx)%>% 
      mutate(site_idx = as.numeric(as.factor(site_idx))) 
    
    nfsp_sc <- nfsp_sc[samp_idx,]
    yrem <- yrem[samp_idx,]
    
    develop <- develop[samp_idx,]
    agri <- agri[samp_idx,]
  }
  
  # fit model ---------------------------
  ## model specification -------------------------------------------------
  ZIbinomcode <- nimbleCode({
    
    for(i in 1:nbeta){
      beta[i] ~ dlogis(0,1)
    }
    
    for(i in 1:2){
      alpha[i] ~ dlogis(0,1)
    }
    sd_pdet ~ dgamma(10,10)
    
    p_n1 ~ dbeta(0.1,7)
    r_n1 ~ dgamma(1,0.1)
    mu_lam ~ dnorm(0,0.15)
    sd_lam ~ dgamma(1,10)
    
    for(i in 1:nea){
      log(lambda[i]) ~ dnorm(mu_lam,sd=sd_lam)
    }
    
    mu_p ~ dnorm(0,1)
    tau_p ~ dgamma(10,10)
    sd_p <- sqrt(1/tau_p)
    
    logit(p_trap) ~ dnorm(mu_p,sd=sd_p)
    logit(p_aerial) ~ dnorm(mu_p,sd=sd_p)
    
    for(k in 1:nsites){ #watershed loop
      for(t in 1:nperiods){
        logit(psi[k,t]) <- beta[1] + beta[2]*nfsp[k,t] + beta[3]*develop[k] + beta[4]*agri[k]
        # logit(psi[k,t]) ~ dnorm(mu_psi[k,t],sd=sd_psi)
        
        #occupancy/detection
        pabs[k,t] <- 1-psi[k,t]
        ea_psi[k,t] <- psi[k,t]*z_ea[ea_site_idx[k],t]
        z_site[k,t] ~ dbern(ea_psi[k,t])
        
        #number of trap-nights required for 95% certainty elim 
        eff_elim[k,t] <- log((pabs[k,t]*(1-elim_prob))/(elim_prob*psi[k,t]))/
                            (log(1-p_sys)) 
        
        #effective detection given efforet
        pstar_site[k,t] <- 1-((1-p_sys)^mn_te)^eff_weeks

        #p(elimination) given effort
        pelim[k,t] <- pabs[k,t]/(pabs[k,t] + psi[k,t]*(1-pstar_site[k,t]))
      }
    }

    for(i in 1:nea){ #elimination area loop
      #initalize N
      N[i,1] ~ dnegbin(size=r_n1,prob=p_n1)
      
      #hierarchical z
      psi_ea[i] ~ dbeta(1,1)
      for(t in 1:nperiods){
        z_ea[i,t] ~ dbern(psi_ea[i])
        N_latent[i,t] <- N[i,t] * z_ea[i,t]
      }
      
      #abundace/removal process
      for(t in 2:nperiods){
        mu_nb[i,t] <- lambda[i] * (N[i,t-1] - yrem_mat[i,t-1])
        # p_nb[i,t] <- mu_nb[i,t]/sig2_nb[i]
        # r_nb[i,t] <- (mu_nb[i,t]^2)/(sig2_nb[i]-mu_nb[i,t])
        # N[i,t] ~ dnegbin(prob=p_nb[i,t],size=r_nb[i,t])
        N[i,t] ~ dpois(mu_nb[i,t])
        # N_change[i,t] <- (N[i,t]-N[i,t-1])/N[i,t-1]
      }
    }
    
    #abundance
    ##aerial
    for(i in 1:nsamp_aerial){
      #effort imputation for helicopter data
      theta_a[pass_idx_a[i],site_idx_a[i],period_idx_a[i]] <-
        1-pow((1-p_aerial),eff_area_a[i])
      
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
      theta_t[pass_idx_t[i],site_idx_t[i],period_idx_t[i]] <-
        1-(1-p_trap)^eff_area_t[i]
      
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
      # pdet[i] <- 1-(1-p_sys)^trap_nights_km[i]
      mu_pdet[i] <- alpha[1] + alpha[2]*trap_nights_km[i]
      logit(pdet[i]) ~ dnorm(mu_pdet[i],sd=sd_pdet)
      
      #occupancy likelihood
      pocc[i] <- pdet[i]*z_site[site_idx_occ[i],period_idx_occ[i]]
      yocc[i] ~ dbinom(prob=pocc[i],size=nweeks[i])
      yocc_pred[i] ~ dbinom(prob=pocc[i],size=nweeks[i])
    }
    
  })

  ## nimble lists  -------------------------------------------------
  modDat <- list(yrem_trap=dat_trap$tot_rem,
                 eff_area_t=dat_trap$eff_area_hrs,
                 gamma_t=gamma_t,
                 yrem_aerial=dat_aerial$tot_rem,
                 eff_area_a=dat_aerial$eff_area_hrs,
                 gamma_a=gamma_a,
                 yocc=dat_occ$detections,
                 nweeks=dat_occ$nweeks,
                 trap_nights_km=dat_occ$trap_nights_km,
                 yrem_mat=yrem)

  ### constants -----------------------------
  const <- list(develop=develop$develop_sc,
                agri=agri$agri_sc,
                nfsp=nfsp_sc,
                nbeta=nbeta,
                nsites=nsites,
                nperiods=nperiods,
                nsamp_occ=nrow(dat_occ),
                site_idx_occ=dat_occ$site_idx,
                period_idx_occ=dat_occ$period_idx,
                ea_site_idx=elim_site_idx$elim_area_idx,
                nea=nea,
                elim_prob=elim_prob,
                mn_te=mn_te,
                eff_weeks=eff_weeks,
                site_idx_t=dat_trap$elim_area_idx,
                site_idx_a=dat_aerial$elim_area_idx,
                period_idx_t=dat_trap$period_idx,
                period_idx_a=dat_aerial$period_idx,
                pass_idx_t=dat_trap$pass_idx,
                pass_idx_a=dat_aerial$pass_idx,
                nsamp_trap=nrow(dat_trap),
                nsamp_aerial=nrow(dat_aerial)
                
  )

  ### initial values -----------------------------
  inits <- list(beta=rnorm(nbeta,c(0,0.5),0.1),
                mu_lam=rnorm(1,0,0.1),
                sd_lam=rnorm(1,0.1,0.01),
                p_sys=rbeta(1,1,3),
                p_aerial=rbeta(1,2,1),
                p_trap=rbeta(1,2,1),
                N = matrix(NA,nea,nperiods),
                z_site = matrix(1,nsites,nperiods),
                p_n1 = rnorm(1,0.001,0.0001),
                r_n1 = rpois(1,7),
                z_ea = matrix(1,nea,nperiods),
                psi_ea = rnorm(nea,0.95,0.01)
  )
  inits$N[,1] <- apply(yrem,1,max)*25

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
                "alpha",
                "sd_pdet",
                "N",
                "pabs",
                "p_sys",
                "pelim",
                "eff_elim",
                "pocc",
                "p_n1",
                "r_n1",
                "N_latent",
                # "N_change",
                "p_trap",
                "p_aerial",
                "pip_a",
                "pip_t",
                "yrem_pred_aerial",
                "yrem_pred_trap",
                "yocc_pred")
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
  
  # Cmod$mod$sd_psi
  # Cmod$mod$N
  # Cmod$mod$lambda
  # Cmod$mod$calculate("cpue")
  # Cmod$mod$calculate("r_n1")
  # Cmod$mod$calculate("p_n1")

  # niter<-100000
  # burnProp<-0.75
  # thin <-5

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
       nsites=nsites,
       nea=nea,
       mn_te=mn_te,
       eff_weeks=eff_weeks,
       dat_occ=dat_occ,
       dat_trap=dat_trap,
       dat_aerial=dat_aerial,
       develop=develop,
       agri=agri
       )
}

