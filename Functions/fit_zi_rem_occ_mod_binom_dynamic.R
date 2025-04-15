### Title: Hurdle model (removal) for southern MO feral swine occupancy/abundance
### Author: Abbey Feuka
### Date: 27NOV2023 edited 22FEB2024
### Notes: 
# 
# monitors = NA
# niter=10000
# burnProp=0.1
# nChains=1
# thin=1
# elim_prob=0.95
# eff_weeks=10
# abund_scale="watersheds"
# dat_occ = sysbait_det_eff
# dat_rem = rem_eff_site


fit_zi_rem_occ <- function(dat_occ,
                           dat_rem, 
                           ssg_formatted, 
                           oak_siteid,
                           nlcd_siteid,
                           nfsp_siteid,
                           # elim_prob=0.95,#threshold for probability of elimination
                           eff_weeks=10,#number of trap nights to determine eliminiation probability
                           abund_scale,#watersheds or ea
                           monitors, #parameters to save mcmc samples for
                           niter, #number of mcmc iterations
                           thin=1,#thinning interval
                           burnProp,
                           nChains
){
  
  source("./R/functions/avail_fun.R")

  ##covariates -----------------------------
  nbeta <- 3
  nbeta_lam <- 4

  
  develop <- data.frame(develop=nlcd_siteid$developed)
  develop$develop_sc <- scale(develop$develop)[,1]
  
  agri <- data.frame(agri=nlcd_siteid$agri)
  agri$agri_sc <- scale(agri$agri)[,1]
  
  ##indices for nimble ---------------------------------
  nsites_occ <- length(unique(dat_occ$site_idx_occ))
  nsites_rem <- length(unique(dat_rem$site_idx_rem))
  nperiods <- length(unique(dat_occ$period_idx))
  
  ## removal data-------------------
  ##split by removal type -------------------
  dat_aerial <- dat_rem %>% filter(method=="Aerial")
  dat_trap <- dat_rem %>% filter(method=="Trap")
  
  gamma_a <- array(0,dim=c(max(dat_aerial$pass_idx),nsites_rem,nperiods))
  gamma_t <- array(0,dim=c(max(dat_trap$pass_idx),nsites_rem,nperiods))
  for(i in 1:nsites_rem){
    for(t in 1:nperiods){
      for(j in 1:max(dat_aerial$pass_idx,dat_trap$pass_idx)){
        x <- dat_aerial$prop_site_impact[dat_aerial$pass_idx==j &
                                           dat_aerial$site_idx_rem==i &
                                           dat_aerial$period_idx==t]
        z <- dat_trap$prop_site_impact[dat_trap$pass_idx==j &
                                         dat_trap$site_idx_rem==i &
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

  
  lat_long_siteid <- lat_long_siteid %>% filter(SiteID%in%site_idx_lookup$SiteID)
  oak_siteid <- oak_siteid %>% filter(SiteID%in%site_idx_lookup$SiteID)
  rem_covs <- rem_covs %>% left_join(oak_siteid)
  
  ##removal matrix ------------------
  remtot <- dat_rem %>% group_by(site_idx_rem,period_idx) %>%
    reframe(removal=sum(tot_rem))
  
  yrem <- matrix(0,nsites_rem,nperiods)
  site_idx <- sort(unique(dat_rem$site_idx_rem))
  for(j in 1:nsites_rem){
    for(t in 1:nperiods){
      x <- remtot$removal[remtot$site_idx_rem==site_idx[j] & remtot$period_idx==t]
      if(length(x)>0){
        yrem[j,t] <- x
      }}}

  #mean trapping effort  -------------------------------
  mn_te <- (((0.5/2.59)*45)-attr(dat_occ$trap_nights_km_sc,"scaled:center"))/
    attr(dat_occ$trap_nights_km_sc,"scaled:scale") 
  
  # fit model ---------------------------
  ## model specification -------------------------------------------------
  ZIbinomcode <- nimbleCode({
    #occupancy regression coefficients
    for(i in 1:nbeta){
      beta[i] ~ dlogis(0,1)
    }
    
    #lambda regression coefficients
    for(i in 1:nbeta_lam){
      beta_lam[i] ~ dnorm(0,5)
    }
    
    for(t in 1:nperiods){
      beta0[t] ~ dnorm(beta[1],sd=sd_beta0)
    }
    sd_beta0 ~ dgamma(10,10)

    for(i in 1:2){
      alpha[i] ~ dlogis(0,1)
      delta_t[i] ~ dlogis(0,1)
      delta_a[i] ~ dlogis(0,1)
    }
  
    sd_pdet ~ dgamma(10,10)
    sd_theta_t ~ dgamma(10,10)
    sd_theta_a ~ dgamma(10,10)
    sd_lam ~ dgamma(1,20)
    # lambda ~ dnorm(1,0.1)

    p_n1 ~ dbeta(0.2,12)
    r_n1 ~ dgamma(6,1) 

    #mean effort detection probability
    mu_pdet_mn <- alpha[1] + alpha[2]*mn_te
    logit(p_sys) ~ dnorm(mu_pdet_mn,sd=sd_pdet)


    for(i in 1:nsites_rem){
      N[i,1] ~ dnegbin(size=r_n1,prob=p_n1)
      
      for(t in 1:nperiods){
        N_latent[i,t] <- N[i,t] * z_site[site_idx_occ_rem[i],t]
      }
      
      #abundace/removal process
      mu_lam[i]<- beta_lam[1] + beta_lam[2]*oak[i] + 
        beta_lam[3]*y[i] + beta_lam[4]*x[i]
      log(lambda[i]) ~ dnorm(mu_lam[i],sd=sd_lam)
      
      for(t in 2:nperiods){
      
        mu_intens[i,t] <- lambda[i] * (N[i,t-1] - yrem_mat[i,t-1])
        N[i,t] ~ dpois(mu_intens[i,t])
      }
    }

    for(k in 1:nsites_occ){ #watershed loop
      for(t in 1:nperiods){
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
      logit(theta_a[pass_idx_a[i],site_idx_a[i],period_idx_a[i]]) ~ 
        dnorm(mu_theta_a[i],sd=sd_theta_a)
      
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
      mu_theta_t[i] <- delta_t[1] + delta_t[2]*eff_area_t[i]
      logit(theta_t[pass_idx_t[i],site_idx_t[i],period_idx_t[i]]) ~
        dnorm(mu_theta_t[i],sd=sd_theta_t)
      
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
      mu_pdet[i] <- alpha[1] + alpha[2]*trap_nights_km[i]
      logit(pdet[i]) ~ dnorm(mu_pdet[i],sd=sd_pdet)
      
      #occupancy likelihood
      pocc[i] <- pdet[i]*z_site[site_idx_occ[i],period_idx_occ[i]]
      yocc[i] ~ dbinom(prob=pocc[i],size=nweeks[i])
      yocc_pred[i] ~ dbinom(prob=pocc[i],size=nweeks[i])
    }
  })#nimblecode end
  
  ## nimble lists  -------------------------------------------------
  modDat <- list(yrem_trap=dat_trap$tot_rem,
                 eff_area_t=dat_trap$effect_area_hrs_sc[,1],
                 gamma_t=gamma_t,
                 yrem_aerial=dat_aerial$tot_rem,
                 eff_area_a=dat_aerial$effect_area_hrs_sc[,1],
                 gamma_a=gamma_a,
                 yocc=dat_occ$detections,
                 nweeks=dat_occ$nweeks,
                 trap_nights_km=dat_occ$trap_nights_km_sc[,1],
                 yrem_mat=yrem,
                 oak=rem_covs$prop_oak_sc[,1],
                 y=rem_covs$lat_sc[,1],
                 x=rem_covs$long_sc[,1]
                 )
  ### constants -----------------------------
  const <- list(develop=develop$develop_sc,
                agri=agri$agri_sc,
                nbeta=nbeta,
                nbeta_lam=nbeta_lam,
                nperiods=nperiods,
                nsamp_occ=nrow(dat_occ),
                nsites_occ=nsites_occ,
                site_idx_occ=dat_occ$site_idx_occ,
                site_idx_occ_rem=site_idx_lookup$site_idx_occ, #indexes to occupancy site idx over rem idx
                period_idx_occ=dat_occ$period_idx,
                nsites_rem=nsites_rem,
                mn_te=mn_te,
                eff_weeks=eff_weeks,
                site_idx_t=dat_trap$site_idx_rem,
                site_idx_a=dat_aerial$site_idx_rem,
                period_idx_t=dat_trap$period_idx,
                period_idx_a=dat_aerial$period_idx,
                pass_idx_t=dat_trap$pass_idx,
                pass_idx_a=dat_aerial$pass_idx,
                nsamp_trap=nrow(dat_trap),
                nsamp_aerial=nrow(dat_aerial)
  )
  
  # const <- const[!sapply(const,is.null)]
  
  # if(abund_scale=="watersheds"){
  #   names(const)[names(const)==""] <- c("site_idx_occ_rem","nsites_rem","site_idx_t","site_idx_a")
  # } else {
  #   names(const)[names(const)==""] <- c("ea_site_idx","elim_idx_rem","nea","nea_rem","site_idx_t","site_idx_a")
  # }

  ### initial values -----------------------------
  inits <- list(#beta=rnorm(nbeta,c(0,0.5),0.1),
                beta=rnorm(nbeta,c(0,0.5),0.1),
                beta_lam=rnorm(nbeta_lam,c(-0.05,0.01,-0.02,0),0.005),
                sd_lam=rnorm(1,0.005,0.0001),
                N = matrix(NA,nsites_rem,nperiods),
                z_site = matrix(1,nsites_occ,nperiods),
                p_n1 = rnorm(1,0.005,0.0001),
                r_n1 = rnorm(1,6,0.1),
                theta_a = array(0,dim=dim(gamma_a)),
                theta_t = array(0,dim=dim(gamma_t))
  )
  
  # inits <- inits[!sapply(inits,is.null)]
  
  # if(abund_scale=="watersheds"){
  #   names(inits)[names(inits)==""] <- c("N","p_n1")
  # } else {
  #   names(inits)[names(inits)==""] <- c("N","p_n1","z_ea","psi_ea")
  # }
  # 
  inits$N[,1] <- apply(yrem,1,max)*20

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
                "sd_lam",
                "beta0",
                "alpha",
                "sd_pdet",
                "delta_t",
                "delta_a",
                "sd_theta_t",
                "sd_theta_a",
                "N",
                "pabs",
                # "pelim",
                "pocc",
                "p_n1",
                "r_n1",
                "N_latent",
                "pip_a",
                "pip_t",
                "yrem_pred_aerial",
                "yrem_pred_trap",
                "yocc_pred"
  )
  # }
  mcmc.conf$setMonitors(monitors)
  
  Cmcmc <- buildMCMC(mcmc.conf,resetFuncitions=T) #uncompiled MCMC
  Cmod <- compileNimble(mod,Cmcmc) #compiled model
  
  ## run mcmc ----------------------------------- 
  # mcmc.out <- runMCMC(Cmod$Cmcmc,
  #                     niter = 0,
  #                     thin= thin,
  #                     nburn = 0,
  #                     nchains = nChains)
  
  samples <- runMCMC(Cmod$Cmcmc,
                     niter = niter,
                     thin= thin,
                     nburn = niter*burnProp,
                     setSeed=1,
                     nchains = nChains)
  
  samples <- as.data.frame(samples)
  
  samples
}

