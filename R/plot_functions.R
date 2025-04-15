require(sf)
require(tidyverse)

plot_trace <- function(samples,
                       dat_rem,
                       out_dir,
                       chain_idx,
                       site_idx_lookup){

  nChains <- max(chain_idx)
  nmcmc <- nrow(samples)/nChains
  # if(class(samples)=="list"){
  #   nChains <- length(samples)
  #   nmcmc <- nrow(samples[[1]])
  # } else {
  #   nChains <- 1
  #   nmcmc <- nrow(samples)
  # }
  
  # if(nChains>1){
  #   beta<- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
  #     cbind.data.frame(samples[[i]][,grepl("beta",colnames(samples[[i]])) &
  #                                     !grepl("beta_lam",colnames(samples[[i]]))],
  #                      chain=i)})),
  #     samp=rep(1:nmcmc,nChains))
  # } else {
    beta<- as.data.frame(samples[,grepl("beta",colnames(samples)) &
                                   !grepl("beta_lam",colnames(samples))])
    beta$chain <- sort(rep(1:nChains,nmcmc))
    beta$samp <- rep(1:nmcmc,nChains)
  # }
  
  beta_names <- c("Intercept","Development","Agriculture")
  nbeta <- length(beta_names)
  colnames(beta)[1:nbeta] <- beta_names
  
  beta_long <- beta %>% 
    pivot_longer(cols=all_of(1:nbeta),names_to="beta",values_to="value") %>% 
    mutate(chain=factor(chain))
  
  g_b<- ggplot(beta_long) + geom_line(aes(x=samp,y=value,col=chain))+
    geom_hline(yintercept=0,col="black",lty=2)+
    facet_wrap(.~beta)
  
  ggsave(g_b, 
         filename = file.path(out_dir,"plots","trace","beta_trace.jpeg"),
         width=7,height=5,units="in",device="jpeg")
  
  beta_sum <- beta_long %>% 
    group_by(beta) %>% 
    summarise(mn=mean(value),
              md=median(value),
              lci=quantile(value,0.025),
              uci=quantile(value,0.975))
  
  ##beta lam coefficients----------------------
  # if(nChains>1){
  #   beta_lam <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
  #     cbind.data.frame(samples[[i]][,grepl("beta_lam",colnames(samples[[i]]))],
  #                      chain=i)})),
  #     samp=rep(1:nmcmc,nChains))
  # } else {
    beta_lam <- as.data.frame(samples[,grepl("beta_lam",colnames(samples))])
    beta_lam$chain <- sort(rep(1:nChains,nmcmc))
    beta_lam$samp <- 1:nmcmc
  # }
  
  beta_lam_names <- c("Intercept","Oak Proportion","Latitiude","Longitude")
  
  colnames(beta_lam)[1:length(beta_lam_names)] <- beta_lam_names
  
  beta_lam_long <- beta_lam %>% 
    pivot_longer(cols=all_of(1:length(beta_lam_names)),
                 names_to="beta_lam",values_to="value") %>% 
    mutate(chain=factor(chain))
  
  g_bl<- ggplot(beta_lam_long) + geom_line(aes(x=samp,y=value,col=chain))+
    geom_hline(yintercept=0,col="black",lty=2)+
    facet_wrap(.~beta_lam)
  
  ggsave(g_bl, 
         filename = file.path(out_dir,"plots","trace","beta_lam_trace.jpeg"),
         width=7,height=5,units="in",device="jpeg")
  
  beta_lam_sum <- beta_lam_long %>% 
    group_by(beta_lam) %>% 
    summarise(mn=mean(value),
              md=median(value),
              lci=quantile(value,0.025),
              uci=quantile(value,0.975))
  
  ##initial n values ------------------
  # if(nChains>1){
  #   r_n1 <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
  #     cbind.data.frame(samples[[i]][,grep("r_n1",colnames(samples[[i]]))],
  #                      chain=i)})),
  #     samp=rep(1:nmcmc,nChains))
  #   
  #   p_n1 <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
  #     cbind.data.frame(samples[[i]][,grep("p_n1",colnames(samples[[i]]))],
  #                      chain=i)})),
  #     samp=rep(1:nmcmc,nChains))
  # } else {
    r_n1 <- as.data.frame(samples[,grep("r_n1",colnames(samples))])
    p_n1 <- as.data.frame(samples[,grep("p_n1",colnames(samples))])
    
    r_n1$chain <- sort(rep(1:nChains,nmcmc))
    p_n1$chain <- sort(rep(1:nChains,nmcmc))
    
    r_n1$samp <- 1:nmcmc
    p_n1$samp <- 1:nmcmc
  # }
  
  colnames(r_n1)[1] <- "r_n1"
  r_n1_long <- r_n1 %>% 
    mutate(chain=factor(chain))
  
  g_r <- ggplot(r_n1_long) + geom_line(aes(x=samp,y=r_n1,col=chain))
  
  colnames(p_n1)[1] <- "p_n1"
  p_n1_long <- p_n1 %>% 
    mutate(chain=factor(chain))
  
  g_p <-ggplot(p_n1_long) + geom_line(aes(x=samp,y=p_n1,col=chain))
  g_n1 <- gridExtra::grid.arrange(g_r,g_p)
  
  ggsave(g_n1, 
         filename = file.path(out_dir,"plots","trace","n1_params_trace.jpeg"),
         width=7,height=5,units="in",device="jpeg")
  
  ## effort removal ----------------------
  # if(nChains>1){
  #   delta <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
  #     cbind.data.frame(samples[[i]][,grepl("delta",colnames(samples[[i]]))],chain=i)})),
  #     samp=rep(1:nmcmc,nChains))
  #   
  #   alpha <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
  #     cbind.data.frame(samples[[i]][,grepl("alpha",colnames(samples[[i]]))],chain=i)})),
  #     samp=rep(1:nmcmc,nChains))
  # } else {
    delta <- as.data.frame(samples[,grepl("delta",colnames(samples))])
    alpha <- as.data.frame(samples[,grepl("alpha",colnames(samples))])
    
    delta$chain <- sort(rep(1:nChains,nmcmc))
    alpha$chain <- sort(rep(1:nChains,nmcmc))
    
    delta$samp <- 1:nmcmc
    alpha$samp <- 1:nmcmc
  # }
  
  cpue_names <- c("Aerial Intercept","Aerial Slope","Trap Intercept","Trap Slope")
  colnames(delta)[1:4] <- cpue_names
  
  delta_long <- delta %>% 
    pivot_longer(cols=all_of(1:4),names_to="delta",values_to="value") %>% 
    mutate(chain=factor(chain))
  
  g_a<- ggplot(delta_long %>% filter(grepl("Aerial",delta))) + 
    geom_line(aes(x=samp,y=value,col=chain))+
    ylab("Coefficient estimate") +
    facet_wrap(.~delta)
  g_t<- ggplot(delta_long %>% filter(grepl("Trap",delta))) + 
    geom_line(aes(x=samp,y=value,col=chain))+
    ylab("Coefficient estimate") +
    facet_wrap(.~delta)
  
  ## effort systematic baiting ------------------
  alpha_names <- c("Systematic Baiting Intercept","Systematic Baiting Slope")
  colnames(alpha)[1:2] <- alpha_names
  
  alpha_long <- alpha %>% 
    pivot_longer(cols=all_of(1:2),names_to="alpha",values_to="value") %>% 
    mutate(chain=factor(chain))
  
  g_s<- ggplot(alpha_long) + 
    geom_line(aes(x=samp,y=value,col=chain))+
    ylab("Coefficient estimate") +
    facet_wrap(.~alpha)
  
  g_alpha<- gridExtra::grid.arrange(g_a,g_t,g_s)
  
  ggsave(g_alpha, 
         filename = file.path(out_dir,"plots","trace","detection_trace.jpeg"),
         width=7,height=5,units="in",device="jpeg")
  
  ## lambda ----------------------
  
  # if(nChains>1){
  #   lambda_N <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
  #     cbind.data.frame(samples[[i]][grepl("lambda",colnames(samples[[i]]))],
  #                      chain=i)})),
  #     samp=rep(1:nmcmc,nChains))
  # } else {
    lambda_N <- samples[grepl("lambda",colnames(samples))]
    lambda_N$chain <- sort(rep(1:nChains,nmcmc))
    lambda_N$samp <- 1:nmcmc
  # }
  
  nsites_rem <- max(dat_rem$site_idx_rem)
  nperiods <- max(dat_rem$period_idx)
  
  lambda_long <- lambda_N %>% 
    pivot_longer(cols=1:(ncol(lambda_N)-2),
                 values_to = "lambda",
                 names_to="idx") %>% 
    select(-idx)
  lambda_long$site_idx_rem <- rep(1:nsites_rem,nmcmc*nChains)
  lambda_long <- lambda_long %>% 
    left_join(site_idx_lookup)
  
  lam_sum <- lambda_long %>% 
    group_by(site_idx_rem) %>% 
    summarise(mn=mean(lambda),
              md=median(lambda),
              lci=quantile(lambda,0.025),
              uci=quantile(lambda,0.975))
  
  list(
    beta_sum=beta_sum,
    beta_lam_sum=beta_lam_sum,
    lam_sum=lam_sum
  )
}

plot_pabs <- function(samples,
                      dat_occ,
                      ssg_formatted,
                      elim_areas,
                      out_dir){
  
  if(class(samples)=="list"){
    nChains <- length(samples)
    nmcmc <- nrow(samples[[1]])
  } else {
    nChains <- 1
    nmcmc <- nrow(samples)
  }
  
  ## map to study site -------------------
  nsites <- nrow(ssg_formatted)
  nperiods <- max(dat_occ$period_idx)
  
  if(nChains>1){
    pabs <- do.call("rbind",lapply(1:nChains,function(i){
      samples[[i]][grep("pabs",colnames(samples[[i]]))]}))
  } else {
    pabs <- samples[grep("pabs",colnames(samples))]
  }
  
  pabs_sum <- data.frame(idx=colnames(pabs),
                         mn=colMeans(pabs),
                         md=sapply(1:ncol(pabs),function(i)quantile(pabs[,i],prob=0.5)),
                         lci=sapply(1:ncol(pabs),function(i)quantile(pabs[,i],prob=0.025)),
                         uci=sapply(1:ncol(pabs),function(i)quantile(pabs[,i],prob=0.975)))
  
  pabs_sum$site_idx <- sapply(1:nrow(pabs_sum),function(x){
    as.numeric(gsub("\\[","",gsub("pabs","",unlist(str_split(pabs_sum[x,"idx"],","))[1])))})
  
  pabs_sum$period_idx <- sapply(1:nrow(pabs_sum),function(x){
    as.numeric(gsub("\\]","",unlist(str_split((pabs_sum[x,"idx"]),","))[2]))})
  
  pabs_sum$per_start <- sapply(1:nrow(pabs_sum),function(i){
    min(dat_occ$per_start[dat_occ$period_idx==pabs_sum$period_idx[i]])})
  pabs_sum$per_start <- as.Date(pabs_sum$per_start,origin='1970-01-01')
  pabs_sum$year <- year(pabs_sum$per_start)
  pabs_sum$fy <- pabs_sum$year
  pabs_sum$fy[month(pabs_sum$per_start)%in%c(10:12)] <- 
    pabs_sum$fy[month(pabs_sum$per_start)%in%c(10:12)] +1
  
  ## elimination area groupings ---------------------------
  pabs_long <- pabs %>% pivot_longer(all_of(1:ncol(pabs)))
  pabs_long$site_idx <- rep(rep(1:nsites,nperiods),nmcmc*nChains)
  pabs_long$period_idx <- rep(sort(rep(1:nperiods,nsites)),nmcmc*nChains)
  pabs_long$samp_idx <- sort(rep(1:(nmcmc*nChains),nsites*nperiods))
  pabs_long <- pabs_long %>%
    left_join(ssg_formatted %>%
                st_drop_geometry() %>%
                rename(site_idx=SiteID) %>%
                select(site_idx,Area_Name))
  
  ##fiscal year ----------------
  ###inidivdual maps------------------
  pabs_sum_fy <- pabs_sum %>% group_by(fy,site_idx) %>% 
    summarise(fy_mn = mean(mn),
              fy_md= median(mn),
              fy_lci = mean(lci),
              fy_uci= mean(uci)) %>% 
    left_join(ssg_formatted %>% rename(site_idx=SiteID)) %>% 
    st_as_sf()
  
  ### facet version ------------------
  ggplot()+
    geom_sf(data=pabs_sum_fy,aes(fill=fy_md))+
    geom_sf(data=elim_areas %>% filter(Area_Name!="0"),aes(col=Area_Name),
            fill="transparent",lwd=1.5)+
    facet_wrap(.~fy,ncol=1)+
    scale_fill_viridis_c(#limits=c(0,1),
      direction=-1,
      name="Feral swine \nelimination probability")+
    scale_color_discrete(name="Elimination Area")+
    theme(panel.background=element_blank(),
          axis.ticks=element_blank(),
          axis.text=element_blank(),
          plot.margin = margin(1,1,1,1, "cm"),
          text=element_text(size=15))
  
  ggsave(filename=file.path(out_dir,"plots","pabs_all_yrs.jpeg"),
         device="jpeg",width=8,height=10,units="in")
  
  ##by elimination area and fy ------------------------------
  pabs_sum_fy <- pabs_sum_fy %>% 
    st_intersection(elim_areas %>% select(Area_Name))
  elim_thresh <- 0.95
  
  ##by elimination area and season --------------------
  pabs_thresh <- pabs_long %>% 
    group_by(Area_Name) %>% 
    mutate(nsites=length(unique(site_idx))) %>% 
    group_by(Area_Name,period_idx,samp_idx) %>% 
    reframe(over_thresh=sum(value>elim_thresh)/nsites,
            nsites=unique(nsites)) %>% 
    group_by(Area_Name,period_idx) %>% 
    summarise(mn=mean(over_thresh),
              md=median(over_thresh),
              lci=quantile(over_thresh,0.025),
              uci=quantile(over_thresh,0.975))
  pabs_thresh$area_name_label <- paste("Elimination Area",pabs_thresh$Area_Name)
  pabs_thresh$area_name_label[pabs_thresh$Area_Name=="0"] <- "Outside EAs"
  
  pabs_thresh <- pabs_thresh %>% 
    left_join(dat_occ %>% ungroup() %>% select(period_idx,per_start) %>% distinct())
  
  ggplot(pabs_thresh%>% filter(Area_Name!="0"))+
    geom_ribbon(aes(x=per_start,ymin=lci,ymax=uci),alpha=0.3)+
    geom_line(aes(x=per_start,y=md,col=Area_Name),lwd=1.5)+
    facet_wrap(.~area_name_label)+
    scale_y_continuous(name=paste0("Proportion of watersheds with p(elim)>",elim_thresh))+
    xlab("Season")+
    scale_color_discrete(name="Elimination area")+
    theme(text=element_text(size=15))+
    guides(color="none")
  
  ggsave(filename=file.path(out_dir,"plots",
                            paste0("pabs_over",elim_thresh,"_season_ea_by_ws.jpeg")),
         device="jpeg",width=10,height=5,units="in")
  
  ## extent of p(abs) threshold ----------------------
  ### entire site ----------------
  pabs_long$above_thresh <- ifelse(pabs_long$value>elim_thresh,1,0)
  pabs_long <- pabs_long %>% left_join(ssg_formatted %>% 
                                         st_drop_geometry() %>% 
                                         rename(site_idx=SiteID) %>%
                                         select(site_idx,area_km))
  tot_area_km <- as.numeric(ssg_formatted %>% st_union() %>% st_area()/1e6)
  
  pabs_thresh_ext <- pabs_long %>% 
    group_by(above_thresh,period_idx,samp_idx) %>% 
    summarise(area_above_thresh=sum(area_km,na.rm = T)) %>% 
    filter(above_thresh==1) %>% 
    group_by(period_idx) %>% 
    summarise(mn=mean(area_above_thresh),
              md=median(area_above_thresh),
              lci=quantile(area_above_thresh,0.025),
              uci=quantile(area_above_thresh,0.975),
              mn_prop=mn/tot_area_km,
              md_prop=md/tot_area_km,
              lci_prop=lci/tot_area_km,
              uci_prop=uci/tot_area_km)%>% 
    left_join(dat_occ %>% select(period_idx,per_start) %>% distinct())
  
  ggplot(pabs_thresh_ext)+ 
    geom_ribbon(aes(x=per_start,ymin=lci,ymax=uci),alpha=0.5)+
    geom_line(aes(x=per_start,y=md),lwd=1)+
    xlab("Season")+
    ylab(expression(paste("K",m^2," with > 0.95 p(elimination)")))+
    theme(text=element_text(size=15))
  
  ggsave(filename=file.path(out_dir,"plots",
                            paste0("km2_eliminated_",elim_thresh,".jpeg")),
         width=7,height=5,units="in",device="jpeg")
  
  ggplot(pabs_thresh_ext)+ 
    geom_ribbon(aes(x=per_start,ymin=lci/2.59,ymax=uci/2.59),alpha=0.5)+
    geom_line(aes(x=per_start,y=md/2.59),lwd=1)+
    xlab("Season")+
    ylab(expression(paste("M",i^2," with > 0.95 p(elimination)")))+
    theme(text=element_text(size=15))
  ggsave(filename=file.path(out_dir,"plots",
                            paste0("mi2_eliminated_",elim_thresh,".jpeg")),
         width=7,height=5,units="in",device="jpeg")
  
  ggplot(pabs_thresh_ext)+ 
    geom_ribbon(aes(x=per_start,ymin=lci_prop,ymax=uci_prop),alpha=0.5)+
    geom_line(aes(x=per_start,y=md_prop),lwd=1)+
    xlab("Season")+
    ylab("Proportion of study area with p(elimination) > 0.95")+
    theme(text=element_text(size=15))
  
  ggsave(filename=file.path(out_dir,"plots",
                            paste0("prop_eliminated_",elim_thresh,".jpeg")),
         width=7,height=5,units="in",device="jpeg")
  
  ### by EA ----------------
  pabs_thresh_ext_samp <- pabs_long %>% 
    group_by(above_thresh,Area_Name,period_idx,samp_idx) %>% 
    summarise(over_thresh=sum(value>elim_thresh)/nsites,
              area_above_thresh=sum(area_km,na.rm=T)) %>% 
    filter(above_thresh==1) %>%
    left_join(elim_areas %>% 
                summarise(ea_area_km=as.numeric(st_area(geometry)/1e6)) %>% 
                st_drop_geometry()) 
  
  max_per <- pabs_thresh_ext_samp %>%
    group_by(Area_Name) %>% 
    summarise(period_idx=max(period_idx))
  
  pabs_thresh_ext_ea <- pabs_thresh_ext_samp %>% 
    group_by(period_idx,Area_Name) %>% 
    summarise(ea_area_km=unique(ea_area_km),
              mn=mean(area_above_thresh,na.rm=T),
              md=median(area_above_thresh,na.rm=T),
              lci=quantile(area_above_thresh,0.025),
              uci=quantile(area_above_thresh,0.975),
              mn_prop=mn/ea_area_km,
              md_prop=md/ea_area_km,
              lci_prop=lci/ea_area_km,
              uci_prop=uci/ea_area_km) 
  
  if(any(max_per$period_idx<nperiods)){
    low_elim <- max_per %>% filter(period_idx<nperiods)
    samp <- list()
    for(i in 1:nrow(low_elim)){
      samp[[i]] <- expand_grid(period_idx=low_elim$period_idx[i]:nperiods,
                               Area_Name=low_elim$Area_Name[i],
                               ea_area_km=unique(pabs_thresh_ext_ea$ea_area_km[pabs_thresh_ext_ea$Area_Name==low_elim$Area_Name[i]]),
                               mn=0,
                               md=0,
                               lci=0,
                               uci=0,
                               mn_prop=0,
                               md_prop=0,
                               lci_prop=0,
                               uci_prop=0
      )
    }
    samp <- do.call("rbind.data.frame",samp)
    
    pabs_thresh_ext_ea <- rbind.data.frame(pabs_thresh_ext_ea,samp)
  }
  
  pabs_thresh_ext_ea <- pabs_thresh_ext_ea %>% 
    left_join(dat_occ %>% select(period_idx,per_start) %>% distinct())
  
  pabs_thresh_ext_ea$area_name_label <- paste("Elimination Area",pabs_thresh_ext_ea$Area_Name)
  
  ggplot(pabs_thresh_ext_ea %>% filter(Area_Name!="0"))+ 
    geom_ribbon(aes(x=per_start,ymin=lci,
                    ymax=uci),alpha=0.3)+
    geom_line(aes(x=per_start,y=md,col=Area_Name),lwd=1)+
    guides(color="none")+
    facet_wrap(.~area_name_label)+
    xlab("Season")+
    ylab(expression(paste("K",m^2," with > 0.95 p(elimination)")))+
    theme(text=element_text(size=15))
  
  ggsave(filename=file.path(out_dir,"plots",
                            paste0("/km2_eliminated_",elim_thresh,"_ea.jpeg")),
         width=7,height=5,units="in",device="jpeg")
  
  ggplot(pabs_thresh_ext_ea %>% filter(Area_Name!="0"))+ 
    geom_ribbon(aes(x=per_start,ymin=lci_prop,
                    ymax=uci_prop),alpha=0.3)+
    geom_line(aes(x=per_start,y=md_prop,col=Area_Name),lwd=1)+
    guides(color="none")+
    facet_wrap(.~area_name_label)+
    xlab("Season")+
    ylab("Proportion of area with p(elimination) > 0.95")+
    theme(text=element_text(size=15))
  
  ggsave(filename=file.path(out_dir,"plots",
                            paste0("prop_area_eliminated_",elim_thresh,"_ea.jpeg")),
         width=7,height=5,units="in",device="jpeg")
  
  ### by year ----------------------
  pabs_long <- pabs_long %>% 
    left_join(dat_occ %>% select(period_idx,per_start) %>% distinct())
  pabs_long$year <- year(pabs_long$per_start)
  pabs_long$fy <- pabs_long$year
  pabs_long$fy[month(pabs_long$per_start)%in%c(10:12)] <- 
    pabs_long$fy[month(pabs_long$per_start)%in%c(10:12)] +1
  
  pabs_thresh_yr <- pabs_long %>% 
    group_by(site_idx,fy,samp_idx) %>% 
    summarise(yr_mn=mean(value),
              area_km=unique(area_km),
              above_thresh=ifelse(yr_mn>elim_thresh,1,0)) %>% 
    filter(above_thresh==1) %>% 
    group_by(fy,samp_idx) %>% 
    summarise(area_above_thresh=sum(area_km,na.rm=T)) %>% 
    group_by(fy) %>% 
    summarise(mn=mean(area_above_thresh,na.rm=T),
              md=median(area_above_thresh,na.rm=T),
              lci=quantile(area_above_thresh,0.025,na.rm=T),
              uci=quantile(area_above_thresh,0.975,na.rm=T),
              mn_prop=mn/tot_area_km,
              md_prop=md/tot_area_km,
              lci_prop=lci/tot_area_km,
              uci_prop=uci/tot_area_km)
  
  write.csv(pabs_thresh_yr,
            file.path(out_dir,"tables",
                      paste0("area_above_,",elim_thresh,"_table.csv")))
  
  ### by year and ea ------------------
  pabs_thresh_yr_ea_samp <- pabs_long %>% 
    group_by(site_idx,Area_Name,fy,samp_idx) %>% 
    summarise(yr_mn=mean(value,na.rm=T),
              area_km=unique(area_km),
              above_thresh=ifelse(yr_mn>elim_thresh,1,0)) %>% 
    filter(above_thresh==1) %>% 
    left_join(elim_areas %>% 
                summarise(ea_area_km=as.numeric(st_area(geometry)/1e6)) %>% 
                st_drop_geometry()) 
  
  pabs_thresh_yr_ea <- pabs_thresh_yr_ea_samp %>% 
    group_by(fy,Area_Name,samp_idx) %>% 
    summarise(area_above_thresh=sum(area_km,na.rm=T),
              ea_area_km=unique(ea_area_km)) %>% 
    group_by(fy,Area_Name,ea_area_km) %>% 
    summarise(mn=mean(area_above_thresh,na.rm=T),
              md=median(area_above_thresh,na.rm=T),
              lci=quantile(area_above_thresh,0.025,na.rm=T),
              uci=quantile(area_above_thresh,0.975,na.rm=T)) %>% 
    mutate(mn_prop=mn/ea_area_km,
           md_prop=md/ea_area_km,
           lci_prop=lci/ea_area_km,
           uci_prop=uci/ea_area_km)
  
  pabs_thresh_yr_ea_prop <- pabs_thresh_yr_ea %>% 
    select(-ea_area_km) %>% 
    pivot_longer(cols=c(mn,md,lci,uci,mn_prop,md_prop,lci_prop,uci_prop),
                 names_to="stat",values_to="value") %>%
    pivot_wider(names_from=Area_Name,values_from=value) %>%
    filter(grepl("prop",stat))
  
  
  pabs_thresh_yr_ea_prop[is.na(pabs_thresh_yr_ea_prop)] <- 0
  
  write.csv(pabs_thresh_yr_ea_prop,
            file.path(out_dir,"tables",
                      paste0("area_above_,",elim_thresh,"_ea_prop_table.csv")))
  
  pabs_thresh_yr_ea_km <- pabs_thresh_yr_ea %>% 
    select(-ea_area_km) %>% 
    pivot_longer(cols=c(mn,md,lci,uci,mn_prop,md_prop,lci_prop,uci_prop),
                 names_to="stat",values_to="value") %>%
    pivot_wider(names_from=Area_Name,values_from=value) %>%
    filter(!grepl("prop",stat))
  
  pabs_thresh_yr_ea_km[is.na(pabs_thresh_yr_ea_km)] <- 0
  
  write.csv(pabs_thresh_yr_ea_km,
            file.path(out_dir,"tables",
                      paste0("area_above_,",elim_thresh,"_ea_km_table.csv")))
  
  ##median p(abs) by season and ea -----------------------
  pabs_sum_ea <- pabs_long %>% 
    group_by(Area_Name,period_idx,samp_idx) %>% 
    reframe(ea_md=median(value)) %>% 
    group_by(Area_Name,period_idx) %>% 
    summarise(mn=mean(ea_md),
              md=median(ea_md),
              lci=quantile(ea_md,0.025),
              uci=quantile(ea_md,0.975))
  pabs_sum_ea$Area_Name <- paste("Elimination Area",pabs_sum_ea$Area_Name)
  pabs_sum_ea$Area_Name[pabs_sum_ea$Area_Name=="Elimination Area 0"] <- "Outside EAs"
  pabs_sum_ea<- pabs_sum_ea %>% 
    left_join(dat_occ %>% ungroup() %>% select(period_idx,per_start) %>% distinct())
  
  ggplot(pabs_sum_ea %>% filter(Area_Name!="Outside EAs"))+
    geom_ribbon(aes(x=per_start,ymin=lci,ymax=uci),alpha=0.3)+
    geom_line(aes(x=per_start,y=md,col=Area_Name),lwd=1.5)+
    facet_wrap(.~Area_Name)+
    xlab("Season")+
    scale_y_continuous(name="Median probability of feral swine absence")+
    scale_color_discrete(name="Elimination area")+
    theme(text=element_text(size=15))+
    guides(color="none")
  
  ggsave(filename=file.path(out_dir,"plots","pabs_md_season_ea.jpeg"),
         device="jpeg",width=10,height=6,units="in")
  
  list(
    pabs_sum=pabs_sum,
    pabs_sum_fy=pabs_sum_fy,
    pabs_sum_ea=pabs_sum_ea,
    pabs_thresh=pabs_thresh, 
    pabs_thresh=pabs_thresh,
    pabs_thresh_yr=pabs_thresh_yr,
    pabs_thresh_ext=pabs_thresh_ext
  )
}

plot_det <- function(samples,
                     dat_occ,
                     out_dir){
  
  if(class(samples)=="list"){
    nChains <- length(samples)
    nmcmc <- nrow(samples[[1]])
  } else {
    nChains <- 1
    nmcmc <- nrow(samples)
  }
  
  if(nChains>1){
    alpha <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
      cbind.data.frame(samples[[i]][,grepl("alpha",colnames(samples[[i]]))],chain=i)})),
      samp=rep(1:nmcmc,nChains))
  } else {
    alpha <- as.data.frame(samples[,grepl("alpha",colnames(samples))])
    alpha$chain <-1 
    alpha$samp <- 1:nmcmc
  }
  
  sys_eff_sc <- seq(min(dat_occ$trap_nights_km_sc),max(dat_occ$trap_nights_km_sc),by=0.5)
  sys_eff <- sys_eff_sc * attr(dat_occ$trap_nights_km_sc,"scaled:scale") +
    attr(dat_occ$trap_nights_km_sc,"scaled:center")
  
  det <- sapply(1:(nmcmc*nChains),function(i){
    boot::inv.logit(alpha[i,1] + alpha[i,2]*sys_eff_sc)})
  
  det_sum <- data.frame(traps_km=sys_eff,
                        mn=rowMeans(det),
                        md=sapply(1:nrow(det),function(i)quantile(det[i,],prob=0.5)),
                        lci=sapply(1:nrow(det),function(i)quantile(det[i,],prob=0.025)),
                        uci=sapply(1:nrow(det),function(i)quantile(det[i,],prob=0.975)))
  
  ggplot(det_sum)+
    geom_ribbon(aes(x=traps_km,ymin=lci,ymax=uci),alpha=0.3)+
    geom_line(aes(x=traps_km,y=md),lwd=1)+
    ylab("Detection probability")+
    xlab(expression(paste("Bait station nights per k",m^2)))+
    theme(text=element_text(size=15))
  
  ggsave(filename=file.path(out_dir,"plots","occ_det_eff_curve.jpeg"),
         device="jpeg",width=7,height=5,units="in")
  
  det_sum
}

plot_elim_samp <- function(samples,
                           dat_occ,
                           out_dir){
  if(class(samples)=="list"){
    nChains <- length(samples)
    nmcmc <- nrow(samples[[1]])
  } else {
    nChains <- 1
    nmcmc <- nrow(samples)
  }
  syseff <- data.frame(trap_nights=c((0.25/2.59)*45,
                                     (0.5/2.59)*45,
                                     (1/2.59)*45))
  syseff$trap_nights_sc <- (syseff$trap_nights-attr(dat_occ$trap_nights_km_sc,"scaled:center"))/
    attr(dat_occ$trap_nights_km_sc,"scaled:scale") 
  
  if(nChains>1){
    sd_pdet <-  do.call("rbind",lapply(1:nChains,function(i){
      samples[[i]][grep("sd_pdet",colnames(samples[[i]]))]}))
    
    alpha <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
      cbind.data.frame(samples[[i]][,grepl("alpha",colnames(samples[[i]]))],chain=i)})),
      samp=rep(1:nmcmc,nChains))
  } else {
    sd_pdet <- samples[grep("sd_pdet",colnames(samples))]
    
    alpha <- as.data.frame(samples[,grepl("alpha",colnames(samples))])
    alpha$chain <-1 
    alpha$samp <- 1:nmcmc
  }
  
  
  pabs_ex <- c(0.25,0.5,0.75)
  elim_prob <- 0.95
  eff_elim <- array(NA,dim=c(nmcmc*nChains,length(pabs_ex),nrow(syseff)))
  eff_elim_sum <- list()
  for(k in 1:nrow(syseff)){
    for(j in 1:length(pabs_ex)){
      for(i in 1:(nmcmc*nChains)){
        p_sys <- boot::inv.logit(
          rnorm(1,alpha[i,1] + alpha[i,2]*syseff$trap_nights_sc[k],sd_pdet[i,1]))
        eff_elim[i,j,k] <- log((pabs_ex[j]*(1-elim_prob))/
                                 (elim_prob*(1-pabs_ex[j])))/(log(1-p_sys)) 
      }
    }
    eff_elim_sum[[k]] <- data.frame(pabs=pabs_ex,
                                    mn=colMeans(eff_elim[,,k]),
                                    md=apply(eff_elim[,,k],2,median),
                                    lci=sapply(1:ncol(eff_elim[,,k]),function(x)quantile(eff_elim[,x,k],0.025)),
                                    uci=sapply(1:ncol(eff_elim[,,k]),function(x)quantile(eff_elim[,x,k],0.975)))
    eff_elim_sum[[k]]$syseff <- syseff$trap_nights[k]
  }
  eff_elim_sum <- do.call("rbind.data.frame",eff_elim_sum)
  eff_elim_sum$syseff_lab<- c(rep("Low",length(pabs_ex)),
                              rep("Medium",length(pabs_ex)),
                              rep("High",length(pabs_ex)))
  eff_elim_sum$syseff_lab <- factor(eff_elim_sum$syseff_lab,levels=c("Low","Medium","High"))
  
  ggplot(eff_elim_sum)+
    geom_bar(aes(y=md,x=factor(pabs),fill=factor(pabs)),
             stat="identity",position="dodge")+
    geom_errorbar(aes(ymin=lci,ymax=uci,x=factor(pabs),width=0))+
    facet_wrap(.~syseff_lab)+
    xlab("Watershed probability of feral swine absence")+
    ylab("Number of 10-day subperiods to determine 95% elimination probability")+
    guides(fill="none")+
    geom_hline(yintercept=max(dat_occ$nweeks),lty=2)
  
  ggsave(filename = file.path(out_dir,"plots","nweeks_to_elim_trapeff.jpeg"),
         width=7,height=5,units="in",device="jpeg")
  
  eff_elim_sum
}

plot_abundance <- function(samples,
                           dat_rem,
                           dat_occ,
                           ssg_formatted,
                           elim_areas,
                           site_idx_lookup,
                           out_dir){
  
  if(class(samples)=="list"){
    nChains <- length(samples)
    nmcmc <- nrow(samples[[1]])
  } else {
    nChains <- 1
    nmcmc <- nrow(samples)
  }
  
  elim_areas$area_km <- as.numeric(st_area(elim_areas)/1e6)
  
  if(nChains>1){
    N <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
      cbind.data.frame(samples[[i]][grepl("N",colnames(samples[[i]])) & 
                                      !grepl("N_latent",colnames(samples[[i]]))],
                       chain=i)})),
      samp=rep(1:nmcmc,nChains))
    
    lambda_N <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
      cbind.data.frame(samples[[i]][grepl("lambda",colnames(samples[[i]]))],
                       chain=i)})),
      samp=rep(1:nmcmc,nChains))
  } else {
    N <- samples[grepl("N",colnames(samples)) & 
                   !grepl("N_latent",colnames(samples))]
    N$chain <-1 
    N$samp <- 1:nmcmc
    
    lambda_N <- samples[grepl("lambda",colnames(samples))]
    lambda_N$chain <-1 
    lambda_N$samp <- 1:nmcmc
  }
  
  nsites_rem <- max(dat_rem$site_idx_rem)
  nperiods <- max(dat_rem$period_idx)
  
  N_long <- N %>% pivot_longer(cols=1:(ncol(N)-2),
                               values_to = "N",
                               names_to="idx") %>% 
    select(-idx)
  
  N_long$site_idx_rem <- rep(rep(1:nsites_rem,nperiods),nmcmc*nChains)
  N_long$period_idx <- rep(sort(rep(1:nperiods,nsites_rem)),nmcmc*nChains)
  N_long <- N_long %>%
    left_join(site_idx_lookup)
  
  per_idx <- dat_occ %>% select(period_idx,per_start) %>% 
    mutate(year=year(per_start)) %>% distinct()
  per_idx$fy <- per_idx$year
  per_idx$fy[month(per_idx$per_start)%in%c(10,11,12)] <- 
    per_idx$fy[month(per_idx$per_start)%in%c(10,11,12)] +1
  
  N_long <- N_long %>% left_join(per_idx)
  
  lambda_long <- lambda_N %>% 
    pivot_longer(cols=1:(ncol(lambda_N)-2),
                 values_to = "lambda",
                 names_to="idx") %>% 
    select(-idx)
  lambda_long$site_idx_rem <- rep(1:nsites_rem,nmcmc*nChains)
  lambda_long <- lambda_long %>% 
    left_join(site_idx_lookup)
  
  N_long <- N_long %>% 
    left_join(lambda_long)%>% 
    left_join(ssg_formatted %>%
                st_drop_geometry() %>% 
                select(SiteID,area_km)) %>% 
    filter(!is.na(N))
  
  #proportion needed to combat pop growth
  N_long$lambda_prop <- rpois(nrow(N_long),N_long$N * (N_long$lambda-1))
  N_long$lambda_prop_dens <- N_long$lambda_prop/N_long$area_km
  N_long$lambda_prop[is.na(N_long$lambda_prop)] <- 0
  N_long$lambda_prop_dens[is.na(N_long$lambda_prop_dens)] <- 0
  
  #standardizes abundance 
  N_ss <- sqrt(sum(N_long$N[!is.na(N_long$N)]^2)/(length(N_long$N[!is.na(N_long$N)])-1))
  N_long$N_std <- N_long$N/N_ss
  N_long$lambda_prop_std <- N_long$lambda_prop/N_ss
  
  #projecting pop growth sans removal 
  N_long$N_no_rem <- NA
  N_long$N_no_rem[N_long$period_idx==1] <- N_long$N[N_long$period_idx==1] 
  for(i in 2:(nperiods)){
    N_long$N_no_rem[N_long$period_idx==i] <-
      rpois(length(N_long$N_no_rem[N_long$period_idx==i]),
            N_long$N_no_rem[N_long$period_idx==(i-1)]*N_long$lambda[N_long$period_idx==(i-1)])
  }
  
  #density 
  N_long <- N_long %>% mutate(dens=N/area_km,
                              dens_no_rem=N_no_rem/area_km,
                              N_ss_no_rem=N_no_rem/N_ss)
  
  N_sum <- N_long %>% 
    pivot_longer(cols=c("N","N_std","dens","N_no_rem","dens_no_rem","N_ss_no_rem",
                        "lambda_prop","lambda_prop_std","lambda_prop_dens"),
                 names_to="metric",values_to="value") %>% 
    group_by(metric,SiteID,period_idx,per_start,year,fy) %>% 
    summarise(mn=mean(value,na.rm=T),
              md=median(value),
              lci=quantile(value,0.025,na.rm=T),
              uci=quantile(value,0.975,na.rm=T)) %>% 
    pivot_wider(names_from=metric,values_from = c(mn,md,lci,uci))
  
  
  dat_rem_sum <- dat_rem %>% 
    group_by(SiteID,method,month=floor_date(subper_start,"month")) %>% 
    summarise(Season=min(subper_start),
              removal=sum(tot_rem),
              tot_hrs=sum(tot_hrs),
              tot_events=sum(num_events)) %>% 
    filter(month<=max(dat_occ$per_start)) %>% 
    left_join(ssg_formatted %>% 
                st_drop_geometry() %>% 
                select(SiteID,area_km)) %>% 
    mutate(tot_rem_km = removal/area_km)
  
  site_rem <- dat_rem %>%
    filter(tot_rem>0) %>%
    group_by(SiteID) %>%
    summarise(n=n()) %>%
    arrange(desc(n))
  
  site_rem_cutoff <- site_rem[1:6,]
  N_sum_top <- N_sum %>% filter(SiteID%in%site_rem_cutoff$SiteID)
  dat_rem_sum_top <- dat_rem_sum %>% filter(SiteID%in%site_rem_cutoff$SiteID)
  
  #map of highglighted sites
  ssg_top <- ssg_formatted %>% filter(SiteID%in%site_rem_cutoff$SiteID)
  
  ggplot()+
    geom_sf(data=ssg_formatted)+
    geom_sf(data=ssg_top,aes(fill=factor(SiteID)))+
    geom_sf_label(data=ssg_top,aes(label=SiteID),size=3,alpha=0.7)+
    guides(fill="none")+
    theme(panel.background = element_blank(),
          text=element_text(size=15),
          axis.ticks=element_blank(),
          axis.title = element_blank(),
          axis.text=element_blank())
  ggsave(filename=file.path(out_dir,"plots","top_removal_sites_map.jpeg"),
         device="jpeg",width=7,height=5,units="in")
  
  #abundance
  axis_trans<- 0.4
  ggplot()+  
    geom_ribbon(data=N_sum_top,
                aes(x=per_start,ymin=lci_N,ymax=uci_N,fill="Estimated Abundance"),alpha=0.2)+
    geom_line(data=N_sum_top,
              aes(x=per_start,y=md_N,col="Estimated Abundance"),lwd=1)+
    geom_ribbon(data=N_sum_top,
                aes(x=per_start,ymin=lci_N_no_rem,ymax=uci_N_no_rem,
                    fill="No Removal"),
                alpha=0.1)+
    geom_line(data=N_sum_top,
              aes(x=per_start,y=md_N_no_rem,col="No Removal"),
              lwd=1,lty=3)+
    facet_wrap(.~SiteID)+
    scale_fill_manual(name="",
                      values=c("Estimated Abundance"="black",
                               "No Removal"="red",
                               "Removal Needed"="blue"))+
    scale_color_manual(name="",
                       values=c("Estimated Abundance"="black",
                                "No Removal"="red",
                                "Removal Needed"="blue"))+
    ylab("Feral swine abundance")+
    xlab("Season")+
    theme(text=element_text(size=15))
  
  ggsave(filename=file.path(out_dir,"plots","abundance_trend_removal.jpeg"),
         device="jpeg",width=10,height=6,units="in")
  
  #density
  axis_trans_d<- 0.05
  ggplot()+  
    geom_ribbon(data=N_sum_top,
                aes(x=per_start,ymin=lci_dens,ymax=uci_dens,
                    fill="Estimated Density"),alpha=0.2)+
    geom_line(data=N_sum_top,
              aes(x=per_start,y=md_dens,col="Estimated Density"),lwd=1)+
    geom_ribbon(data=N_sum_top,
                aes(x=per_start,ymin=lci_dens_no_rem,ymax=uci_dens_no_rem,fill="No Removal"),
                alpha=0.1)+
    geom_line(data=N_sum_top,
              aes(x=per_start,y=md_dens_no_rem,col="No Removal"),lwd=1,lty=3)+
    facet_wrap(.~SiteID)+
    scale_fill_manual(name="",
                      values=c("Estimated Density"="black",
                               "No Removal"="red"))+
    
    scale_color_manual(name="",
                       values=c("Estimated Density"="black",
                                "No Removal"="red"))+
    ylab(expression(paste("Feral swine density (swine/k",m^2,")")))+
    xlab("Season")+
    theme(text=element_text(size=15))
  
  ggsave(filename=file.path(out_dir,"plots","density_trend_removal.jpeg"),
         device="jpeg",width=10,height=6,units="in")
  
  #removal needed
  axis_trans_d<- 0.03
  ggplot()+  
    geom_bar(data=dat_rem_sum_top,
             aes(x=month,y=tot_rem_km/axis_trans_d),stat="identity",
             alpha=0.6)+
    geom_ribbon(data=N_sum_top,
                aes(x=per_start,
                    ymin=lci_lambda_prop_dens/axis_trans_d,
                    ymax=uci_lambda_prop_dens/axis_trans_d,
                    fill="Removal Needed"),
                alpha=0.1)+
    geom_line(data=N_sum_top,aes(x=per_start,
                                 y=md_lambda_prop_dens/axis_trans_d,group=SiteID,
                                 col="Removal Needed"),lwd=1,lty=2)+
    
    facet_wrap(.~SiteID)+
    scale_fill_manual(name="",
                      values=c("Estimated Density"="black",
                               "No Removal"="red",
                               "Removal Needed"="blue"))+
    scale_color_manual(name="",
                       values=c("Estimated Density"="black",
                                "No Removal"="red",
                                "Removal Needed"="blue"))+
    ylab(expression(paste("Feral swine density (swine/k",m^2,")")))+
    xlab("Season")+
    theme(text=element_text(size=15))
  
  ggsave(filename=file.path(out_dir,"plots","density_trend_removal_needed.jpeg"),
         device="jpeg",width=10,height=6,units="in")
  
  
  #standardized abundance
  axis_scale <- 3
  ggplot()+  
    geom_ribbon(data=N_sum_top,
                aes(x=per_start,ymin=lci_N_std,ymax=uci_N_std,
                    fill="Estimated \nStandardized Abundance"),alpha=0.2)+
    geom_line(data=N_sum_top,
              aes(x=per_start,y=md_N_std,col="Estimated \nStandardized Abundance"),lwd=1)+
    geom_ribbon(data=N_sum_top,
                aes(x=per_start,ymin=lci_N_ss_no_rem,ymax=uci_N_ss_no_rem,
                    fill="No Removal"),
                alpha=0.1)+
    geom_line(data=N_sum_top,
              aes(x=per_start,y=md_N_ss_no_rem,col="No Removal"),
              lty=3,lwd=1)+
    facet_wrap(.~SiteID)+
    scale_fill_manual(name="",
                      values=c("Estimated \nStandardized Abundance"="black",
                               "No Removal"="red",
                               "Removal Needed"="blue"))+
    scale_color_manual(name="",
                       values=c("Estimated \nStandardized Abundance"="black",
                                "No Removal"="red",
                                "Removal Needed"="blue"))+
    ylab("Standardized feral swine abundance")+
    xlab("Season")+
    theme(text=element_text(size=15))
  
  ggsave(filename=file.path(out_dir,"plots","std_abundance_trend_removal.jpeg"),
         device="jpeg",width=10,height=6,units="in")
  
  #annual change in abundance table
  N_yr <- N_long %>% 
    # filter(SiteID%in%site_rem_cutoff$SiteID) %>% 
    filter(month(per_start)==10)
  
  N_yr$per_temp <- as.numeric(as.factor(N_yr$period_idx))
  # N_yr$N_prev<- NA
  # for(i in 1:nrow(N_yr)){
  #   if(N_yr$per_temp[i]>1){
  #     N_yr$N_prev[i] <- N_yr$N[N_yr$per_temp==(N_yr$per_temp[i]-1) &
  #                                N_yr$SiteID==N_yr$SiteID[i] &
  #                                N_yr$samp==N_yr$samp[i] &
  #                                N_yr$chain==N_yr$chain[i]]
  #   }
  # }
  N_yr$N_prev <- sapply(1:nrow(N_yr),function(i){
    if(N_yr$per_temp[i]>1){
      N_yr$N[N_yr$per_temp==(N_yr$per_temp[i]-1) &
               N_yr$SiteID==N_yr$SiteID[i] &
               N_yr$samp==N_yr$samp[i] &
               N_yr$chain==N_yr$chain[i]]
    } else {
      NA
    }
  })
  
  N_yr$per_change <- (N_yr$N-N_yr$N_prev)/N_yr$N_prev
  
  N_yr_sum <- N_yr %>% 
    filter(year>2020) %>% 
    group_by(SiteID,period_idx,per_start) %>% 
    summarise(mn=mean(per_change,na.rm=T),
              md=median(per_change,na.rm=T),
              lci=quantile(per_change,0.025,na.rm=T),
              uci=quantile(per_change,0.975,na.rm=T)) 
  
  write.csv(N_yr_sum,file.path(out_dir,"tables","N_change_table.csv"))
  
  list(
    N_sum=N_sum,
    N_yr_sum=N_yr_sum,
    dat_rem_sum=dat_rem_sum
  )
}

plot_removal_prob <- function(samples,
                              dat_rem,
                              out_dir
){
  if(class(samples)=="list"){
    nChains <- length(samples)
    nmcmc <- nrow(samples[[1]])
  } else {
    nChains <- 1
    nmcmc <- nrow(samples)
  }
  
  if(nChains>1){
    delta <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
      cbind.data.frame(samples[[i]][,grepl("delta",colnames(samples[[i]]))],chain=i)})),
      samp=rep(1:nmcmc,nChains))
  } else {
    delta <- as.data.frame(samples[,grepl("delta",colnames(samples))])
    delta$chain <-1 
    delta$samp <- 1:nmcmc
  }
  
  cpue_names <- c("Aerial Intercept","Aerial Slope","Trap Intercept","Trap Slope")
  colnames(delta)[1:4] <- cpue_names
  
  dat_aerial <- dat_rem %>% filter(method=="Aerial")
  dat_trap <- dat_rem %>% filter(method=="Trap")
  
  #aerial
  aerial_eff_sc <- seq(min(dat_aerial$effect_area_hrs),max(dat_aerial$effect_area_hrs),by=0.005)
  aerial_eff <- aerial_eff_sc * attr(scale(dat_aerial$effect_area_hrs),"scaled:scale") +
    attr(scale(dat_aerial$effect_area_hrs),"scaled:center")
  aerial_area_mn <- mean(dat_aerial$prop_site_impact[dat_aerial$prop_site_impact!=0])
  
  det_a <- sapply(1:(nmcmc*nChains),function(i){aerial_area_mn*
      boot::inv.logit(delta$`Aerial Intercept`[i] + delta$`Aerial Slope`[i]*aerial_eff_sc)})
  
  det_a_sum <- data.frame(effect_area_hrs=aerial_eff,
                          mn=rowMeans(det_a),
                          md=sapply(1:nrow(det_a),function(i)quantile(det_a[i,],prob=0.5)),
                          lci=sapply(1:nrow(det_a),function(i)quantile(det_a[i,],prob=0.025)),
                          uci=sapply(1:nrow(det_a),function(i)quantile(det_a[i,],prob=0.975)))
  
  g_a <- ggplot(det_a_sum)+
    geom_ribbon(aes(x=effect_area_hrs,ymin=lci,ymax=uci),alpha=0.3)+
    geom_line(aes(x=effect_area_hrs,y=md))+
    ylab("Detection probability")+
    xlab(expression(paste("Hours flown per k",m^2)))+
    theme(text=element_text(size=15))
  
  #trap
  trap_eff_sc <- seq(min(dat_trap$effect_area_hrs),max(dat_trap$effect_area_hrs),by=0.05)
  trap_eff <- trap_eff_sc * attr(scale(dat_trap$effect_area_hrs),"scaled:scale") +
    attr(scale(dat_trap$effect_area_hrs),"scaled:center")
  trap_area_mn <- mean(dat_trap$prop_site_impact[dat_trap$prop_site_impact!=0])
  
  det_t <- sapply(1:(nmcmc*nChains),function(i){trap_area_mn*
      boot::inv.logit(delta$`Trap Intercept`[i] + delta$`Trap Slope`[i]*trap_eff_sc)})
  
  det_t_sum <- data.frame(effect_area_hrs=trap_eff,
                          mn=rowMeans(det_t),
                          md=sapply(1:nrow(det_t),function(i)quantile(det_t[i,],prob=0.5)),
                          lci=sapply(1:nrow(det_t),function(i)quantile(det_t[i,],prob=0.025)),
                          uci=sapply(1:nrow(det_t),function(i)quantile(det_t[i,],prob=0.975)))
  
  g_t <- ggplot(det_t_sum)+
    geom_ribbon(aes(x=effect_area_hrs,ymin=lci,ymax=uci),alpha=0.3)+
    geom_line(aes(x=effect_area_hrs,y=md))+
    ylab("Proportion of population removed in 10 days")+
    xlab(expression(paste("Hours spent trapping per k",m^2)))+
    theme(text=element_text(size=15))
  
  ##all together --------------------
  g_all <- gridExtra::arrangeGrob(g_a, g_t,nrow=1) #generates g
  
  ggsave(g_all,file=file.path(out_dir,"plots","rem_det_curve_all.jpeg"),
         device="jpeg",width=13,height=6,units="in")
  
  ##removal comparisons --------------------
  a_mn <- max(dat_aerial$prop_site_impact)*
    boot::inv.logit(delta[,"Aerial Intercept"] + delta[,"Aerial Slope"]* 
                      max(scale(dat_aerial$effect_area_hrs)[dat_aerial$effect_area_hrs!=0]))
  
  max(dat_aerial$effect_area_hrs)
  max(dat_aerial$prop_site_impact)
  
  t_mn <- max(dat_trap$prop_site_impact)*
    boot::inv.logit(delta[,"Trap Intercept"] + delta[,"Trap Slope"]* 
                      max(scale(dat_trap$effect_area_hrs)[dat_trap$effect_area_hrs!=0]))
  max(dat_trap$effect_area_hrs)
  max(dat_trap$prop_site_impact)
  
  rem_mn <- cbind.data.frame(aerial=a_mn,trap=t_mn)
  
  rem_df <- rem_mn %>% pivot_longer(cols=1:2,names_to="rem_typ",values_to="value") %>% 
    group_by(rem_typ) %>% 
    summarise(mn=mean(value),
              md=median(value),
              lci=quantile(value,0.025),
              uci=quantile(value,0.975))
  rem_df$effect_area_hrs <- c(max(dat_aerial$effect_area_hrs),max(dat_trap$effect_area_hrs))
  rem_df$prop_site_impact <- c(max(dat_aerial$prop_site_impact),max(dat_trap$prop_site_impact))
  
  list(
    rem_df=rem_df,
    det_a_sum=det_a_sum,
    det_t_sum=det_t_sum
  )
}