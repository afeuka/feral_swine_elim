### Title: Plotting hurdle model results
### Author: Abbey Feuka
### Date: 29NOV2023
### Notes: uses output from fit_zi_rem_mod.R

library(tidyverse)
library(lubridate)
library(sf)
library(tigris)

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
# load("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/ziBinMod_area_04SEP24_logit_det.Rdata")

# if(nbeta==3){
#   subfolder<-"No NFSP"
# } else {
#   subfolder<-"NFSP2"
# }
subfolder <- "No NFSP 2024"

## removal data -----------
dat_rem <- rbind.data.frame(dat_aerial,dat_trap)
dat_rem_sum <- dat_rem %>% 
  group_by(Area_Name,method,month=floor_date(Date,"month")) %>% 
  summarise(Season=min(Date),
            removal=sum(tot_rem),
            tot_hrs=sum(tot_hrs),
            tot_events=sum(num_events)) %>% 
  filter(month<=max(dat_occ$per_start))

dat_rem_sum_all <- rem_eff_ea %>% 
  group_by(Area_Name,method,month=floor_date(Date,"month")) %>% 
  summarise(Season=min(Date),
            removal=sum(tot_rem),
            tot_hrs=sum(tot_hrs),
            tot_events=sum(num_events)) %>% 
  filter(month<=max(dat_occ$per_start)) %>% 
  filter(method!="Ground")

nmcmc <- nrow(samples[[1]])
nChains <- length(samples)

#study site map -----------------------------
# us_states <- states()
# mo <- us_states %>% filter(STUSPS=="MO")
# huc10 <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/huc10_mo.shp")
# 
# ggplot()+geom_sf(data=mo)+
#   geom_sf(data=huc10,fill="grey88")+
#   geom_sf(data=study_site_grid,fill="grey55")+
#   geom_sf(data=elim_areas %>% filter(Area_Name!="0"),
#           aes(col=Area_Name),fill="transparent",
#           lwd=1.5)+
#   scale_color_discrete(name="Elimination Area")+
#   theme(panel.background = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         text=element_text(size=15))
# ggsave(filename = "./Model outputs/Plots/",subfolder,"/mo_map.jpeg",
#        width=7,height=5,units="in",device="jpeg")

#trace plots -------------------------------
##beta coefficients----------------------
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
if(subfolder=="No NFSP 2024"){
  beta <- beta %>% select(-c(`Feral swine range`))
  nbeta <- ncol(beta %>% select(-c(chain,samp)))
}

beta_long <- beta %>% 
  pivot_longer(cols=all_of(1:nbeta),names_to="beta",values_to="value") %>% 
  mutate(chain=factor(chain))

g_b<- ggplot(beta_long) + geom_line(aes(x=samp,y=value,col=chain))+
  geom_hline(yintercept=0,col="black",lty=2)+
  facet_wrap(.~beta)
g_b

ggsave(g_b, 
       filename = paste0("./Model outputs/Plots/",subfolder,"/Trace/beta_trace.jpeg"),
       width=7,height=5,units="in",device="jpeg")

##initial n values ------------------
r_n1 <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
  cbind.data.frame(samples[[i]][,grep("r_n1",colnames(samples[[i]]))],
                   chain=i)})),
  samp=rep(1:nmcmc,nChains))
colnames(r_n1)[1] <- "r_n1"
r_n1_long <- r_n1 %>% 
  mutate(chain=factor(chain))

g_r <- ggplot(r_n1_long) + geom_line(aes(x=samp,y=r_n1,col=chain))

p_n1 <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
  cbind.data.frame(samples[[i]][,grep("p_n1",colnames(samples[[i]]))],
                   chain=i)})),
  samp=rep(1:nmcmc,nChains))
colnames(p_n1)[1] <- "p_n1"
p_n1_long <- p_n1 %>% 
  mutate(chain=factor(chain))

g_p <-ggplot(p_n1_long) + geom_line(aes(x=samp,y=p_n1,col=chain))
g_n1 <- gridExtra::grid.arrange(g_r,g_p)

ggsave(g_n1, 
       filename = paste0("./Model outputs/Plots/",subfolder,"/Trace/n1_params_trace.jpeg"),
       width=7,height=5,units="in",device="jpeg")

## effort removal ----------------------
delta <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
  cbind.data.frame(samples[[i]][,grepl("delta",colnames(samples[[i]]))],chain=i)})),
  samp=rep(1:nmcmc,nChains))

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

## effort systematic baiting -------------------
alpha <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
  cbind.data.frame(samples[[i]][,grepl("alpha",colnames(samples[[i]]))],chain=i)})),
  samp=rep(1:nmcmc,nChains))

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
       filename = paste0("./Model outputs/Plots/",subfolder,"/Trace/detection_trace.jpeg"),
       width=7,height=5,units="in",device="jpeg")

##lambda -----------------
lambda <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
  cbind.data.frame(samples[[i]][,grep("lambda",colnames(samples[[i]]))],
                   chain=i)})),
  samp=rep(1:nmcmc,nChains))

colnames(lambda)[1] <- "Outside Elimination Areas"
colnames(lambda)[2:nea] <- paste("Elimination Area",1:(nea-1))

lambda <- lambda %>% select(c("Elimination Area 4","Elimination Area 6","chain","samp"))

lambda_long <- lambda %>% 
  pivot_longer(cols=all_of(1:2),names_to="elim_area",values_to="value") %>% 
  mutate(chain=factor(chain))

g_l <- ggplot(lambda_long) + geom_line(aes(x=samp,y=value,col=chain))+
  facet_wrap(.~elim_area)+
  geom_hline(yintercept=1,col="black",lty=2)
g_l

ggsave(g_l, 
       filename = paste0("./Model outputs/Plots/",subfolder,"/Trace/lambda_trace.jpeg"),
       width=7,height=5,units="in",device="jpeg")

rm(g_l,g_b,g_a,g_t,g_n1,g_r,g_p)

#occupancy regression---------------------
beta_long %>% 
  group_by(beta) %>% 
  summarise(mn=mean(value),
            md=median(value),
            lci=quantile(value,0.025),
            uci=quantile(value,0.975))

#absence probability --------------------
## map to study site -------------------
pabs <- do.call("rbind",lapply(1:nChains,function(i){samples[[i]][grep("pabs",colnames(samples[[i]]))]}))

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

pabs_sum_sf <- pabs_sum %>% 
  left_join(study_site_grid %>% rename(site_idx=SiteID)) %>% 
  st_as_sf()

## elimination area groupings ---------------------------
pabs_long <- pabs %>% pivot_longer(all_of(1:ncol(pabs)))
pabs_long$site_idx <- rep(rep(1:nsites,nperiods),nmcmc*nChains)
pabs_long$period_idx <- rep(sort(rep(1:nperiods,nsites)),nmcmc*nChains)
pabs_long$samp_idx <- sort(rep(1:(nmcmc*nChains),nsites*nperiods))
pabs_long <- pabs_long %>%
  left_join(study_site_grid %>%
              st_drop_geometry() %>%
              rename(site_idx=SiteID) %>%
              select(site_idx,Area_Name))
# 
# pabs_long %>% 
#   group_by(Area_Name,period_idx) %>% 
#   summarise(md=median(value))

# mean(apply(pabs,1,median))
median(apply(pabs,1,median))
quantile(apply(pabs,1,median),c(0.025,0.975))

pabs_sum[which.max(pabs_sum$md),]
pabs_sum[which.min(pabs_sum$md),]

#elimination area medians (of)
# pabs_sum <- pabs_sum %>% 
#   left_join(study_site_grid %>% 
#               st_drop_geometry() %>% 
#               rename(site_idx=SiteID) %>% 
#               select(site_idx,Area_Name))
# pabs_ea_sum <- pabs_sum %>% 
#   group_by(Area_Name,per_start) %>% 
#   summarise(ea_md=median(mn))
# 
# pabs_ea_sum[which.max(pabs_ea_sum$ea_md),]
# pabs_ea_sum[which.min(pabs_ea_sum$ea_md),]

##fiscal year ----------------
###inidivdual maps------------------
pabs_sum_fy <- pabs_sum %>% group_by(fy,site_idx) %>% 
  summarise(fy_mn = mean(mn),
            fy_md= median(mn),
            fy_lci = mean(lci),
            fy_uci= mean(uci)) %>% 
  left_join(study_site_grid %>% rename(site_idx=SiteID)) %>% 
  st_as_sf()

# fys <- unique(pabs_sum_fy$fy)
# g<-list()
# for(i in 1:length(fys)){
#   g[[i]] <-
#     ggplot()+
#     geom_sf(data=pabs_sum_fy %>%
#               filter(fy==fys[i]),aes(fill=fy_mn))+
#     geom_sf(data=elim_areas %>% filter(Area_Name!="0"),aes(col=Area_Name),
#             fill="transparent",lwd=1.5)+
#     scale_fill_viridis_c(#limits=c(0,1),
#       direction=-1,
#       name="Feral swine \nelimination probability")+
#     scale_color_discrete(name="Elimination Area")+
#     theme(panel.background=element_blank(),
#           axis.ticks=element_blank(),
#           axis.text=element_blank(),
#           plot.margin = margin(1,1,1,1, "cm"),
#           text=element_text(size=15))+
#     ggtitle(paste("FY",fys[i]))
#   ggsave(g[[i]],filename=paste0("./Model Outputs/Plots/Maps/pabs_fy",fys[i],".jpeg"),
#          device="jpeg",width=10,height=5,units="in")
# }
# g[[3]]

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

ggsave(filename=paste0("./Model Outputs/Plots/",subfolder,"/pabs_all_yrs.jpeg"),
       device="jpeg",width=8,height=10,units="in")

##by elimination area and fy ------------------------------
pabs_sum_fy <- pabs_sum_fy %>% st_intersection(elim_areas %>% select(Area_Name))

elim_thresh <- 0.95
  
# pabs_prop_thresh <- pabs_sum_fy %>% st_drop_geometry() %>% 
#   group_by(Area_Name,fy) %>% 
#   summarise(mn_thresh=sum(fy_mn>elim_thresh),
#             lci_thresh=sum(fy_lci>elim_thresh),
#             uci_thresh=sum(fy_uci>elim_thresh),
#             prop_mn_thresh=mn_thresh/n(),
#             prop_lci_thresh=lci_thresh/n(),
#             prop_uci_thresh=uci_thresh/n())
# pabs_prop_thresh$Area_Name[pabs_prop_thresh$Area_Name=="0"] <- "Outside EAs"

# ggplot()+
#   geom_line(data=pabs_prop_thresh,
#               aes(x=factor(fy),y=prop_mn_thresh,col=Area_Name,group=Area_Name),lwd=1)+
#   geom_errorbar(data=pabs_prop_thresh,
#                 aes(x=factor(fy),ymin=prop_lci_thresh,ymax=prop_uci_thresh,col=Area_Name),
#                 width=0)+
#   facet_wrap(.~Area_Name)+
#   xlab("FY")+
#   ylab(paste0("Prop. of watersheds with p(elim) >",elim_thresh))+
#   scale_color_discrete(name="Elimination Area")+
#   theme(text=element_text(size=15))+
#   guides(color="none")

# ggsave(filename=paste0("./Model Outputs/Plots/",subfolder,"/pelim_elim_areas_fy_thresh.jpeg"),
#        device="jpeg",width=10,height=5,units="in")

##by elimination area and season --------------------
# pabs_sum_sf <- pabs_sum_sf %>% st_intersection(elim_areas)

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

# pabs_prop_thresh_season <- pabs_sum_sf %>% st_drop_geometry() %>% 
#   group_by(Area_Name,period_idx) %>% 
#   summarise(Season=min(per_start),
#             mn_thresh=sum(mn>elim_thresh),
#             lci_thresh=sum(lci>elim_thresh),
#             uci_thresh=sum(uci>elim_thresh),
#             prop_mn_thresh=mn_thresh/n(),
#             prop_lci_thresh=lci_thresh/n(),
#             prop_uci_thresh=uci_thresh/n())
# pabs_prop_thresh_season$Area_Name[is.na(pabs_prop_thresh_season$Area_Name)] <- "Outside EAs"

# axis_scale <- 15000

ggplot(pabs_thresh%>% filter(Area_Name!="0"))+
  geom_ribbon(aes(x=per_start,ymin=lci,ymax=uci),alpha=0.3)+
  geom_line(aes(x=per_start,y=md,col=Area_Name),lwd=1.5)+
  # geom_bar(data=dat_rem_sum,aes(x=month,y=removal/axis_scale),
  #          stat="identity",alpha=0.7)+
  facet_wrap(.~area_name_label)+
  scale_y_continuous(name=paste0("Proportion of watersheds with p(elim)>",elim_thresh)#,
                     # sec.axis = sec_axis(~.*axis_scale,name="No. feral swine removed",
                     #                     breaks=c(0,2500,5000))
                     )+
  xlab("Season")+
  scale_color_discrete(name="Elimination area")+
  # scale_fill_discrete(name="Removal type",labels=c("Aerial","Ground","Trap"))+
  theme(text=element_text(size=15))+
  guides(color="none")

ggsave(filename=paste0("./Model Outputs/Plots/",subfolder,"/pabs_over",elim_thresh,"_season_ea.jpeg"),
       device="jpeg",width=10,height=5,units="in")

## extent of p(abs) threshold ----------------------
### entire site ----------------
pabs_long$above_thresh <- ifelse(pabs_long$value>elim_thresh,1,0)
pabs_long <- pabs_long %>% left_join(study_site_grid %>% 
                          st_drop_geometry() %>% 
                          rename(site_idx=SiteID) %>%
                          select(site_idx,area_km))
tot_area_km <- as.numeric(study_site_grid %>% st_union() %>% st_area()/1e6)

pabs_thresh_ext <- pabs_long %>% 
  group_by(above_thresh,period_idx,samp_idx) %>% 
  summarise(area_above_thresh=sum(area_km)) %>% 
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

ggsave(filename=paste0("./Model outputs/Plots/",subfolder,"/km2_eliminated_",elim_thresh,".jpeg"),
       width=7,height=5,units="in",device="jpeg")

ggplot(pabs_thresh_ext)+ 
  geom_ribbon(aes(x=per_start,ymin=lci_prop,ymax=uci_prop),alpha=0.5)+
  geom_line(aes(x=per_start,y=md_prop),lwd=1)+
  xlab("Season")+
  ylab("Proportion of study area with p(elimination) > 0.95")+
  theme(text=element_text(size=15))

ggsave(filename=paste0("./Model outputs/Plots/",subfolder,"/prop_eliminated_",elim_thresh,".jpeg"),
       width=7,height=5,units="in",device="jpeg")

### by EA ----------------
pabs_thresh_ext_ea <- pabs_long %>% 
  group_by(above_thresh,Area_Name,period_idx,samp_idx) %>% 
  summarise(area_above_thresh=sum(area_km)) %>% 
  filter(above_thresh==1) %>% 
  group_by(period_idx,Area_Name) %>% 
  summarise(mn=mean(area_above_thresh),
            md=median(area_above_thresh),
            lci=quantile(area_above_thresh,0.025),
            uci=quantile(area_above_thresh,0.975),
            mn_prop=mn/tot_area_km,
            md_prop=md/tot_area_km,
            lci_prop=lci/tot_area_km,
            uci_prop=uci/tot_area_km)%>% 
  left_join(dat_occ %>% select(period_idx,per_start) %>% distinct())
pabs_thresh_ext_ea$area_name_label <- paste("Elimination Area",pabs_thresh_ext_ea$Area_Name)

# ggplot(pabs_thresh_ext_ea)+ 
#   geom_ribbon(aes(x=per_start,ymin=lci,ymax=uci),alpha=0.5)+
#   geom_line(aes(x=per_start,y=mn),lwd=1)+
#   facet_wrap(.~Area_Name)+
#   xlab("Season")+
#   ylab(expression(paste("K",m^2," with > 0.95 p(elimination)")))
# 
# ggsave(filename=paste0("./Model outputs/Plots/",subfolder,"/km2_eliminated_",elim_thresh,"_ea.jpeg"),
#        width=7,height=5,units="in",device="jpeg")

ggplot(pabs_thresh_ext_ea %>% filter(Area_Name!="0"))+ 
  geom_ribbon(aes(x=per_start,ymin=lci_prop,
                  ymax=uci_prop),alpha=0.3)+
  geom_line(aes(x=per_start,y=md_prop,col=Area_Name),lwd=1)+
  guides(color="none")+
  facet_wrap(.~area_name_label)+
  xlab("Season")+
  ylab("Proportion of area with p(elimination) > 0.95")+
  theme(text=element_text(size=15))

ggsave(filename=paste0("./Model outputs/Plots/",subfolder,"/prop_eliminated_",elim_thresh,"_ea.jpeg"),
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
  summarise(area_above_thresh=sum(area_km)) %>% 
  group_by(fy) %>% 
  summarise(mn=mean(area_above_thresh),
            md=median(area_above_thresh),
            lci=quantile(area_above_thresh,0.025),
            uci=quantile(area_above_thresh,0.975),
            mn_prop=mn/tot_area_km,
            md_prop=md/tot_area_km,
            lci_prop=lci/tot_area_km,
            uci_prop=uci/tot_area_km)
write.csv(pabs_thresh_yr,paste0("./Model outputs/Plots/",subfolder,"/area_above_,",elim_thresh,"_table.csv"))

### by year and ea ------------------
pabs_thresh_yr_ea <- pabs_long %>% 
  group_by(site_idx,Area_Name,fy,samp_idx) %>% 
  summarise(yr_mn=mean(value),
            area_km=unique(area_km),
            above_thresh=ifelse(yr_mn>elim_thresh,1,0)) %>% 
  filter(above_thresh==1) %>% 
  group_by(fy,Area_Name,samp_idx) %>% 
  summarise(area_above_thresh=sum(area_km)) %>% 
  group_by(fy,Area_Name) %>% 
  summarise(mn=mean(area_above_thresh),
            md=median(area_above_thresh),
            lci=quantile(area_above_thresh,0.025),
            uci=quantile(area_above_thresh,0.975),
            mn_prop=mn/tot_area_km,
            md_prop=md/tot_area_km,
            lci_prop=lci/tot_area_km,
            uci_prop=uci/tot_area_km)
pabs_thresh_yr_ea_prop <- pabs_thresh_yr_ea %>% 
  pivot_longer(cols=c(mn,md,lci,uci,mn_prop,md_prop,lci_prop,uci_prop),
               names_to="stat",values_to="value") %>%
  pivot_wider(names_from=Area_Name,values_from=value) %>%
  filter(grepl("prop",stat))

write.csv(pabs_thresh_yr_ea_prop,paste0("./Model outputs/Plots/",subfolder,"/area_above_,",elim_thresh,"_ea_prop_table.csv"))

pabs_thresh_yr_ea_km <- pabs_thresh_yr_ea %>% 
  pivot_longer(cols=c(mn,md,lci,uci,mn_prop,md_prop,lci_prop,uci_prop),
               names_to="stat",values_to="value") %>%
  pivot_wider(names_from=Area_Name,values_from=value) %>%
  filter(!grepl("prop",stat))

write.csv(pabs_thresh_yr_ea_km,paste0("./Model outputs/Plots/",subfolder,"/area_above_,",elim_thresh,"_ea_km_table.csv"))


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

ggsave(filename=paste0("./Model Outputs/Plots/",subfolder,"/pabs_md_season_ea.jpeg"),
       device="jpeg",width=10,height=6,units="in")

pabs_sum_ea %>% filter(grepl("4",Area_Name) | 
                         grepl("6",Area_Name))
pabs_sum_ea[which.min(pabs_sum_ea$md),]
pabs_sum_ea[which.max(pabs_sum_ea$md),]
pabs_sum_ea %>% filter(Area_Name=="Elimination Area 1")

#calculate by fy all sites
pabs_sum_fy %>% 
  # filter(!is.na(Area_Name)) %>% 
  st_drop_geometry() %>% 
  group_by(fy) %>% 
  summarise(median(fy_md),
            mean(fy_mn))

#calculate by fy and elim area
pabs_sum_fy <- pabs_sum_fy %>% st_drop_geometry() %>% 
  group_by(fy,Area_Name) %>% 
  summarise(mn=mean(fy_mn),
            md=median(fy_md))
pabs_sum %>% 
  summarise(mn=mean(mn))

#detection probability with trap effort -------------------
sys_eff_sc <- seq(min(dat_occ$trap_nights_km_sc),max(dat_occ$trap_nights_km_sc),by=0.5)
sys_eff <- sys_eff_sc * attr(dat_occ$trap_nights_km_sc,"scaled:scale") +
  attr(dat_occ$trap_nights_km_sc,"scaled:center")

det <- sapply(1:(nmcmc*nChains),function(i)boot::inv.logit(alpha[i,1] + alpha[i,2]*sys_eff_sc))

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

ggsave(filename=paste0("./Model Outputs/Plots/",subfolder,"/occ_det_eff_curve.jpeg"),
       device="jpeg",width=7,height=5,units="in")

#conditional absence probability maps --------------------
## map to study site -------------------
pelim <- do.call("rbind",lapply(1:nChains,function(i){samples[[i]][grep("pelim",colnames(samples[[i]]))]}))

pelim_sum <- data.frame(idx=colnames(pelim),
                       mn=colMeans(pelim),
                       md=sapply(1:ncol(pelim),function(i)quantile(pelim[,i],prob=0.5)),
                       lci=sapply(1:ncol(pelim),function(i)quantile(pelim[,i],prob=0.025)),
                       uci=sapply(1:ncol(pelim),function(i)quantile(pelim[,i],prob=0.975)))

pelim_sum$site_idx <- sapply(1:nrow(pelim_sum),function(x){
  as.numeric(gsub("\\[","",gsub("pelim","",unlist(str_split(pelim_sum[x,"idx"],","))[1])))})

pelim_sum$period_idx <- sapply(1:nrow(pelim_sum),function(x){
  as.numeric(gsub("\\]","",unlist(str_split((pelim_sum[x,"idx"]),","))[2]))})

pelim_sum$per_start <- sapply(1:nrow(pelim_sum),function(i){
  min(dat_occ$per_start[dat_occ$period_idx==pelim_sum$period_idx[i]])})
pelim_sum$per_start <- as.Date(pelim_sum$per_start,origin='1970-01-01')
pelim_sum$year <- year(pelim_sum$per_start)
pelim_sum$fy <- pelim_sum$year
pelim_sum$fy[month(pelim_sum$per_start)%in%c(10:12)] <- 
  pelim_sum$fy[month(pelim_sum$per_start)%in%c(10:12)] +1

pelim_sum_sf <- pelim_sum %>% 
  left_join(study_site_grid %>% rename(site_idx=SiteID)) %>% 
  st_as_sf()

mean(rowMeans(pelim))
quantile(rowMeans(pelim),c(0.025,0.975))

pelim_sum[which.max(pelim_sum$md),]
pelim_sum[which.min(pelim_sum$md),]

pelim_sum <- pelim_sum %>% 
  left_join(study_site_grid %>% 
              st_drop_geometry() %>% 
              rename(site_idx=SiteID) %>% 
              select(site_idx,Area_Name))
pabs_ea_sum <- pelim_sum %>% 
  group_by(Area_Name,per_start) %>% 
  summarise(ea_md=median(md))

pabs_ea_sum[which.max(pabs_ea_sum$ea_md),]
pabs_ea_sum[which.min(pabs_ea_sum$ea_md),]

## fiscal year summary ------------------
###individiual plots -------------------
pelim_sum_fy <- pelim_sum %>% group_by(fy,site_idx) %>% 
  summarise(fy_mn = mean(mn),
            fy_md= median(md)) %>% 
  left_join(study_site_grid %>% rename(site_idx=SiteID)) %>% 
  st_as_sf()

# fys <- unique(pelim_sum_fy$fy)
# g<-list()
# for(i in 1:length(fys)){
#   g[[i]] <- ggplot()+
#     geom_sf(data=pelim_sum_fy %>% 
#               filter(fy==fys[i]),aes(fill=fy_mn))+
#     geom_sf(data=elim_areas %>% filter(Area_Name!="0"),aes(col=Area_Name),
#             fill="transparent",lwd=1.5)+
#     scale_fill_viridis_c(#limits=c(0,1),
#                          direction=-1,
#                          name="Feral swine \nelimination probability \nafter 10 negative \ntrap nights")+
#     scale_color_discrete(name="Elimination Area")+
#     theme(panel.background=element_blank(),
#           axis.ticks=element_blank(),
#           axis.text=element_blank(),
#           plot.margin = margin(1,1,1,1, "cm"),
#           text=element_text(size=15))+
#     ggtitle(paste("FY",fys[i]))
#   ggsave(g[[i]],filename=paste0("./Model Outputs/Plots/",subfolder,"/pelim_fy",fys[i],".jpeg"),
#          device="jpeg",width=10,height=5,units="in")
# }
# g[[i]]

### facet version ------------------
ggplot()+
  geom_sf(data=pelim_sum_fy,aes(fill=fy_md))+
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

ggsave(filename=paste0("./Model Outputs/Plots/",subfolder,"/pelim_all_yrs.jpeg"),
       device="jpeg",width=8,height=10,units="in")

#calculate by fy all sites
pelim_sum_fy %>% 
  st_drop_geometry() %>% 
  group_by(fy) %>% 
  summarise(median(fy_md),
            mean(fy_mn))

##calculate by fy and elim area --------------------
pelim_sum_fy <- pelim_sum_fy %>% st_intersection(elim_areas %>% select(Area_Name))
pelim_ea_fy <- pelim_sum_fy %>% st_drop_geometry() %>% 
  group_by(fy,Area_Name) %>% 
  summarise(mn=mean(fy_mn),
            md=median(fy_md))

ggplot(pelim_ea_fy)+
  geom_line(aes(x=factor(fy),y=mn,col=Area_Name,group=Area_Name),lwd=1)+
  ylab("Mean p(elim)")+xlab("FY")+
  scale_color_discrete(name="Elimination Area")+
  theme(text=element_text(size=15))

ggplot(pelim_ea_fy)+
  geom_line(aes(x=factor(fy),y=md,col=Area_Name,group=Area_Name),lwd=1)+
  ylab("Median p(elim)")+xlab("FY")+
  scale_color_discrete(name="Elimination Area")+
  theme(text=element_text(size=15))

##calculate by season and elim area------------------------
pelim_sum_sf <- pelim_sum_sf %>% st_intersection(elim_areas %>% select(Area_Name))
pelim_sum_sf <- pelim_sum_sf %>% left_join(dat_occ %>% select(period_idx,per_start) %>% distinct())
pelim_ea <- pelim_sum_sf %>% st_drop_geometry() %>% 
  group_by(period_idx,Area_Name) %>% 
  summarise(season=min(per_start),
            mn=mean(mn),
            md=median(md))

# ggplot(pelim_ea)+geom_line(aes(x=season,y=mn,col=Area_Name),lwd=1)+
#   ylab("Mean p(elim)")+xlab("Season")+
#   scale_color_discrete(name="Elimination Area")
# 
# ggplot(pelim_ea)+geom_line(aes(x=season,y=md,col=Area_Name),lwd=1)+
#   ylab("Median p(elim)")+xlab("Season")+
#   scale_color_discrete(name="Elimination Area")

#number of samples to elimination --------------------
syseff <- data.frame(trap_nights=c((0.25/2.59)*45,
                                   (0.5/2.59)*45,
                                   (1/2.59)*45))
syseff$trap_nights_sc <- (syseff$trap_nights-attr(dat_occ$trap_nights_km_sc,"scaled:center"))/
  attr(dat_occ$trap_nights_km_sc,"scaled:scale") 

sd_pdet <-  do.call("rbind",lapply(1:nChains,function(i){samples[[i]][grep("sd_pdet",colnames(samples[[i]]))]}))

pabs_ex <- c(0.25,0.5,0.75)
elim_prob <- 0.95
eff_elim <- array(NA,dim=c(nmcmc*nChains,length(pabs_ex),nrow(syseff)))
eff_elim_sum <- list()
for(k in 1:nrow(syseff)){
  for(j in 1:length(pabs_ex)){
    for(i in 1:(nmcmc*nChains)){
      p_sys <- boot::inv.logit(rnorm(1,alpha[i,1] + alpha[i,2]*syseff$trap_nights_sc[k],sd_pdet[i,1]))
      eff_elim[i,j,k] <- log((pabs_ex[j]*(1-elim_prob))/(elim_prob*(1-pabs_ex[j])))/(log(1-p_sys)) 
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
  
ggsave(filename = paste0("./Model outputs/Plots/",subfolder,"/nweeks_to_elim_trapeff.jpeg"),
       width=7,height=5,units="in",device="jpeg")

#abundance  -----------------------------
elim_areas$area_km <- as.numeric(st_area(elim_areas)/1e6)
N <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
  cbind.data.frame(samples[[i]][grepl("N",colnames(samples[[i]])) & 
                                  !grepl("N_latent",colnames(samples[[i]]))],
                   chain=i)})),
  samp=rep(1:nmcmc,nChains))

# N <- do.call("rbind",lapply(1:nChains,function(i){
#   samples[[i]][grepl("N",colnames(samples[[i]])) & !grepl("N_latent",colnames(samples[[i]]))]}))

N <- N[,grepl("\\[5,",colnames(N)) | grepl("\\[7,",colnames(N)) |
         grepl("samp",colnames(N)) | grepl("chain",colnames(N))]
colnames(N)
N_long <- N %>% pivot_longer(cols=1:(ncol(N)-2),
                             values_to = "N",
                             names_to="idx") %>% 
  select(-idx)
N_long$elim_area_idx <- rep(rep(c(5,7),nperiods),nmcmc*nChains)
N_long$period_idx <- rep(sort(rep(1:nperiods,2)),nmcmc*nChains)
N_long <- N_long %>%
  left_join(elim_areas %>% st_drop_geometry())

per_idx <- dat_occ %>% select(period_idx,per_start) %>% 
  mutate(year=year(per_start)) %>% distinct()
per_idx$fy <- per_idx$year
per_idx$fy[month(per_idx$per_start)%in%c(10,11,12)] <- per_idx$fy[month(per_idx$per_start)%in%c(10,11,12)] +1

N_long <- N_long %>% left_join(per_idx)
lambda_N <- lambda
colnames(lambda_N) <- c("4","6","chain","samp")
lambda_N <- lambda_N %>% pivot_longer(cols=c("4","6"),
                                      names_to="Area_Name",
                                      values_to="lambda")
N_long <- N_long %>% left_join(lambda_N)

#proportion needed to combat pop growth
N_long$lambda_prop <- rpois(nrow(N_long),N_long$N * (N_long$lambda-1))
N_long$lambda_prop_dens <- N_long$lambda_prop/N_long$area_km
N_long$lambda_prop[is.na(N_long$lambda_prop)] <- 0
N_long$lambda_prop_dens[is.na(N_long$lambda_prop_dens)] <- 0

#standardizes abundance 
N_ss <- sqrt(sum(N_long$N^2)/(length(N_long$N)-1))
N_long$N_std <- N_long$N/N_ss
N_long$lambda_prop_std <- N_long$lambda_prop/N_ss

#projecting pop growth sans removal 
N_long$N_no_rem <- NA
N_long$N_no_rem[N_long$period_idx==1] <- N_long$N[N_long$period_idx==1] 
for(i in 2:nperiods){
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
  group_by(metric,Area_Name,period_idx,per_start,year,fy) %>% 
  summarise(mn=mean(value),
            md=median(value),
            lci=quantile(value,0.025),
            uci=quantile(value,0.975)) %>% 
  pivot_wider(names_from=metric,values_from = c(mn,md,lci,uci))

N_sum$Area_Name_label <- paste("Elimination Area",N_sum$Area_Name)
dat_rem_sum$Area_Name_label <- paste("Elimination Area",dat_rem_sum$Area_Name)

# N_sum_sf <- N_sum %>% 
#   left_join(elim_areas) %>% 
#   st_as_sf() %>% 
  # left_join(dat_rem %>% ungroup() %>% 
  #             group_by(period_idx) %>% 
  #             dplyr::select(period_idx,Date) %>% 
  #             summarise(per_start=min(Date)))

N_sum_4 <- N_sum %>% filter(Area_Name==4)
N_sum_4[which.min(N_sum_4$md_dens),"per_start"]
round(N_sum_4[which.min(N_sum_4$md_dens),c("md_dens","lci_dens","uci_dens")],2)
N_sum_4[which.max(N_sum_4$md_dens),"per_start"]
round(N_sum_4[which.max(N_sum_4$md_dens),c("md_dens","lci_dens","uci_dens")],2)

N_sum_6 <- N_sum %>% filter(Area_Name==6)
N_sum_6[which.min(N_sum_6$md_dens),"per_start"]
round(N_sum_6[which.min(N_sum_6$md_dens),c("md_dens","lci_dens","uci_dens")],2)
N_sum_6[which.max(N_sum_6$md_dens),"per_start"]
round(N_sum_6[which.max(N_sum_6$md_dens),c("md_dens","lci_dens","uci_dens")],2)

####raw abundance EA 4 and 6 -------------------------
dat_rem_sum <- dat_rem_sum %>% left_join(elim_areas %>% st_drop_geometry())
dat_rem_sum$tot_rem_km <- dat_rem_sum$removal/dat_rem_sum$area_km

axis_trans<- 0.1
ggplot()+  
  geom_ribbon(data=N_sum,
              aes(x=per_start,ymin=lci_N,ymax=uci_N,fill="Estimated Abundance"),alpha=0.2)+
  geom_line(data=N_sum,
            aes(x=per_start,y=md_N,col="Estimated Abundance"),lwd=1)+
  geom_bar(data=dat_rem_sum,
           aes(x=month,y=removal/axis_trans),stat="identity",
           alpha=0.6)+
  geom_ribbon(data=N_sum,
              aes(x=per_start,
                  ymin=lci_lambda_prop/axis_trans,ymax=uci_lambda_prop/axis_trans,
                  fill="Removal Needed"),
              alpha=0.1)+
  geom_line(data=N_sum,aes(x=per_start,y=md_lambda_prop/axis_trans,group=Area_Name,
                           col="Removal Needed"),
            lwd=1,lty=2)+
  geom_ribbon(data=N_sum,
              aes(x=per_start,ymin=lci_N_no_rem,ymax=uci_N_no_rem,
                  fill="No Removal"),
              alpha=0.1)+
  geom_line(data=N_sum,
            aes(x=per_start,y=md_N_no_rem,col="No Removal"),
            lwd=1,lty=3)+
  facet_wrap(.~Area_Name_label)+
  scale_y_continuous(sec.axis = sec_axis(transform=~ . * axis_trans,
                        name = "No. feral swine removed"
                        ))+
  scale_fill_manual(name="",
                    values=c("Estimated Abundance"="black",
                             "No Removal"="red",
                             "Removal Needed"="blue"))+
  scale_color_manual(name="",
                     values=c("Estimated Abundance"="black",
                              "No Removal"="red",
                              "Removal Needed"="blue"))+
  # scale_fill_manual(name="Removal method",values=c("lightskyblue","mediumblue"))+
  ylab("Feral swine abundance")+
  xlab("Season")+
  # scale_color_manual(name="Elimination area",
  #                    values=scales::hue_pal()(6)[c(4,6)])+
  # guides(color="none")+
  theme(text=element_text(size=15))
ggsave(filename=paste0("./Model Outputs/Plots/",subfolder,"/abundance_trend_removal.jpeg"),
       device="jpeg",width=10,height=6,units="in")

#### density EA 4 and 6 -------------------------
axis_trans_d<- 0.1
ggplot()+  
  geom_ribbon(data=N_sum,
              aes(x=per_start,ymin=lci_dens,ymax=uci_dens,
                  fill="Estimated Density"),alpha=0.2)+
  geom_line(data=N_sum,
            aes(x=per_start,y=md_dens,col="Estimated Density"),lwd=1)+
  geom_bar(data=dat_rem_sum,
           aes(x=month,y=tot_rem_km/axis_trans_d),stat="identity",
           alpha=0.6)+
  geom_ribbon(data=N_sum,
              aes(x=per_start,
                  ymin=lci_lambda_prop_dens/axis_trans_d,
                  ymax=uci_lambda_prop_dens/axis_trans_d,
                  fill="Removal Needed"),
              alpha=0.1)+
  geom_line(data=N_sum,aes(x=per_start,
                           y=md_lambda_prop_dens/axis_trans_d,group=Area_Name,
                           col="Removal Needed"),lwd=1,lty=2)+
  geom_ribbon(data=N_sum,
              aes(x=per_start,ymin=lci_dens_no_rem,ymax=uci_dens_no_rem,fill="No Removal"),
              alpha=0.1)+
  geom_line(data=N_sum,
            aes(x=per_start,y=md_dens_no_rem,col="No Removal"),lwd=1,lty=3)+
  facet_wrap(.~Area_Name_label)+
  scale_y_continuous(
    sec.axis = sec_axis(transform=~.*axis_trans_d,
                        name=expression(paste("Feral swine removed (swine/k",m^2,")"))
    ))+
  scale_fill_manual(name="",
                     values=c("Estimated Density"="black",
                              "No Removal"="red",
                              "Removal Needed"="blue"))+
  scale_color_manual(name="",
                    values=c("Estimated Density"="black",
                             "No Removal"="red",
                             "Removal Needed"="blue"))+
  # scale_fill_manual(name="Removal method",
  #                   values=c("lightskyblue","mediumblue"))+
  ylab(expression(paste("Feral swine density (swine/k",m^2,")")))+
  xlab("Season")+
  # scale_color_manual(name="Elimination area",
  #                    values=scales::hue_pal()(6)[c(4,6)])+
  # guides(color="none")+
  theme(text=element_text(size=15))
ggsave(filename=paste0("./Model Outputs/Plots/",subfolder,"/density_trend_removal.jpeg"),
       device="jpeg",width=10,height=6,units="in")

####standardized abundance EA 4 and 6 -------------------------
axis_scale <- 6000
ggplot()+  
  geom_ribbon(data=N_sum,
              aes(x=per_start,ymin=lci_N_std,ymax=uci_N_std,
                  fill="Estimated \nStandardized Abundance"),alpha=0.2)+
  geom_line(data=N_sum,
            aes(x=per_start,y=md_N_std,col="Estimated \nStandardized Abundance"),lwd=1)+
  geom_bar(data=dat_rem_sum,
           aes(x=month,y=removal/axis_scale),stat="identity",
           alpha=0.6)+
  geom_ribbon(data=N_sum,
              aes(x=per_start,
                  ymin=lci_lambda_prop/axis_scale,
                  ymax=uci_lambda_prop/axis_scale,
                  fill="Removal Needed"),
              alpha=0.2)+
  geom_line(data=N_sum,aes(x=per_start,
                           y=md_lambda_prop/axis_scale,
                           group=Area_Name,col="Removal Needed"),
            lwd=1,lty=2)+
  geom_ribbon(data=N_sum,
              aes(x=per_start,ymin=lci_N_ss_no_rem,ymax=uci_N_ss_no_rem,
                  fill="No Removal"),
              alpha=0.1)+
  geom_line(data=N_sum,
            aes(x=per_start,y=md_N_ss_no_rem,col="No Removal"),lty=3,lwd=1)+
  facet_wrap(.~Area_Name)+
  scale_y_continuous(sec.axis = sec_axis(transform=~.*axis_scale, 
                                         name = "No. feral swine removed"
                                         ))+
  scale_fill_manual(name="",
                    values=c("Estimated \nStandardized Abundance"="black",
                             "No Removal"="red",
                             "Removal Needed"="blue"))+
  scale_color_manual(name="",
                     values=c("Estimated \nStandardized Abundance"="black",
                              "No Removal"="red",
                              "Removal Needed"="blue"))+
  # scale_fill_manual(name="Removal method",values=c("lightskyblue","mediumblue"))+
  ylab("Standardized feral swine abundance")+
  xlab("Season")+
  # scale_color_manual(name="Elimination area",
  #                    values=scales::hue_pal()(6)[c(4,6)])+
  # guides(color="none")+
  theme(text=element_text(size=15))

ggsave(filename=paste0("./Model Outputs/Plots/",subfolder,"/std_abundance_trend_removal.jpeg"),
       device="jpeg",width=10,height=6,units="in")

### annual ------------------------
N_yr <- N_long %>% filter(month(per_start)==10 | (month(per_start)==7 & year==2023))
N_yr$per_temp <- as.numeric(as.factor(N_yr$period_idx))

N_yr <- N_yr %>% filter(elim_area_idx%in%c(5,7))
N_yr$N_prev<- NA
for(i in 1:nrow(N_yr)){
  if(N_yr$per_temp[i]>1){
    N_yr$N_prev[i] <- N_yr$N[N_yr$per_temp==(N_yr$per_temp[i]-1) &
                                       N_yr$elim_area_idx==N_yr$elim_area_idx[i] &
                                       N_yr$samp==N_yr$samp[i] &
                                       N_yr$chain==N_yr$chain[i]]
  }
}

N_yr$per_change <- (N_yr$N-N_yr$N_prev)/N_yr$N_prev

N_yr_sum <- N_yr %>% 
  filter(year>2020) %>% 
  group_by(Area_Name,period_idx,per_start) %>% 
  summarise(mn=mean(per_change,na.rm=T),
            md=median(per_change,na.rm=T),
            lci=quantile(per_change,0.025,na.rm=T),
            uci=quantile(per_change,0.975,na.rm=T)) 

write.csv(N_yr_sum,paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/Plots/",subfolder,"/N_change_table.csv"))


# removal probability -----------------
## aerial ----------------
aerial_eff_sc <- seq(min(dat_aerial$eff_area_hrs),max(dat_aerial$eff_area_hrs),by=0.005)
aerial_eff <- aerial_eff_sc * attr(dat_aerial$eff_area_hrs_sc,"scaled:scale") +
  attr(dat_aerial$eff_area_hrs_sc,"scaled:center")
aerial_area_mn <- mean(dat_aerial$prop_ea_impact[dat_aerial$prop_ea_impact!=0])

det_a <- sapply(1:(nmcmc*nChains),function(i){aerial_area_mn*
    boot::inv.logit(delta$`Aerial Intercept`[i] + delta$`Aerial Slope`[i]*aerial_eff_sc)})

det_a_sum <- data.frame(eff_area_hrs=aerial_eff,
                        mn=rowMeans(det_a),
                        md=sapply(1:nrow(det_a),function(i)quantile(det_a[i,],prob=0.5)),
                        lci=sapply(1:nrow(det_a),function(i)quantile(det_a[i,],prob=0.025)),
                        uci=sapply(1:nrow(det_a),function(i)quantile(det_a[i,],prob=0.975)))

g_a <- ggplot(det_a_sum)+
  geom_ribbon(aes(x=eff_area_hrs,ymin=lci,ymax=uci),alpha=0.3)+
  geom_line(aes(x=eff_area_hrs,y=md))+
  ylab("Detection probability")+
  xlab(expression(paste("Hours flown per k",m^2)))+
  theme(text=element_text(size=15))

## trap ----------------------
trap_eff_sc <- seq(min(dat_trap$eff_area_hrs),max(dat_trap$eff_area_hrs),by=0.05)
trap_eff <- trap_eff_sc * attr(dat_trap$eff_area_hrs_sc,"scaled:scale") +
  attr(dat_trap$eff_area_hrs_sc,"scaled:center")
trap_area_mn <- mean(dat_trap$prop_ea_impact[dat_trap$prop_ea_impact!=0])

det_t <- sapply(1:(nmcmc*nChains),function(i){trap_area_mn*
    boot::inv.logit(delta$`Trap Intercept`[i] + delta$`Trap Slope`[i]*trap_eff_sc)})

det_t_sum <- data.frame(eff_area_hrs=trap_eff,
                        mn=rowMeans(det_t),
                        md=sapply(1:nrow(det_t),function(i)quantile(det_t[i,],prob=0.5)),
                        lci=sapply(1:nrow(det_t),function(i)quantile(det_t[i,],prob=0.025)),
                        uci=sapply(1:nrow(det_t),function(i)quantile(det_t[i,],prob=0.975)))

g_t <- ggplot(det_t_sum)+
  geom_ribbon(aes(x=eff_area_hrs,ymin=lci,ymax=uci),alpha=0.3)+
  geom_line(aes(x=eff_area_hrs,y=md))+
  ylab("Detection probability")+
  xlab(expression(paste("Hours spent trapping per k",m^2)))+
  theme(text=element_text(size=15))

##all together --------------------
g_all <- gridExtra::arrangeGrob(g_a+ylim(0,0.001), g_t+ylim(0,0.001),nrow=1) #generates g

ggsave(g_all,file=paste0("./Model Outputs/Plots/",subfolder,"/rem_det_curve_all.jpeg"),
       device="jpeg",width=13,height=6,units="in")

##removal comparisons --------------------
a_mn <- max(dat_aerial$prop_ea_impact)*
  boot::inv.logit(delta[,"Aerial Intercept"] + delta[,"Aerial Slope"]* 
                    max(dat_aerial$eff_area_hrs_sc[dat_aerial$eff_area_hrs!=0]))

max(dat_aerial$eff_area_hrs)
max(dat_aerial$prop_ea_impact)

t_mn <- max(dat_trap$prop_ea_impact)*
  boot::inv.logit(delta[,"Trap Intercept"] + delta[,"Trap Slope"]* 
                    max(dat_trap$eff_area_hrs_sc[dat_trap$eff_area_hrs!=0]))
max(dat_trap$eff_area_hrs)
max(dat_trap$prop_ea_impact)

rem_mn <- cbind.data.frame(aerial=a_mn,trap=t_mn)

rem_df <- rem_mn %>% pivot_longer(cols=1:2,names_to="rem_typ",values_to="value") %>% 
  group_by(rem_typ) %>% 
  summarise(mn=mean(value),
            md=median(value),
            lci=quantile(value,0.025),
            uci=quantile(value,0.975))
rem_df$mn[1]/rem_df$mn[2]

range(elim_areas %>% filter(Area_Name!="0") %>% st_area()/1e6)

#population growth rate--------------
if(!exists("lambda")){
  lambda <- do.call("rbind",lapply(1:nChains,function(i){
    samples[[i]][grepl("lambda",colnames(samples[[i]]))]}))
  colnames(lambda)[1] <- "Outside Elimination Areas"
  colnames(lambda)[2:nea] <- paste("Elimination Area",1:(nea-1))
  lambda <- lambda %>% select(c("Elimination Area 4","Elimination Area 6"))
}

lam_sum <- data.frame(mn=colMeans(lambda[,1:2]),
           lci=sapply(1:2,function(i)quantile(lambda[,i],probs=c(0.025))),
           uci=sapply(1:2,function(i)quantile(lambda[,i],probs=c(0.975))))

round(min(N_sum_4$md_lambda_prop),0)
round(N_sum_4$lci_lambda_prop[which.min(N_sum_4$md_lambda_prop)],0)
round(N_sum_4$uci_lambda_prop[which.min(N_sum_4$md_lambda_prop)],0)
round(max(N_sum_4$md_lambda_prop),0)
round(N_sum_4$lci_lambda_prop[which.max(N_sum_4$md_lambda_prop)],0)
round(N_sum_4$uci_lambda_prop[which.max(N_sum_4$md_lambda_prop)],0)

round(min(N_sum_4$md_lambda_prop_dens),2)
round(N_sum_4$lci_lambda_prop_dens[which.min(N_sum_4$md_lambda_prop_dens)],2)
round(N_sum_4$uci_lambda_prop_dens[which.min(N_sum_4$md_lambda_prop_dens)],2)
round(max(N_sum_4$md_lambda_prop_dens),2)
round(N_sum_4$lci_lambda_prop_dens[which.max(N_sum_4$md_lambda_prop_dens)],2)
round(N_sum_4$uci_lambda_prop_dens[which.max(N_sum_4$md_lambda_prop_dens)],2)

round(min(N_sum_6$md_lambda_prop),0)
round(N_sum_6$lci_lambda_prop[which.min(N_sum_6$md_lambda_prop)],0)
round(N_sum_6$uci_lambda_prop[which.min(N_sum_6$md_lambda_prop)],0)
round(max(N_sum_6$md_lambda_prop),0)
round(N_sum_6$lci_lambda_prop[which.max(N_sum_6$md_lambda_prop)],0)
round(N_sum_6$uci_lambda_prop[which.max(N_sum_6$md_lambda_prop)],0)

round(min(N_sum_6$md_lambda_prop_dens),2)
round(N_sum_6$lci_lambda_prop_dens[which.min(N_sum_6$md_lambda_prop_dens)],2)
round(N_sum_6$uci_lambda_prop_dens[which.min(N_sum_6$md_lambda_prop_dens)],2)
round(max(N_sum_6$md_lambda_prop_dens),2)
round(N_sum_6$lci_lambda_prop_dens[which.max(N_sum_6$md_lambda_prop_dens)],2)
round(N_sum_6$uci_lambda_prop_dens[which.max(N_sum_6$md_lambda_prop_dens)],2)


#pig removals ---------------------
##by EA--------------------
# dat_rem_sum %>% 
#   filter(Area_Name!="0") %>% 
#   ggplot()+geom_bar(aes(y=removal,x=month,fill=Area_Name),
#                     stat="identity")+
#   scale_fill_discrete(name="Elimination Area")+
#   xlab("Season")+ylab("No. feral swine removed")+
#   theme(text=element_text(size=15))
dat_rem_sum$Area_Name_label <- paste("Eliminaton Area",dat_rem_sum$Area_Name)
dat_rem_sum_all$Area_Name_label <- paste("Eliminaton Area",dat_rem_sum_all$Area_Name)

dat_rem_sum_all %>% 
  filter(Area_Name!="0") %>% 
  ggplot()+geom_bar(aes(y=removal,x=month,fill=method),
                    stat="identity")+
  facet_wrap(.~Area_Name_label)+
  scale_fill_discrete(name="Removal method")+
  xlab("Season")+ylab("No. feral swine removed")+
  theme(text=element_text(size=15))
ggsave(filename = "./Model outputs/Plots/Raw Data/raw_removal_ea.jpeg",
       device="jpeg",height=5,width=7,units="in")

#effort ---------------------
dat_rem_sum_all %>% 
  filter(Area_Name!="0") %>% 
  ggplot()+geom_bar(aes(y=tot_hrs,x=month,fill=method),
                    stat="identity")+
  facet_wrap(.~Area_Name_label)+
  scale_fill_discrete(name="Removal method")+
  xlab("Season")+ylab("No. hours spent removing feral swine")+
  theme(text=element_text(size=15))
ggsave(filename = "./Model outputs/Plots/Raw Data/raw_effort_hours.jpeg",
       device="jpeg",height=5,width=7,units="in")

#raw cpue -------------------
dat_rem_sum_all %>% 
  filter(Area_Name!="0") %>% 
  ggplot()+geom_bar(aes(y=removal/tot_hrs,x=month,fill=method),
                    stat="identity",position="dodge")+
  facet_wrap(.~Area_Name_label)+
  scale_fill_discrete(name="Removal method")+
  xlab("Season")+ylab("No. of pigs removed per hour of effort")+
  theme(text=element_text(size=15))

ggsave(filename = "./Model outputs/Plots/Raw Data/raw_cpue.jpeg",
       device="jpeg",height=5,width=7,units="in")


# summary stats ------------------
##raw data ------------------------
dat_rem %>% 
  group_by(method) %>% 
  summarise(sum(tot_rem),
            sum(tot_hrs))

sysbait_det_eff %>% 
  ungroup() %>% 
  summarise(dets=sum(detection),
            trap_nights=sum(trap_nights))

# subfolder <- "NFSP"
#save posterior summaries ---------------------------------
# st_write(eff_sum_sf,paste0("./Model outputs/Plots/",subfolder,"/eff_sum_sf.shp"),append=F)
# st_write(N_sum_sf,paste0("./Model outputs/Plots/",subfolder,"/N_sum_sf.shp"),append=F)
st_write(pabs_sum_sf,paste0("./Model outputs/Plots/",subfolder,"/pabs_sum_sf.shp"),append=F)
save(eff_elim_sum,N_sum,N_yr_sum,pabs_sum_ea,pabs_sum,
     pabs_sum_fy,pabs_thresh,pabs_thresh,pabs_thresh_ext,
     pelim_ea,pelim_ea_fy,pelim_sum,pelim_sum_fy,
     det_a_sum,det_t_sum,
     pabs_thresh_yr,det_sum,nea,nsites,nperiods,nmcmc,nChains,
     dat_occ,dat_rem,dat_rem_sum,elim_thresh,rem_df,
     file=paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/Plots/",subfolder,"/posterior_summaries_07OCT24.RData"))


          