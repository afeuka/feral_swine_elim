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
load("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/ziBinMod_area_25JUL24_multi_lambda_21_24.Rdata")

## removal data -----------
dat_rem <- rbind.data.frame(dat_aerial,dat_trap)
dat_rem_sum <- dat_rem %>% 
  group_by(Area_Name,method,month=floor_date(Date,"month")) %>% 
  summarise(Season=min(Date),
            removal=sum(tot_rem),
            tot_hrs=sum(tot_hrs),
            tot_events=sum(num_events)) %>% 
  filter(month<=max(dat_occ$per_start))

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
# ggsave(filename = "./Model outputs/Plots/Manuscript/mo_map.jpeg",
#        width=7,height=5,units="in",device="jpeg")

#trace plots -------------------------------
##beta coefficients----------------------
beta<- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
  cbind.data.frame(samples[[i]][,grep("beta",colnames(samples[[i]]))],
                   chain=i)})),
  samp=rep(1:nmcmc,nChains))

beta_names <- c("Intercept","Feral swine range","Development","Agriculture")
colnames(beta)[1:nbeta] <- beta_names

beta_long <- beta %>% 
  pivot_longer(cols=all_of(1:nbeta),names_to="beta",values_to="value") %>% 
  mutate(chain=factor(chain))

g_b<- ggplot(beta_long) + geom_line(aes(x=samp,y=value,col=chain))+
  geom_hline(yintercept=0,col="black",lty=2)+
  facet_wrap(.~beta)
g_b

ggsave(g_b, 
       filename = "./Model outputs/Plots/Manuscript/Trace/beta_trace.jpeg",
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
       filename = "./Model outputs/Plots/Manuscript/Trace/n1_params_trace.jpeg",
       width=7,height=5,units="in",device="jpeg")

##catch per unit effort ----------------------
cpue <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
  cbind.data.frame(samples[[i]][,grepl("p_",colnames(samples[[i]])) &
                                  !grepl("p_n1",colnames(samples[[i]]))],
                   chain=i)})),
  samp=rep(1:nmcmc,nChains))

cpue_names <- c("Aerial","Ground","Systematic Baiting","Trap")
colnames(cpue)[1:4] <- cpue_names

cpue_long <- cpue %>% 
  pivot_longer(cols=all_of(1:4),names_to="cpue",values_to="value") %>% 
  mutate(chain=factor(chain))

g_a<- ggplot(cpue_long %>% filter(cpue=="Aerial")) + 
  geom_line(aes(x=samp,y=value,col=chain))+
  ylab("Aerial")
g_g<-ggplot(cpue_long %>% filter(cpue=="Ground")) + 
  geom_line(aes(x=samp,y=value,col=chain))+
  ylab("Ground")
g_s<-ggplot(cpue_long %>% filter(cpue=="Systematic Baiting")) + 
  geom_line(aes(x=samp,y=value,col=chain))+
  ylab("Systematic Baiting")
g_t<-ggplot(cpue_long %>% filter(cpue=="Trap")) + 
  geom_line(aes(x=samp,y=value,col=chain))+
 ylab("Trap")
g_cpue<- gridExtra::grid.arrange(g_a,g_g,g_t,g_s)

ggsave(g_cpue, 
       filename = "./Model outputs/Plots/Manuscript/Trace/cpue_trace.jpeg",
       width=7,height=5,units="in",device="jpeg")

##lambda -----------------
lambda <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
  cbind.data.frame(samples[[i]][,grep("lambda",colnames(samples[[i]]))],
                   chain=i)})),
  samp=rep(1:nmcmc,nChains))

colnames(lambda)[1] <- "Outside Elimination Areas"
colnames(lambda)[2:nea] <- paste("Elimination Area",1:(nea-1))
lambda_long <- lambda %>% 
  pivot_longer(cols=all_of(1:nea),names_to="elim_area",values_to="value") %>% 
  mutate(chain=factor(chain))

g_l <- ggplot(lambda_long) + geom_line(aes(x=samp,y=value,col=chain))+
  facet_wrap(.~elim_area)+
  geom_hline(yintercept=1,col="black",lty=2)
g_l

ggsave(g_l, 
       filename = "./Model outputs/Plots/Manuscript/Trace/lambda_trace.jpeg",
       width=7,height=5,units="in",device="jpeg")

rm(g_l,g_b,g_cpue,g_a,g_t,g_g,g_n1,g_r,g_p,g_s)

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

mean(rowMeans(pabs))
quantile(rowMeans(pabs),c(0.025,0.975))

pabs_sum[which.max(pabs_sum$mn),]
pabs_sum[which.min(pabs_sum$mn),]

#elimination area medians (of)
pabs_sum <- pabs_sum %>% 
  left_join(study_site_grid %>% 
              st_drop_geometry() %>% 
              rename(site_idx=SiteID) %>% 
              select(site_idx,Area_Name))
pabs_ea_sum <- pabs_sum %>% 
  group_by(Area_Name,per_start) %>% 
  summarise(ea_md=median(mn))
pabs_ea_sum[which.max(pabs_ea_sum$ea_md),]
pabs_ea_sum[which.min(pabs_ea_sum$ea_md),]

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
  geom_sf(data=pabs_sum_fy,aes(fill=fy_mn))+
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

ggsave(filename=paste0("./Model Outputs/Plots/Manuscript/pabs_all_yrs.jpeg"),
       device="jpeg",width=8,height=10,units="in")

##by elimination area and fy ------------------------------
pabs_sum_fy <- pabs_sum_fy %>% st_intersection(elim_areas %>% select(Area_Name))

elim_thresh <- 0.5
  
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

# ggsave(filename=paste0("./Model Outputs/Plots/Manuscript/pelim_elim_areas_fy_thresh.jpeg"),
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
pabs_thresh$Area_Name[pabs_thresh$Area_Name=="0"] <- "Outside EAs"

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
ggplot(pabs_thresh)+
  geom_ribbon(aes(x=per_start,ymin=lci,ymax=uci),alpha=0.3)+
  geom_line(aes(x=per_start,y=mn,col=Area_Name),lwd=1.5)+
  # geom_bar(data=dat_rem_sum,aes(x=month,y=removal/axis_scale),
  #          stat="identity",alpha=0.7)+
  facet_wrap(.~Area_Name)+
  scale_y_continuous(name=paste0("Prop. of watersheds with p(elim)>",elim_thresh)#,
                     # sec.axis = sec_axis(~.*axis_scale,name="No. feral swine removed",
                     #                     breaks=c(0,2500,5000))
                     )+
  xlab("Season")+
  scale_color_discrete(name="Elimination area")+
  # scale_fill_discrete(name="Removal type",labels=c("Aerial","Ground","Trap"))+
  theme(text=element_text(size=15))+
  guides(color="none")

ggsave(filename=paste0("./Model Outputs/Plots/Manuscript/pabs_over50_season_ea.jpeg"),
       device="jpeg",width=10,height=5,units="in")

## extent of p(abs) threshold ----------------------
elim_thresh<-0.5
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

ggplot(pabs_thresh_ext)+ geom_ribbon(aes(x=per_start,ymin=lci,ymax=uci),alpha=0.5)+
  geom_line(aes(x=per_start,y=mn),lwd=1)+
  xlab("Season")+
  ylab(expression(paste("K",m^2," with > 0.75 p(elimination)")))+
  ylim(0,tot_area_km/4)

ggsave(filename="./Model outputs/Plots/km2_eliminated_75.jpeg",
       width=7,height=5,units="in",device="jpeg")

ggplot(pabs_thresh_ext)+ geom_ribbon(aes(x=per_start,ymin=lci_prop,ymax=uci_prop),alpha=0.5)+
  geom_line(aes(x=per_start,y=mn_prop),lwd=1)+
  xlab("Season")+
  ylab("Proportion of study area with p(elimination) > 0.75")+
  ylim(0,0.25)

ggsave(filename="./Model outputs/Plots/prop_eliminated_75.jpeg",
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
write.csv(pabs_thresh_yr,paste0("./Model outputs/area_above_,",elim_thresh,"_table.csv"))

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
  geom_line(aes(x=per_start,y=mn,col=Area_Name),lwd=1.5)+
  facet_wrap(.~Area_Name)+
  xlab("Season")+
  scale_y_continuous(name="Median elimination probability")+
  scale_color_discrete(name="Elimination area")+
  theme(text=element_text(size=15))+
  guides(color="none")

ggsave(filename=paste0("./Model Outputs/Plots/Manuscript/pabs_md_season_ea.jpeg"),
       device="jpeg",width=10,height=6,units="in")

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
if(!exists("cpue")){
  cpue <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
    cbind.data.frame(samples[[i]][,grepl("p_",colnames(samples[[i]])) &
                                    !grepl("p_n1",colnames(samples[[i]]))],
                     chain=i)})),
    samp=rep(1:nmcmc,nChains))
  
  cpue_names <- c("Aerial","Ground","Systematic Baiting","Trap")
  colnames(cpue)[1:4] <- cpue_names
}
p_sys <- cpue %>% select(`Systematic Baiting`,chain,samp)

sys_eff <- seq(min(dat_occ$trap_nights_km),max(dat_occ$trap_nights_km),by=0.5)

det <- sapply(1:(nmcmc*nChains),function(i)1-(1-p_sys[i,1])^sys_eff)

det_sum <- data.frame(traps_km=sys_eff,
                      mn=rowMeans(det),
                      md=sapply(1:nrow(det),function(i)quantile(det[i,],prob=0.5)),
                      lci=sapply(1:nrow(det),function(i)quantile(det[i,],prob=0.025)),
                      uci=sapply(1:nrow(det),function(i)quantile(det[i,],prob=0.975)))

ggplot(det_sum)+
  geom_ribbon(aes(x=traps_km,ymin=lci,ymax=uci),alpha=0.3)+
  geom_line(aes(x=traps_km,y=mn))+
  ylab("Detection probability")+
  xlab(expression(paste("Bait stations per k",m^2)))+
  theme(text=element_text(size=15))

ggsave(filename=paste0("./Model Outputs/Plots/Manuscript/occ_det_eff_curve.jpeg"),
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

pelim_sum[which.max(pelim_sum$mn),]
pelim_sum[which.min(pelim_sum$mn),]

pelim_sum <- pelim_sum %>% 
  left_join(study_site_grid %>% 
              st_drop_geometry() %>% 
              rename(site_idx=SiteID) %>% 
              select(site_idx,Area_Name))
pabs_ea_sum <- pelim_sum %>% 
  group_by(Area_Name,per_start) %>% 
  summarise(ea_md=median(mn))

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
#   ggsave(g[[i]],filename=paste0("./Model Outputs/Plots/Manuscript/pelim_fy",fys[i],".jpeg"),
#          device="jpeg",width=10,height=5,units="in")
# }
# g[[i]]

### facet version ------------------
ggplot()+
  geom_sf(data=pelim_sum_fy,aes(fill=fy_mn))+
  geom_sf(data=elim_areas %>% filter(Area_Name!="0"),aes(col=Area_Name),
          fill="transparent",lwd=1.5)+
  facet_wrap(.~fy,ncol=1)+
  scale_fill_viridis_c(#limits=c(0,1),
    direction=-1,
    name="Feral swine \nelimination probability \nafter 10 negative \ntrap nights")+
  scale_color_discrete(name="Elimination Area")+
  theme(panel.background=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        text=element_text(size=15))

ggsave(filename=paste0("./Model Outputs/Plots/Manuscript/pelim_all_yrs.jpeg"),
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

#number of samples maps --------------------
eff_elim <- do.call("rbind",lapply(1:nChains,function(i){samples[[i]][grep("eff_elim",colnames(samples[[i]]))]}))

eff_sum <- data.frame(idx=colnames(eff_elim),
           mn=colMeans(eff_elim,na.rm=T),
           md=sapply(1:ncol(eff_elim),function(i)quantile(eff_elim[,i],prob=0.5,na.rm=T)),
           lci=sapply(1:ncol(eff_elim),function(i)quantile(eff_elim[,i],prob=0.025,na.rm=T)),
           uci=sapply(1:ncol(eff_elim),function(i)quantile(eff_elim[,i],prob=0.975,na.rm=T)))

eff_sum$site_idx <- sapply(1:nrow(eff_sum),function(x){
  as.numeric(gsub("\\[","",gsub("eff_elim","",unlist(str_split((eff_sum[x,"idx"]),","))[1])))})

eff_sum$period_idx <- sapply(1:nrow(eff_sum),function(x){
  as.numeric(gsub("\\]","",unlist(str_split((eff_sum[x,"idx"]),","))[2]))})
  
eff_sum_sf <- eff_sum %>% 
  left_join(study_site_grid %>% rename(site_idx=SiteID)) %>% 
  st_as_sf() 

i<-nperiods
ggplot()+
  geom_sf(data=eff_sum_sf %>% filter(period_idx==i),
                 aes(fill=mn))+
  geom_sf(data=elim_areas %>% filter(Area_Name!="0"),aes(col=Area_Name),
          fill="transparent",lwd=1.5)+
  scale_fill_viridis_c(name="No. of negative \nbait nights required to \ndetermine 95% \nprobability of elimination",
                       option="B")+
  scale_color_discrete(name="Elimination Area")+
  theme_void()+
  theme(plot.margin = margin(1,1,1,1, "cm"),
        text=element_text(size=15))

ggsave(filename="./Model Outputs/Plots/Manuscript/no_samp_ws_summer2023.jpeg",
       device="jpeg",width=10,height=5,units="in")

eff_sum[which.min(eff_sum$mn),]
eff_sum[which.max(eff_sum$mn),]

#abundance  -----------------------------
N <- do.call("rbind",lapply(1:nChains,function(i){
  samples[[i]][grepl("N",colnames(samples[[i]])) & !grepl("N_latent",colnames(samples[[i]]))]}))

N_sum=data.frame(idx=colnames(N),
           mn=colMeans(N),
           md=sapply(1:ncol(N),function(i)quantile(N[,i],prob=0.5,na.rm=T)),
           lci=sapply(1:ncol(N),function(i)quantile(N[,i],prob=0.025,na.rm=T)),
           uci=sapply(1:ncol(N),function(i)quantile(N[,i],prob=0.975,na.rm=T)))

# N_sum <- data.frame(idx=names(N_mn),mn=N_mn,md=N_md,lci=N_lci,uci=N_uci)
N_sum$elim_area_idx <- sapply(1:nrow(N_sum),function(x){
  as.numeric(gsub("\\[","",gsub("N","",unlist(str_split(N_sum[x,"idx"],","))[1])))})

N_sum$period_idx <- sapply(1:nrow(N_sum),function(x){
  as.numeric(gsub("\\]","",unlist(str_split((N_sum[x,"idx"]),","))[2]))})

N_sum$per_start <- sapply(1:nrow(N_sum),function(i){
  min(dat_occ$per_start[dat_occ$period_idx==N_sum$period_idx[i]])})
N_sum$per_start <- as.Date(N_sum$per_start,origin='1970-01-01')
N_sum$year <- year(N_sum$per_start)
N_sum$fy <- N_sum$year
N_sum$fy[month(N_sum$per_start)%in%c(10:12)] <- 
  N_sum$fy[month(N_sum$per_start)%in%c(10:12)] +1

#standardized abundance
all_N <- unlist(c(N_sum[,c("mn","md","lci","uci")]))
N_ss <- sqrt(sum(all_N^2)/(length(all_N)-1))

N_sum$mn_std <- N_sum$mn/N_ss
N_sum$md_std <- N_sum$md/N_ss
N_sum$lci_std <- N_sum$lci/N_ss
N_sum$uci_std <- N_sum$uci/N_ss

N_sum_sf <- N_sum %>% 
  left_join(elim_areas) %>% 
  st_as_sf() %>% 
  left_join(dat_rem %>% ungroup() %>% 
              group_by(period_idx) %>% 
              dplyr::select(period_idx,Date) %>% 
              summarise(per_start=min(Date)))

N_sum_sf <- N_sum_sf %>% 
  left_join(elim_areas %>%
              summarise(ea_area_km=as.numeric(st_area(geometry)/1e6)) %>% 
              st_drop_geometry())
N_sum_sf <- N_sum_sf %>% mutate(mn_dens=mn/ea_area_km,
                                md_dens=md/ea_area_km,
                                lci_dens=lci/ea_area_km,
                                uci_dens=uci/ea_area_km)

N_sum_sf$Area_Name_label <- paste("Elimination Area ",N_sum_sf$Area_Name)
dat_rem_sum$Area_Name_label <- paste("Elimination Area ",dat_rem_sum$Area_Name)

N_sum_4 <- N_sum_sf %>% st_drop_geometry() %>% filter(Area_Name==4)
N_sum_4[which.min(N_sum_4$mn_dens),]
N_sum_4[which.max(N_sum_4$mn_dens),]

N_sum_6 <- N_sum_sf %>% st_drop_geometry() %>% filter(Area_Name==6)
N_sum_6[which.min(N_sum_6$mn_dens),]
N_sum_6[which.max(N_sum_6$mn_dens),]

####raw abundance EA 4 and 6 -------------------------
axis_trans<- 0.015
ggplot()+  
  geom_ribbon(data=N_sum_sf %>% filter(Area_Name%in%c(4,6) & period_idx!=1),
              aes(x=per_start,ymin=lci,ymax=uci),alpha=0.2)+
  geom_line(data=N_sum_sf %>% filter(Area_Name%in%c(4,6)& period_idx!=1),
            aes(x=per_start,y=mn,col=Area_Name),lwd=0.75)+
  geom_bar(data=dat_rem_sum %>% filter(Area_Name%in%c(4,6)),
           aes(x=month,y=removal/axis_trans,fill=method),stat="identity",
           alpha=0.7)+
  facet_wrap(.~Area_Name_label)+
  scale_y_continuous(#limits=c(0,0.15),
    sec.axis = sec_axis(transform=~ . * axis_trans,
                        name = "No. feral swine removed"
                        # breaks=c(0,1000,2000)
                        ))+
  scale_fill_discrete(name="Removal method")+
  ylab("Feral swine abundance")+
  xlab("Season")+
  scale_color_manual(name="Elimination area",
                     values=scales::hue_pal()(6)[c(4,6)])+
  guides(color="none")+
  theme(text=element_text(size=15))
ggsave(filename="./Model Outputs/Plots/Manuscript/abundance_trend_removal.jpeg",
       device="jpeg",width=10,height=6,units="in")

#### density EA 4 and 6 -------------------------
axis_trans_d<- 200
ggplot()+  
  geom_ribbon(data=N_sum_sf %>% filter(Area_Name%in%c(4,6)),
              aes(x=per_start,ymin=lci_dens,ymax=uci_dens),alpha=0.2)+
  geom_line(data=N_sum_sf %>% filter(Area_Name%in%c(4,6)),
            aes(x=per_start,y=mn_dens,col=Area_Name),lwd=0.75)+
  geom_bar(data=dat_rem_sum %>% filter(Area_Name%in%c(4,6)),
           aes(x=month,y=removal/axis_trans_d,fill=method),stat="identity",
           alpha=0.7)+
  facet_wrap(.~Area_Name_label)+
  scale_y_continuous(#limits=c(0,0.15),
    sec.axis = sec_axis(transform=~ . * axis_trans_d,
                        name = "No. feral swine removed"
                        # breaks=c(0,1000,2000)
    ))+
  scale_fill_discrete(name="Removal method")+
  ylab(expression(paste("Feral swine density (swine/k",m^2,")")))+
  xlab("Season")+
  scale_color_manual(name="Elimination area",
                     values=scales::hue_pal()(6)[c(4,6)])+
  guides(color="none")+
  theme(text=element_text(size=15))
ggsave(filename="./Model Outputs/Plots/Manuscript/density_trend_removal.jpeg",
       device="jpeg",width=10,height=6,units="in")

####standardized abundance EA 4 and 6 -------------------------
axis_scale <- 700
ggplot()+  
  geom_ribbon(data=N_sum_sf %>% filter(Area_Name%in%c(4,6)& period_idx!=1),
              aes(x=per_start,ymin=lci_std,ymax=uci_std),alpha=0.2)+
  geom_line(data=N_sum_sf %>% filter(Area_Name%in%c(4,6)& period_idx!=1),
            aes(x=per_start,y=mn_std,col=Area_Name),lwd=1)+
  geom_bar(data=dat_rem_sum %>% filter(Area_Name%in%c(4,6)),
           aes(x=month,y=removal/axis_scale,fill=method),stat="identity",
           alpha=0.7)+
  facet_wrap(.~Area_Name)+
  scale_y_continuous(#limits=c(0,0.15),
                     sec.axis = sec_axis(transform=~.*axis_scale, 
                                         name = "No. feral swine removed"
                                         # breaks=c(0,1000,2000)
                                         ))+
  scale_fill_discrete(name="Removal method")+
  ylab("Standardized feral swine abundance")+
  xlab("Season")+
  scale_color_manual(name="Elimination area",
                     values=scales::hue_pal()(6)[c(4,6)])+
  guides(color="none")+
  theme(text=element_text(size=15))

ggsave(filename="./Model Outputs/Plots/Manuscript/std_abundance_trend_removal.jpeg",
       device="jpeg",width=10,height=6,units="in")

## abundance trends ------------------------------
# yrs <- unique(N_sum_sf$fy)
# gN <- N_trend <- list()
# for(i in 1:length(yrs)){
#   if(i==1){
#     N_trend[[i]] <- N_sum_sf %>% 
#       filter(fy==yrs[i]) %>% 
#       filter(period_idx!=1) %>% 
#       group_by(elim_area_idx) %>% 
#       summarise(trend = lm(md_std~period_idx)$coefficients[2],
#                 fy=min(fy))
#   } else {
#     N_trend[[i]] <- N_sum_sf %>% 
#       filter(fy==yrs[i]) %>% 
#       group_by(elim_area_idx) %>% 
#       summarise(trend = lm(md_std~period_idx)$coefficients[2],
#                 fy=min(fy))
#   }
# }
# 
# trend_max <- max(unlist(lapply(1:length(yrs),function(i)max(N_trend[[i]]$trend))))
# trend_min <- min(unlist(lapply(1:length(yrs),function(i)min(N_trend[[i]]$trend))))

# for(i in 1:length(yrs)){
  # yr_idx <- N_sum_sf %>% st_drop_geometry() %>% 
  #   filter(fy==yrs[i]) %>% summarise(period=min(period_idx)) %>% 
  #   as.numeric()
  
  # gN[[i]] <- ggplot()+
  #   geom_sf(data=elim_areas,fill="grey88")+
  #   geom_sf(data=N_trend[[i]]%>% filter(elim_area_idx!=1),aes(fill=trend))+
  #   geom_sf(data=elim_areas %>% filter(Area_Name!="0"),aes(col=Area_Name),
  #           fill="transparent",lwd=1) +
  #   scale_color_discrete(name="Elimination Area")+
  #   scale_fill_gradient2(name=paste("Stanardized \nabundance trend"),
  #                        high="red",mid="white",low="blue",na.value="grey50",
  #                        limits=c(-0.35,0.1)
  #                        )+
  #   theme(panel.background = element_blank(),
  #         text=element_text(size=15),
  #         axis.text=element_blank(),
  #         axis.ticks=element_blank())+
  #   ggtitle(paste("FY",yrs[i]))
  # ggsave(gN[[i]],
  #        filename = paste0("./Model outputs/Plots/Maps/abund_trend_fy",yrs[i],".jpeg"),
  #        device="jpeg",width=6,height=5,units="in")
# }
# gN[[1]]
# gN[[2]]
# gN[[3]]

# N_trend <- do.call("rbind",N_trend)

##abundance trends wi elim areas fy -----------------
# N_trend <- N_trend %>% left_join(elim_areas %>% st_drop_geometry())
# N_trend$Area_Name[N_trend$Area_Name==0] <- "Outside EAs"
# N_trend_elim_fy <- N_trend %>% st_drop_geometry() %>% 
#   group_by(fy,Area_Name) %>% 
#   summarise(mn=mean(trend,na.rm=T),
#             md=median(trend,na.rm=T))

# ggplot(N_trend_elim_fy)+
#   geom_hline(yintercept = 0,col="grey33",lty=2)+
#   geom_line(aes(x=factor(fy),y=mn,col=factor(Area_Name),group=Area_Name),lwd=1)+
#   scale_color_discrete(name="Elimination Area")+
#   xlab("FY")+ylab("Mean trend in feral swine abundance")+
#   theme(text=element_text(size=15))
# 
# ggplot(N_trend_elim_fy)+
#   geom_hline(yintercept = 0,col="grey33",lty=2)+
#   geom_line(aes(x=factor(fy),y=md,col=factor(Area_Name),group=Area_Name),lwd=1)+
#   scale_color_discrete(name="Elimination Area")+
#   xlab("FY")+ylab("Median trend in feral swine abundance")+
#   theme(text=element_text(size=15))

##overall trends by fy -------------------
# N_trend %>% 
#   st_drop_geometry() %>% 
#   group_by(fy) %>% 
#   summarise(md=median(trend,na.rm=T))
# 
# N_sum %>% 
#   group_by(period_idx,fy) %>% 
#   mutate(grndmn=mean(md_std)) %>% 
#   group_by(fy) %>% 
#   # ungroup() %>% 
#   summarise(trend = lm(grndmn~period_idx)$coefficients[2])

## proportional change ---------------------------
### season ------------------------
# N_change <- do.call("rbind",lapply(1:nChains,function(i){
#   samples[[i]][grepl("N_change",colnames(samples[[i]]))]}))
# 
# N_change_sum=data.frame(idx=colnames(N_change),
#                         mn=colMeans(N_change),
#                         md=sapply(1:ncol(N_change),function(i)quantile(N_change[,i],prob=0.5,na.rm=T)),
#                         lci=sapply(1:ncol(N_change),function(i)quantile(N_change[,i],prob=0.025,na.rm=T)),
#                         uci=sapply(1:ncol(N_change),function(i)quantile(N_change[,i],prob=0.975,na.rm=T)))
# 
# N_change_sum$elim_area_idx <- sapply(1:nrow(N_change_sum),function(x){
#   as.numeric(gsub("\\[","",gsub("N_change","",unlist(str_split(N_change_sum[x,"idx"],","))[1])))})
# 
# N_change_sum$period_idx <- sapply(1:nrow(N_change_sum),function(x){
#   as.numeric(gsub("\\]","",unlist(str_split((N_change_sum[x,"idx"]),","))[2]))})
# 
# N_change_sum$per_start <- sapply(1:nrow(N_change_sum),function(i){
#   min(dat_occ$per_start[dat_occ$period_idx==N_change_sum$period_idx[i]])})
# N_change_sum$per_start <- as.Date(N_change_sum$per_start,origin='1970-01-01')
# N_change_sum$year <- year(N_change_sum$per_start)
# N_change_sum$fy <- N_change_sum$year
# N_change_sum$fy[month(N_change_sum$per_start)%in%c(10:12)] <- 
#   N_change_sum$fy[month(N_change_sum$per_start)%in%c(10:12)]+1
# 
# ggplot(N_change_sum)+geom_line(aes(x=per_start,y=mn,col=factor(elim_area_idx)))

### annual ------------------------
N_long <- N %>% pivot_longer(cols=all_of(1:ncol(N)),names_to="idx",values_to="value")
N_long$elim_area_idx <- rep(rep(1:nea,nperiods),nmcmc*nChains)
N_long$period_idx <- rep(sort(rep(1:nperiods,nea)),nmcmc*nChains)
N_long$samp_idx <- sort(rep(1:(nmcmc*nChains),nea*nperiods))
N_long <- N_long %>%
  left_join(elim_areas %>% st_drop_geometry())

per_idx <- dat_occ %>% select(period_idx,per_start) %>% 
  mutate(year=year(per_start)) %>% distinct()
per_idx$fy <- per_idx$year
per_idx$fy[month(per_idx$per_start)%in%c(10,11,12)] <- per_idx$fy[month(per_idx$per_start)%in%c(10,11,12)] +1

N_long <- N_long %>% left_join(per_idx)

N_yr <- N_long %>% filter(month(per_start)==10 | (month(per_start)==7 & year==2023))
N_yr$per_temp <- as.numeric(as.factor(N_yr$period_idx))

N_yr$value_prev<- NA
for(i in 1:nrow(N_yr)){
  if(N_yr$per_temp[i]>1){
    N_yr$value_prev[i] <- N_yr$value[N_yr$per_temp==(N_yr$per_temp[i]-1) &
                                       N_yr$elim_area_idx==N_yr$elim_area_idx[i] &
                                       N_yr$samp_idx==N_yr$samp_idx[i]]
  }
}

N_yr$per_change <- (N_yr$value-N_yr$value_prev)/N_yr$value_prev

N_yr_sum <- N_yr %>% 
  filter(year>2020) %>% 
  group_by(Area_Name,period_idx,per_start) %>% 
  summarise(mn=mean(per_change,na.rm=T)*100,
            md=median(per_change,na.rm=T)*100,
            lci=quantile(per_change,0.025,na.rm=T)*100,
            uci=quantile(per_change,0.975,na.rm=T)*100) 

write.csv(N_yr_sum,"./Model outputs/N_change_table.csv")


# removal probability -----------------
if(!exists("cpue")){
  cpue <- cbind.data.frame(do.call("rbind",lapply(1:nChains,function(i){
    cbind.data.frame(samples[[i]][,grepl("p_",colnames(samples[[i]])) &
                                    !grepl("p_n1",colnames(samples[[i]]))],
                     chain=i)})),
    samp=rep(1:nmcmc,nChains))
  
  cpue_names <- c("Aerial","Ground","Systematic Baiting","Trap")
  colnames(cpue)[1:4] <- cpue_names
}

p_rem <- cpue %>% select(Aerial,Ground,Trap,chain,samp)
p_rem_sum <- p_rem %>%
  pivot_longer(cols=all_of(1:3),names_to="rem_typ",values_to="value") %>%
  group_by(rem_typ) %>%
  summarise(mn=mean(value),
            md=median(value),
            lci=quantile(value,0.025),
            uci=quantile(value,0.975))
p_rem_sum

## aerial --------------------------
aerial_eff <- seq(min(dat_aerial$eff_area_hrs),max(dat_aerial$eff_area_hrs),by=0.01)
aerial_area_mn <- mean(dat_aerial$prop_ea_impact[dat_aerial$prop_ea_impact!=0])
det_aerial <- sapply(1:(nmcmc*nChains),function(i){
  aerial_area_mn*(1-(1-p_rem[i,"Aerial"])^aerial_eff)})


det_aerial_sum <- data.frame(flight_hrs_km=aerial_eff,
                      mn=rowMeans(det_aerial),
                      md=sapply(1:nrow(det_aerial),function(i)quantile(det_aerial[i,],prob=0.5)),
                      lci=sapply(1:nrow(det_aerial),function(i)quantile(det_aerial[i,],prob=0.025)),
                      uci=sapply(1:nrow(det_aerial),function(i)quantile(det_aerial[i,],prob=0.975)))
g_a <- ggplot(det_aerial_sum)+
  geom_ribbon(aes(x=flight_hrs_km,ymin=lci,ymax=uci),alpha=0.3)+
  geom_line(aes(x=flight_hrs_km,y=mn))+
  ylab("Removal rate")+
  xlab(expression(paste("Flight hours per k",m^2)))+
  ggtitle("Aerial Operations")+
  theme(text=element_text(size=15))

# ggsave(g_a,filename=paste0("./Model Outputs/Plots/Manuscript/rem_det_curve_aerial.jpeg"),
#        device="jpeg",width=7,height=5,units="in")

## trapping --------------------------
trap_eff <- seq(min(dat_trap$eff_area_hrs),max(dat_trap$eff_area_hrs),by=0.01)
trap_area_mn <- mean(dat_trap$prop_ea_impact[dat_trap$prop_ea_impact!=0])
det_trap <- sapply(1:(nmcmc*nChains),function(i){
  trap_area_mn*(1-(1-p_rem[i,"Trap"])^trap_eff)})

det_trap_sum <- data.frame(trap_hrs_km=trap_eff,
                             mn=rowMeans(det_trap),
                             md=sapply(1:nrow(det_trap),function(i)quantile(det_trap[i,],prob=0.5)),
                             lci=sapply(1:nrow(det_trap),function(i)quantile(det_trap[i,],prob=0.025)),
                             uci=sapply(1:nrow(det_trap),function(i)quantile(det_trap[i,],prob=0.975)))
g_t <- ggplot(det_trap_sum)+
  geom_ribbon(aes(x=trap_hrs_km,ymin=lci,ymax=uci),alpha=0.3)+
  geom_line(aes(x=trap_hrs_km,y=mn))+
  xlim(0,2.55)+
  ylab("Removal rate")+
  xlab(expression(paste("Trapping hours per k",m^2)))+
  ggtitle("Trapping")+
  theme(text=element_text(size=15))

# ggsave(g_t,filename=paste0("./Model Outputs/Plots/Manuscript/rem_det_curve_trap.jpeg"),
       # device="jpeg",width=7,height=5,units="in")

## ground -------------------------
ground_eff <- seq(min(dat_ground$eff_area_events),max(dat_ground$eff_area_events),by=0.01)
ground_area_mn <- mean(dat_ground$prop_ea_impact[dat_ground$prop_ea_impact!=0])
det_ground <- sapply(1:(nmcmc*nChains),function(i){
  ground_area_mn*(1-(1-p_rem[i,"Ground"])^ground_eff)})

det_ground_sum <- data.frame(events_km=ground_eff,
                             mn=rowMeans(det_ground),
                             md=sapply(1:nrow(det_ground),function(i)quantile(det_ground[i,],prob=0.5)),
                             lci=sapply(1:nrow(det_ground),function(i)quantile(det_ground[i,],prob=0.025)),
                             uci=sapply(1:nrow(det_ground),function(i)quantile(det_ground[i,],prob=0.975)))
g_g<- ggplot(det_ground_sum)+
  geom_ribbon(aes(x=events_km,ymin=lci,ymax=uci),alpha=0.3)+
  geom_line(aes(x=events_km,y=mn))+
  ylab("Removal rate")+
  xlab(expression(paste("No. of ground shooting events per k",m^2)))+
  ggtitle("Ground Shooting")+
  theme(text=element_text(size=15))

# ggsave(g_g,filename=paste0("./Model Outputs/Plots/Manuscript/rem_det_curve_ground.jpeg"),
       # device="jpeg",width=7,height=5,units="in")

## all together --------------------
g_all <- gridExtra::arrangeGrob(g_a+ylim(0,0.035), g_t+ylim(0,0.035), g_g+ylim(0,0.035),nrow=1) #generates g

ggsave(g_all,file=paste0("./Model Outputs/Plots/Manuscript/rem_det_curve_all.jpeg"),
       device="jpeg",width=13,height=6,units="in")

##removal comparisons --------------------
a_mn <- max(dat_aerial$prop_ea_impact)*(1-(1-p_rem[,"Aerial"])^max(dat_aerial$eff_area_hrs[dat_aerial$eff_area_hrs!=0]))
t_mn <- max(dat_trap$prop_ea_impact)*(1-(1-p_rem[,"Trap"])^max(dat_trap$eff_area_hrs[dat_trap$eff_area_hrs!=0]))
g_mn <- max(dat_ground$prop_ea_impact)*(1-(1-p_rem[,"Ground"])^max(dat_ground$eff_area_events[dat_ground$eff_area_events!=0]))

rem_mn <- cbind.data.frame(aerial=a_mn,trap=t_mn,ground=g_mn)
rem_df <- rem_mn %>% pivot_longer(cols=1:3,names_to="rem_typ",values_to="value") %>% 
  group_by(rem_typ) %>% 
  summarise(mn=mean(value),
            lci=quantile(value,0.025),
            uci=quantile(value,0.975))
rem_df$mn[1]/rem_df$mn[2]
rem_df$mn[1]/rem_df$mn[3]
range(elim_areas %>% filter(Area_Name!="0") %>% st_area()/1e6)

#population growth rate--------------
if(!exists("lambda")){
  lambda <- do.call("rbind",lapply(1:nChains,function(i){
  samples[[i]][grepl("lambda",colnames(samples[[i]]))]}))
}

data.frame(mn=colMeans(lambda[,1:nea]),
           lci=sapply(1:nea,function(i)quantile(lambda[,i],probs=c(0.025))),
           uci=sapply(1:nea,function(i)quantile(lambda[,i],probs=c(0.975))))

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

dat_rem_sum %>% 
  filter(Area_Name!="0") %>% 
  ggplot()+geom_bar(aes(y=removal,x=month,fill=method),
                    stat="identity")+
  facet_wrap(.~Area_Name_label)+
  scale_fill_discrete(name="Removal method")+
  xlab("Season")+ylab("No. feral swine removed")+
  theme(text=element_text(size=15))
ggsave(filename = "./Model outputs/Plots/Manuscript/raw_removal_ea.jpeg",
       device="jpeg",height=5,width=7,units="in")

#effort ---------------------
dat_rem_sum %>% 
  filter(Area_Name!="0") %>% 
  ggplot()+geom_bar(aes(y=tot_hrs,x=month,fill=method),
                    stat="identity")+
  facet_wrap(.~Area_Name_label)+
  scale_fill_discrete(name="Removal method")+
  xlab("Season")+ylab("No. hours spent removing feral swine")+
  theme(text=element_text(size=15))
ggsave(filename = "./Model outputs/Plots/raw_effort_hours.jpeg",
       device="jpeg",height=5,width=7,units="in")

#raw cpue -------------------
dat_rem_sum %>% 
  filter(Area_Name!="0") %>% 
  ggplot()+geom_bar(aes(y=removal/tot_hrs,x=month,fill=method),
                    stat="identity",position="dodge")+
  facet_wrap(.~Area_Name_label)+
  scale_fill_discrete(name="Removal method")+
  xlab("Season")+ylab("No. of pigs removed per hour of effort")+
  theme(text=element_text(size=15))

ggsave(filename = "./Model outputs/Plots/raw_cpue.jpeg",
       device="jpeg",height=5,width=7,units="in")

# summary stats ------------------
##raw data ------------------------
dat_rem_sum %>% 
  group_by(method) %>% 
  summarise(tot_rem=sum(removal),
            tot_hrs=sum(tot_hrs),
            num_events=sum(tot_events))

dat_rem_sum %>% 
  group_by(Area_Name) %>% 
  summarise(tot_rem=sum(removal),
            tot_hrs=sum(tot_hrs),
            num_events=sum(tot_events))

sysbait_det_eff %>% 
  ungroup() %>% 
  summarise(dets=sum(detection),
            trap_nights=sum(trap_nights))

#save posterior summaries ---------------------------------
save(eff_sum_sf,N_sum_sf,N_yr_sum,p_rem_sum,pabs_sum_ea,pabs_sum,
     pabs_sum_fy,pabs_thresh,pabs_thresh,pabs_thresh_ext,
     pelim_ea,pelim_ea_fy,pelim_sum,pelim_sum_fy,rem_df,
     det_aerial_sum,det_trap_sum,det_ground_sum,
     pabs_thresh_yr,det_sum,nea,nsites,nperiods,nmcmc,nChains,
     dat_occ,dat_rem,dat_rem_sum,
     file="./Model outputs/posterior_summaries_25JUL24.RData")


          