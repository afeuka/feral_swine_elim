### Title: Plotting hurdle model results - posterior summaries only
### Author: Abbey Feuka
### Date: 07AUG24
### Notes: uses posterior_summaries_25JUL24.RData

library(tidyverse)
library(sf)

setwd("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble")

##study site--------------------
study_site_grid <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/huc10_siteIDs_cond_EA.shp")
study_site_grid <- study_site_grid %>% rename(elim_area_idx=elm_r_d,
                                              Area_Name=Area_Nm)
## elimination areas -----------------
elim_areas <- study_site_grid %>%
  group_by(Area_Name,elim_area_idx) %>% 
  summarise(geometry=st_union(geometry))

# load posterior summaries ---------------------
load("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/Plots/NFSP/posterior_summaries_04SEP24.Rdata")


#absence probability --------------------
## text -----------------
pabs_sum[which.max(pabs_sum$md),]
pabs_sum[which.min(pabs_sum$md),]

## map to study site -------------------
ggplot()+
  geom_sf(data=pabs_sum_fy,aes(fill=fy_mn))+
  geom_sf(data=elim_areas %>% filter(Area_Name!="0"),aes(col=Area_Name),
          fill="transparent",lwd=1.5)+
  facet_wrap(.~fy,ncol=1)+
  scale_fill_viridis_c(direction=-1,name="Feral swine \nelimination probability")+
  scale_color_discrete(name="Elimination Area")+
  theme(panel.background=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        plot.margin = margin(1,1,1,1, "cm"),
        text=element_text(size=15))

ggsave(filename=paste0("./Model Outputs/Plots/Manuscript/pabs_all_yrs.jpeg"),
       device="jpeg",width=8,height=10,units="in")

##by elimination area and fy ------------------------------
ggplot(pabs_thresh)+
  geom_ribbon(aes(x=per_start,ymin=lci,ymax=uci),alpha=0.3)+
  geom_line(aes(x=per_start,y=mn,col=Area_Name),lwd=1.5)+
  facet_wrap(.~Area_Name)+
  scale_y_continuous(name=paste0("Prop. of watersheds with p(elim)>",elim_thresh))+
  xlab("Season")+
  scale_color_discrete(name="Elimination area")+
  theme(text=element_text(size=15))+
  guides(color="none")

ggsave(filename=paste0("./Model Outputs/Plots/Manuscript/pabs_over",elim_thresh,"_season_ea.jpeg"),
       device="jpeg",width=10,height=5,units="in")

## extent of p(abs) threshold ----------------------
ggplot(pabs_thresh_ext)+ 
  geom_ribbon(aes(x=per_start,ymin=lci,ymax=uci),alpha=0.5)+
  geom_line(aes(x=per_start,y=mn),lwd=1)+
  xlab("Season")+
  ylab(expression(paste("K",m^2," with > 0.5 p(elimination)")))

ggsave(filename="./Model outputs/Plots/km2_eliminated_50.jpeg",
       width=7,height=5,units="in",device="jpeg")

ggplot(pabs_thresh_ext)+ 
  geom_ribbon(aes(x=per_start,ymin=lci_prop,ymax=uci_prop),alpha=0.5)+
  geom_line(aes(x=per_start,y=mn_prop),lwd=1)+
  xlab("Season")+
  ylab("Proportion of study area with p(elimination) > 0.95")

ggsave(filename="./Model outputs/Plots/prop_eliminated_",elim_thresh,".jpeg",
       width=7,height=5,units="in",device="jpeg")

## elimination threshold extent tables ----------------
pabs_thesh_ext_tab <- read.csv(paste0("./Model outputs/area_above_,",elim_thresh,"_table.csv"))
pabs_thesh_ext_tab

pabs_thresh_ext_ea_prop_tab <- read.csv(paste0("./Model outputs/Plots/",subfolder,"/area_above_,",elim_thresh,"_ea_prop_table.csv"))
pabs_thresh_ext_ea_prop_tab

pabs_thresh_ext_ea_km_tab <- read.csv(paste0("./Model outputs/Plots/",subfolder,"/area_above_,",elim_thresh,"_ea_km_table.csv"))
pabs_thresh_ext_ea_km_tab

##median p(abs) by season and ea -----------------------
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

## text ------------
pabs_sum_ea %>% filter(grepl("4",Area_Name) | 
                         grepl("6",Area_Name))
pabs_sum_ea[which.min(pabs_sum_ea$md),]
pabs_sum_ea[which.max(pabs_sum_ea$md),]
pabs_sum_ea %>% filter(Area_Name=="Elimination Area 1")

#calculate by fy all sites
pabs_sum_fy %>% 
  st_drop_geometry() %>% 
  group_by(fy) %>% 
  summarise(median(fy_md),
            mean(fy_mn))

#calculate by fy and elim area
pabs_sum_fy %>% st_drop_geometry() %>% 
  group_by(fy,Area_Name) %>% 
  summarise(mn=mean(fy_mn),
            md=median(fy_md))

pabs_sum %>% 
  summarise(mn=mean(mn))

#detection probability with trap effort -------------------
ggplot(det_sum)+
  geom_ribbon(aes(x=traps_km,ymin=lci,ymax=uci),alpha=0.3)+
  geom_line(aes(x=traps_km,y=mn))+
  ylab("Detection probability")+
  xlab(expression(paste("Bait stations per k",m^2)))+
  theme(text=element_text(size=15))

ggsave(filename=paste0("./Model Outputs/Plots/Manuscript/occ_det_eff_curve.jpeg"),
       device="jpeg",width=7,height=5,units="in")

#conditional absence probability maps --------------------
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

## text -------------
pelim_sum[which.max(pelim_sum$md),]
pelim_sum[which.min(pelim_sum$md),]

pabs_ea_sum[which.max(pabs_ea_sum$ea_md),]
pabs_ea_sum[which.min(pabs_ea_sum$ea_md),]

#calculate by fy all sites
pelim_sum_fy %>% 
  st_drop_geometry() %>% 
  group_by(fy) %>% 
  summarise(median(fy_md),
            mean(fy_mn))

#number of samples map --------------------
ggplot()+
  geom_sf(data=eff_sum_sf %>% filter(period_idx==nperiods),
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

#abundance  -----------------------------
## text --------------
N_yr_sum <- read.csv(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/Plots/",subfolder,"/N_change_table.csv"))
N_yr_sum

## density ----------
dat_rem <- rbind.data.frame(dat_aerial,dat_trap)
dat_rem_sum <- dat_rem %>% 
  group_by(SiteID,method,month=floor_date(Date,"month")) %>% 
  summarise(Season=min(Date),
            removal=sum(tot_rem),
            tot_hrs=sum(tot_hrs),
            tot_events=sum(num_events)) %>% 
  filter(month<=max(dat_occ$per_start))

dat_rem_sum <- dat_rem_sum %>% 
  left_join(study_site_grid %>% 
              st_drop_geometry() %>% 
              select(SiteID,area_km)) %>% 
  mutate(tot_rem_km = removal/area_km)


axis_trans_d<- 0.03
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
              aes(x=per_start,ymin=lci_dens_no_rem,ymax=uci_dens_no_rem,fill="No Removal"),
              alpha=0.1)+
  geom_line(data=N_sum,
            aes(x=per_start,y=md_dens_no_rem,col="No Removal"),lwd=1,lty=3)+
  facet_wrap(.~SiteID)+
  scale_y_continuous(sec.axis = sec_axis(transform=~.*axis_trans_d,
                                         name="No. feral swine removed per km2"))+
  scale_fill_manual(name="",
                    values=c("Estimated Density"="black",
                             "No Removal"="red"
                    ))+
  scale_color_manual(name="",
                     values=c("Estimated Density"="black",
                              "No Removal"="red"
                     ))+
  ylab(expression(paste("Feral swine density (swine/k",m^2,")")))+
  xlab("Season")+
  theme(text=element_text(size=15))
ggsave(filename=paste0("./Model Outputs/Plots/",subfolder,"/density_trend_removal.jpeg"),
       device="jpeg",width=10,height=6,units="in")

## removal needed -----------
axis_trans_d<- 0.03
ggplot()+  
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
ggsave(filename=paste0("./Model Outputs/Plots/",subfolder,"/density_trend_removal_needed.jpeg"),
       device="jpeg",width=10,height=6,units="in")

# removal probability -----------------
## text ------
rem_df 
rem_df$mn[1]/rem_df$mn[2]

## aerial --------------------------
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

## all together --------------------
g_all <- gridExtra::arrangeGrob(g_a+ylim(0,0.035), g_t+ylim(0,0.035), g_g+ylim(0,0.035),nrow=1) #generates g

ggsave(g_all,file=paste0("./Model Outputs/Plots/Manuscript/rem_det_curve_all.jpeg"),
       device="jpeg",width=13,height=6,units="in")

#pig removals ---------------------
##by EA--------------------
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

dat_rem %>% 
  group_by(method) %>% 
  summarise(max(prop_ea_impact))

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


