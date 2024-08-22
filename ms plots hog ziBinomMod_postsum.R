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
load("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/posterior_summaries_25JUL24.Rdata")

#absence probability --------------------
## map to study site -------------------
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
elim_thresh <- 0.95

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

## elimination threshold extent table ----------------
write.csv(pabs_thresh_yr,paste0("./Model outputs/area_above_,",elim_thresh,"_table.csv"))

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
## map to study site -------------------
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
## density EA 4 and 6 -------------------------
axis_trans_d<- 300
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
axis_scale <- 600
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

# removal probability -----------------
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

## ground -------------------------
# g_g<- ggplot(det_ground_sum)+
#   geom_ribbon(aes(x=events_km,ymin=lci,ymax=uci),alpha=0.3)+
#   geom_line(aes(x=events_km,y=mn))+
#   ylab("Removal rate")+
#   xlab(expression(paste("No. of ground shooting events per k",m^2)))+
#   ggtitle("Ground Shooting")+
#   theme(text=element_text(size=15))

# ggsave(g_g,filename=paste0("./Model Outputs/Plots/Manuscript/rem_det_curve_ground.jpeg"),
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


