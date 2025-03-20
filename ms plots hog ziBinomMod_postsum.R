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
# load("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/Plots/NFSP/posterior_summaries_04SEP24.Rdata")
load("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/Plots/Single Lambda/posterior_summaries_17JAN25.Rdata")

subfolder <- "Single Lambda"

#absence probability --------------------
## text -----------------
pabs_sum[which.min(pabs_sum$md),]
pabs_sum[which.max(pabs_sum$md),]

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
pabs_sum_ea %>% filter(Area_Name=="Elimination Area 1") %>% 
  summarise(min(md),
            max(md))

pabs_sum_ea %>% filter(Area_Name=="Elimination Area 4") %>% 
  summarise(min(md),
            max(md))

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

## text ----
det_sum[which.max(det_sum$md),] %>% select(traps_km,md,lci,uci)

#conditional absence probability maps --------------------
# ggplot()+
#   geom_sf(data=pelim_sum_fy,aes(fill=fy_mn))+
#   geom_sf(data=elim_areas %>% filter(Area_Name!="0"),aes(col=Area_Name),
#           fill="transparent",lwd=1.5)+
#   facet_wrap(.~fy,ncol=1)+
#   scale_fill_viridis_c(#limits=c(0,1),
#     direction=-1,
#     name="Feral swine \nelimination probability \nafter 10 negative \ntrap nights")+
#   scale_color_discrete(name="Elimination Area")+
#   theme(panel.background=element_blank(),
#         axis.ticks=element_blank(),
#         axis.text=element_blank(),
#         plot.margin = margin(1,1,1,1, "cm"),
#         text=element_text(size=15))
# 
# ggsave(filename=paste0("./Model Outputs/Plots/Manuscript/pelim_all_yrs.jpeg"),
#        device="jpeg",width=8,height=10,units="in")

# ## text -------------
# pelim_sum[which.max(pelim_sum$md),]
# pelim_sum[which.min(pelim_sum$md),]
# 
# pabs_ea_sum[which.max(pabs_ea_sum$ea_md),]
# pabs_ea_sum[which.min(pabs_ea_sum$ea_md),]
# 
# #calculate by fy all sites
# pelim_sum_fy %>% 
#   st_drop_geometry() %>% 
#   group_by(fy) %>% 
#   summarise(median(fy_md),
            # mean(fy_mn))

# beta summary -----------
beta_sum %>% select(-mn)

#lambda summary --------------
lam_sum

#abundance  -----------------------------
## text --------------
N_yr_sum <- read.csv(paste0("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/Plots/",subfolder,"/N_change_table.csv"))
N_yr_sum

N_sum[which.min(N_sum$md_dens),] %>% select(SiteID,per_start,md_dens,lci_dens,uci_dens)
N_sum[which.max(N_sum$md_dens),] %>% select(SiteID,per_start,md_dens,lci_dens,uci_dens)

N_sum[which.min(N_sum$md_lambda_prop_dens),] %>% select(SiteID,per_start,md_lambda_prop_dens,lci_lambda_prop_dens,uci_lambda_prop_dens)
N_sum[which.max(N_sum$md_lambda_prop_dens),] %>% select(SiteID,per_start,md_lambda_prop_dens,lci_lambda_prop_dens,uci_lambda_prop_dens)


#removal needed ------------
N_sum[which.min(N_sum$md_lambda_prop),] %>% select(SiteID,per_start,md_lambda_prop,lci_lambda_prop,uci_lambda_prop)
N_sum[which.max(N_sum$md_lambda_prop),] %>% select(SiteID,per_start,md_lambda_prop,lci_lambda_prop,uci_lambda_prop)

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

##denstiy --------
axis_trans_d<- 0.03
dens_plot <- ggplot()+  
  geom_ribbon(data=N_sum,
              aes(x=per_start,ymin=lci_dens,ymax=uci_dens,
                  fill="Estimated Density"),alpha=0.2)+
  geom_line(data=N_sum,
            aes(x=per_start,y=md_dens,col="Estimated Density"),lwd=1)+
  # geom_bar(data=dat_rem_sum,
  #          aes(x=month,y=tot_rem_km/axis_trans_d),stat="identity",
  #          alpha=0.6)+
  geom_ribbon(data=N_sum,
              aes(x=per_start,ymin=lci_dens_no_rem,ymax=uci_dens_no_rem,fill="No Removal"),
              alpha=0.1)+
  geom_line(data=N_sum,
            aes(x=per_start,y=md_dens_no_rem,col="No Removal"),lwd=1,lty=3)+
  facet_wrap(.~SiteID)+
  # scale_y_continuous(sec.axis = sec_axis(transform=~.*axis_trans_d,
  #                                        name="No. feral swine removed per km2"))+
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
  theme(text=element_text(size=18),
        plot.margin = margin(c(0.5,0.5,0.5,1),unit="cm"))
# ggsave(filename=paste0("./Model Outputs/Plots/",subfolder,"/density_trend_removal.jpeg"),
#        device="jpeg",width=10,height=6,units="in")

#standardized abundnace -------
text_size<- 18
std_abund <- ggplot()+  
  geom_ribbon(data=N_sum,
              aes(x=per_start,ymin=lci_N_std,ymax=uci_N_std,
                  fill="Estimated \nstandardized \nabundance"),alpha=0.2)+
  geom_line(data=N_sum,
            aes(x=per_start,y=md_N_std,col="Estimated \nstandardized \nabundance"),lwd=1)+
  geom_ribbon(data=N_sum,
              aes(x=per_start,ymin=lci_N_ss_no_rem,ymax=uci_N_ss_no_rem,
                  fill="No removal"),
              alpha=0.1)+
  geom_line(data=N_sum,
            aes(x=per_start,y=md_N_ss_no_rem,col="No removal"),
            lty=3,lwd=1)+
  facet_wrap(.~SiteID)+
  scale_fill_manual(name="",
                    values=c("Estimated \nstandardized \nabundance"="black",
                             "No removal"="red"))+
  scale_color_manual(name="",
                     values=c("Estimated \nstandardized \nabundance"="black",
                              "No removal"="red"))+
  ylab("Standardized feral swine abundance")+
  xlab("Season")+
  theme(text=element_text(size=text_size),
        plot.margin = margin(c(0.5,0.5,0.5,1),unit="cm"))

## removal needed -----------
lambda_plot <- ggplot()+  
  geom_bar(data=dat_rem_sum,
           aes(x=month,y=tot_rem_km),stat="identity",
           alpha=0.6)+
  geom_ribbon(data=N_sum,
              aes(x=per_start,
                  ymin=lci_lambda_prop_dens,
                  ymax=uci_lambda_prop_dens,
                  fill="Removal needed \nto combat growth"),
              alpha=0.1)+
  geom_line(data=N_sum,aes(x=per_start,
                           y=md_lambda_prop_dens,group=SiteID,
                           col="Removal needed \nto combat growth"),lwd=1,lty=2)+
  
  facet_wrap(.~SiteID)+
  scale_fill_manual(name="",
                    values=c("Removal needed \nto combat growth"="blue"))+
  scale_color_manual(name="",
                     values=c("Removal needed \nto combat growth"="blue"))+
  # guides(col="none",fill="none")+
  ylab(expression(paste("Feral swine removal required to combat growth (swine/k",m^2,")")))+
  xlab("Season")+
  theme(text=element_text(size=text_size),
        plot.margin = margin(c(0.5,0.5,0.5,1),unit="cm"))
# ggsave(filename=paste0("./Model Outputs/Plots/",subfolder,"/density_trend_removal_needed.jpeg"),
#        device="jpeg",width=10,height=6,units="in")

cowplot::plot_grid(std_abund,lambda_plot,ncol=1,
                   labels = c("A","B"))
ggsave(filename=paste0("./Model Outputs/Plots/",subfolder,"/std_trend_removal_needed.jpeg"),
       device="jpeg",width=12,height=14,units="in")

# removal probability -----------------
## text ------
rem_df %>% select(-mn)

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
dat_rem_sum %>% 
  group_by(method) %>% 
  summarise(tot_rem=sum(removal),
            tot_hrs=sum(tot_hrs))

## removals plot --------
ggplot(dat_rem_sum)+
  geom_bar(aes(x=month,y=removal,fill=method),stat='identity')+
  facet_wrap(.~SiteID)+
  ylab("Number of feral swine removed")+
  xlab("Month")+
  scale_fill_discrete(name="Method",labels=c("Aerial operation","Trapping"))+
  theme(text=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1))
ggsave(filename="C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/Plots/Raw Data/raw_removal_ws.jpeg",
       width=8,height=6,units="in")

## effort plot --------
ggplot(dat_rem_sum)+
  geom_bar(aes(x=month,y=tot_hrs,fill=method),stat='identity')+
  facet_wrap(.~SiteID)+
  ylab("Number of effort hours spent removing feral swine")+
  xlab("Month")+
  scale_fill_discrete(name="Method",labels=c("Aerial operation","Trapping"))+
  theme(text=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1))
ggsave(filename="C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/Plots/Raw Data/raw_effort_ws.jpeg",
       width=8,height=6,units="in")

## cpue plot --------
ggplot(dat_rem_sum)+
  geom_bar(aes(x=month,y=removal/tot_hrs,fill=method),stat='identity')+
  facet_wrap(.~SiteID)+
  ylab("Number of feral swine removed per effort hour (CPUE)")+
  xlab("Month")+
  scale_fill_discrete(name="Method",labels=c("Aerial operation","Trapping"))+
  theme(text=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1))
ggsave(filename="C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model outputs/Plots/Raw Data/raw_cpue_ws.jpeg",
       width=8,height=6,units="in")

#occupancy data -----------
sysbait_det_eff %>% 
  summarise(ndets=sum(detection),
            ntraps=sum(trap_nights))

# validation ------------
load("C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/Model validation/Within sample/Plots/Single Lambda/modval_wi_occ_rem_17JAN25.Rdata")
pVal_occ
pVal_rem

pVal_occ_mn
pVal_rem_mn

pVal_occ_var
pVal_rem_var

auc_occ
auc_prev_pocc
