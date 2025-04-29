### Title: Map of study study for MO ocucpancy
### Author: Abbey Feuka
### Date: 15JAN25
### Notes: 
library(tigris)
library(sf)
library(tidyverse)

#study site map -----------------------------
us_states <- states()
mo <- us_states %>% filter(STUSPS=="MO")
huc10 <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/huc10_mo.shp")
study_site_grid <- st_read('C:/Users/Abigail.Feuka/OneDrive - USDA/Feral Hogs/Missouri/nimble/feral_swine_elim/data/input/HUC10_MO/huc10_siteIDs_cond_EA.shp')
elim_areas <- study_site_grid %>% 
  rename(Area_Name=Area_Nm) %>% 
  group_by(Area_Name) %>% 
  summarise(geometry=st_union(geometry))

ggplot()+geom_sf(data=mo)+
  geom_sf(data=huc10,fill="grey88")+
  geom_sf(data=study_site_grid,fill="grey55")+
  geom_sf(data=elim_areas %>% filter(Area_Name!="0"),
          aes(col=Area_Name),fill="transparent",
          lwd=1.5)+
  scale_color_discrete(name="Elimination Area")+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        text=element_text(size=15))
ggsave(filename = file.path('output','plots',"mo_map.jpeg"),
       width=7,height=5,units="in",device="jpeg")

#removal sites
rem_sites <- study_site_grid %>% filter(SiteID%in%c(115,161,160,163,43,44))

ggplot()+
  geom_sf(data=study_site_grid %>% 
            filter(Area_Nm%in%c(2,4)),fill="grey88")+
  geom_sf(data=elim_areas %>% filter(Area_Name%in%c(2,4)),
          aes(col=Area_Name),fill="transparent",
          lwd=1.5)+
  scale_color_discrete(name="Elimination Area")+
  # scale_fill_discrete(name="Removal site")+
  geom_sf(data=rem_sites,aes(fill=factor(SiteID))) +
  geom_sf_label(data=rem_sites,aes(label=SiteID),size=6,alpha=0.6,
                nudge_x = c(0,0,0,10000,-1000,-10000),
                nudge_y = c(-10000,0,0,5000,1000,0))+
  xlab("")+ylab("")+
  guides(fill='none')+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        text=element_text(size=15))

ggsave(filename = file.path('output','plots',"rem_map.jpeg"),
       width=7,height=5,units="in",device="jpeg")
