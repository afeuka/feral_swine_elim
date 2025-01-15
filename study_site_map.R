### Title: Map of study study for MO ocucpancy
### Author: Abbey Feuka
### Date: 15JAN25
### Notes: 

subfolder <- "Single Lambda"

#study site map -----------------------------
us_states <- states()
mo <- us_states %>% filter(STUSPS=="MO")
huc10 <- st_read("C:/Users/Abigail.Feuka/OneDrive - USDA/GIS Data/Missouri/huc10_mo.shp")

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
ggsave(filename = "./Model outputs/Plots/",subfolder,"/mo_map.jpeg",
       width=7,height=5,units="in",device="jpeg")