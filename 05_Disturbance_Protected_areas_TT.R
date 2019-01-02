
## Script to read in ESRI data tables and combine outputs

## this script calculates the areas of protection based on the consolidated protected layer (no-overlaps) 
## you will need to pre-clip the area of protection with the boundary of the widest extent in Arcmap. 


##  help files referenced: 
##https://gis.stackexchange.com/questions/265863/how-to-read-a-feature-class-in-an-esri-personal-geodatabase-using-r
##http://rstudio-pubs-static.s3.amazonaws.com/255550_1c983c8a6a5a45a0aee5a923206f4818.html
#http://www.rspatial.org/spatial/rst/7-vectmanip.html#spatial-queries

# read in libraries 
library(dplyr)
library(rgdal)
library(sp)
library(raster)
library(rgeos)
library(maptools)
library(magrittr)
library(tibble)
library(tidyr)
library(sf)
library(lwgeom)


# set work drive 
#out.dir = "X:/projects/Desktop_Analysis/data/output1/"
#temp.dir = "X:/projects/Desktop_Analysis/data/temp/"

# to run analysis on C drive: 
out.dir = "C:/Temp/TweedTelkwa/Temp/Perkins/Outputs/"
temp.dir = "C:/Temp/TweedTelkwa/Temp/Perkins/Data/"

Base= "C:/Temp/TweedTelkwa/Temp/Perkins/Data/Base_data.gdb"

# List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name))
base_list <- ogrListLayers(Base); print(base_list)

#####################################################################
# Read in herd boundary layers 
all.range <- st_read(dsn=Base,layer="TT_boundary"); plot(st_geometry(all.range))
all.range <- st_zm(all.range,drop = TRUE)
all.range <- st_cast(all.range,"MULTIPOLYGON")
all.range <- st_make_valid(all.range)
all.range.out <- data.frame(all.range)%>%
  dplyr::select(SiteName,V17_CH )

# read in the herd Key and detailed herd key 
Herd_key <- read.csv(paste(out.dir,"Herd_key.csv",sep = ""))
Herd_key_detail <- read.csv(paste(out.dir,"Herd_key_detail.csv",sep = ""))


##############################################################################################
# Read in protected layers 
prot <- st_read(dsn=Base,layer="Protect_clip") # ; st_crs(prot)# check CRS
prot <- st_zm(prot,drop = TRUE)
#prot <- st_transform(prot,crs = 3005) # ; plot(st_geometry(prot), add = T)
#prot = st_cast(prot,"POLYGON")
st_is_valid(prot)
prot <- st_make_valid(prot)
#unique(prot$desig)
st_is_valid(all.range)
all.range<- st_make_valid(all.range)

## Range boundary # intersect protection layer with ranges 
    r.prot <-  st_intersection(all.range,prot)
    r.prot$area.m = as.numeric(st_area(r.prot)) 
  
    ### Convert to tables and then calculate the areas 
    r.prot.df = data.frame(r.prot)        # calculate the length per range 
    r.prot.all.summary <- r.prot.df %>% 
      group_by(SiteName,V17_CH) %>% 
      summarise(R_protect = sum(area.m))

    # calculate per type and range 
    r.prot.summary <-  r.prot.df %>% 
      group_by(SiteName,V17_CH,desig) %>% 
      summarise(R_Pro.Type.area.m2 = sum(area.m))
   
 # Add the disturbance data to the herd key information
prot.key <- left_join(Herd_key_detail,r.prot.summary)  

write.csv(prot.key,paste(out.dir,"Draft_protection.csv",sep= "")) 
    
######################################

# Format into table to match other formats     
tout =prot.key

# detailed breakdown by protection type 
tout.by.protect <- tout %>% 
  mutate(Area_ha =  R_Pro.Type.area.m2/10000) %>%# convert to ha for all values and calculate the % of protected areas.
  mutate(Percentage = (Area_ha/R_area_ha) * 100) %>%
  dplyr::select(c(SiteName,V17_CH,desig,Area_ha,Percentage))

tout.area <- tout.by.protect %>% dplyr::select(-c(Percentage)) %>%  spread(., desig, value = Area_ha)
tout.perc<- tout.by.protect %>% dplyr::select(-c(Area_ha)) %>%  spread(., desig, value =Percentage)

# summary for all protection combined
tout.all <- tout%>% 
  group_by(SiteName, V17_CH,R_area_ha)%>%
  summarise(Area.m = sum(R_Pro.Type.area.m2)) %>%
  mutate(Area.ha = Area.m /10000) %>%
  mutate(Percentage = (Area.ha/R_area_ha) * 100) %>%
  dplyr::select(c(SiteName,V17_CH,R_area_ha,Area.ha,Percentage))

# join the tables together

protect.summary.ha = left_join(tout.all,tout.area, by = c('SiteName','V17_CH'))
protect.summary.pc = left_join(tout.all,tout.perc, by = c('SiteName','V17_CH'))

## may also want to reshape to convert to habitat catergories as columns rather than row. 
out.ha =protect.summary.ha %>% reshape2:::melt(.) %>% spread(.,V17_CH,value)
out.ha[is.na(out.ha)]<-0

out.pc=protect.summary.pc %>% reshape2:::melt(.) %>% spread(.,V17_CH,value)
out.pc[is.na(out.pc)]<-0

## Perhaps fix this part also to 
write.csv(out.ha,paste(out.dir,"Final_protect_summary_ha.csv",sep= ""))  
write.csv(out.pc,paste(out.dir,"Final_protect_summary_pc.csv",sep= ""))

