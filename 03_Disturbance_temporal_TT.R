# Copyright 2019 Province of British Columbia
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and limitations under the License.





##
## Caribou disturbance analysis 2018
##
## August 14th 2018
##
## written by genevieve perkins (genevieve.perkins@gov.bc.ca)
##
## This requires initial preparation of layers within arcmap.
## Step 1)
##
## Assemble layers as per datasheet in arcmap mxd
## Clip the layers to the range with the largest extent (for example range boundary for Boreal (not core))
## Create a filegeodatabdase and output these to the given data base.
##
## Step 2)
## Run through the script 01_disturbance layers first!
##
## Step 3)
## Run this script
## You may need to adjust
## - the names of files to match your arcmap exports
## - the directory/folder sructure.
##
## General notes:
## This script is used to look at more detailed break-down of cutblocks, fire and pests and to create plots


## Read in packages and libraries required:

#install.packages(c("rgdal","sp","dplyr","raster","rgeos","maptools","magrittr","tibble",
#				"tidyr","sf","lwgeom","mapview"),dep = T )
library(ggplot2)
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
library(mapview)

# set work drive
#out.dir = "X:/projects/Desktop_Analysis/data/output1/"
#temp.dir = "X:/projects/Desktop_Analysis/data/temp/"

# to run analysis on C drive:
out.dir = "C:/Temp/TweedTelkwa/Temp/Perkins/Outputs/"
temp.dir = "C:/Temp/TweedTelkwa/Temp/Perkins/Data/"
Base= "C:/Temp/TweedTelkwa/Temp/Perkins/Data/Base_data.gdb"

# to run on D drive
out.dir = "D:/Temp/TweedTelkwa/Outputs/"
temp.dir = "D:/Temp/TweedTelkwa/Data/"    #Z:\01.Projects\Wildlife\Caribou\02.Disturbance\TweedTelkwa\Temp\Perkins\Data
shape.output.dir = "D:/Temp/TweedTelkwa/Outputs/disturb_layers/" #Z:\01.Projects\Wildlife\Caribou\02.Disturbance\TweedTelkwa\Temp\Perkins\Outputs\disturb_layers
Base= "D:/Temp/TweedTelkwa/Data/Base_data.gdb" #"Z:\01.Projects\Wildlife\Caribou\02.Disturbance\TweedTelkwa\Temp\Perkins\Data\Base_data.gdb"



# List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name))
base_list <- ogrListLayers(Base); print(base_list)

#####################################################################
# Read in herd boundary layers
all.range <- st_read(dsn=Base,layer="TT_boundary") #; plot(st_geometry(all.range))
all.range <- st_zm(all.range,drop = TRUE)
all.range <- st_cast(all.range,"MULTIPOLYGON")
all.range <- st_make_valid(all.range)
all.range.out <- data.frame(all.range)%>%
  dplyr::select(SiteName,V17_CH )

# read in the herd Key and detailed herd key
Herd_key <- read.csv(paste(out.dir,"Herd_key.csv",sep = ""))
Herd_key_detail <- read.csv(paste(out.dir,"Herd_key_detail.csv",sep = ""))

#dha <- read.csv(paste(out.dir,"Final_TT_summary_ha.csv",sep = ""))
#dpc <- read.csv(paste(out.dir,"Final_TT_summary_pc.csv",sep = ""))

##############################################################################################

#  TEMPORAL DISTURBANCE  - read in data generated in 01_Disturbance_layers_TT.R

# Cutblocks,
# Burns
# Pests

############################################################
# Cutblocks
############################################################

# HERD 1)  ## Telkwa
b.r.c= st_read(dsn = Base , layer = "cutblock_union_Te")      # read in file
#b.r.c <- st_cast(b.r.c,"POLYGON")                             # confirm the geometry is polygon
b.r.c <- st_zm(b.r.c ,drop = TRUE)                            # drop the z co-ordinate
# get both values of HARVEST YEAR (as this this two cutblock layers unioned together - we need to get a single value for harvest year )
b.r.c$HARVEST_YEAR_ALL = ifelse(b.r.c$HARVEST_YEAR > 0,b.r.c$HARVEST_YEAR,b.r.c$HARVEST_YEAR_1 ) #sort(unique(b.r.c$HARVEST_YEAR_ALL)) # error check
b.r.c$TimeSinceCut = 2018-b.r.c$HARVEST_YEAR_ALL; # create new column with age since cut convert this to years
#sort(unique(b.r.c$TimeSinceCut)) # check the range in years since cut #note in the boreal the oldest age cut is 78 years
b.r.c <- st_intersection(all.range,b.r.c)
r.cut.df <- b.r.c

# add a column to differentiate the age brackets of cutblocks by decades (this is using the spatial data)
r.cut.df <- mutate(r.cut.df,dec.period = ifelse(HARVEST_YEAR_ALL >= 1940 & HARVEST_YEAR_ALL <= 1949,1940,0))
r.cut.df <- mutate(r.cut.df,dec.period = ifelse(HARVEST_YEAR_ALL >= 1950 & HARVEST_YEAR_ALL <= 1959,1950,dec.period))
r.cut.df <- mutate(r.cut.df,dec.period = ifelse(HARVEST_YEAR_ALL >= 1960 & HARVEST_YEAR_ALL <= 1969,1960,dec.period))
r.cut.df <- mutate(r.cut.df,dec.period = ifelse(HARVEST_YEAR_ALL >= 1970 & HARVEST_YEAR_ALL <= 1979,1970,dec.period))
r.cut.df <- mutate(r.cut.df,dec.period = ifelse(HARVEST_YEAR_ALL >= 1980 & HARVEST_YEAR_ALL <= 1989,1980,dec.period))
r.cut.df <- mutate(r.cut.df,dec.period = ifelse(HARVEST_YEAR_ALL >= 1990 & HARVEST_YEAR_ALL <= 1999,1990,dec.period))
r.cut.df <- mutate(r.cut.df,dec.period = ifelse(HARVEST_YEAR_ALL >= 2000 & HARVEST_YEAR_ALL <= 2009,2000,dec.period))
r.cut.df <- mutate(r.cut.df,dec.period = ifelse(HARVEST_YEAR_ALL >= 2010 & HARVEST_YEAR_ALL <= 2019,2010,dec.period))
#r.cut.df[r.cut.df$dec.period == 0,] ; unique(r.cut.df$dec.period)  # error check

## generate table output the amount of cutblock by range (all years (0-80))
#r.cut.df.df = as.data.frame(r.cut.df)

################################################################################################################
# Figure 1: make a plot of temporal and core.p.range for cutblocks
  b.range.temp <-  all.range[all.range$SiteName == "Telkwa",]
  p1 = ggplot() +
          geom_sf(data = b.range.temp) +
          #geom_sf(data = r.cut.df, col = "red") + facet_grid(.~ dec.period)+
          geom_sf(data = r.cut.df, col = "red") + facet_wrap(~dec.period) #+
   #theme_()
  p1
  ggsave(paste(out.dir,"Telkwa_cutblock_decades.png",sep = ""))

############################################################################

  # HERD 2)  ## Tweedsmuir herds
  # cutblock data from Canfor 2016 - 2018 data set
  #b.r.c1 = st_read(paste("C:\\Temp\\TweedTelkwa\\Temp\\Perkins\\Data\\whitesail_shapefiles_2018_12_10","Whitesail_harvested_blocks_cfp_2018_12_11.shp",sep = "\\"))
  b.r.c1 = st_read(paste(temp.dir,"\\whitesail_shapefiles_2018_12_10\\Whitesail_harvested_blocks_cfp_2018_12_11.shp",sep = ""))
  b.r.c1  = st_transform(b.r.c1 ,3005)
  #unique(b.r.c1$HS_DATE)
  b.r.c1<- st_intersection(all.range,b.r.c1) ; plot(st_geometry(b.r.c1 ))
  ########################
  b.r.c2= st_read(dsn = Base , layer = "cutblock_union_Tw")
  b.r.c2 <- st_zm(b.r.c2 ,drop = TRUE)
  # get both values of HARVEST YEAR
  b.r.c2$HARVEST_YEAR_ALL = ifelse(b.r.c2$HARVEST_YEAR > 0,b.r.c2$HARVEST_YEAR,b.r.c2$HARVEST_YEAR_1 ) #sort(unique(b.r.c$HARVEST_YEAR_ALL)) # error check b.r.c2$TimeSinceCut = 2018-b.r.c2$HARVEST_YEAR_ALL; # create new column with age since cut
  b.r.c2$TimeSinceCut = 2018-b.r.c2$HARVEST_YEAR_ALL; # create new column with age since cut
  #sort(unique(b.r.c$TimeSinceCut)) # check the range in years since cut #note in the boreal the oldest age cut is 78 years
  #b.r.c2 <- st_cast(b.r.c2,"POLYGON")
  b.r.c2 <- st_make_valid(b.r.c2)

  #up to here

  r.cut.df2 = b.r.c2
  r.cut.df3 = b.r.c1
  # add a column to differentiate the age brackets of cutblocks
  r.cut.df2 <- mutate(r.cut.df2,dec.period = ifelse(HARVEST_YEAR_ALL >= 1940 & HARVEST_YEAR_ALL <= 1949,1940,0))
  r.cut.df2 <- mutate(r.cut.df2,dec.period = ifelse(HARVEST_YEAR_ALL >= 1950 & HARVEST_YEAR_ALL <= 1959,1950,dec.period))
  r.cut.df2 <- mutate(r.cut.df2,dec.period = ifelse(HARVEST_YEAR_ALL >= 1960 & HARVEST_YEAR_ALL <= 1969,1960,dec.period))
  r.cut.df2 <- mutate(r.cut.df2,dec.period = ifelse(HARVEST_YEAR_ALL >= 1970 & HARVEST_YEAR_ALL <= 1979,1970,dec.period))
  r.cut.df2 <- mutate(r.cut.df2,dec.period = ifelse(HARVEST_YEAR_ALL >= 1980 & HARVEST_YEAR_ALL <= 1989,1980,dec.period))
  r.cut.df2 <- mutate(r.cut.df2,dec.period = ifelse(HARVEST_YEAR_ALL >= 1990 & HARVEST_YEAR_ALL <= 1999,1990,dec.period))
  r.cut.df2 <- mutate(r.cut.df2,dec.period = ifelse(HARVEST_YEAR_ALL >= 2000 & HARVEST_YEAR_ALL <= 2009,2000,dec.period))
  r.cut.df2 <- mutate(r.cut.df2,dec.period = ifelse(HARVEST_YEAR_ALL >= 2010 & HARVEST_YEAR_ALL <= 2019,2010,dec.period))

  # add a column to differentiate the age brackets of cutblocks
  r.cut.df3 <- mutate(r.cut.df3,dec.period = 2010)

  ################################################################################################################
  # Figure 1: make a plot of temporal and core.p.range for cutblocks
  b.range.temp <-  all.range[all.range$SiteName == "Tweedsmuir",]
  p1 = ggplot() +
    geom_sf(data = b.range.temp) +
    #geom_sf(data = r.cut.df, col = "red") + facet_grid(.~ dec.period)+
    geom_sf(data = r.cut.df2 , col = "red") + facet_wrap(~dec.period) +
    geom_sf(data = r.cut.df3,col = "red") + facet_wrap(~dec.period)
  #theme_()
  p1
  ggsave(paste(out.dir,"Tweeds_cutblock_decades.png",sep = ""))


  ############################################################################
  ## Figure 2: Percent Area of cutblock per habitat type.

  out.data = read.csv(paste(out.dir,"Total_disturb_summary_unformatted.csv",sep =""))
  out.ha =out.data %>% reshape2:::melt(.) %>% spread(.,V17_CH,value)
  out.ha[is.na(out.ha)]<-0

  dha <- read.csv(paste(out.dir,"Final_TT_unbuf_summary_ha.csv",sep = ""))
  dpc <- read.csv(paste(out.dir,"Final_TT_unbuf_summary_pc.csv",sep = ""))

# Figure 2: make a bar graph of the habitat types per herd.
# format data together


  ## STILL TO FIX UP ####################################


  ### NOT RUN

## select only the pc rows.
#dpc=dpc[-(grep("X",dpc$variable)),]
#dpc=dpc[(grep("_pc",dpc$variable)),]
#dpc = dpc[-(grep("R_area_ha_pc",dpc$variable)),]

# for cut data # generate a year for each time period.
dpcc = dpc[(grep("_cut_",dpc$variable)),]
dpcc = dpcc %>%
  mutate(dec.period  = gsub("R_cut_","",variable))%>%
  mutate(dec.period = gsub("_m2_pc","",dec.period))


# output the details per core area:
herds = c('Telkwa','Tweedsmuir')

for (i in herds) {
  #i = herds[1]
  # write out core area details.
  t.c = dpcc %>% filter(SiteName == paste(i)) %>%
    dplyr:: select(SiteName,dec.period,HWSR,LSR,LWR,MR)

  if(i == 'Telkwa'){
  t.c = t.c %>% dplyr::select(-c('LSR','MR'))
  t.c.c = t.c %>% gather(key = var_name, value = value, 3:ncol(t.c))
  # plot 1
  p = ggplot(t.c.c, aes(dec.period,value)) + try(facet_wrap(~var_name))+ geom_bar(stat = 'identity')
  p = p + ylab("Percentage") + xlab("Decade") + ggtitle(paste("Percent Cutblock (per area) from 1960 - 2017 (",i," herd)",sep = "" ))
  p = p + theme_bw()
  ggsave(paste(out.dir,i,"_cutblock_barchart.png",sep = ""))
  } else {

  t.c.c = t.c %>% gather(key = var_name, value = value, 3:ncol(t.c))
  # plot 1
  p = ggplot(t.c.c, aes(dec.period,value)) + try(facet_wrap(~var_name))+ geom_bar(stat = 'identity')
  p = p + ylab("Percentage") + xlab("Decade") + ggtitle(paste("Percent Cutblock (per area) from 1960 - 2017 (",i," herd)",sep = "" ))
  p = p + theme_bw()
  p
  ggsave(paste(out.dir,i,"_cutblock_barchart.png",sep = ""))
  }

#}

  ## NOT RUN UP TO HERE




#################################################################
# Burns
#################################################################
  b.r.0 = st_read(dsn = Base, layer ="fire_clip")
  b.r.0<- st_zm(b.r.0 ,drop = TRUE) # this is a linear feature so need to buffer to estimate area calcs
  #st_is_valid(b.r.0)
  b.r.0$TimeSinceBurn = 2018-b.r.0$FIRE_YEAR #; plot(b.r.0$Shape)
  b.r.0 <- st_intersection(all.range,b.r.0) #; st_is_valid(b.r.0)
  #b.r.0 <- st_cast(b.r.0,"POLYGON")

# split burns into decades (using all the entire data set)
b.r.00.df <- b.r.0[b.r.0$TimeSinceBurn <81,];

# add a column to differentiate the age brackets of cutblocks
b.r.00.df <- mutate(b.r.00.df,dec.period = ifelse(FIRE_YEAR >= 1940 & FIRE_YEAR <= 1949,1940,0))
b.r.00.df <- mutate(b.r.00.df,dec.period = ifelse(FIRE_YEAR >= 1950 & FIRE_YEAR <= 1959,1950,dec.period))
b.r.00.df<- mutate(b.r.00.df,dec.period = ifelse(FIRE_YEAR >= 1960 & FIRE_YEAR <= 1969,1960,dec.period))
b.r.00.df<- mutate(b.r.00.df,dec.period = ifelse(FIRE_YEAR >= 1970 & FIRE_YEAR <= 1979,1970,dec.period))
b.r.00.df<- mutate(b.r.00.df,dec.period = ifelse(FIRE_YEAR >= 1980 & FIRE_YEAR <= 1989,1980,dec.period))
b.r.00.df<- mutate(b.r.00.df,dec.period = ifelse(FIRE_YEAR >= 1990 & FIRE_YEAR <= 1999,1990,dec.period))
b.r.00.df<- mutate(b.r.00.df,dec.period = ifelse(FIRE_YEAR >= 2000 & FIRE_YEAR <= 2009,2000,dec.period))
b.r.00.df<- mutate(b.r.00.df,dec.period = ifelse(FIRE_YEAR >= 2010 & FIRE_YEAR <= 2019,2010,dec.period))

#b.r.00.df[b.r.00.df$dec.period == 0,] # error check to see what is over 80 years old
b.r.00.df = b.r.00.df %>% filter(dec.period>0)
#unique(b.r.00.df$dec.period) # error check

plot(st_geometry(b.r.00.df))
plot(st_geometry(all.range),add = TRUE)

# generate cumulative burn disturbance shapefiles to be added sequentially to "static Disturbance"
Burn.dec = b.r.00.df


###############

## Figure 3:

# output the details per core area:
herds = c('Telkwa','Tweedsmuir')

for (i in herds) {
  #i = herds[1]
  # write out core area details.
  t.c = b.r.00.df %>% filter(SiteName == paste(i))

  b.range.temp <-  all.range[all.range$SiteName ==  paste(i),]
  p1 = ggplot() +
    geom_sf(data = b.range.temp) +
    #geom_sf(data = r.cut.df, col = "red") + facet_grid(.~ dec.period)+
    geom_sf(data = t.c , col = "red") + facet_wrap(~dec.period) #+
  #theme_()
  p1
  ggsave(paste(out.dir,paste(i),"_burn_decades.png",sep = ""))

}



##################################################################################
# DATA SET :  PEST DATA
##################################################################################

r.pest.te <-  st_read(dsn = Base, layer ="Pest_clip_Te_IBMIBS") #; plot(st_geometry(r.pest.te))

# Telkwa
r.pest <- r.pest.te
r.pest<- st_zm(r.pest ,drop = TRUE) # this is a linear feature so need to buffer to estimate area calcs
#r.pest = st_cast(r.pest,"POLYGON")
r.pest <- st_make_valid(r.pest)
r.pest$TimeSincePest = 2018-r.pest$CAPTURE_YEAR
r.pest<- st_intersection(all.range,r.pest) #; st_is_valid(b.r.0)
#r.pest <- st_cast(r.pest,"POLYGON")
#plot(st_geometry(r.pest), col = 'red')
#plot(st_geometry(all.range),add = T)

r.pest.df = r.pest

# add decade data   # add a column to differentiate the age brackets of pest capture
r.pest.df<- r.pest.df %>% mutate(dec.period = ifelse(CAPTURE_YEAR >= 1960 & CAPTURE_YEAR <= 1969,1960,0))
r.pest.df<- mutate(r.pest.df,dec.period = ifelse(CAPTURE_YEAR >= 1970 & CAPTURE_YEAR <= 1979,1970,dec.period))
r.pest.df<- mutate(r.pest.df,dec.period = ifelse(CAPTURE_YEAR >= 1980 & CAPTURE_YEAR <= 1989,1980,dec.period))
r.pest.df<- mutate(r.pest.df,dec.period = ifelse(CAPTURE_YEAR >= 1990 & CAPTURE_YEAR <= 1999,1990,dec.period))
r.pest.df<- mutate(r.pest.df,dec.period = ifelse(CAPTURE_YEAR >= 2000 & CAPTURE_YEAR <= 2009,2000,dec.period))
r.pest.df<- mutate(r.pest.df,dec.period = ifelse(CAPTURE_YEAR >= 2010 & CAPTURE_YEAR <= 2018,2010,dec.period))

# output the details per core area:
  b.range.temp <-  all.range[all.range$SiteName == "Telkwa",]
  r.pest.df.temp <- r.pest.df[r.pest.df$SiteName == "Telkwa",]


  p1 = ggplot() +
    geom_sf(data = b.range.temp) +
    #geom_sf(data = r.cut.df, col = "red") + facet_grid(.~ dec.period)+
    geom_sf(data = r.pest.df.temp , col = "red") + facet_wrap(~dec.period) #+
  #theme_()
  p1
  ggsave(paste(out.dir,"Telkwa_pest_decades.png",sep = ""))

################################

# Tweedsmuir
r.pest.tw <-  st_read(dsn = Base, layer ="Pest_clip_Tw_IBMIBS") #; plot(st_geometry(r.pest.tw))
r.pest2 <- r.pest.tw
r.pest2<- st_zm(r.pest2 ,drop = TRUE) # this is a linear feature so need to buffer to estimate area calcs
#r.pest2 = st_cast(r.pest2,"POLYGON")#; st_is_valid(r.pest2)
r.pest2 = st_make_valid(r.pest2) # ; st_is_valid(r.pest2)
r.pest2$TimeSincePest = 2018-r.pest2$CAPTURE_YEAR

## Up to here
r.pest2<- st_intersection(all.range,r.pest2) #; st_is_valid(all.range)
#plot(st_geometry(r.pest2), col = 'red') ; plot(st_geometry(all.range),add = T)

r.pest.df2 = r.pest2

# add decade data   # add a column to differentiate the age brackets of pest capture
r.pest.df2<- r.pest.df2 %>% mutate(dec.period = ifelse(CAPTURE_YEAR >= 1960 & CAPTURE_YEAR <= 1969,1960,0))
r.pest.df2<- mutate(r.pest.df2,dec.period = ifelse(CAPTURE_YEAR >= 1970 & CAPTURE_YEAR <= 1979,1970,dec.period))
r.pest.df2<- mutate(r.pest.df2,dec.period = ifelse(CAPTURE_YEAR >= 1980 & CAPTURE_YEAR <= 1989,1980,dec.period))
r.pest.df2<- mutate(r.pest.df2,dec.period = ifelse(CAPTURE_YEAR >= 1990 & CAPTURE_YEAR <= 1999,1990,dec.period))
r.pest.df2<- mutate(r.pest.df2,dec.period = ifelse(CAPTURE_YEAR >= 2000 & CAPTURE_YEAR <= 2009,2000,dec.period))
r.pest.df2<- mutate(r.pest.df2,dec.period = ifelse(CAPTURE_YEAR >= 2010 & CAPTURE_YEAR <= 2018,2010,dec.period))


# output the details per core area:
b.range.temp <-  all.range[all.range$SiteName == "Tweedsmuir",]
r.pest.df.temp <- r.pest.df2[r.pest.df2$SiteName =="Tweedsmuir",]

p1 = ggplot() +
  geom_sf(data = b.range.temp) +
  #geom_sf(data = r.cut.df, col = "red") + facet_grid(.~ dec.period)+
  geom_sf(data = r.pest.df.temp , col = "red") + facet_wrap(~dec.period) #+
#theme_()
p1
ggsave(paste(out.dir,"Tweedsmuir_pest_decades.png",sep = ""))


################################################################################

