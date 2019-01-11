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



## Caribou disturbance analysis 2018
##
## August 14th 2018
##
## written by genevieve perkins (genevieve.perkins@gov.bc.ca)
##
## This requires initial preparation of layers within arcmap.
## Step 1)
## Assemble layers as per datasheet in arcmap mxd
## Clip the layers to the range with the largest extent (for example range boundary for Boreal (not core))
## Create a filegeodatabdase and output these to the given data base.
##
## Step 2)
## Run through the script below. You may need to adjust
## - the names of files to match your arcmap exports
## - the directory/folder sructure.
##
## General notes:
## For each disturbance layers the script will read in, intersect with range and core areas and calculate the area and or length. The peripery area will be calculated for each herd as well.
## With each layer the compiled disturbance will also be calculated.

## Associated help files for reference:
##https://gis.stackexchange.com/questions/265863/how-to-read-a-feature-class-in-an-esri-personal-geodatabase-using-r
##http://rstudio-pubs-static.s3.amazonaws.com/255550_1c983c8a6a5a45a0aee5a923206f4818.html
#http://www.rspatial.org/spatial/rst/7-vectmanip.html#spatial-queries
#https://r-spatial.github.io/sf/articles/sf1.html
#https://github.com/r-spatial/sf/issues/290

## Read in packages and libraries required:

#install.packages(c("rgdal","sp","dplyr","raster","rgeos","maptools","magrittr","tibble",
#			"tidyr","sf","lwgeom","mapview"),dep = T )

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

## set your output directory

# to run analysis on C drive:
#out.dir = "C:/Temp/TweedTelkwa/Temp/Perkins/Outputs/"
#temp.dir = "C:/Temp/TweedTelkwa/Temp/Perkins/Data/"
#shape.output.dir = "C:/Temp/TweedTelkwa/Temp/Perkins/Outputs/disturb_layers/"

# to run on D drive
out.dir = "D:/Temp/TweedTelkwa/Outputs/"
temp.dir = "D:/Temp/TweedTelkwa/Data/"    #Z:\01.Projects\Wildlife\Caribou\02.Disturbance\TweedTelkwa\Temp\Perkins\Data
shape.output.dir = "D:/Temp/TweedTelkwa/Outputs/disturb_layers/" #Z:\01.Projects\Wildlife\Caribou\02.Disturbance\TweedTelkwa\Temp\Perkins\Outputs\disturb_layers
Base= "D:/Temp/TweedTelkwa/Data/Base_data.gdb" #"Z:\01.Projects\Wildlife\Caribou\02.Disturbance\TweedTelkwa\Temp\Perkins\Data\Base_data.gdb"


## Set your input geodatabases (this will be where you saved your arcmap exports)
## edit these to your filepath and name of gdb
#Base= "C:/Temp/TweedTelkwa/Temp/Perkins/Data/Base_data.gdb" # contains

## List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name))
base_list <- ogrListLayers(Base); print(base_list)

##############################################################################################
# Read in herd boundary layers

all.range <- st_read(dsn=Base,layer="TT_boundary"); plot(st_geometry(all.range))
all.range <- st_zm(all.range,drop = TRUE)
all.range <- st_cast(all.range,"MULTIPOLYGON")
all.range <- st_make_valid(all.range)
all.range.out <- data.frame(all.range)%>%
  dplyr::select(SiteName,V17_CH )

# Read in summary of area and habitat created in script 1.
Herd_key = read.csv(paste(out.dir,"Herd_key.csv",sep = ""))
Herd_key_detail= read.csv(paste(out.dir,"Herd_key_detail.csv",sep = ""))

##############################################################################################
# Read in individual Disturbance layers:

# 1) pipeline Area - 500 m buffer
      r.pipe <- st_read(dsn=Base,layer="Pipeline_clip") # multistring
      r.pipe.int = st_intersection(all.range,r.pipe)   # intersect with single all ranges
      # RANGE: intersect with range
      # buffer to the 500m distance
      r.pipe.int <- st_buffer(r.pipe.int,5)  # buffer to 10m as per the other disturbance analysis
      r.pipe.int <- st_buffer(r.pipe.int,250)  ; all.pipe = sum(st_area(r.pipe)) ; plot(st_geometry(r.pipe))
      r.pipe.int.u <- st_union(r.pipe.int) # union to
      st_write(r.pipe.int.u,paste(shape.output.dir,"Dist_pipe_buf.shp",sep = "")) #write out individual dist_layer for Range
      #st_area(r.pipe.int.u )
      #st_is_valid(r.pipe.int.u)                       # check valid geometr     # fix overlaps
      r.pipe.int = st_intersection(all.range,r.pipe.int.u)
      r.pipe.int$Area.m <- as.numeric(st_area(r.pipe.int)) #; plot(st_geometry(r.pipe.int))

      r.pipe.int.df = data.frame(r.pipe.int)        # calcaulte area
      r.pipe.int.df.out  = r.pipe.int.df %>%
        group_by(SiteName,V17_CH) %>%
        summarise(R_Pipe_area_m = sum(Area.m))

      out.pipe = r.pipe.int
      all.range.out<- left_join(all.range.out, r.pipe.int.df.out)
      all.range.out[is.na(all.range.out)]<-0

# 2) transmission (length and area)
      r.tran.sf <- st_read(dsn=Base,layer="Trans_clip") # multistring   # read in data
      # 1) RANGE: calculate the range extent to use range extent
      r.tran <- st_intersection(all.range,r.tran.sf)# intersect with ranges
      r.tran <- st_buffer(r.tran,5) #;all.tran = sum(st_area(r.tran)) # buffer to 1m to convert from a line to a polygon
      r.tran <- st_buffer(r.tran,250) # add 500m buffer
      r.tran <- st_cast(r.tran,"POLYGON") # check geometry
      r.tran$area_m = as.numeric(st_area(r.tran)) # calculate the area of each polygon
      #st_is_valid(r.tran)                   # check valid geometr
      r.tran.df = data.frame(r.tran)        # Convert to dataframe
      r.tran.df.out  = r.tran.df %>%        # calculate the area per range
        group_by(SiteName,V17_CH) %>%
        summarise(R_Trans_area_m = sum(area_m))

      ##plot(st_geometry(r.tran))
      r.tran.u <- st_union(r.tran)
      st_write(r.tran,paste(shape.output.dir,"Dist_tran_buf.shp",sep = ""))       #write out individual dist_layer for Range
      #out.trans <- r.tran.df.out
      all.range.out<- left_join(all.range.out, r.tran.df.out)
      all.range.out[is.na(all.range.out)]<-0

      ## 3) ALL DISTURBANCE UNION 1
      out = st_union(r.pipe.int,r.tran) ; plot(st_geometry(out))
      out = st_union(out); plot(st_geometry(out),add = T)
      out = st_cast(out,"POLYGON")
      ##x.area = sum(st_area(out)) ; x.area #298446.6

# 3) mine
      r.mine.sf <- st_read(dsn=Base,layer="Mining_clip") # multipoly
      r.mine <- st_zm(r.mine.sf,drop = TRUE) # drop the z portion of the shapefile.
      r.mine <- st_buffer(r.mine,250) # add 500m buffer
      r.mine <- st_union(r.mine)
      # 1) RANGE: calculate the range extent to use range extent
      r.mine <- st_intersection(all.range,r.mine)

      #st_is_valid(r.mine)                   # check valid geometr
      r.mine$Area.m <- as.numeric(st_area(r.mine))
      r.mine.u = st_union(r.mine)
      st_write(r.mine.u,paste(shape.output.dir,"Dist_mine_buf.shp",sep = ""))

      r.mine.df = data.frame(r.mine)        # calculate the length per range
      r.mine.df.out  = r.mine.df %>%
        group_by(SiteName,V17_CH) %>%
        summarise(R_mine_m2 = sum(Area.m))

      # combine into disturbance by layer
      all.range.out <- left_join(all.range.out,r.mine.df.out) ;
      all.range.out[is.na(all.range.out)]<-0

      ## 3) ALL DISTURBANCE UNION 2 # join together the spatial data for the successive disturbance layers
      out = st_make_valid(out)
      r.mine = st_make_valid(r.mine)
      out1 = st_union(out,r.mine)  ; plot(st_geometry(out1)) ; rm(out)
      #  #x.area = sum(st_area(out1)) ; x.area #20570863 m2 # error check


# BC oil and gas layer tested and no area found within the TT boundaries


# 4) agriculture
      r.agr.sf <- st_read(dsn=Base,layer="Agri_clip") # multipoly
      r.agr.sf <- st_zm(r.agr.sf,drop = TRUE)
      r.agr.sf <- st_buffer(r.agr.sf,250)
      r.agr.sf <- st_union(r.agr.sf) # need to union before intersect to avoid overlapping polys in area calcs
      # 1) RANGE: calculate the range extent to use range extent
      r.agr <- st_intersection(all.range,r.agr.sf)# intersect with ranges
      #r.agr <- st_cast( r.agr,"POLYGON")
      #st_is_valid(r.agr)                   # check valid geometr
      r.agr$Area.m <- as.numeric(st_area(r.agr))
      r.agr.u <- st_union(r.agr)
      st_write(r.agr.u,paste(shape.output.dir,"Dist_agri_buf.shp",sep = ""))

      r.agr.df = data.frame(r.agr)        # calculate the length per range
      r.agr.df.out  = r.agr.df %>%
        group_by(SiteName,V17_CH) %>%
        summarise(R_agr_m2 = sum(Area.m))
      #plot(st_geometry(r.agr))

      # combine into disturbance by layer
      all.range.out <- left_join(all.range.out,r.agr.df.out) ;
      all.range.out[is.na(all.range.out)]<-0

      ## ALL DISTURBANCE COMBINED: UNION 3
      out2 = st_union(out1,r.agr.sf) ; plot(st_geometry(out2))
      out2 = st_union(out2); plot(st_geometry(out2)) ; rm (out1)
      #x.area = sum(st_area(out2)); x.area # 144954516  #out2 = st_cast(out2,"POLYGON")

# 5) air strips
      r.air.sf <- st_read(dsn=Base,layer="Airstrip_clip") # multipoly
      r.air.sf <- st_zm(r.air.sf,drop = TRUE)
      r.air.sf <- st_buffer( r.air.sf ,250)
      r.air.sf  <- st_union( r.air.sf )

      #r.air.sf <- st_cast(r.air.sf,"POLYGON")

      # 1) RANGE: calculate the range extent to use range extent
      r.air <- st_intersection(all.range,r.air.sf)# intersect with ranges
      #r.air <- st_cast( r.air,"POLYGON")
      #st_is_valid(r.air)                   # check valid geometr
      r.air$Area.m <- as.numeric(st_area(r.air))
      r.air.u <- st_union(r.air)
      st_write(r.air.u,paste(shape.output.dir,"Dist_air_buf.shp",sep = ""))

      r.air.df = data.frame(r.air)        # calculate the length per range
      r.air.df.out  = r.air.df %>%
        group_by(SiteName,V17_CH) %>%
        summarise(R_air_m2 = sum(Area.m))
      #plot(st_geometry(r.air))

      # combine into disturbance by layer
      all.range.out <- left_join(all.range.out,r.air.df.out) ;
      all.range.out[is.na(all.range.out)]<-0

      ## UNION 4
      out3 = st_union(out2,r.air.sf) #; plot(st_geometry(out3))
      out3 = st_union(out3); plot(st_geometry(out3)) ; rm(out2)
      ##x.area = sum(st_area(out3)) ;x.area #145235941 [m^2]

# 6) dams
      r.dam.sf <- st_read(dsn=Base,layer="dams_clip" ) # multipoly
      r.dam.sf  <- st_zm(r.dam.sf ,drop = TRUE)
      r.dam.sf<- st_buffer(r.dam.sf ,250)
      r.dam.sf  <- st_union(r.dam.sf )

      ## RANGE: intersect with range and calculate length per range
      r.dams = st_intersection(all.range,r.dam.sf)   # intersect with ranges
      r.dams$Area.m <- as.numeric(st_area(r.dams))
      ##plot(st_geometry(r.dams))

      r.dams.df = data.frame(r.dams)        # calculate the length per range
      r.dams.df.out  = r.dams.df %>%
        group_by(SiteName,V17_CH) %>%
        summarise(R_Dams_m2 = sum(Area.m))
      r.dams.u <- st_union(r.dams)
      st_write(r.dams.u,paste(shape.output.dir,"Dist_dams_buf.shp",sep = ""))       #write out individual dist_layer for Range

      # combine into disturbance by layer
      all.range.out <- left_join(all.range.out,r.dams.df.out)
      all.range.out[is.na(all.range.out)]<-0

      # # UNION 5
      out4 = st_union(out3,r.dams.u) ; plot(st_geometry(out4))
      out4 = st_union(out4); plot(st_geometry(out4)) ; rm(out3)
      #x.area = sum(st_area(out4))  ;x.area    #145240406 [m^2]

# 7) reservoirs
##r.res.sf <- st_read(dsn=Base,layer="Resevoir_clip") ;  plot(st_geometry(r.res.sf))
# only one feature found (Nechako Reservoir) This is not included as it is not calculated as part of the area of the herds.

# 9) wells
      r.wells.sf <- st_read(dsn=Base,layer="Wells_clip" ) # multipoly
      # no wells present in these herd

      #r.wells.sf <- st_union(r.wells.sf)
      #all.wells = sum(st_area(r.wells.sf)) ; plot(st_geometry(r.wells.sf))

      ## All Disturbance : UNION 5 # may need to run this in stand alone R rather than R -studio
      #out6 = st_union(out4,r.wells.sf) ; plot(st_geometry(out6))
      ##out6 = st_union(out6); plot(st_geometry(out6))
      ##x.area = sum(st_area(out6)) ; x.area #392383839

      #out6 = out5

      ## Intersect with RANGE: intersect with range and calculate length per range
      r.wells = st_intersection(all.range,r.wells.sf)   # intersect with range
      r.wells <- st_cast(r.wells,"POLYGON")
      r.wells$area <- as.numeric(st_area(r.wells))

      ##plot(st_geometry(r.wells))
      #st_write(r.dams,"Dist_R_wells.shp")       #write out individual dist_layer for Range
      r.wells.df = data.frame(r.wells)        # calculate the length per range
      r.wells.df.out  = r.wells.df %>%
        group_by( SiteName,V17_CH ) %>%
        summarise(R_Wells_m2 = sum(area))

      # combine into disturbance by layer
      all.range.out <- left_join(all.range.out,r.wells.df.out)
      all.range.out[is.na(all.range.out)]<-0

 # 10) urban
r.urban.sf <- st_read(dsn=Base,layer="urban_clip" ) # multipoly
      r.urban.sf <- st_zm(r.urban.sf ,drop = TRUE)
      r.urban.sf <- st_buffer(r.urban.sf,250);
      r.urban.sf <- st_union(r.urban.sf);
     # r.urban.sf <- st_cast(r.urban.sf,"POLYGON")
      #all.urban = sum(st_area(r.urban.sf)) ; plot(st_geometry(r.urban.sf))

      ## All disturbance: UNION 6 # may need to run this in stand alone R rather than R -studio
      out7 = st_union(out4,r.urban.sf)  # or line below if there is wells present ; rm(out4)
      #out7 = st_union(out6,r.urban.sf)
      out7 = st_union(out7);
      #x.area = sum(st_area(out7)); x.area #396359128


      ## Intersect with RANGE: intersect with range and calculate length per range
      r.urban = st_intersection(all.range,r.urban.sf)   # intersect with range
      #r.urban = st_cast(r.urban,"POLYGON")
      r.urban$area <- as.numeric(st_area(r.urban))
      r.urban.u <- st_union(r.urban)
      ##plot(st_geometry(r.urban))
      st_write(r.urban.u,paste(shape.output.dir,"Dist_urban_buf.shp", sep = ""))      #write out individual dist_layer for Range
      r.urban.df = data.frame(r.urban)        # calculate the length per range
      r.urban.df.out  = r.urban.df %>%
        group_by( SiteName,V17_CH) %>%
        summarise(R_Urban_m2 = sum(area))
      # combine into disturbance by layer
      all.range.out <- left_join(all.range.out,r.urban.df.out)
      all.range.out[is.na(all.range.out)]<-0


# 11) Rail  # no rail in the herds
      r.rail.sf <- st_read(dsn=Base,layer="Rail_clip") # multipoly
      r.rail.sf <- st_intersection(all.range,r.rail.sf)
      #r.rail.sf <- st_cast(r.rail.sf,"POLYGON")
      #st_is_valid(r.rail.sf)
      #r.rail.sf <- st_make_valid((r.rail.sf))
      r.rail <- st_buffer(r.rail.sf,1)
      r.rail <- st_cast(r.rail,"POLYGON") ; all.rail= sum(st_area(r.rail)) ; plot(st_geometry(r.rail))

      # This union is not tested
      ## UNION 7 # may need to run this in stand alone R rather than R -studio
      #out8 = st_union(out7,r.rail) ;# plot(st_geometry(out7))
      #  out8 = st_union(out8);# plot(st_geometry(out8))
      #  #x.area = sum(st_area(out7))

      #out8 = out7 # need to comment out this line if you are running the union about (union7)

      r.rail$area.m <- as.numeric(st_area(r.rail))
      r.rail.df = data.frame(r.rail)        # calculate the length per range
      r.rail.df.out  = r.rail.df%>%
        group_by(SiteName,V17_CH) %>%
        summarise(R_Rail_m2 = sum(area.m))
      ##plot(st_geometry(r.rail))
      #st_write(r.rail,"Dist_rail_buf.shp")       #write out individual dist_layer for Range

      # combine into disturbance by layer
      all.range.out <- left_join(all.range.out,r.rail.df.out)
      all.range.out[is.na(all.range.out)]<-0

# 12)  Recreation sites
r.rec.sf <- st_read(dsn=Base,layer="Rec_clip") # multipoly
      r.rec.sf <- st_zm(r.rec.sf ,drop = TRUE)
      r.rec.sf <- st_buffer(r.rec.sf,250)
      r.rec.sf <- st_union(r.rec.sf)
      r.rec.sf <- st_intersection(all.range,r.rec.sf)

      #all.rec = sum(st_area(r.rec.sf)) ; plot(st_geometry(r.rec.sf))

      ## ALL DISTURBANCE: UNION 8 # may need to run this in stand alone R rather than R -studio
      out8 = out7; rm(out7)
      out9 = st_union(out8,r.rec.sf) ; plot(st_geometry(out9)); rm(out8)
      out9 = st_union(out9); plot(st_geometry(out9))
      #x.area = sum(st_area(out8))     ; x.area  # 396359128 [m^2]

      # 1) RANGE: calculate the range extent to use range extent
      #r.rec <- st_intersection(all.range,r.rec.sf)# intersect with range
      #r.rec <- st_buffer(r.rec,0)
      #r.rec <- st_cast(r.rec,"POLYGON")
      r.rec.sf$area = as.numeric(st_area(r.rec.sf))

      #all.rec = sum(st_area(r.rec))
      #st_is_valid(r.tran)                   # check valid geometr
      r.rec.df = data.frame(r.rec.sf)        # calculate the length per range
      r.rec.df.out  = r.rec.df%>%
        group_by(SiteName,V17_CH) %>%
        summarise(R_Rec_m2 = sum(area))

      r.rec.u <- st_union(r.rec.sf)
      ##plot(st_geometry(r.rec))
      st_write(r.rec.u,paste(shape.output.dir,"Dist_rec_buf.shp",sep = ""))       #write out individual dist_layer for Range

      # combine into disturbance by layer
      all.range.out <- left_join(all.range.out,r.rec.df.out)
      all.range.out[is.na(all.range.out)]<-0

# 13) Seismic Lines
      b.s1 = sf::st_read(dsn = Base , layer ="Seismic_clip")
      b.s1 <- st_zm(b.s1,drop = TRUE)
      b.s1 <- st_make_valid(b.s1)
      b.s1<- st_buffer(b.s1, 250)
      b.s1 <- st_union(b.s1)

      ## ALL DISTURBANCE: UNION 10 # may need to run this in stand alone R rather than R -studio
      out10 = st_union(out9,b.s1 ) ; plot(st_geometry(out10)) ; rm(out9)
      out10 = st_union(out10); plot(st_geometry(out10))
      ##x.area = sum(st_area(out8))     ; x.area  # 420881436 [m^2]

      b.s1 <- st_intersection(all.range,b.s1)# intersect with range
      #b.s1 <- st_cast(b.s1 , "POLYGON")
      #st_is_valid(b.s1)

      b.s1$area.m <- as.numeric(st_area(b.s1))

      b.s1.u <- st_union(b.s1)
      st_write(b.s1.u,paste(shape.output.dir,"Dist_seismic_buf.shp",sep = ""))  # write out individual layer

      b.s1.df = data.frame(b.s1)        # calculate the length per range
      b.s1.df.out  =  b.s1.df%>%
        group_by(SiteName,V17_CH) %>%
        summarise(R_Seis_m2 = sum(area.m))

      # combine into disturbance by layer
      all.range.out <- left_join(all.range.out,b.s1.df.out)
      all.range.out[is.na(all.range.out)]<-0


## 14) Roads
# this requires pre-processing in Arcmap map. Spli
# Buffered to total width of 15m (based on previous disturbance analysis and large number of smaller forestry roads in the tweedsmuir herd)
# Preparation includes split into herds and buffer to 15m width the dissolve into single polgons.
# due to the size of the herds this was split into smaller areas for the tweedsmuir herds due to the size

## Telkwa herd
      b.r1.sf = sf::st_read(dsn = Base , layer ="Te_road_buf_diss" )
      b.r1.sf = st_make_valid(b.r1.sf)
      b.r1.sf <- st_buffer(b.r1.sf,250)
      b.r1.sf <- st_union(b.r1.sf)
      b.r1 = st_intersection(all.range,b.r1.sf )  # intersect with single all ranges
      #st_is_valid(b.r1)
      b.r1$Area.m <- as.numeric(st_area(b.r1 ))
      #plot(st_geometry(b.r1))
      b.r1.u <- st_union(b.r1)
      st_write(b.r1.u,paste(shape.output.dir,"Dist_Te_road_buf.shp",sep = ""))      #write out individual dist_layer for Range
      b.r1.df = data.frame(b.r1)        # calcaulte area
      b.r1.df.out  = b.r1.df %>%
        group_by(SiteName,V17_CH) %>%
        summarise(R_Road_area_m = sum(Area.m))

## Tweedsmuir herd # Need to break this up into smaller sections then buffer and disolve within ArcMap
# HWSR
      b.r2.sf = sf::st_read(dsn = Base , layer ="Tw_HWSR_Rd_bu" )
      #st_is_valid(b.r2.sf)
      b.r2.sf = st_make_valid(b.r2.sf)
      b.r2.sf <- st_buffer(b.r2.sf,250)
      b.r2.sf <- st_union(b.r2.sf)

      b.r2 = st_intersection(all.range,b.r2.sf)
      b.r2$Area.m <- as.numeric(st_area(b.r2))
      b.r2.df = as.data.frame(b.r2)
      b.r2.df.out  = b.r2.df %>%
        group_by(SiteName,V17_CH) %>%
        filter(V17_CH == "High Elevation Winter/Summer Range") %>%
        summarise(R_Road_area_m = sum(Area.m))
      # add to the Telkwa data table and the spatial file.
      b.df.out <- rbind(b.r1.df.out,b.r2.df.out) # add to the roads summary table
      roads.union = st_union(b.r1,b.r2); rm(b.r1); rm(b.r2);  plot(st_geometry(roads.union))
      roads.union = st_union(roads.union)

# LWR
      b.r2.sf = sf::st_read(dsn = Base , layer ="Tw_LWR_Rd_bu"  )
      b.r2.sf = st_make_valid(b.r2.sf)
      b.r2.sf <- st_zm( b.r2.sf ,drop = TRUE)
      b.r2.sf <- st_buffer(b.r2.sf,250)
      b.r2.sf <- st_union(b.r2.sf)
      b.r2 = st_intersection(all.range,b.r2.sf)
      b.r2$Area.m <- as.numeric(st_area(b.r2))
      b.r2.df = as.data.frame(b.r2)
      b.r2.df.out  = b.r2.df %>%
        group_by(SiteName,V17_CH) %>%
        filter(V17_CH == "Low Elevation Winter Range") %>%
        summarise(R_Road_area_m = sum(Area.m))
      # add to the Telkwa data table and the spatial file.
      b.df.out <- rbind( b.df.out,b.r2.df.out) # add to the roads summary table
      roads.union = st_union(roads.union,b.r2); rm(b.r2);  plot(st_geometry(roads.union))
      #roads.union = st_cast(roads.union,"POLYGON")
      roads.union = st_union(roads.union)
     # roads.union = st_cast(roads.union,"POLYGON")

# LSR
      b.r2.sf = sf::st_read(dsn = Base , layer ="Tw_LSR_Rd_bd"  )
      b.r2.sf = st_make_valid(b.r2.sf)
      b.r2.sf <- st_buffer( b.r2.sf ,250)
      #b.r2.sf <- st_cast(b.r2.sf,"POLYGON")
      b.r2.sf <- st_union(b.r2.sf)
      b.r2 = st_intersection(all.range,b.r2.sf)
      b.r2$Area.m <- as.numeric(st_area(b.r2))
      b.r2.df = as.data.frame(b.r2)
      b.r2.df.out  = b.r2.df %>%
        group_by(SiteName,V17_CH) %>%
        filter(V17_CH == "Low Elevation Summer Range") %>%
        summarise(R_Road_area_m = sum(Area.m))
      # add to the Telkwa data table and the spatial file.
      b.df.out <- rbind( b.df.out,b.r2.df.out) # add to the roads summary table
      roads.union = st_union( roads.union,b.r2); rm(b.r2);  plot(st_geometry(roads.union))# takes some time
      #roads.union = st_cast(roads.union,"POLYGON")
      roads.union = st_union(roads.union)

# Matrix
      b.r2.sf = sf::st_read(dsn = Base , layer ="Tw_MR_Rd_bd"  )
      b.r2.sf = st_make_valid(b.r2.sf)
      b.r2.sf <- st_buffer( b.r2.sf ,250)
      #b.r2.sf <- st_zm( b.r2.sf ,drop = TRUE)
      #b.r2.sf <- st_cast(b.r2.sf,"POLYGON")
      b.r2.sf <- st_union(b.r2.sf)
      b.r2 = st_intersection(all.range,b.r2.sf)
      #b.r2 <- st_cast(b.r2,"POLYGON")
      b.r2$Area.m <- as.numeric(st_area(b.r2))
      b.r2.u = st_union(b.r2)
      b.r2.df = as.data.frame(b.r2)
      b.r2.df.out  = b.r2.df %>%
        group_by(SiteName,V17_CH) %>%
        filter(V17_CH == "Matrix Range") %>%
        summarise(R_Road_area_m = sum(Area.m))

      # add to the Telkwa data table and the spatial file.
      b.df.out <- rbind(b.df.out,b.r2.df.out) # add to the roads summary table

      roads.union = st_union(roads.union)
      roads.union = st_union(roads.union,b.r2.u); rm(b.r2.u);  plot(st_geometry(roads.union))
      #roads.union = st_cast(roads.union,"POLYGON")
      roads.union = st_union(roads.union)
      st_write(roads.union,paste(shape.output.dir,"Dist_road_buf.shp",sep = ""))


      ##Join the Roads layers back to the "combined static disturbance
      # combine into disturbance by layer
      all.range.out <- left_join(all.range.out,b.df.out)
      all.range.out[is.na(all.range.out)]<-0

      ## ALL DISTURBANCE: UNION 11 # may need to run this in stand alone R rather than R -studio
      #out10 = st_cast(out10,"POLYGON")
      out11 = st_union(out10,roads.union ) ; plot(st_geometry(out11))
      out11 = st_union(out11)
      rm(out10); rm(b.r1.u); rm(roads.union); rm(b.r2.sf)
      #plot(st_geometry(out11))


#######################################################################

# end of part 1: write out all the static disturbance types;

# write out all the layers combined into a single disturbance layer ( note this excludes Roads and Seismic )
st_write(out11,paste(shape.output.dir,"Static_disturb_TT_buf.shp",sep = "")) # this writes out as single layer

# Write out the datasheet onto a temp file
write.csv(all.range.out,paste(temp.dir,"Static_dist_TT_buf.csv",sep = "") )

all.range.out = read.csv(paste(temp.dir,"Static_dist_TT_buf.csv",sep = ""), header = TRUE)


################################################################################
################################################################################

## PART TWO : TEMPORAL DISTURBANCE

# Cutblocks,
# Burns
# Pests

##################################################################################

## 8) Cutblock ## this is all years of consolidated cutblock layer and Blairs error checked cutblocks

# HERD 1)  ## Telkwa
b.r.c= st_read(dsn = Base , layer = "cutblock_union_Te")      # read in file
    b.r.c <- st_zm(b.r.c ,drop = TRUE)                            # drop the z co-ordinate
    b.r.c <- st_buffer(b.r.c ,250)
    # get both values of HARVEST YEAR (as this this two cutblock layers unioned together - we need to get a single value for harvest year )
    b.r.c$HARVEST_YEAR_ALL = ifelse(b.r.c$HARVEST_YEAR > 0,b.r.c$HARVEST_YEAR,b.r.c$HARVEST_YEAR_1 ) #sort(unique(b.r.c$HARVEST_YEAR_ALL)) # error check
    b.r.c$TimeSinceCut = 2018-b.r.c$HARVEST_YEAR_ALL; # create new column with age since cut convert this to years
    #sort(unique(b.r.c$TimeSinceCut)) # check the range in years since cut #note in the boreal the oldest age cut is 78 years

    # Calculate the disturbance for 0-80 and 0-40 years
    # cutblocks 0-80 years
    b.r.c0.80 = b.r.c[b.r.c$TimeSinceCut < 81,] ; unique(b.r.c0.80$TimeSinceCut); #plot(b.r.c0.80$Shape) # subset all cutblocks from 0 - 40 years old
    b.r.c0.80 = st_union(b.r.c0.80)     # union (equivalent to dissolve)
    b.r.c0.80 = st_cast(b.r.c0.80,"POLYGON"); # check it is polygon and geometry is complete
    b.r.c0.80 = st_make_valid( b.r.c0.80)     # fix geometry if not complete
    b.r.c0.80 = st_intersection(all.range, b.r.c0.80) # intersect with herd boundaries
    b.r.c0.80$area.m = as.numeric(st_area(b.r.c0.80)) # calculate the area

      b.r.c0.80.df <- as.data.frame(b.r.c0.80 ) # create a new object from the "atribute" table of the spatial data
      r.cut.df.out  = b.r.c0.80.df%>%  # create a summary table of the area for each herd (sitename) and each habitat/ boundary type
        group_by(SiteName,V17_CH ) %>%
        summarise(R_cut0_80_m2 = sum(area.m))

    # cutblocks 0-40 years # repeat the process above for the 0 - 80 years cutblock
    b.r.c0.40 = b.r.c[b.r.c$TimeSinceCut < 41,] ; unique(b.r.c0.40$TimeSinceCut); #plot(st_geometry(b.r.c0.40))
    b.r.c0.40 = st_union(b.r.c0.40)
    b.r.c0.40 = st_cast(b.r.c0.40,"POLYGON")
    b.r.c0.40 = st_make_valid(b.r.c0.40) #; st_is_valid(b.r.c0.40)
    b.r.c0.40 = st_intersection(all.range, b.r.c0.40)
    b.r.c0.40$area.m = as.numeric(st_area(b.r.c0.40))

    b.r.c0.40.df <- as.data.frame(b.r.c0.40 )
    r.cut.df.out40  = b.r.c0.40.df%>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_cut0_40_m2 = sum(area.m))

    # calculate the disturbance for telkwa (0-40 and 0 - 80 in table form - use this later to add to Tweeds herds)
    r.cut.out = merge(r.cut.df.out,r.cut.df.out40  ) # merge the two summary tables together

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

    # generate cumulative burn disturbance shapefiles to be added sequentially to "static disturbance"
    head(r.cut.df) ; unique(r.cut.df$dec.period)
    Cut.dec.1950 <- r.cut.df %>% filter(dec.period == 1950) # empty
    Cut.dec.1960 <- r.cut.df %>% filter(dec.period < 1961 )
    Cut.dec.1970 <-r.cut.df %>% filter(dec.period < 1971 )
    Cut.dec.1980 <-r.cut.df %>% filter(dec.period < 1981 )
    Cut.dec.1990 <- r.cut.df %>% filter(dec.period < 1991 )
    Cut.dec.2000 <- r.cut.df%>% filter(dec.period < 2001 )
    Cut.dec.2010 <- r.cut.df %>% filter(dec.period < 2011 )

    ## write out the shapefiles to Data.drive
    st_write(Cut.dec.1950,paste(shape.output.dir,"Cut.te.dec.1950.buf.shp",sep = "")) # this writes out as single layer
    st_write(Cut.dec.1960,paste(shape.output.dir,"Cut.te.dec.1960.buf.shp",sep = "")) # this writes out as single layer
    st_write(Cut.dec.1970,paste(shape.output.dir,"Cut.te.dec.1970.buf.shp",sep = "")) # this writes out as single layer
    st_write(Cut.dec.1980,paste(shape.output.dir,"Cut.te.dec.1980.buf.shp",sep = "")) # this writes out as single layer
    st_write(Cut.dec.1990,paste(shape.output.dir,"Cut.te.dec.1990.buf.shp",sep = "")) # this writes out as single layer
    st_write(Cut.dec.2000,paste(shape.output.dir,"Cut.te.dec.2000.buf.shp",sep = "")) # this writes out as single layer
    st_write(Cut.dec.2010,paste(shape.output.dir,"Cut.te.dec.2010.buf.shp",sep = "")) # this writes out as single layer


          # generate decade burn disturbance shapefiles to be added sequentially to "static Disturbance"
          #head(r.cut.df2) ; unique(r.cut.df2$dec.period)
          Cut.d1950 <- r.cut.df %>% filter(dec.period == 1950)
          Cut.d1960<- r.cut.df %>% filter(dec.period == 1960 )
          Cut.d1970<-r.cut.df %>% filter(dec.period== 1970 )
          Cut.d1980<-r.cut.df %>% filter(dec.period == 1980 )
          Cut.d1990 <- r.cut.df %>% filter(dec.period == 1990 )
          Cut.d2000 <- r.cut.df%>% filter(dec.period == 2000 )
          Cut.d2010<- r.cut.df %>% filter(dec.period == 2010 )

          st_write(Cut.d1960,paste(shape.output.dir,"Cut.te.d1960.buf.shp",sep = "")) # this writes out as single layer
          st_write(Cut.d1970,paste(shape.output.dir,"Cut.te.d1970.buf.shp",sep = "")) # this writes out as single layer
          st_write(Cut.d1980,paste(shape.output.dir,"Cut.te.d1980.buf.shp",sep = "")) # this writes out as single layer
          st_write(Cut.d1990,paste(shape.output.dir,"Cut.te.d1990.buf.shp",sep = "")) # this writes out as single layer
          st_write(Cut.d2000,paste(shape.output.dir,"Cut.te.d2000.buf.shp",sep = "")) # this writes out as single layer
          st_write(Cut.d2010,paste(shape.output.dir,"Cut.te.d2010.buf.shp",sep = "")) # this writes out as single layer


    # generate consolidated output summary tables for each decade (same as the 0-40 and 0-80 as above) overtime per decage
    c.1950 <- st_union(Cut.dec.1950) # empty

    c.1960 <- st_union(Cut.dec.1960)
    c.1960 <- st_cast(c.1960,"POLYGON") #; st_is_valid(c.1960)
    c.1960 = st_intersection(all.range, c.1960) ; c.1960 <- st_make_valid(c.1960)
    #c.1960 <- st_cast(c.1960,"POLYGON")
    c.1960$area.m = as.numeric(st_area(c.1960))
    c.1960.df <- as.data.frame(c.1960)
    c.1960.df.out <-  c.1960.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_cut_1960_m2 = sum(area.m))
    r.cut.out = merge(r.cut.out,c.1960.df.out) # add to the data summary

    c.19601 <- c.1960

    c.1970 <- st_union(Cut.dec.1970)
    c.1970 <- st_cast(c.1970,"POLYGON") #; st_is_valid(c.1960)
    c.1970 = st_intersection(all.range, c.1970) ; c.1970 <- st_make_valid(c.1970)
    #c.1970 <- st_cast(c.1970,"POLYGON")
    c.1970$area.m = as.numeric(st_area(c.1970))
    c.1970.df <- as.data.frame(c.1970)
    c.1970.df.out <-  c.1970.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_cut_1970_m2 = sum(area.m))
    r.cut.out = merge(r.cut.out,c.1970.df.out) # add to the data summary
    c.19701 <- c.1970

    c.1980 <- st_union(Cut.dec.1980)
    c.1980 <- st_cast(c.1980,"POLYGON") #; st_is_valid(c.1960)
    c.1980 = st_intersection(all.range, c.1980) ; c.1980 <- st_make_valid(c.1980)
    #c.1980 <- st_cast(c.1980,"POLYGON")
    c.1980$area.m = as.numeric(st_area(c.1980))
    c.1980.df <- as.data.frame(c.1980)
    c.1980.df.out <-  c.1980.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_cut_1980_m2 = sum(area.m))
    r.cut.out = merge(r.cut.out,c.1980.df.out) # add to the data summary
    c.19801 <- c.1980

    c.1990 <- st_union(Cut.dec.1990)
    c.1990 <- st_cast(c.1990,"POLYGON") #; st_is_valid(c.1960)
    c.1990 = st_intersection(all.range, c.1990) ; c.1990 <- st_make_valid(c.1990)
    #c.1990 <- st_cast(c.1990,"POLYGON")
    c.1990$area.m = as.numeric(st_area(c.1990))
    c.1990.df <- as.data.frame(c.1990)
    c.1990.df.out <-  c.1990.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_cut_1990_m2 = sum(area.m))
    r.cut.out = merge(r.cut.out,c.1990.df.out) # add to the data summary
    c.19901 <- c.1990

    c.2000 <- st_union(Cut.dec.2000)
    c.2000 <- st_cast(c.2000,"POLYGON") #; st_is_valid(c.2000)
    c.2000 = st_intersection(all.range, c.2000) ; c.2000 <- st_make_valid(c.2000)
    #c.2000 <- st_cast(c.2000,"POLYGON")
    c.2000$area.m = as.numeric(st_area(c.2000))
    c.2000.df <- as.data.frame(c.2000)
    c.2000.df.out <-  c.2000.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_cut_2000_m2 = sum(area.m))
    r.cut.out = merge(r.cut.out,c.2000.df.out) # add to the data summary
    c.20001 <- c.2000

    c.2010 <- st_union(Cut.dec.2010)
    c.2010 <- st_cast(c.2010,"POLYGON") ; st_is_valid(c.2010); c.2010 <- st_make_valid(c.2010)
    c.2010 = st_intersection(all.range, c.2010)
    #c.2010 <- st_cast(c.2010,"POLYGON")
    c.2010$area.m = as.numeric(st_area(c.2010))
    c.2010.df <- as.data.frame(c.2010)
    c.2010.df.out <-  c.2010.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_cut_2010_m2 = sum(area.m))
    r.cut.out = merge(r.cut.out,c.2010.df.out) # add to the data summary
    r.cut.out.te = r.cut.out
    c.20101 <- c.2010

    # outputs   = data table with all 80-0 and 40-0 + each decade (cumulative)
    #           = spatial data at each decade.

# HERD 2)  ## Tweedsmuir

# cutblock data from Canfor 2016 - 2018 data set
#b.r.c1 = st_read(paste("C:\\Temp\\TweedTelkwa\\Temp\\Perkins\\Data\\whitesail_shapefiles_2018_12_10","Whitesail_harvested_blocks_cfp_2018_12_11.shp",sep = "\\"))
b.r.c1 = st_read(paste(temp.dir,"\\whitesail_shapefiles_2018_12_10\\Whitesail_harvested_blocks_cfp_2018_12_11.shp",sep = ""))

    b.r.c1  = st_transform(b.r.c1 ,3005)
    b.r.c1 = st_buffer(b.r.c1,250)
    #unique(b.r.c1$HS_DATE)
    b.r.c1  = st_union(b.r.c1) #;plot(st_geometry(b.r.c1))
    b.r.c1<- st_intersection(all.range,b.r.c1) ; plot(st_geometry(b.r.c1 ))
    #b.r.c1 = st_cast(b.r.c1,"POLYGON")
    b.r.c1  = st_union(b.r.c1)

########################
b.r.c2= st_read(dsn = Base , layer = "cutblock_union_Tw")
    b.r.c2 <- st_zm(b.r.c2 ,drop = TRUE)
    b.r.c2 = st_buffer(b.r.c2,250)
    head(b.r.c2)
    # get both values of HARVEST YEAR
    b.r.c2$HARVEST_YEAR_ALL = ifelse(b.r.c2$HARVEST_YEAR > 0,b.r.c2$HARVEST_YEAR,b.r.c2$HARVEST_YEAR_1 ) #sort(unique(b.r.c$HARVEST_YEAR_ALL)) # error check b.r.c2$TimeSinceCut = 2018-b.r.c2$HARVEST_YEAR_ALL; # create new column with age since cut
    b.r.c2$TimeSinceCut = 2018-b.r.c2$HARVEST_YEAR_ALL; # create new column with age since cut
    #sort(unique(b.r.c$TimeSinceCut)) # check the range in years since cut #note in the boreal the oldest age cut is 78 years
    b.r.c2 <- st_cast(b.r.c2,"POLYGON")
    b.r.c2 <- st_make_valid(b.r.c2)

    # cutblocks 0-80 years
    b.r.c0.802 = b.r.c2[b.r.c2$TimeSinceCut < 81,] ; unique(b.r.c0.802$TimeSinceCut); #plot(b.r.c0.80$Shape)
    b.r.c0.802 = st_make_valid( b.r.c0.802)
    b.r.c0.802 = st_union(b.r.c0.802)
    b.r.c0.802 = st_union(b.r.c0.802,b.r.c1)     # add the white sale additional data
    #b.r.c0.802 = st_cast(b.r.c0.802,"POLYGON"); #
    b.r.c0.802 = st_make_valid( b.r.c0.802)
    b.r.c0.802 = st_intersection(all.range, b.r.c0.802)
    b.r.c0.802$area.m = as.numeric(st_area(b.r.c0.802))
    b.r.c0.802.df <- as.data.frame(b.r.c0.802 )
    r.cut.df2.out  = b.r.c0.802.df%>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_cut0_80_m2 = sum(area.m))

    # cutblocks 0-40 years
    b.r.c0.402 = b.r.c2[b.r.c2$TimeSinceCut < 41,] ; unique(b.r.c0.402$TimeSinceCut); #plot(b.r.c0.40$Shape)
    b.r.c0.402 = st_union(b.r.c0.402)
    b.r.c1 <- st_union(b.r.c1)
    b.r.c0.402 = st_union(b.r.c0.402,b.r.c1) # add the white sale additional data
    b.r.c0.402 = st_cast(b.r.c0.402,"POLYGON"); #
    b.r.c0.402 = st_make_valid( b.r.c0.402)
    b.r.c0.402 = st_intersection(all.range, b.r.c0.402)
    b.r.c0.402$area.m = as.numeric(st_area(b.r.c0.402))
    b.r.c0.402.df = as.data.frame(b.r.c0.402)

    r.cut.df2.out40  = b.r.c0.402.df%>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_cut0_40_m2 = sum(area.m))

    # calculate the disturbance for telkwa (0-40 and 0 - 80 in table form - use this later to add to Tweeds herds)
    r.cut.out2 = merge(r.cut.df2.out,r.cut.df2.out40)

    r.cut.df2 = b.r.c2

    # add a column to differentiate the age brackets of cutblocks
    r.cut.df2 <- mutate(r.cut.df2,dec.period = ifelse(HARVEST_YEAR_ALL >= 1940 & HARVEST_YEAR_ALL <= 1949,1940,0))
    r.cut.df2 <- mutate(r.cut.df2,dec.period = ifelse(HARVEST_YEAR_ALL >= 1950 & HARVEST_YEAR_ALL <= 1959,1950,dec.period))
    r.cut.df2 <- mutate(r.cut.df2,dec.period = ifelse(HARVEST_YEAR_ALL >= 1960 & HARVEST_YEAR_ALL <= 1969,1960,dec.period))
    r.cut.df2 <- mutate(r.cut.df2,dec.period = ifelse(HARVEST_YEAR_ALL >= 1970 & HARVEST_YEAR_ALL <= 1979,1970,dec.period))
    r.cut.df2 <- mutate(r.cut.df2,dec.period = ifelse(HARVEST_YEAR_ALL >= 1980 & HARVEST_YEAR_ALL <= 1989,1980,dec.period))
    r.cut.df2 <- mutate(r.cut.df2,dec.period = ifelse(HARVEST_YEAR_ALL >= 1990 & HARVEST_YEAR_ALL <= 1999,1990,dec.period))
    r.cut.df2 <- mutate(r.cut.df2,dec.period = ifelse(HARVEST_YEAR_ALL >= 2000 & HARVEST_YEAR_ALL <= 2009,2000,dec.period))
    r.cut.df2 <- mutate(r.cut.df2,dec.period = ifelse(HARVEST_YEAR_ALL >= 2010 & HARVEST_YEAR_ALL <= 2019,2010,dec.period))
    #r.cut.df[r.cut.df$dec.period == 0,] ; unique(r.cut.df$dec.period)  # error check

    # generate cumulative burn disturbance shapefiles to be added sequentially to "static Disturbance"
    head(r.cut.df2) ; unique(r.cut.df2$dec.period)
    Cut.dec.19502 <- r.cut.df2 %>% filter(dec.period == 1950)
    Cut.dec.19602 <- r.cut.df2 %>% filter(dec.period < 1961 )
    Cut.dec.19702 <-r.cut.df2 %>% filter(dec.period < 1971 )
    Cut.dec.19802 <-r.cut.df2 %>% filter(dec.period < 1981 )
    Cut.dec.19902 <- r.cut.df2 %>% filter(dec.period < 1991 )
    Cut.dec.20002 <- r.cut.df2%>% filter(dec.period < 2001 )
    Cut.dec.20102 <- r.cut.df2 %>% filter(dec.period < 2011 )
    Cut.dec.20102 <- st_union(Cut.dec.20102 ,b.r.c1)

    ## write out the shapefiles to Data.drive
    st_write(Cut.dec.19502,paste(shape.output.dir,"Cut.tw.dec.1950.buf.shp",sep = "")) # this writes out as single layer
    st_write(Cut.dec.19602,paste(shape.output.dir,"Cut.tw.dec.1960.buf.shp",sep = "")) # this writes out as single layer
    st_write(Cut.dec.19702,paste(shape.output.dir,"Cut.tw.dec.1970.buf.shp",sep = "")) # this writes out as single layer
    st_write(Cut.dec.19802,paste(shape.output.dir,"Cut.tw.dec.1980.buf.shp",sep = "")) # this writes out as single layer
    st_write(Cut.dec.19902,paste(shape.output.dir,"Cut.tw.dec.1990.buf.shp",sep = "")) # this writes out as single layer
    st_write(Cut.dec.20002,paste(shape.output.dir,"Cut.tw.dec.2000.buf.shp",sep = "")) # this writes out as single layer
    st_write(Cut.dec.20102,paste(shape.output.dir,"Cut.tw.dec.2010.buf.shp",sep = "")) # this writes out as single layer

      # generate decade burn disturbance shapefiles to be added sequentially to "static Disturbance"
          head(r.cut.df2) ; unique(r.cut.df2$dec.period)
          Cut.d19502 <- r.cut.df2 %>% filter(dec.period == 1950)
          Cut.d19602 <- r.cut.df2 %>% filter(dec.period == 1960 )
          Cut.d19702 <-r.cut.df2 %>% filter(dec.period== 1970 )
          Cut.d19802 <-r.cut.df2 %>% filter(dec.period == 1980 )
          Cut.d19902 <- r.cut.df2 %>% filter(dec.period == 1990 )
          Cut.d20002 <- r.cut.df2%>% filter(dec.period == 2000 )
          Cut.d20102 <- r.cut.df2 %>% filter(dec.period == 2010 )
          ## write out the shapefiles to Data.drive
          st_write(Cut.d19502,paste(shape.output.dir,"Cut.tw.d1950.buf.shp",sep = "")) # this writes out as single layer
          st_write(Cut.d19602,paste(shape.output.dir,"Cut.tw.d1960.buf.shp",sep = "")) # this writes out as single layer
          st_write(Cut.d19702,paste(shape.output.dir,"Cut.tw.d1970.buf.shp",sep = "")) # this writes out as single layer
          st_write(Cut.d19802,paste(shape.output.dir,"Cut.tw.d1980.buf.shp",sep = "")) # this writes out as single layer
          st_write(Cut.d19902,paste(shape.output.dir,"Cut.tw.d1990.buf.shp",sep = "")) # this writes out as single layer
          st_write(Cut.d20002,paste(shape.output.dir,"Cut.tw.d2000.buf.shp",sep = "")) # this writes out as single layer
          st_write(Cut.d20102,paste(shape.output.dir,"Cut.tw.d2010.buf.shp",sep = "")) # this writes out as single layer


    ## add the section for each decade (cumulative )
    # generate consolidated outputs overtime per decage
    c.1950 <- st_union(Cut.dec.19502) # empty

    c.1960 <- st_union(Cut.dec.19602)
    c.1960 <- st_cast(c.1960,"POLYGON") #; st_is_valid(c.1960)
    c.1960 = st_intersection(all.range, c.1960) ; c.1960 <- st_make_valid(c.1960)
    #c.1960 <- st_cast(c.1960,"POLYGON")
    c.1960$area.m = as.numeric(st_area(c.1960))
    c.1960.df <- as.data.frame(c.1960)
    c.1960.df.out <-  c.1960.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_cut_1960_m2 = sum(area.m))
    r.cut.out = merge(r.cut.out2,c.1960.df.out,all.x= T) # add to the data summary

    c.1970 <- st_union(Cut.dec.19702)
    c.1970 <- st_cast(c.1970,"POLYGON") #; st_is_valid(c.1960)
    c.1970 = st_intersection(all.range, c.1970) ; c.1970 <- st_make_valid(c.1970)
    #c.1970 <- st_cast(c.1970,"POLYGON")
    c.1970$area.m = as.numeric(st_area(c.1970))
    c.1970.df <- as.data.frame(c.1970)
    c.1970.df.out <-  c.1970.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_cut_1970_m2 = sum(area.m))
    r.cut.out = merge(r.cut.out,c.1970.df.out) # add to the data summary

    c.1980 <- st_union(Cut.dec.19802)
    c.1980 <- st_cast(c.1980,"POLYGON") #; st_is_valid(c.1960)
    c.1980 = st_intersection(all.range, c.1980) ; c.1980 <- st_make_valid(c.1980)
    #c.1980 <- st_cast(c.1980,"POLYGON")
    c.1980$area.m = as.numeric(st_area(c.1980))
    c.1980.df <- as.data.frame(c.1980)
    c.1980.df.out <-  c.1980.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_cut_1980_m2 = sum(area.m))
    r.cut.out = merge(r.cut.out,c.1980.df.out) # add to the data summary

    c.1990 <- st_union(Cut.dec.19902)
    c.1990 <- st_cast(c.1990,"POLYGON") #; st_is_valid(c.1960)
    c.1990 = st_intersection(all.range, c.1990) ; c.1990 <- st_make_valid(c.1990)
    #c.1990 <- st_cast(c.1990,"POLYGON")
    c.1990$area.m = as.numeric(st_area(c.1990))
    c.1990.df <- as.data.frame(c.1990)
    c.1990.df.out <-  c.1990.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_cut_1990_m2 = sum(area.m))
    r.cut.out = merge(r.cut.out,c.1990.df.out) # add to the data summary

    c.2000 <- st_union(Cut.dec.20002)
    c.2000 <- st_cast(c.2000,"POLYGON") #; st_is_valid(c.2000)
    c.2000 = st_intersection(all.range, c.2000) ; c.2000 <- st_make_valid(c.2000)
    #c.2000 <- st_cast(c.2000,"POLYGON")
    c.2000$area.m = as.numeric(st_area(c.2000))
    c.2000.df <- as.data.frame(c.2000)
    c.2000.df.out <-  c.2000.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_cut_2000_m2 = sum(area.m))
    r.cut.out = merge(r.cut.out,c.2000.df.out) # add to the data summary

    c.2010 <- st_union(Cut.dec.20102)
    c.2010 <- st_cast(c.2010,"POLYGON") ; st_is_valid(c.2010); c.2010 <- st_make_valid(c.2010)
    c.2010 = st_intersection(all.range, c.2010)
    #c.2010 <- st_cast(c.2010,"POLYGON")
    c.2010$area.m = as.numeric(st_area(c.2010))
    c.2010.df <- as.data.frame(c.2010)
    c.2010.df.out <-  c.2010.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_cut_2010_m2 = sum(area.m))
    r.cut.out = merge(r.cut.out,c.2010.df.out) # add to the data summary

    r.cut.out.tw =  r.cut.out

    ###############
    # Join the telkwa and Tweedsmuir herd info together
    r.cut.out.all = rbind(r.cut.out.te,r.cut.out.tw)

    # combine into disturbance by layer
    all.range.out <- left_join(all.range.out,r.cut.out.all)
    all.range.out[is.na(all.range.out)]<-0

    # Write out the datasheet onto a temp file
    write.csv(all.range.out,paste(temp.dir,"Static_dist_TT_buf_temp.csv",sep = "") )
    all.range.out = read.csv(paste(temp.dir,"Static_dist_TT_buf_temp.csv",sep = "") )

########################################################################
## BURN:  NO BUffer for Burns section
# Add to Burn information Break down of burns into 0-40 and 0-80 years (and total)
# Break down burns into fire years (same as cutblocks)
# this generates the outputs for the 0-80 years, 0-40 years and for cumulative decades (1950 - 2010)

b.r.0 = st_read(dsn = Base, layer ="fire_clip")
    b.r.0<- st_zm(b.r.0 ,drop = TRUE) # this is a linear feature so need to buffer to estimate area calcs
    #st_is_valid(b.r.0)
    b.r.0$TimeSinceBurn = 2018-b.r.0$FIRE_YEAR #; plot(b.r.0$Shape)
    b.r.0 <- st_intersection(all.range,b.r.0) #; st_is_valid(b.r.0)
    #b.r.0 <- st_cast(b.r.0,"POLYGON")

    #plot(st_geometry(b.r.0), col = 'red')
    #plot(st_geometry(all.range),add = T)

    # burns 0-40 years
    b.r.0.40 = b.r.0[b.r.0$TimeSinceBurn <41,]; #sort(unique(b.r.0$TimeSinceBurn))
    b.r.0.40 = st_union(b.r.0.40)
    b.r.0.40 <- st_cast(b.r.0.40,"POLYGON")
    b.r.0.40 = st_make_valid( b.r.0.40 )
    b.r.0.40 = st_intersection(all.range,b.r.0.40)
    b.r.0.40$area.m = as.numeric(st_area( b.r.0.40 ))
    b.r.0.40.df = data.frame(b.r.0.40)        # calculate the length per range
    b.r.0.40.out  =  b.r.0.40.df%>%
      group_by(SiteName,V17_CH) %>%
      summarise(R_Burn040_m2 = sum(area.m))

    # burns 0-80 years
    b.r.0.80 =b.r.0[b.r.0$TimeSinceBurn <81,];
    b.r.0.80 = st_union(b.r.0.80)
    b.r.0.80 <- st_cast(b.r.0.80,"POLYGON")
    b.r.0.80 = st_make_valid( b.r.0.80 )
    b.r.0.80 = st_intersection(all.range,b.r.0.80)
    b.r.0.80$area.m = as.numeric(st_area(b.r.0.80))
    b.r.0.80.df = data.frame(b.r.0.80)        # calculate the length per range
    b.r.0.80.out  =  b.r.0.80.df%>%
      group_by(SiteName,V17_CH) %>%
      summarise(R_Burn080_m2 = sum(area.m))

    r.burn.out = merge(b.r.0.40.out,b.r.0.80.out )  # keep this to add to other layers
    r.burn.out = left_join(Herd_key_detail, r.burn.out) # add to the total of all range habitats
    # combine into disturbance by layer (from 0 - 80 years)

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

    # generate cumulative burn disturbance shapefiles to be added sequentially to "static Disturbance"
    Burn.dec = b.r.00.df

    Burn.dec.1950 <- Burn.dec %>% filter(dec.period == 1950)
    Burn.dec.1960 <- Burn.dec %>% filter(dec.period < 1961 )
    Burn.dec.1970 <- Burn.dec %>% filter(dec.period < 1971 )
    Burn.dec.1980 <- Burn.dec %>% filter(dec.period < 1981 )
    Burn.dec.1990 <- Burn.dec %>% filter(dec.period < 1991 )
    Burn.dec.2000 <- Burn.dec %>% filter(dec.period < 2001 )
    Burn.dec.2010 <- Burn.dec %>% filter(dec.period < 2011 )
    # Burn.dec.2010 <- Burn.dec %>% filter(dec.period < 2011 )

    ## write out the shapefiles to Data.drive
    st_write(Burn.dec.1950,paste(shape.output.dir,"Burn.dec.1950.buf.shp",sep = "")) # this writes out as single layer
    st_write(Burn.dec.1960,paste(shape.output.dir,"Burn.dec.1960.buf.shp",sep = "")) # this writes out as single layer
    st_write(Burn.dec.1970,paste(shape.output.dir,"Burn.dec.1970.buf.shp",sep = "")) # this writes out as single layer
    st_write(Burn.dec.1980,paste(shape.output.dir,"Burn.dec.1980.buf.shp",sep = "")) # this writes out as single layer
    st_write(Burn.dec.1990,paste(shape.output.dir,"Burn.dec.1990.buf.shp",sep = "")) # this writes out as single layer
    st_write(Burn.dec.2000,paste(shape.output.dir,"Burn.dec.2000.buf.shp",sep = "")) # this writes out as single layer
    st_write(Burn.dec.2010,paste(shape.output.dir,"Burn.dec.2010.buf.shp",sep = "")) # this writes out as single layer

                ## Create Burns per decade
                Burn.d1950 <- Burn.dec %>% filter(dec.period == 1950)
                Burn.d1960 <- Burn.dec %>% filter(dec.period == 1960 )
                Burn.d1970 <- Burn.dec %>% filter(dec.period == 1970 )
                Burn.d1980 <- Burn.dec %>% filter(dec.period == 1980 )
                Burn.d1990 <- Burn.dec %>% filter(dec.period == 1990 )
                Burn.d2000 <- Burn.dec %>% filter(dec.period == 2000 )
                Burn.d2010 <- Burn.dec %>% filter(dec.period == 2010 )

                ## write out the shapefiles to Data.drive
                st_write(Burn.d1950,paste(shape.output.dir,"Burn.d1950.buf.shp",sep = "")) # this writes out as single layer
                st_write(Burn.d1960,paste(shape.output.dir,"Burn.d1960.buf.shp",sep = "")) # this writes out as single layer
                st_write(Burn.d1970,paste(shape.output.dir,"Burn.d1970.buf.shp",sep = "")) # this writes out as single layer
                st_write(Burn.d1980,paste(shape.output.dir,"Burn.d1980.buf.shp",sep = "")) # this writes out as single layer
                st_write(Burn.d1990,paste(shape.output.dir,"Burn.d1990.buf.shp",sep = "")) # this writes out as single layer
                st_write(Burn.d2000,paste(shape.output.dir,"Burn.d2000.buf.shp",sep = "")) # this writes out as single layer
                st_write(Burn.d2010,paste(shape.output.dir,"Burn.d2010.buf.shp",sep = "")) # this writes out as single layer

    ## add the section for each decade (cumulative )
    # generate consolidated outputs overtime per decage
    b.1950 <- st_union(Burn.dec.1950) # empty
    b.1950 <- st_cast(b.1950,"POLYGON") #; st_is_valid(b.1950)
    b.1950 = st_intersection(all.range, b.1950) ; b.1950 <- st_make_valid(b.1950)
    #b.1950 <- st_cast(b.1950,"POLYGON")
    b.1950$area.m = as.numeric(st_area(b.1950))
    b.1950.df <- as.data.frame(b.1950)
    b.1950.df.out <-  b.1950.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_burn_1950_m2 = sum(area.m))
    r.burn.out.total = left_join(r.burn.out, b.1950.df.out) # add to the data summary

    b.1960 <- st_union(Burn.dec.1960)
    b.1960 <- st_cast(b.1960,"POLYGON") #; st_is_valid(c.1960)
    b.1960 = st_intersection(all.range, b.1960) ; b.1960 <- st_make_valid(b.1960)
    # b.1960 <- st_cast(b.1960,"POLYGON")
    b.1960$area.m = as.numeric(st_area(b.1960))
    b.1960.df <- as.data.frame(b.1960)
    b.1960.df.out <-  b.1960.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_burn_1960_m2 = sum(area.m))
    r.burn.out.total = left_join(r.burn.out.total,b.1960.df.out) # add to the data summary

    b.1970 <- st_union(Burn.dec.1970)
    b.1970 <- st_cast(b.1970,"POLYGON") #; st_is_valid(c.1960)
    b.1970 = st_intersection(all.range, b.1970) ; b.1970 <- st_make_valid(b.1970)
    #b.1970 <- st_cast(b.1970,"POLYGON")
    b.1970$area.m = as.numeric(st_area(b.1970))
    b.1970.df <- as.data.frame(b.1970)
    b.1970.df.out <-  b.1970.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_burn_1970_m2 = sum(area.m))
    r.burn.out.total = left_join(r.burn.out.total,b.1970.df.out) # add to the data summary

    b.1980 <- st_union(Burn.dec.1980)
    b.1980 <- st_cast(b.1980,"POLYGON") #; st_is_valid(c.1960)
    b.1980 = st_intersection(all.range, b.1980) #; b.1980 <- st_make_valid(b.1980)
    #b.1980 <- st_cast(b.1980,"POLYGON") #; st_is_valid(b.1980)
    b.1980$area.m = as.numeric(st_area(b.1980))
    b.1980.df <- as.data.frame(b.1980)
    b.1980.df.out <-  b.1980.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_burn_1980_m2 = sum(area.m))
    r.burn.out.total = left_join(r.burn.out.total,b.1980.df.out) # add to the data summary

    b.1990 <- st_union(Burn.dec.1990)
    b.1990 <- st_cast(b.1990,"POLYGON") #; st_is_valid(c.1960)
    b.1990 = st_intersection(all.range, b.1990) #; b.1990 <- st_make_valid(c.1990)
    #b.1990 <- st_cast(b.1990,"POLYGON")
    b.1990$area.m = as.numeric(st_area(b.1990))
    b.1990.df <- as.data.frame(b.1990)
    b.1990.df.out <-  b.1990.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_burn_1990_m2 = sum(area.m))
    r.burn.out.total = left_join(r.burn.out.total,b.1990.df.out) # add to the data summary

    b.2000 <- st_union(Burn.dec.2000)
    b.2000 <- st_cast(b.2000,"POLYGON") #; st_is_valid(c.2000)
    b.2000 = st_intersection(all.range, b.2000) ; b.2000 <- st_make_valid(b.2000)
    #b.2000 <- st_cast(c.2000,"POLYGON")
    b.2000$area.m = as.numeric(st_area(b.2000))
    b.2000.df <- as.data.frame(b.2000)
    b.2000.df.out <-  b.2000.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_burn_2000_m2 = sum(area.m))
    r.burn.out.total = left_join(r.burn.out.total,b.2000.df.out) # add to the data summary

    b.2010 <- st_union(Burn.dec.2010)
    b.2010 <- st_cast(b.2010,"POLYGON") #; st_is_valid(b.2010); b.2010 <- st_make_valid(b.2010)
    b.2010 = st_intersection(all.range, b.2010)
    #head(b.2010)
    #b.2010 <- st_cast(b.2010,"POLYGON")
    b.2010$area.m = as.numeric(st_area(b.2010))
    b.2010.df <- as.data.frame(b.2010)
    b.2010.df.out <-  b.2010.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_burn_2010_m2 = sum(area.m))
    r.burn.out.total = left_join(r.burn.out.total,b.2010.df.out) # add to the data summary


    ###############

    # combine into disturbance by layer
    all.range.out <- left_join(all.range.out, r.burn.out.total)
    all.range.out[is.na(all.range.out)]<-0

    # Write out the datasheet onto a temp file
    write.csv(all.range.out,paste(temp.dir,"Static_dist_TT_buf_temp.csv",sep = "") )


############################################
### PEST
# this needs pre-processing in arcmap select ibm and ibs species and clip to boundary
# split into two sections as this is too large

r.pest.te <-  st_read(dsn = Base, layer ="Pest_clip_Te_IBMIBS") #; plot(st_geometry(r.pest.te))

    # Telkwa
    r.pest <- r.pest.te
    r.pest<- st_zm(r.pest ,drop = TRUE) # this is a linear feature so need to buffer to estimate area calcs
    r.pest = st_buffer(r.pest,250)
    r.pest <- st_make_valid(r.pest)
    r.pest$TimeSincePest = 2018-r.pest$CAPTURE_YEAR
    r.pest<- st_intersection(all.range,r.pest) #; st_is_valid(b.r.0)
    r.pest <- st_cast(r.pest,"POLYGON")
    #plot(st_geometry(r.pest), col = 'red')
    #plot(st_geometry(all.range),add = T)

    # pest 0-40 years
    r.p.0.40 = r.pest[r.pest$TimeSincePest <41,]; #sort(unique(b.r.0$TimeSinceBurn))
    r.p.0.40 = st_union(r.p.0.40)
    r.p.0.40 <- st_cast(r.p.0.40,"POLYGON") #; st_is_valid(r.p.0.40)
    #r.p.0.40 = st_make_valid( r.p.0.40)
    r.p.0.40= st_intersection(all.range,r.p.0.40)
    r.p.0.40$area.m = as.numeric(st_area(r.p.0.40))
    r.p.0.40.df = data.frame(r.p.0.40)        # calculate the length per range
    r.p.0.40.out  =  r.p.0.40.df%>%
      group_by(SiteName,V17_CH) %>%
      summarise(R_Pest040_m2 = sum(area.m))

    r.pest.df = r.pest

    # add decade data   # add a column to differentiate the age brackets of pest capture
    r.pest.df<- r.pest.df %>% mutate(dec.period = ifelse(CAPTURE_YEAR >= 1950 & CAPTURE_YEAR <= 1965,1950,0))
    r.pest.df<- mutate(r.pest.df,dec.period = ifelse(CAPTURE_YEAR >= 1960 & CAPTURE_YEAR <= 1969,1960,dec.period))
    r.pest.df<- mutate(r.pest.df,dec.period = ifelse(CAPTURE_YEAR >= 1970 & CAPTURE_YEAR <= 1979,1970,dec.period))
    r.pest.df<- mutate(r.pest.df,dec.period = ifelse(CAPTURE_YEAR >= 1980 & CAPTURE_YEAR <= 1989,1980,dec.period))
    r.pest.df<- mutate(r.pest.df,dec.period = ifelse(CAPTURE_YEAR >= 1990 & CAPTURE_YEAR <= 1999,1990,dec.period))
    r.pest.df<- mutate(r.pest.df,dec.period = ifelse(CAPTURE_YEAR >= 2000 & CAPTURE_YEAR <= 2009,2000,dec.period))
    r.pest.df<- mutate(r.pest.df,dec.period = ifelse(CAPTURE_YEAR >= 2010 & CAPTURE_YEAR <= 2018,2010,dec.period))

    # Generate the cumulative pest damage into decades to add to "static Disturbance"
    Pest.dec.1950 <- r.pest.df %>% filter(dec.period == 1950)
    Pest.dec.1960 <- r.pest.df %>% filter(dec.period < 1961 ); Pest.dec.1960 <-st_union(Pest.dec.1960 )
    Pest.dec.1970 <- r.pest.df %>% filter(dec.period < 1971 ); Pest.dec.1970 <-st_union(Pest.dec.1970 )
    Pest.dec.1980 <- r.pest.df %>% filter(dec.period < 1981 ); Pest.dec.1980 <-st_union(Pest.dec.1980 )
    Pest.dec.1990 <- r.pest.df %>% filter(dec.period < 1991 ); Pest.dec.1990 <-st_union(Pest.dec.1990 )
    Pest.dec.2000 <- r.pest.df %>% filter(dec.period < 2001 ); Pest.dec.2000 <-st_union(Pest.dec.2000 )
    Pest.dec.2010 <- r.pest.df %>% filter(dec.period < 2011 ); Pest.dec.2010 <-st_union(Pest.dec.2010 )
    # Burn.dec.2010 <- Burn.dec %>% filter(dec.period < 2011 )

    ## output shapefiles.
    ## write out the shapefiles to Data.drive
    st_write(Pest.dec.1950 ,paste(shape.output.dir,"Pest.te.dec.1950.buf.shp",sep = "")) # this writes out as single layer
    st_write(Pest.dec.1960 ,paste(shape.output.dir,"Pest.te.dec.1960.buf.shp",sep = "")) # this writes out as single layer
    st_write(Pest.dec.1970 ,paste(shape.output.dir,"Pest.te.dec.1970.buf.shp",sep = "")) # this writes out as single layer
    st_write(Pest.dec.1980 ,paste(shape.output.dir,"Pest.te.dec.1980.buf.shp",sep = "")) # this writes out as single layer
    st_write(Pest.dec.1990 ,paste(shape.output.dir,"Pest.te.dec.1990.buf.shp",sep = "")) # this writes out as single layer
    st_write(Pest.dec.2000 ,paste(shape.output.dir,"Pest.te.dec.2000.buf.shp",sep = "")) # this writes out as single layer
    st_write(Pest.dec.2010 ,paste(shape.output.dir,"Pest.te.dec.2010.buf.shp",sep = "")) # this writes out as single layer

    ## add the section for each decade (cumulative )
    # generate consolidated outputs overtime per decage
    #p.1950 <- st_union(Pest.dec.1950) # empty

    p.1960 <- st_union(Pest.dec.1960)
    p.1960 <- st_cast(p.1960,"POLYGON") #; st_is_valid(c.1960)
    p.1960 = st_intersection(all.range, p.1960) ; p.1960 <- st_make_valid(p.1960)
    p.1960$area.m = as.numeric(st_area(p.1960))
    p.1960.df <- as.data.frame(p.1960)
    p.1960.df.out <-  p.1960.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_pest_1960_m2 = sum(area.m))
    r.pest.out.total = left_join(r.p.0.40.out,p.1960.df.out) # add to the data summary

    p.1970 <- st_union(Pest.dec.1970)
    p.1970 <- st_cast(p.1970,"POLYGON") #; st_is_valid(c.1960)
    p.1970 = st_intersection(all.range, p.1970) ; p.1970 <- st_make_valid(p.1970)
    #b.1970 <- st_cast(b.1970,"POLYGON")
    p.1970$area.m = as.numeric(st_area(p.1970))
    p.1970.df <- as.data.frame(p.1970)
    p.1970.df.out <-  p.1970.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_pest_1970_m2 = sum(area.m))
    r.pest.out.total = left_join( r.pest.out.total,p.1970.df.out) # add to the data summary

    p.1980 <- st_union(Pest.dec.1980)
    p.1980 <- st_cast(p.1980,"POLYGON") #; st_is_valid(c.1960)
    p.1980 = st_intersection(all.range, p.1980) #; b.1980 <- st_make_valid(b.1980)
    #b.1980 <- st_cast(b.1980,"POLYGON") #; st_is_valid(b.1980)
    p.1980$area.m = as.numeric(st_area(p.1980))
    p.1980.df <- as.data.frame(p.1980)
    p.1980.df.out <-  p.1980.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_pest_1980_m2 = sum(area.m))
    r.pest.out.total = left_join(r.pest.out.total,p.1980.df.out) # add to the data summary

    p.1990 <- st_union(Pest.dec.1990)
    p.1990 <- st_cast(p.1990,"POLYGON") #; st_is_valid(c.1960)
    p.1990 = st_intersection(all.range, p.1990) #; b.1990 <- st_make_valid(c.1990)
    #p.1990 <- st_cast(b.1990,"POLYGON")
    p.1990$area.m = as.numeric(st_area(p.1990))
    p.1990.df <- as.data.frame(p.1990)
    p.1990.df.out <-  p.1990.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_pest_1990_m2 = sum(area.m))
    r.pest.out.total = left_join(r.pest.out.total,p.1990.df.out) # add to the data summary

    p.2000 <- st_union(Pest.dec.2000)
    p.2000 <- st_cast(p.2000,"POLYGON") #; st_is_valid(c.2000)
    p.2000 = st_intersection(all.range, p.2000) ; p.2000 <- st_make_valid(p.2000)
    #p.2000 <- st_cast(c.2000,"POLYGON")
    p.2000$area.m = as.numeric(st_area(p.2000))
    p.2000.df <- as.data.frame(p.2000)
    p.2000.df.out <-  p.2000.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_pest_2000_m2 = sum(area.m))
    r.pest.out.total = left_join(r.pest.out.total,p.2000.df.out) # add to the data summary

    p.2010 <- st_union(Pest.dec.2010)
    p.2010 <- st_cast(p.2010,"POLYGON") #; st_is_valid(b.2010); b.2010 <- st_make_valid(b.2010)
    p.2010 = st_intersection(all.range, p.2010)
    #head(b.2010)
    #p.2010 <- st_cast(b.2010,"POLYGON")
    p.2010$area.m = as.numeric(st_area(p.2010))
    p.2010.df <- as.data.frame(p.2010)
    p.2010.df.out <-  p.2010.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_pest_2010_m2 = sum(area.m))
    r.pest.out.total = left_join(r.pest.out.total,p.2010.df.out) # add to the data summary

    r.pest.out.total.te = r.pest.out.total


# Tweedsmuir
r.pest.tw <-  st_read(dsn = Base, layer ="Pest_clip_Tw_IBMIBS") #; plot(st_geometry(r.pest.tw))
    r.pest2 <- r.pest.tw
    r.pest2<- st_zm(r.pest2 ,drop = TRUE) # this is a linear feature so need to buffer to estimate area calcs
    r.pest2 = st_buffer(r.pest2,250)#; st_is_valid(r.pest2)
    #r.pest2 = st_make_valid(r.pest2)
    r.pest2$TimeSincePest = 2018-r.pest2$CAPTURE_YEAR
    #plot(st_geometry(r.pest2), col = 'red') ; plot(st_geometry(all.range),add = T)

    # burns 0-40 years only one
    r.p.0.402 = r.pest2[r.pest2$TimeSincePest <41,]; #sort(unique(b.r.0$TimeSinceBurn))
    r.p.0.402 = st_union(r.p.0.402)
    r.p.0.402 <- st_cast(r.p.0.402,"POLYGON") #; st_is_valid(r.p.0.40)
    r.p.0.402= st_intersection(all.range,r.p.0.402)

    r.p.0.402$area.m = as.numeric(st_area(r.p.0.402))
    r.p.0.402.df = data.frame(r.p.0.402)        # calculate the length per range
    r.p.0.402.out  =  r.p.0.402.df%>%
      group_by(SiteName,V17_CH) %>%
      summarise(R_Pest040_m2 = sum(area.m))

    r.pest.df2 = r.pest2

    # add decade data   # add a column to differentiate the age brackets of pest capture
    r.pest.df2<- r.pest.df2 %>% mutate(dec.period = ifelse(CAPTURE_YEAR >= 1950 & CAPTURE_YEAR <= 1959,1950,0))
    r.pest.df2<-mutate(r.pest.df2,dec.period = ifelse(CAPTURE_YEAR  >= 1960 & CAPTURE_YEAR <= 1969,1960,dec.period))
    r.pest.df2<- mutate(r.pest.df2,dec.period = ifelse(CAPTURE_YEAR >= 1970 & CAPTURE_YEAR <= 1979,1970,dec.period))
    r.pest.df2<- mutate(r.pest.df2,dec.period = ifelse(CAPTURE_YEAR >= 1980 & CAPTURE_YEAR <= 1989,1980,dec.period))
    r.pest.df2<- mutate(r.pest.df2,dec.period = ifelse(CAPTURE_YEAR >= 1990 & CAPTURE_YEAR <= 1999,1990,dec.period))
    r.pest.df2<- mutate(r.pest.df2,dec.period = ifelse(CAPTURE_YEAR >= 2000 & CAPTURE_YEAR <= 2009,2000,dec.period))
    r.pest.df2<- mutate(r.pest.df2,dec.period = ifelse(CAPTURE_YEAR >= 2010 & CAPTURE_YEAR <= 2018,2010,dec.period))

    # Generate the cumulative pest damage into decades to add to "static Disturbance"
    Pest.dec.19502 <- r.pest.df2 %>% filter(dec.period == 1950)
    Pest.dec.19602 <- r.pest.df2 %>% filter(dec.period < 1961 ); Pest.dec.19602 <-st_union(Pest.dec.19602 )
    Pest.dec.19702 <- r.pest.df2 %>% filter(dec.period < 1971 ); Pest.dec.19702 <-st_union(Pest.dec.19702 )
    Pest.dec.19802 <- r.pest.df2 %>% filter(dec.period < 1981 ); Pest.dec.19802 <-st_union(Pest.dec.19802 )
    Pest.dec.19902 <- r.pest.df2 %>% filter(dec.period < 1991 ); Pest.dec.19902 <-st_union(Pest.dec.19902 )
    Pest.dec.20002 <- r.pest.df2 %>% filter(dec.period < 2001 );
    Pest.dec.20102 <- r.pest.df2 %>% filter(dec.period < 2011 )#; Pest.dec.20102 <-st_union(Pest.dec.20102 )
    # Burn.dec.2010 <- Burn.dec %>% filter(dec.period < 2011 )

    #  split this into seperate pest files for 2000 and 2010 as these files are too large to handle in R
    Pest.dec.20002.only <- r.pest.df2 %>% filter(dec.period == 2000)
    Pest.dec.20102.only  <- r.pest.df2 %>% filter(dec.period == 2011 )
    st_write(Pest.dec.20002.only,paste(shape.output.dir,"Pest.tw.dec.2000.only.buf.shp",sep = "")) # this writes out as single layer
    st_write(Pest.dec.20102.only,paste(shape.output.dir,"Pest.tw.dec.2010.only.buf.shp",sep = "")) # this writes out as single layer


    ## write out the shapefiles to Data.drive
    st_write(Pest.dec.19502 ,paste(shape.output.dir,"Pest.tw.dec.1950.buf.shp",sep = "")) # this writes out as single layer
    st_write(Pest.dec.19602 ,paste(shape.output.dir,"Pest.tw.dec.1960.buf.shp",sep = "")) # this writes out as single layer
    st_write(Pest.dec.19702 ,paste(shape.output.dir,"Pest.tw.dec.1970.buf.shp",sep = "")) # this writes out as single layer
    st_write(Pest.dec.19802 ,paste(shape.output.dir,"Pest.tw.dec.1980.buf.shp",sep = "")) # this writes out as single layer
    st_write(Pest.dec.19902,paste(shape.output.dir,"Pest.tw.dec.1990.buf.shp",sep = "")) # this writes out as single layer
    #st_write(Pest.dec.20002,paste(shape.output.dir,"Pest.tw.dec.2000.buf.shp",sep = "")) # this writes out as single layer
    #st_write(Pest.dec.20102,paste(shape.output.dir,"Pest.tw.dec.2010.buf.shp",sep = "")) # this writes out as single layer

    ## add the section for each decade (cumulative )
    # generate consolidated outputs overtime per decage

    #p.1950 <- st_union(Pest.dec.1950) # empty

    # check if there is a 1950;s
    p.19602 <- st_union(Pest.dec.19602)
    p.19602 <- st_cast(p.19602,"POLYGON") #; st_is_valid(p.19602)
    p.19602 = st_intersection(all.range, p.19602) #; p.19602 <- st_make_valid(p.19602)
    p.19602$area.m = as.numeric(st_area(p.19602))
    p.19602.df <- as.data.frame(p.19602)
    p.19602.df.out <-  p.19602.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_pest_1960_m2 = sum(area.m))
    r.pest.out.total2 = left_join(r.p.0.402.out,p.19602.df.out) # add to the data summary

    p.19702 <- st_union(Pest.dec.19702)
    p.19702 <- st_cast(p.19702,"POLYGON") #; st_is_valid(p.19702)
    p.19702 = st_intersection(all.range, p.19702) ; p.19702 <- st_make_valid(p.19702)
    #b.19702 <- st_cast(b.19702,"POLYGON")
    p.19702$area.m = as.numeric(st_area(p.19702))
    p.19702.df <- as.data.frame(p.19702)
    p.19702.df.out <-  p.19702.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_pest_1970_m2 = sum(area.m))
    r.pest.out.total2 = left_join( r.pest.out.total2,p.19702.df.out) # add to the data summary

    p.19802 <- st_union(Pest.dec.19802)
    p.19802 <- st_cast(p.19802,"POLYGON") #; st_is_valid(p.19802)
    p.19802 = st_intersection(all.range, p.19802) #; b.1980 <- st_make_valid(b.1980)
    #b.19802 <- st_cast(b.1980,"POLYGON") #; st_is_valid(b.1980)
    p.19802$area.m = as.numeric(st_area(p.19802))
    p.19802.df <- as.data.frame(p.19802)
    p.19802.df.out <-  p.19802.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_pest_1980_m2 = sum(area.m))
    r.pest.out.total2 = left_join(r.pest.out.total2,p.19802.df.out) # add to the data summary

    p.19902 <- st_union(Pest.dec.19902)
    p.19902 <- st_cast(p.19902,"POLYGON") #; st_is_valid(p.19902)
    p.19902 = st_intersection(all.range, p.19902) #; b.1990 <- st_make_valid(c.1990) ; head(p.19902)
    #p.1990 <- st_cast(b.1990,"POLYGON")
    p.19902$area.m = as.numeric(st_area(p.19902))
    p.19902.df <- as.data.frame(p.19902)
    p.19902.df.out <-  p.19902.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_pest_1990_m2 = sum(area.m))
    r.pest.out.total2 = left_join(r.pest.out.total2,p.19902.df.out) # add to the data summary

          #r.pest.out.total.tw = r.pest.out.total2
        ###########################################
        ## Join the telkwa and Tweedsmuir herd info together
        #r.pest.out.all = rbind(r.pest.out.total.te,r.pest.out.total.tw)
        #all.range.out <- left_join(all.range.out,r.pest.out.all)
        #all.range.out[is.na(all.range.out)]<-0
        #write.csv(all.range.out,paste(temp.dir,"All_data_summary_buffer_temp.csv",sep =""))

        all.range.out<- read.csv(paste(temp.dir,"All_data_summary_buffer_temp.csv",sep =""))

    ##################################
    # need to dissolve these files in arcmap then read back in here:
    processed.files = "Z:/01.Projects/Wildlife/Caribou/02.Disturbance/TweedTelkwa/Outputs/disturb_layers/DistLayer_buffered/"
    p.20002 <- st_read(paste(processed.files,"Pesttw2000bufd.shp",sep = ""))
    p.20002 <- st_union(p.20002)
    p.20002 <- st_cast(p.20002,"POLYGON") #; st_is_valid(p.20002)
    p.20002 = st_intersection(all.range, p.20002) #; p.20002 <- st_make_valid(p.20002)
    p.20002$area.m = as.numeric(st_area(p.20002))
    p.20002.df <- as.data.frame(p.20002)
    p.20002.df.out <-  p.20002.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_pest_2000_m2 = sum(area.m))
    r.pest.out.total2 = left_join(r.pest.out.total2,p.20002.df.out) # add to the data summary
    #r.pest.out.total2 = p.20002.df.out

    p.20102 <- st_read(paste(processed.files,"Pesttw2010bufd.shp",sep = ""))
    p.20102 <- st_union(p.20102)
    p.20102 <- st_cast(p.20102,"POLYGON") #; st_is_valid(p.20102); p.20102 <- st_make_valid(p.20102)
    p.20102 <- st_make_valid(p.20102)
    p.20102 = st_intersection(all.range, p.20102)
    p.20102$area.m = as.numeric(st_area(p.20102))
    p.20102.df <- as.data.frame(p.20102)
    p.20102.df.out <-  p.20102.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(R_pest_2010_m2 = sum(area.m))
    r.pest.out.total2 = left_join(r.pest.out.total2,p.20102.df.out) # add to the data summary

    r.pest.out.total.tw = r.pest.out.total2

    # or
    # temporary format fix start:
   # r.pest.temp = all.range.out[,c("SiteName","V17_CH", "R_Pest040_m2","R_pest_1960_m2", "R_pest_1970_m2", "R_pest_1980_m2",
  #                                 "R_pest_1990_m2","R_pest_2000_m2","R_pest_2010_m2")]
#
#    r.pest.temp.te = r.pest.temp %>% dplyr::filter(SiteName == "Telkwa")
#    r.pest.temp.tw = r.pest.temp %>% dplyr::filter(SiteName == "Tweedsmuir")
#    r.pest.temp.tw = r.pest.temp.tw[,1:7]
#    r.pest.temp.tw = left_join(r.pest.temp.tw,r.pest.out.total.tw)
#    r.pest.out.all = rbind( r.pest.temp.te,r.pest.temp.tw)#
#
#    all.range.out =  all.range.out[,c(4:36)]
#    all.range.out <- left_join(all.range.out,r.pest.out.all)
    # Temp fix end.


    ###########################################
# Join the telkwa and Tweedsmuir herd info together
r.pest.out.all = rbind(r.pest.out.total.te,r.pest.out.total.tw)

all.range.out <- left_join(all.range.out,r.pest.out.all)
all.range.out[is.na(all.range.out)]<-0

write.csv(all.range.out,paste(temp.dir,"All_data_summary_buffer_temp.csv",sep =""))

all.range.out <- read.csv(paste(temp.dir,"All_data_summary_buffer_temp.csv",sep = ""))
############################################################################

## agregate all the temporal datasets with the static data sets.
## This was done in Arcmap using union and dissolve due to the large size of the files.
# output the files into a single geodatabase, this can then be read into the following script

## Step 1:

arc.gdb= "Z:/01.Projects/Wildlife/Caribou/02.Disturbance/TweedTelkwa/Outputs/disturb_layers/DistFinal_buff.gdb" # contains
#arc.gdb= "C:\Temp\TweedTelkwa\Temp\Perkins\Outputs\disturb_layers#
# check on the C drive


## List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name))
final_list <- ogrListLayers(arc.gdb); print(final_list)

# read in the union and dissolved data and calculate the areas to add to the table.

#1950
  dist.1950<- st_read(dsn=arc.gdb,layer="dist1950buf_d") # multipoly
  dist.1950<- st_set_crs(dist.1950,3005)
  dist.1950 <- st_union(dist.1950)
  #head(dist.1950); plot(dist.1950)
  dist.1950 <-st_cast(dist.1950,"POLYGON")
  dist.1950 <- st_intersection(all.range, dist.1950)
  dist.1950$area.m = as.numeric(st_area(dist.1950))
  dist.1950.df <- as.data.frame(dist.1950)
  dist.1950.df.out <-  dist.1950.df %>%
    group_by(SiteName,V17_CH ) %>%
    summarise(Total_1950_m2 = sum(area.m))
  all.dist.tally_buf = dist.1950.df.out
  rm("dist.1950")

#1960
  dist.1960<- st_read(dsn=arc.gdb,layer="dist1960buf1_d") # multipoly # read in the buffered polygon
  dist.1960<- st_set_crs(dist.1960,3005) # set the CRS
  dist.1960 <- st_union(dist.1960)# dissolve
  #head(dist.1960); plot(dist.1960)
  dist.1960 <-st_cast(dist.1960,"POLYGON") # fix mulitpolygons
  dist.1960 <- st_intersection(all.range, dist.1960)
  dist.1960$area.m = as.numeric(st_area(dist.1960))
  dist.1960.df <- as.data.frame(dist.1960)
  dist.1960.df.out <-  dist.1960.df %>%
    group_by(SiteName,V17_CH ) %>%
    summarise(Total_1960_m2 = sum(area.m))# calculate the total area
  all.dist.tally_buf = left_join(all.dist.tally_buf,dist.1960.df.out)# add the total area to the previous decade
  rm("dist.1960")

#1970
  dist.1970<- st_read(dsn=arc.gdb,layer="dist1970buf_d") # multipoly
  dist.1970<- st_set_crs(dist.1970,3005)
  dist.1970 <- st_union(dist.1970)
  #head(dist.1970); plot(dist.1970)
  dist.1970 <-st_cast(dist.1970,"POLYGON")
  dist.1970 <- st_intersection(all.range, dist.1970)
  dist.1970$area.m = as.numeric(st_area(dist.1970))
  dist.1970.df <- as.data.frame(dist.1970)
  dist.1970.df.out <-  dist.1970.df %>%
    group_by(SiteName,V17_CH ) %>%
    summarise(Total_1970_m2 = sum(area.m))# calculate the total area
  all.dist.tally_buf =left_join(all.dist.tally_buf, dist.1970.df.out)
  rm("dist.1970")

  #1980
  dist.1980<- st_read(dsn=arc.gdb,layer="dist1980buf_d") # multipoly
  dist.1980<- st_set_crs(dist.1980,3005)
  dist.1980 <- st_union(dist.1980)
  #head(dist.1980); plot(dist.1980)
  dist.1980 <-st_cast(dist.1980,"POLYGON")
  dist.1980 <- st_intersection(all.range, dist.1980)
  dist.1980$area.m = as.numeric(st_area(dist.1980))
  dist.1980.df <- as.data.frame(dist.1980)
  dist.1980.df.out <-  dist.1980.df %>%
    group_by(SiteName,V17_CH ) %>%
    summarise(Total_1980_m2 = sum(area.m))# calculate the total area
  all.dist.tally_buf = left_join(all.dist.tally_buf, dist.1980.df.out)
  rm("dist.1980")

  #1990
  dist.1990<- st_read(dsn=arc.gdb,layer="dist1990buf_d") # multipoly
  dist.1990<- st_set_crs(dist.1990,3005)
  dist.1990 <- st_union(dist.1990)
  #head(dist.1990); plot(dist.1990)
  dist.1990 <-st_cast(dist.1990,"POLYGON")
  dist.1990 <- st_intersection(all.range, dist.1990)
  dist.1990$area.m = as.numeric(st_area(dist.1990))
  dist.1990.df <- as.data.frame(dist.1990)
  dist.1990.df.out <-  dist.1990.df %>%
    group_by(SiteName,V17_CH ) %>%
    summarise(Total_1990_m2 = sum(area.m))# calculate the total area
  all.dist.tally_buf = left_join(all.dist.tally_buf, dist.1990.df.out)
  rm("dist.1990")

#2000
  dist.2000<- st_read(dsn=arc.gdb,layer="dist2000buf_d") # multipoly
  dist.2000<- st_set_crs(dist.2000,3005)
  dist.2000 <- st_union(dist.2000)
  #head(dist.2000); plot(dist.2000)
  dist.2000 <-st_cast(dist.2000,"POLYGON")
  dist.2000 <- st_intersection(all.range, dist.2000)
  dist.2000$area.m = as.numeric(st_area(dist.2000))
  dist.2000.df <- as.data.frame(dist.2000)
  dist.2000.df.out <- dist.2000.df %>%
    group_by(SiteName,V17_CH ) %>%
    summarise(Total_2000_m2 = sum(area.m))# calculate the total area
  all.dist.tally_buf = left_join(all.dist.tally_buf, dist.2000.df.out)
  rm("dist.2000")

  # still to create the 2000_buf_dissolve in arcmap
 #2010
  dist.2010<- st_read(dsn=arc.gdb,layer="dist2010buf_d") # multipoly
  dist.2010<- st_set_crs( dist.2010,3005)
  dist.2010 <- st_union( dist.2010)
  head( dist.2010); plot( dist.2010)

  dist.2010 <-st_cast( dist.2010,"POLYGON")
  dist.2010 <- st_intersection(all.range, dist.2010)
  dist.2010$area.m = as.numeric(st_area(dist.2010))
  dist.2010.df <- as.data.frame(dist.2010)
  dist.2010.df.out <- dist.2010.df %>%
    group_by(SiteName,V17_CH ) %>%
    summarise(Total_2010_m2 = sum(area.m))# calculate the total area

  all.dist.tally_buf = left_join(all.dist.tally_buf, dist.2010.df.out)
  rm("dist.2010")

write.csv(all.dist.tally_buf,paste(temp.dir,"Total_disturb_summary_buffer_temp.csv",sep =""))
all.dist.tally_buf = read.csv(paste(temp.dir,"Total_disturb_summary_buffer_temp.csv",sep =""))


# create disturbance per decade without pest (only use the 40 years rather than consolidated)
# These files are created in Arcmap by joining together burns and cuts (for each decade - 1950,1960,1970,etc) then combined
# with the disturbance for the last 40 yeare (1950-1989, 1960 - 1999, 1970 - 2009,1980 - 2018) and static disturbnace.
# these are unioned and dissolved then output to the gdb

    #1950 - 1989
    dist.195089<- st_read(dsn=arc.gdb,layer="dist19501989buf_d") # multipoly
    dist.195089<- st_set_crs(dist.195089,3005)
    dist.195089 <- st_union(dist.195089)
    plot(st_geometry(dist.195089))
    dist.195089 <-st_cast(dist.195089,"POLYGON")
    dist.195089 <- st_intersection(all.range, dist.195089)
    dist.195089$area.m = as.numeric(st_area(dist.195089))
    dist.195089.df <- as.data.frame(dist.195089)
    dist.195089.df.out <-  dist.195089.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(Total_195089_m2 = sum(area.m))
    all.dist.tally.buf.dec = dist.195089.df.out
    rm("dist.195089")

    #1960 - 1999
    dist.196099<- st_read(dsn=arc.gdb,layer="dist196099buf_dd") # multipoly
    dist.196099<- st_set_crs(dist.196099,3005)
    dist.196099 <- st_union(dist.196099)
    #head(dist.1950); plot(dist.1950)
    dist.196099 <-st_cast(dist.196099,"POLYGON")
    dist.196099 <- st_intersection(all.range, dist.196099)
    dist.196099$area.m = as.numeric(st_area(dist.196099))
    dist.196099.df <- as.data.frame(dist.196099)
    dist.196099.df.out <-  dist.196099.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(Total_196099_m2 = sum(area.m))

    all.dist.tally.buf.dec = left_join(all.dist.tally.buf.dec,dist.196099.df.out)
    rm("dist.196099")

    #1970 - 2009
    dist.197009<- st_read(dsn=arc.gdb,layer="dist19702009buf_d") # multipoly
    dist.197009<- st_set_crs(dist.197009,3005)
    dist.197009 <- st_union(dist.197009)
    #head(dist.1950); plot(dist.1970)
    dist.197009 <-st_cast(dist.197009,"POLYGON")
    dist.197009 <- st_intersection(all.range, dist.197009)
    dist.197009$area.m = as.numeric(st_area(dist.197009))
    dist.197009.df <- as.data.frame(dist.197009)
    dist.197009.df.out <-  dist.197009.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(Total_197009_m2 = sum(area.m))

    all.dist.tally.buf.dec = left_join(all.dist.tally.buf.dec,dist.197009.df.out)
    rm("dist.197009")

    #1980 - 2019
    dist.198018<- st_read(dsn=arc.gdb,layer="dist198018buf_d") # multipoly
    dist.198018<- st_set_crs(dist.198018,3005)
    dist.198018 <- st_union(dist.198018)
    #head(dist.1950); plot(dist.1970)
    dist.198018 <-st_cast(dist.198018,"POLYGON")
    dist.198018 <- st_intersection(all.range, dist.198018)
    dist.198018$area.m = as.numeric(st_area(dist.198018))
    dist.198018.df <- as.data.frame(dist.198018)
    dist.198018.df.out <-  dist.198018.df %>%
      group_by(SiteName,V17_CH ) %>%
      summarise(Total_198018_m2 = sum(area.m))

    all.dist.tally.buf.dec = left_join(all.dist.tally.buf.dec,dist.198018.df.out)
    rm("dist.198018")

    write.csv(all.dist.tally.buf.dec,paste(temp.dir,"Total_disturb_nopest_summary_dec_buf.csv",sep =""))
    all.dist.tally.buf.dec = read.csv(paste(temp.dir,"Total_disturb_nopest_summary_dec_buf.csv",sep =""))

#############################################################################

# Step 2 aggregate table and output and write out

# Combine the static data with the temporal data (all.temp.out)
all.range.out = read.csv(paste(temp.dir,"All_data_summary_buffer_temp.csv",sep = ""))
all.range.out = left_join(all.range.out,all.dist.tally_buf)
all.range.out = left_join(all.range.out,all.dist.tally.buf.dec)
all.range.out = all.range.out[,-c(1,2)]
all.range.out = all.range.out[,-c(1)]

# Read in summary of area and habitat created in script 1.
Herd_key_detail= read.csv(paste(out.dir,"Herd_key_detail.csv",sep = ""))
all.range.out <- left_join(Herd_key_detail,all.range.out)

# split the data into ha and percent cover tables
tout.key <-all.range.out %>%                     # grab the non-numeric cols only to add back to converted table
  dplyr::select(X,SiteName, V17_CH)

# convert all values to ha.
tout.ha <- all.range.out %>%                     # convert meters squared to ha
  dplyr::select(-c(SiteName, V17_CH)) %>% # remove non-numeric columns
  mutate_all(funs(ha =  ./10000))



### STILL TO DO
### ADJUST THE VALUES HEER


tout.ha = tout.ha[,c(1,2,52:98)]
tout.ha <- left_join(tout.key,tout.ha)  # join back together key and numeric data

# calculate the percentage
tout.pc <- all.range.out%>%                     # convert meters squared to ha
  dplyr::select(-c(SiteName, V17_CH)) %>% # remove non-numeric columns
  mutate_all(funs(pc = (( ./10000)/R_area_ha *100)))
tout.pc = tout.pc[,c(1,2,52:98)]
#tout.pc = tout.pc[,c(1,54:100)] # check these

tout.pc <- left_join(tout.key,tout.pc)  # join back together key and numeric data

# reformat the ha table
## may also want to reshape to convert to habitat catergories as columns rather than row.
out.ha =tout.ha %>% reshape2:::melt(.) %>% spread(.,V17_CH,value)
out.ha[is.na(out.ha)]<-0

# change the columns names to something short and manageble
colnames(out.ha)[colnames(out.ha)=='High Elevation Winter/Summer Range'] <- "HWSR"
colnames(out.ha)[colnames(out.ha)=='Low Elevation Summer Range'] <- "LSR"
colnames(out.ha)[colnames(out.ha)=='Low Elevation Winter Range'] <- "LWR"
colnames(out.ha)[colnames(out.ha)=='Matrix Range'] <- "MR"

write.csv(out.ha,paste(out.dir,"Final_TT_buf_summary_ha.csv",sep= ""))  # write out the ha table

# reformat the pc table
## may also want to reshape to convert to habitat catergories as columns rather than row.
out.pc =tout.pc %>% reshape2:::melt(.) %>% spread(.,V17_CH,value)
out.pc[is.na(out.pc)]<-0

# change the columns names to something short and manageble
colnames(out.pc)[colnames(out.pc)=='High Elevation Winter/Summer Range'] <- "HWSR"
colnames(out.pc)[colnames(out.pc)=='Low Elevation Summer Range'] <- "LSR"
colnames(out.pc)[colnames(out.pc)=='Low Elevation Winter Range'] <- "LWR"
colnames(out.pc)[colnames(out.pc)=='Matrix Range'] <- "MR"

# round the percent values to 1 digit
out.pc = out.pc%>%
  mutate_if(is.numeric,round,digits = 1)
out.pc [is.na(out.pc)]<- 0

write.csv(out.pc,paste(out.dir,"Final_TT_buf_summary_pc.csv",sep= ""))  # write out





