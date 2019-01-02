
## SCript to aggregate all the outputs of script 1 and 2. 
## This script will generate the "total disturbance % or ha values" 


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

## set your output directory 

# to run analysis on C drive: 
out.dir = "C:/Temp/TweedTelkwa/Temp/Perkins/Outputs/"
temp.dir = "C:/Temp/TweedTelkwa/Temp/Perkins/Data/"

##########################################################################################
# Read in data tables 
# Table 1 : all disturbance 
# Table 2 : Protection at the Range and the core levels

dis1 <- read.csv(paste(temp.dir,"Temp_output_summary_all.csv",sep = ""),header = TRUE)     
herd_key <- read.csv(paste(out.dir,"Herd_key_detail.csv",sep = ""))                # import Key 

# clean up the data table 
dis1 = dis1 %>% dplyr::select(-c(X,X.1)) 
dis1 = dis1[-5,] # remove the duplicate row for low elevation (only for Telkwa)

# join all the disturbances together 
tout = left_join(herd_key,dis1) 

# split the data into ha and percent cover tables 
tout.key <-tout %>%                     # grab the non-numeric cols only to add back to converted table
  dplyr::select(X,SiteName, V17_CH)

# convert all values to ha. 
tout.ha <- tout %>%                     # convert meters squared to ha 
  dplyr::select(-c(SiteName, V17_CH)) %>% # remove non-numeric columns
  mutate_all(funs(ha =  ./10000)) 
tout.ha <- left_join(tout.key,tout.ha)  # join back together key and numeric data 

# calculate the percentage   
tout.pc <- tout %>%                     # convert meters squared to ha 
  dplyr::select(-c(SiteName, V17_CH)) %>% # remove non-numeric columns
  mutate_all(funs(pc = (( ./10000)/R_area_ha *100))) 
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
      
      # write out the ha table 
      write.csv(out.ha,paste(out.dir,"Final_TT_summary_ha.csv",sep= ""))  

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
      
      # write out 
      write.csv(out.pc,paste(out.dir,"Final_TT_summary_pc.csv",sep= ""))  


      