###
### 9/17/2025
###

### Notes from previous versions of this script

{

#### Create Data Frame For Rocky Ford 
#### 1) Use All instrument types (i.e., all imagery avail), but keep separate as independent daily obs.
#### 2) Get all met from RFD01
#### 3) Get all discharge/river data from CDWR gage
#### 4) Vegetation Types are: (1) Cottonwoods, (2) Tam, (3) Russian Olive, (4) Prairie Grass
####                                (5) Wetland grasses, (6) Willow
####
####


#
# Best examples of grasses are:
#
# 1) Certain type of grass all over near well A1 and A2 (see pics), lots of points
# 
# 2) Mixed forb, prairie grass at well c3 and c4 (see pics), lots of points
#
# 3) Certain type of grass community near well A3 (see pics), 6-ish points, similar to A1 and A2
#
# 4) Wetland like grasses near well C1, 4-ish points
# 
# 5) Prairie grass community near well B4, few points
#
# 6) Grass near well A4 is like C3 and C4, only one point
#

}



library(dplyr)
library(tidyr) # only need the "gather" function
library(sf)
library(terra)
library(ggplot2)




# Bring in GROUNDWATER data for Rocky Ford (all 12 wells)
# This chunk is the same script on GitHub repository
# "arkansasR-gw-data"

# Both PhD and postdoc datasets have NA's in them.
# The PhD data has NA's on 11/3/2018 for approx 12 hours worth of data.
# PhD data is from 5/14/2018 (3:30:00PM) to 10/17/2020 (9:30:00AM)
# Postdoc data is 10/25/2023 (1:00:00PM) to 6/26/2025 (9:30:00AM)

{
  
  
  library(dplyr)
  
  # Bring in GROUNDWATER data for Rocky Ford (all 12 wells)
  
  # Both datasets have NA's in them. The PhD data has NA's on 11/3/2018 for approx 12 hours worth of data.
  # PhD data is from 5/14/2018 (3:30:00PM) to 10/17/2020 (9:30:00AM)
  # Postdoc data is 10/25/2023 (1:00:00PM) to 6/26/2025 (9:30:00AM)
  
  tranA_phd <- read.csv(file = "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/To_GitHub/raw_gw_data/transectA_phd.csv",header = TRUE)
  tranB_phd <- read.csv(file = "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/To_GitHub/raw_gw_data/transectB_phd.csv",header = TRUE)
  tranC_phd <- read.csv(file = "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/To_GitHub/raw_gw_data/transectC_phd.csv",header = TRUE)
  
  tranA_postdoc <- read.csv(file = "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/To_GitHub/raw_gw_data/transectA_postdoc.csv",header = TRUE)
  tranB_postdoc <- read.csv(file = "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/To_GitHub/raw_gw_data/transectB_postdoc.csv",header = TRUE)
  tranC_postdoc <- read.csv(file = "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/To_GitHub/raw_gw_data/transectC_postdoc.csv",header = TRUE)
  
  
  # Convert date column to actual date recognized by R
  
  tranA_phd$datez <- as.Date(tranA_phd$Date_corr, format = "%m/%d/%Y")
  tranB_phd$datez <- as.Date(tranB_phd$Date_corr, format = "%m/%d/%Y")
  tranC_phd$datez <- as.Date(tranC_phd$Date_corr, format = "%m/%d/%Y")
  
  tranA_postdoc$datez <- as.Date(tranA_postdoc$Date_corr, format = "%m/%d/%Y")
  tranB_postdoc$datez <- as.Date(tranB_postdoc$Date_corr, format = "%m/%d/%Y")
  tranC_postdoc$datez <- as.Date(tranC_postdoc$Date_corr, format = "%m/%d/%Y")
  
  
  # Drop the last day of data as it is filled with error data
  
  tranA_phd <- tranA_phd[tranA_phd$datez < "2020-10-17", ]
  tranB_phd <- tranB_phd[tranB_phd$datez < "2020-10-17", ]
  tranC_phd <- tranC_phd[tranC_phd$datez < "2020-10-17", ]
  
  tranA_postdoc <- tranA_postdoc[tranA_postdoc$datez < "2025-06-26", ]
  tranB_postdoc <- tranB_postdoc[tranB_postdoc$datez < "2025-06-26", ]
  tranC_postdoc <- tranC_postdoc[tranC_postdoc$datez < "2025-06-26", ]
  
  
  ### PhD Data
  
  
  
  # Replace outliers with NA instead of removing rows
  
  tranA_phd_list <- list()
  
  for(i in 3:(ncol(tranA_phd)-1)){
    
    a <- tranA_phd[,c(i,17)]
    colnames(a) <- c("varz","datez")
    
    data_with_na_outliers <- a %>%
      group_by(datez) %>%
      mutate(
        Q1 = quantile(varz, 0.25, na.rm = TRUE),
        Q3 = quantile(varz, 0.75, na.rm = TRUE),
        IQR = Q3 - Q1,
        lower_bound = Q1 - 10.0 * IQR,
        upper_bound = Q3 + 10.0 * IQR,
        varz = ifelse(varz < lower_bound | varz > upper_bound, NA, varz)
      )%>%
      ungroup() %>%
      select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound)
    
    tranA_phd_list[[i]] <- data_with_na_outliers
    
  }
  
  tranA_phd_clean <- bind_cols(tranA_phd_list)
  tranA_phd_clean <- tranA_phd_clean[,-c(2,4,6,8,10,12,14,16,18,20,22,24,26)]
  colnames(tranA_phd_clean) <- colnames(tranA_phd[,3:17])
  
  tranA_phd_clean <- as.data.frame(tranA_phd_clean)
  
  # Take average of values by day
  
  tranA_phd_avg <- aggregate(tranA_phd_clean[1:14],
                             by = list(datez = tranA_phd_clean$datez), FUN = function(x) mean(x, na.rm = TRUE))
  
  # Convert all instances of NaN's to NA's
  
  tranA_phd_avg <- tranA_phd_avg %>%
    mutate(across(everything(), ~ifelse(is.nan(.), NA, .)))
  
  
  
  # Replace outliers with NA instead of removing rows
  
  tranB_phd_list <- list()
  
  for(i in 3:(ncol(tranB_phd)-1)){
    
    a <- tranB_phd[,c(i,17)]
    colnames(a) <- c("varz","datez")
    
    data_with_na_outliers <- a %>%
      group_by(datez) %>%
      mutate(
        Q1 = quantile(varz, 0.25, na.rm = TRUE),
        Q3 = quantile(varz, 0.75, na.rm = TRUE),
        IQR = Q3 - Q1,
        lower_bound = Q1 - 10.0 * IQR,
        upper_bound = Q3 + 10.0 * IQR,
        varz = ifelse(varz < lower_bound | varz > upper_bound, NA, varz)
      )%>%
      ungroup() %>%
      select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound)
    
    tranB_phd_list[[i]] <- data_with_na_outliers
    
  }
  
  tranB_phd_clean <- bind_cols(tranB_phd_list)
  tranB_phd_clean <- tranB_phd_clean[,-c(2,4,6,8,10,12,14,16,18,20,22,24,26)]
  colnames(tranB_phd_clean) <- colnames(tranB_phd[,3:17])
  
  tranB_phd_clean <- as.data.frame(tranB_phd_clean)
  
  # Take average of values by day
  
  tranB_phd_avg <- aggregate(tranB_phd_clean[1:14],
                             by = list(datez = tranB_phd_clean$datez), FUN = function(x) mean(x, na.rm = TRUE))
  
  # Convert all instances of NaN's to NA's
  
  tranB_phd_avg <- tranB_phd_avg %>%
    mutate(across(everything(), ~ifelse(is.nan(.), NA, .)))
  
  
  # Replace outliers with NA instead of removing rows
  
  tranC_phd_list <- list()
  
  for(i in 3:(ncol(tranC_phd)-1)){
    
    a <- tranC_phd[,c(i,17)]
    colnames(a) <- c("varz","datez")
    
    data_with_na_outliers <- a %>%
      group_by(datez) %>%
      mutate(
        Q1 = quantile(varz, 0.25, na.rm = TRUE),
        Q3 = quantile(varz, 0.75, na.rm = TRUE),
        IQR = Q3 - Q1,
        lower_bound = Q1 - 10.0 * IQR,
        upper_bound = Q3 + 10.0 * IQR,
        varz = ifelse(varz < lower_bound | varz > upper_bound, NA, varz)
      )%>%
      ungroup() %>%
      select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound)
    
    tranC_phd_list[[i]] <- data_with_na_outliers
    
  }
  
  tranC_phd_clean <- bind_cols(tranC_phd_list)
  tranC_phd_clean <- tranC_phd_clean[,-c(2,4,6,8,10,12,14,16,18,20,22,24,26)]
  colnames(tranC_phd_clean) <- colnames(tranC_phd[,3:17])
  
  tranC_phd_clean <- as.data.frame(tranC_phd_clean)
  
  # Take average of values by day
  
  tranC_phd_avg <- aggregate(tranC_phd_clean[1:14],
                             by = list(datez = tranC_phd_clean$datez), FUN = function(x) mean(x, na.rm = TRUE))
  
  # Convert all instances of NaN's to NA's
  
  tranC_phd_avg <- tranC_phd_avg %>%
    mutate(across(everything(), ~ifelse(is.nan(.), NA, .)))
  
  
  
  ### Postdoc Data
  
  
  
  # Replace outliers with NA instead of removing rows
  
  tranA_postdoc_list <- list()
  
  for(i in 3:(ncol(tranA_postdoc)-1)){
    
    a <- tranA_postdoc[,c(i,17)]
    colnames(a) <- c("varz","datez")
    
    data_with_na_outliers <- a %>%
      group_by(datez) %>%
      mutate(
        Q1 = quantile(varz, 0.25, na.rm = TRUE),
        Q3 = quantile(varz, 0.75, na.rm = TRUE),
        IQR = Q3 - Q1,
        lower_bound = Q1 - 10.0 * IQR,
        upper_bound = Q3 + 10.0 * IQR,
        varz = ifelse(varz < lower_bound | varz > upper_bound, NA, varz)
      )%>%
      ungroup() %>%
      select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound)
    
    tranA_postdoc_list[[i]] <- data_with_na_outliers
    
  }
  
  tranA_postdoc_clean <- bind_cols(tranA_postdoc_list)
  tranA_postdoc_clean <- tranA_postdoc_clean[,-c(2,4,6,8,10,12,14,16,18,20,22,24,26)]
  colnames(tranA_postdoc_clean) <- colnames(tranA_postdoc[,3:17])
  
  tranA_postdoc_clean <- as.data.frame(tranA_postdoc_clean)
  
  # Take average of values by day
  
  tranA_postdoc_avg <- aggregate(tranA_postdoc_clean[1:14],
                                 by = list(datez = tranA_postdoc_clean$datez), FUN = function(x) mean(x, na.rm = TRUE))
  
  # Convert all instances of NaN's to NA's
  
  tranA_postdoc_avg <- tranA_postdoc_avg %>%
    mutate(across(everything(), ~ifelse(is.nan(.), NA, .)))
  
  
  
  # Replace outliers with NA instead of removing rows
  
  tranB_postdoc_list <- list()
  
  for(i in 3:(ncol(tranB_postdoc)-1)){
    
    a <- tranB_postdoc[,c(i,17)]
    colnames(a) <- c("varz","datez")
    
    data_with_na_outliers <- a %>%
      group_by(datez) %>%
      mutate(
        Q1 = quantile(varz, 0.25, na.rm = TRUE),
        Q3 = quantile(varz, 0.75, na.rm = TRUE),
        IQR = Q3 - Q1,
        lower_bound = Q1 - 10.0 * IQR,
        upper_bound = Q3 + 10.0 * IQR,
        varz = ifelse(varz < lower_bound | varz > upper_bound, NA, varz)
      )%>%
      ungroup() %>%
      select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound)
    
    tranB_postdoc_list[[i]] <- data_with_na_outliers
    
  }
  
  tranB_postdoc_clean <- bind_cols(tranB_postdoc_list)
  tranB_postdoc_clean <- tranB_postdoc_clean[,-c(2,4,6,8,10,12,14,16,18,20,22,24,26)]
  colnames(tranB_postdoc_clean) <- colnames(tranB_postdoc[,3:17])
  
  tranB_postdoc_clean <- as.data.frame(tranB_postdoc_clean)
  
  # Take average of values by day
  
  tranB_postdoc_avg <- aggregate(tranB_postdoc_clean[1:14],
                                 by = list(datez = tranB_postdoc_clean$datez), FUN = function(x) mean(x, na.rm = TRUE))
  
  # Convert all instances of NaN's to NA's
  
  tranB_postdoc_avg <- tranB_postdoc_avg %>%
    mutate(across(everything(), ~ifelse(is.nan(.), NA, .)))
  
  
  # Replace outliers with NA instead of removing rows
  
  tranC_postdoc_list <- list()
  
  for(i in 3:(ncol(tranC_postdoc)-1)){
    
    a <- tranC_postdoc[,c(i,17)]
    colnames(a) <- c("varz","datez")
    
    data_with_na_outliers <- a %>%
      group_by(datez) %>%
      mutate(
        Q1 = quantile(varz, 0.25, na.rm = TRUE),
        Q3 = quantile(varz, 0.75, na.rm = TRUE),
        IQR = Q3 - Q1,
        lower_bound = Q1 - 10.0 * IQR,
        upper_bound = Q3 + 10.0 * IQR,
        varz = ifelse(varz < lower_bound | varz > upper_bound, NA, varz)
      )%>%
      ungroup() %>%
      select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound)
    
    tranC_postdoc_list[[i]] <- data_with_na_outliers
    
  }
  
  tranC_postdoc_clean <- bind_cols(tranC_postdoc_list)
  tranC_postdoc_clean <- tranC_postdoc_clean[,-c(2,4,6,8,10,12,14,16,18,20,22,24,26)]
  colnames(tranC_postdoc_clean) <- colnames(tranC_postdoc[,3:17])
  
  tranC_postdoc_clean <- as.data.frame(tranC_postdoc_clean)
  
  # Take average of values by day
  
  tranC_postdoc_avg <- aggregate(tranC_postdoc_clean[1:14],
                                 by = list(datez = tranC_postdoc_clean$datez), FUN = function(x) mean(x, na.rm = TRUE))
  
  # Convert all instances of NaN's to NA's
  
  tranC_postdoc_avg <- tranC_postdoc_avg %>%
    mutate(across(everything(), ~ifelse(is.nan(.), NA, .)))
  
  
  
  
  
  ###
  ### Plot data to make sure everything looks good
  ###
  
  # plot(tranA_phd$A1_Temp_F)
  # plot(tranA_phd_clean$A1_Temp_F)
  # plot(tranA_phd_avg$A2_Temp_F)
  # 
  # plot(tranA_phd$A1_WSEL_masl)
  # plot(tranA_phd_clean$A1_WSEL_masl)
  # plot(tranA_phd_avg$A1_WSEL_masl)
  # 
  # plot(tranA_phd$A1_DTGW_m)
  # plot(tranA_phd_clean$A1_DTGW_m)
  # plot(tranA_phd_avg$A1_DTGW_m)
  # 
  # 
  # plot(tranA_postdoc$A1_Temp_C)
  # plot(tranA_postdoc_clean$A1_Temp_C)
  # plot(tranA_postdoc_avg$A2_Temp_C)
  # 
  # plot(tranA_postdoc$A1_WSEL_masl)
  # plot(tranA_postdoc_clean$A1_WSEL_masl)
  # plot(tranA_postdoc_avg$A1_WSEL_masl)
  # 
  # plot(tranA_postdoc$A1_DTGW_m)
  # plot(tranA_postdoc_clean$A1_DTGW_m)
  # plot(tranA_postdoc_avg$A1_DTGW_m)
  
  
  
  # Drop unnecessary variables that were created along the way....
  
  rm(a, data_with_na_outliers,
     tranA_phd_list,tranB_phd_list,tranC_phd_list,
     tranA_postdoc_list,tranB_postdoc_list,tranC_postdoc_list)
  
  
  
  
  
}


# Turn 12 ArkansasR monitoring well locations into sf object

{
  
  phd.site.loc <- matrix(0,ncol=3,nrow=12)
  
  phd.site.loc[1,] = c("A1",-103.6871574,38.0588224)
  phd.site.loc[2,] = c("A2",-103.6850763,38.05829884)
  phd.site.loc[3,] = c("A3",-103.685835,38.05563578)
  phd.site.loc[4,] = c("A4",-103.6845869,38.0548223)
  phd.site.loc[5,] = c("B1",-103.6834083,38.04795136)
  phd.site.loc[6,] = c("B2",-103.6834615,38.04925915)
  phd.site.loc[7,] = c("B3",-103.6829593,38.05026509)
  phd.site.loc[8,] = c("B4",-103.6830133,38.05199954)
  phd.site.loc[9,] = c("C1",-103.6749107,38.05060644)
  phd.site.loc[10,] = c("C2",-103.6725324,38.0522361)
  phd.site.loc[11,] = c("C3",-103.6720018,38.05290651)
  phd.site.loc[12,] = c("C4",-103.6710846,38.05374287)
  
  phd.site.loc <- as.data.frame(phd.site.loc)
  colnames(phd.site.loc) <- c("WELL_SITE_ID","long","lat")
  phd.site.loc$WELL_SITE_ID <- format(phd.site.loc$WELL_SITE_ID,scientific = FALSE)
  
  phd.site.loc$long <- as.numeric(phd.site.loc$long)
  phd.site.loc$lat <- as.numeric(phd.site.loc$lat)
  
  
  rf.gw.locs <- st_as_sf(x = phd.site.loc, 
                         coords = c("long", "lat"),
                         crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
  
}


# Bring in CoAgMET from RFD01 (3288 rows), 01012016 to 12312024
# This data was pre-downloaded from CoAgMET's website

{
  
  
  rfd01.raw <- read.csv(
    "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/from_CoAgMET/arkansas/raw_rfd01_01012016_12312024.csv",
    header = FALSE)
  
  
  rfd01.namez <- paste(unlist(rfd01.raw[1,]),unlist(rfd01.raw[2,]),sep="_")
  
  colnames(rfd01.raw) <- rfd01.namez
  
  rfd01.raw <- rfd01.raw[-c(1,2),]
  
  rfd01.raw$datez <- as.Date(rfd01.raw$Date_date, format = "%m/%d/%Y")
  
  met.data <- rfd01.raw
  
  
}


# Bring in Discharge from CO-DWR (3288 rows)
# This data was pre-downloaded

{
  
  arcrocco.raw <- read.csv(
    "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/from_ColoradoDWR/arkansas/Selected_Station_Observations_Daily_Xtab_202505252012.csv",
    header = TRUE) 
  
  q.data <- arcrocco.raw
  
  # Determine datez from 2016-01-01 to 2024-12-31
  
  datez2 <- seq(as.Date("2016-01-01"), as.Date("2024-12-31"),by="days")
  datez2 <- as.data.frame(datez2)
  colnames(datez2) <- c("datez")
  
  q.data$datez <- datez2$datez
  
  
  
}


# Bring in List of Vegetation Points (Rocky Ford Only)
# Collected during field campaigns using GPSMAP67i model

{
  
  rf.raw.veg.pnt <- read.csv(
    "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/from_GPSMAP67i/rockyford/09272024_points_supp2.csv",header=TRUE)
  
}


# Bring in NDVI for Rocky Ford (ndvi_"instrument type" from 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024)
# NDVI has already been computed, per instrument type, before this script

{

  
  
#
# For Instrument type ps2
# 

{
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2016\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  datez.2016 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2016 <- unique(datez.2016)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2017\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  datez.2017 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2017 <- unique(datez.2017)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2018\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  datez.2018 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2018 <- unique(datez.2018)
  
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2019\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  datez.2019 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2019 <- unique(datez.2019)
  
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2020\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  datez.2020 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2020 <- unique(datez.2020)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2021\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  datez.2021 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2021 <- unique(datez.2021)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2022\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  datez.2022 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2022 <- unique(datez.2022)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2023\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  datez.2023 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2023 <- unique(datez.2023)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2024\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  datez.2024 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2024 <- unique(datez.2024)
  
}


datez.u.all.ps2 <- c(datez.u.2016, datez.u.2017, datez.u.2018, datez.u.2019, datez.u.2020,
                     datez.u.2021, datez.u.2022, datez.u.2023, datez.u.2024)


#
# For Instrument type psbsd
# 

{
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2016\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  datez.2016 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2016 <- unique(datez.2016)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2017\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  datez.2017 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2017 <- unique(datez.2017)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2018\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  datez.2018 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2018 <- unique(datez.2018)
  
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2019\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  datez.2019 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2019 <- unique(datez.2019)
  
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2020\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  datez.2020 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2020 <- unique(datez.2020)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2021\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  datez.2021 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2021 <- unique(datez.2021)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2022\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  datez.2022 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2022 <- unique(datez.2022)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2023\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  datez.2023 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2023 <- unique(datez.2023)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2024\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  datez.2024 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2024 <- unique(datez.2024)
  
}


datez.u.all.psbsd <- c(datez.u.2016, datez.u.2017, datez.u.2018, datez.u.2019, datez.u.2020,
                       datez.u.2021, datez.u.2022, datez.u.2023, datez.u.2024)


#
# For Instrument type ps2sd
# 

{
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2016\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  datez.2016 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2016 <- unique(datez.2016)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2017\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  datez.2017 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2017 <- unique(datez.2017)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2018\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  datez.2018 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2018 <- unique(datez.2018)
  
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2019\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  datez.2019 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2019 <- unique(datez.2019)
  
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2020\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  datez.2020 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2020 <- unique(datez.2020)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2021\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  datez.2021 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2021 <- unique(datez.2021)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2022\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  datez.2022 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2022 <- unique(datez.2022)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2023\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  datez.2023 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2023 <- unique(datez.2023)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2024\\unzip\\PSScene\\"
  setwd(workdir)
  
  #bring in all tif files recursively into R.
  
  namez <- list.files(pattern="ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  datez.2024 <- gsub(".*([0-9]{8}).*", "\\1",namez)
  datez.u.2024 <- unique(datez.2024)
  
}


datez.u.all.ps2sd <- c(datez.u.2016, datez.u.2017, datez.u.2018, datez.u.2019, datez.u.2020,
                       datez.u.2021, datez.u.2022, datez.u.2023, datez.u.2024)





###
### Now bring in all _ndvi_ps2.tif to R environment
###

{
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2016\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  tiflist.16 <- lapply(list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2017\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  tiflist.17 <- lapply(list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2018\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  tiflist.18 <- lapply(list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2019\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  tiflist.19 <- lapply(list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2020\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  tiflist.20 <- lapply(list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2021\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  tiflist.21 <- lapply(list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2022\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  tiflist.22 <- lapply(list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2023\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  tiflist.23 <- lapply(list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2024\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE)
  tiflist.24 <- lapply(list.files(pattern="_ndvi_ps2.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  tiflist.all.ndvi.ps2 <- c(tiflist.16,tiflist.17,tiflist.18,tiflist.19,tiflist.20,
                            tiflist.21,tiflist.22,tiflist.23,tiflist.24)
  
}

###
### Now bring in all _ndvi_psbsd.tif to R environment 
###

{
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2016\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.16 <- lapply(list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2017\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.17 <- lapply(list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2018\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.18 <- lapply(list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2019\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.19 <- lapply(list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2020\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.20 <- lapply(list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2021\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.21 <- lapply(list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2022\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.22 <- lapply(list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2023\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.23 <- lapply(list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2024\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.24 <- lapply(list.files(pattern="_ndvi_psbsd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  tiflist.all.ndvi.psbsd <- c(tiflist.16,tiflist.17,tiflist.18,tiflist.19,tiflist.20,
                              tiflist.21,tiflist.22,tiflist.23,tiflist.24)
  
}

###
### Now bring in all _ndvi_ps2sd.tif to R environment
###

{
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2016\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.16 <- lapply(list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2017\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.17 <- lapply(list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2018\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.18 <- lapply(list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2019\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.19 <- lapply(list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2020\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.20 <- lapply(list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2021\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.21 <- lapply(list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2022\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.22 <- lapply(list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2023\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.23 <- lapply(list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  workdir <- "C:\\Users\\19139\\Desktop\\NSF_Hypoth2_repeat\\Analysis_planet\\rockyford\\2024\\unzip\\PSScene\\"
  setwd(workdir)
  namez <- list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE)
  tiflist.24 <- lapply(list.files(pattern="_ndvi_ps2sd.tif", full.names=TRUE, recursive=TRUE),FUN = rast)
  
  tiflist.all.ndvi.ps2sd <- c(tiflist.16,tiflist.17,tiflist.18,tiflist.19,tiflist.20,
                              tiflist.21,tiflist.22,tiflist.23,tiflist.24)
  
}






}

# Keep track of crs that is given to each raster of a given aoi

rast1 <- rast(tiflist.22[1])
new.crs <- crs(rast1)



# FOR WELL B1 ONLY:
# Create line shapefile across field where well B1 (fallow field) is located in geojson.io
# Import that feature...
b1_line <- st_read("C:/Users/19139/Desktop/NSF_Hypoth2_repeat/from_Geojsonio/B1_transect/POLYLINE.shp")
b1_line_prj <- st_transform(b1_line, crs = crs(new.crs))

# Generate points along the line
# Option A: Specify the number of points (n)
b1_pnts <- st_line_sample(b1_line_prj, n = 10) # Generates 5 points regularly spaced

# Convert back to lon-lat and make it look look other vegetation sample df's
b1_pnts <- st_transform(b1_pnts, crs = crs(b1_line))
b1_pnts_mat <- st_coordinates(b1_pnts)

rf.fl.b1 <- as.data.frame(b1_pnts_mat[,c(1,2)])
colnames(rf.fl.b1) <- c("lon","lat")

rf.fl.b1$gw_well_id <- 5



# Spatially relate my veg points (rf.raw.veg.pnt) to my groundwater locs (rf.gw.locs)

{
  
  rf.veg.locs <- st_as_sf(x = rf.raw.veg.pnt, 
                          coords = c("lon", "lat"),
                          crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
  # Check crs' of each variable
  same.crs(rf.veg.locs,rf.gw.locs) #TRUE
  
  gw.veg.near <- st_nearest_feature(x=rf.veg.locs, y=rf.gw.locs, pairwise=TRUE)
  
}


# Make subsetted groups of points that are associated with a given GW well (based on proximity)

rf.raw.veg.pnt$gw_well_id <- gw.veg.near

{
  
  
  
  
  ###
  ### For lists Cottonwoods at various wells in Rocky Ford
  ###
  
  {
    
    
    # ***********
    #For cottonwood only, near well C1
    # ***********
    
    rf.pf.c1 <- dplyr::filter(rf.raw.veg.pnt, cmt == "", gw_well_id == 9)
    rf.pf.c1 <- rbind(rf.pf.c1, rf.raw.veg.pnt[124,])
    
    #For cottonwood only, near well A4
    
    rf.pf.a4 <- dplyr::filter(rf.raw.veg.pnt, cmt == "", gw_well_id == 4)
    
    #For cottonwood only, near well A3
    
    rf.pf.a3 <- dplyr::filter(rf.raw.veg.pnt, cmt == "", gw_well_id == 3)
    
    #For cottonwood only, near well B4
    
    rf.pf.b4 <- dplyr::filter(rf.raw.veg.pnt, cmt == "", gw_well_id == 8)
    
    #For cottonwood only, near well C4
    
    rf.pf.c4 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Pf_c4", gw_well_id == 12)
    
    #For cottonwood only, near well C2
    
    rf.pf.c2 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Pf_c1c2", gw_well_id == 10)
    
  }
  
  
  ###
  ### For lists of Tamarisk at various wells in Rocky Ford
  ###
  
  
  {
    
    
    # ***********
    #For tamarisk only, near well A2
    # ***********
    
    rf.tm.a2 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Tam", gw_well_id == 2)
    
    #For tamarisk only, near well B2
    
    rf.tm.b2 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Tam", gw_well_id == 6)
    
    #For tamarisk only, near well A1 
    
    rf.tm.a1 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Tam", gw_well_id == 1)
    
    #For tamarisk only, near well A3 
    
    rf.tm.a3 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Tam", gw_well_id == 3)
    
    #For tamarisk only, near well B4 
    
    rf.tm.b4 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Tam", gw_well_id == 8)
    
  }
  
  
  ###
  ### For lists of Russian Olive at various wells in Rocky Ford
  ###
  
  
  {
    
    # ***********
    #For RO only, near well A1
    # ***********
    
    rf.ro.a1 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Ro", gw_well_id == 1)
    
    #For RO only, near well A4
    
    rf.ro.a4 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Ro", gw_well_id == 4)
    
    #For RO only, near well B4
    
    rf.ro.b4 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Ro", gw_well_id == 8)
    
    
  } 
  
  
  ###
  ### For lists of Grasses at various wells in Rocky Ford
  ###
  
  
  {
    
    # ***********
    #For Grasses only, near well A1
    # ***********
    
    rf.pg.a1 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Grs_a1a2", gw_well_id == 1)
    rf.pg.a1 <- rbind(rf.pg.a1,dplyr::filter(rf.raw.veg.pnt, cmt == "Grs_a1a2", gw_well_id == 2))
    rf.pg.a1 <- rbind(rf.pg.a1,dplyr::filter(rf.raw.veg.pnt, cmt == "Grs", gw_well_id == 2))
    
    # ***********
    #For Grasses only, near well C3
    # ***********
    
    rf.pg.c3 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Grs_c3c4", gw_well_id == 11)
    
    
    # ***********
    #For Grasses only, near well C4
    # ***********
    
    rf.pg.c4 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Grs_c3c4", gw_well_id == 12) 
    
    
    # ***********
    #For Grasses only, near well A3
    # ***********
    
    rf.pg.a3 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Grs", gw_well_id == 3) 
    
    
    #For Grasses only, near well B4
    
    rf.pg.b4 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Grs", gw_well_id == 8)
    
    
    #For Wetland-like Grasses only, near well C1
    
    rf.wg.c1 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Wtl grs", gw_well_id == 9) 
    
    
    #For Mixed Forbs only, near well B4
    
    rf.mf.b4 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Mxd frb", gw_well_id == 8)
    
    
    #For Mixed Forbs only, near well A3 (only one point)
    
    rf.mf.a3 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Mxd forb", gw_well_id == 3)
    
    
  }
  
  
  ###
  ### For lists of Willow at various wells in Rocky Ford
  ###
  
  
  {
    
    # ***********
    #For Willow only, near well B3
    # ***********
    
    rf.wl.b3 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Wil", gw_well_id == 7)
    
    
    #For Willow only, near well B2
    
    rf.wl.b2 <- dplyr::filter(rf.raw.veg.pnt, cmt == "Wil", gw_well_id == 6)
    
    
    
  }
  
  
  
  
  
  
  
}




###########################################################################
###########################################################################

# Best examples of each veg type, daily and monthly

###########################################################################
###########################################################################


# rf.pf.c1 for cottonwoods near well c1
# rf.tm.a2 for tamarisk near well A2
# rf.ro.a1 for Russian Olive near well A1
# rf.wl.b3 for Willow near well B3

### GRASSES 
# rf.pg.a1 ***
# rf.pg.a3 *** on Bird Farm Side

# rf.pg.c3 
# rf.pg.c4 *** on Bird Farm Side
# rf.mf.b4

# rf.wg.c1



#  FOR C1, COTTONWOODS

{
  
  # Compute NDVI at Well C1, Cottonwoods, Instrument Type PS2
  
  tiflist.all.ndvi <- tiflist.all.ndvi.ps2
  datez.u.all <- datez.u.all.ps2
  
  # Makes df: "pf.c1.df.rep.ps2", keep repetitions in observations
  
  # rf.pf.c1 
  
  {
    
    
    pf.rf.lst.ps2 <- list()
    
    for(i in 1:nrow(rf.pf.c1)){
      
      a.pnt <- cbind(rf.pf.c1[i,5],rf.pf.c1[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      pf.rf.lst.ps2[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    

    pf.c1.mat.rep.ps2 <- matrix(0, nrow=nrow(pf.rf.lst.ps2[[1]]), ncol = (length(pf.rf.lst.ps2)+1))
    
    for (i in 1:length(pf.rf.lst.ps2)){
      
      a <- as.matrix(pf.rf.lst.ps2[[i]]$value, ncol=1)
      pf.c1.mat.rep.ps2[,i] <- a
      
      if(i == length(pf.rf.lst.ps2)) {
        
        b <- as.matrix(pf.rf.lst.ps2[[i]]$datez.char, ncol=1)
        pf.c1.mat.rep.ps2[,i+1] <- b
        
      }
      
      
    }
    
    
    
    pf.c1.df.rep.ps2 <- as.data.frame(pf.c1.mat.rep.ps2)
    pf.c1.df.rep.ps2$it <- "ps2"
    
    
    
    
  }
  
  
  
  # Compute NDVI at Well C1, Cottonwoods, Instrument Type psbsd
  
  tiflist.all.ndvi <- tiflist.all.ndvi.psbsd
  datez.u.all <- datez.u.all.psbsd
  

  # Makes df: "pf.c1.df.rep.psbsd", keep repetitions in observations
  
  # rf.pf.c1 
  
  {
    
    
    pf.rf.lst.psbsd <- list()
    
    for(i in 1:nrow(rf.pf.c1)){
      
      a.pnt <- cbind(rf.pf.c1[i,5],rf.pf.c1[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      pf.rf.lst.psbsd[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    pf.c1.mat.rep.psbsd <- matrix(0, nrow=nrow(pf.rf.lst.psbsd[[1]]), ncol = (length(pf.rf.lst.psbsd)+1))
    
    for (i in 1:length(pf.rf.lst.psbsd)){
      
      a <- as.matrix(pf.rf.lst.psbsd[[i]]$value, ncol=1)
      pf.c1.mat.rep.psbsd[,i] <- a
      
      if(i == length(pf.rf.lst.psbsd)) {
        
        b <- as.matrix(pf.rf.lst.psbsd[[i]]$datez.char, ncol=1)
        pf.c1.mat.rep.psbsd[,i+1] <- b
        
      }
      
      
    }
    
    
    
    pf.c1.df.rep.psbsd <- as.data.frame(pf.c1.mat.rep.psbsd)
    pf.c1.df.rep.psbsd$it <- "psbsd"
    
    
    
    
  }
  
  
  
  
  # Compute NDVI at Well C1, Cottonwoods, Instrument Type ps2sd
  
  tiflist.all.ndvi <- tiflist.all.ndvi.ps2sd
  datez.u.all <- datez.u.all.ps2sd
  

  # Makes df: "pf.c1.df.rep.ps2sd", keep repetitions in observations
  
  # rf.pf.c1 
  
  {
    
    
    pf.rf.lst.ps2sd <- list()
    
    for(i in 1:nrow(rf.pf.c1)){
      
      a.pnt <- cbind(rf.pf.c1[i,5],rf.pf.c1[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      pf.rf.lst.ps2sd[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    pf.c1.mat.rep.ps2sd <- matrix(0, nrow=nrow(pf.rf.lst.ps2sd[[1]]), ncol = (length(pf.rf.lst.ps2sd)+1))
    
    for (i in 1:length(pf.rf.lst.ps2sd)){
      
      a <- as.matrix(pf.rf.lst.ps2sd[[i]]$value, ncol=1)
      pf.c1.mat.rep.ps2sd[,i] <- a
      
      if(i == length(pf.rf.lst.ps2sd)) {
        
        b <- as.matrix(pf.rf.lst.ps2sd[[i]]$datez.char, ncol=1)
        pf.c1.mat.rep.ps2sd[,i+1] <- b
        
      }
      
      
    }
    
    
    
    pf.c1.df.rep.ps2sd <- as.data.frame(pf.c1.mat.rep.ps2sd)
    pf.c1.df.rep.ps2sd$it <- "ps2sd"
    
    
    
    
  }
  
  
 
  
   
}


#  FOR A2, Tamarisk

{
  
  tiflist.all.ndvi <- tiflist.all.ndvi.ps2
  datez.u.all <- datez.u.all.ps2
  
  # Makes df: "tm.a2.df.rep.ps2", keep repetitions in observations
  
  # rf.tm.a2 
  
  {
    
    
    tm.rf.lst.ps2 <- list()
    
    for(i in 1:nrow(rf.tm.a2)){
      
      a.pnt <- cbind(rf.tm.a2[i,5],rf.tm.a2[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      tm.rf.lst.ps2[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    tm.a2.mat.rep.ps2 <- matrix(0, nrow=nrow(tm.rf.lst.ps2[[1]]), ncol = (length(tm.rf.lst.ps2)+1))
    
    for (i in 1:length(tm.rf.lst.ps2)){
      
      a <- as.matrix(tm.rf.lst.ps2[[i]]$value, ncol=1)
      tm.a2.mat.rep.ps2[,i] <- a
      
      if(i == length(tm.rf.lst.ps2)) {
        
        b <- as.matrix(tm.rf.lst.ps2[[i]]$datez.char, ncol=1)
        tm.a2.mat.rep.ps2[,i+1] <- b
        
      }
      
      
    }
    
    
    
    tm.a2.df.rep.ps2 <- as.data.frame(tm.a2.mat.rep.ps2)
    tm.a2.df.rep.ps2$it <- "ps2"
    
    
    
    
  }
  
  
  

  tiflist.all.ndvi <- tiflist.all.ndvi.psbsd
  datez.u.all <- datez.u.all.psbsd
  
  
  # Makes df: "tm.a2.df.rep.psbsd", keep repetitions in observations
  
  # rf.tm.a2 
  
  {
    
    
    tm.rf.lst.psbsd <- list()
    
    for(i in 1:nrow(rf.tm.a2)){
      
      a.pnt <- cbind(rf.tm.a2[i,5],rf.tm.a2[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      tm.rf.lst.psbsd[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    tm.a2.mat.rep.psbsd <- matrix(0, nrow=nrow(tm.rf.lst.psbsd[[1]]), ncol = (length(tm.rf.lst.psbsd)+1))
    
    for (i in 1:length(tm.rf.lst.psbsd)){
      
      a <- as.matrix(tm.rf.lst.psbsd[[i]]$value, ncol=1)
      tm.a2.mat.rep.psbsd[,i] <- a
      
      if(i == length(tm.rf.lst.psbsd)) {
        
        b <- as.matrix(tm.rf.lst.psbsd[[i]]$datez.char, ncol=1)
        tm.a2.mat.rep.psbsd[,i+1] <- b
        
      }
      
      
    }
    
    
    
    tm.a2.df.rep.psbsd <- as.data.frame(tm.a2.mat.rep.psbsd)
    tm.a2.df.rep.psbsd$it <- "psbsd"
    
    
    
    
  }
  
  
  
  
  tiflist.all.ndvi <- tiflist.all.ndvi.ps2sd
  datez.u.all <- datez.u.all.ps2sd
  
  
  # Makes df: "tm.a2.df.rep.ps2sd", keep repetitions in observations
  
  # rf.tm.a2 
  
  {
    
    
    tm.rf.lst.ps2sd <- list()
    
    for(i in 1:nrow(rf.tm.a2)){
      
      a.pnt <- cbind(rf.tm.a2[i,5],rf.tm.a2[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      tm.rf.lst.ps2sd[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    tm.a2.mat.rep.ps2sd <- matrix(0, nrow=nrow(tm.rf.lst.ps2sd[[1]]), ncol = (length(tm.rf.lst.ps2sd)+1))
    
    for (i in 1:length(tm.rf.lst.ps2sd)){
      
      a <- as.matrix(tm.rf.lst.ps2sd[[i]]$value, ncol=1)
      tm.a2.mat.rep.ps2sd[,i] <- a
      
      if(i == length(tm.rf.lst.ps2sd)) {
        
        b <- as.matrix(tm.rf.lst.ps2sd[[i]]$datez.char, ncol=1)
        tm.a2.mat.rep.ps2sd[,i+1] <- b
        
      }
      
      
    }
    
    
    
    tm.a2.df.rep.ps2sd <- as.data.frame(tm.a2.mat.rep.ps2sd)
    tm.a2.df.rep.ps2sd$it <- "ps2sd"
    
    
    
    
  }
  
  
  
  
  
}


#  FOR A1, Russian Olive

{
  
  tiflist.all.ndvi <- tiflist.all.ndvi.ps2
  datez.u.all <- datez.u.all.ps2
  
  # Makes df: "ro.a1.df.rep.ps2", keep repetitions in observations
  
  # rf.ro.a1 
  
  {
    
    
    ro.rf.lst.ps2 <- list()
    
    for(i in 1:nrow(rf.ro.a1)){
      
      a.pnt <- cbind(rf.ro.a1[i,5],rf.ro.a1[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      ro.rf.lst.ps2[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    ro.a1.mat.rep.ps2 <- matrix(0, nrow=nrow(ro.rf.lst.ps2[[1]]), ncol = (length(ro.rf.lst.ps2)+1))
    
    for (i in 1:length(ro.rf.lst.ps2)){
      
      a <- as.matrix(ro.rf.lst.ps2[[i]]$value, ncol=1)
      ro.a1.mat.rep.ps2[,i] <- a
      
      if(i == length(ro.rf.lst.ps2)) {
        
        b <- as.matrix(ro.rf.lst.ps2[[i]]$datez.char, ncol=1)
        ro.a1.mat.rep.ps2[,i+1] <- b
        
      }
      
      
    }
    
    
    
    ro.a1.df.rep.ps2 <- as.data.frame(ro.a1.mat.rep.ps2)
    ro.a1.df.rep.ps2$it <- "ps2"
    
    
    
    
  }
  
  
  
  
  tiflist.all.ndvi <- tiflist.all.ndvi.psbsd
  datez.u.all <- datez.u.all.psbsd
  
  
  # Makes df: "ro.a1.df.rep.psbsd", keep repetitions in observations
  
  # rf.ro.a1 
  
  {
    
    
    ro.rf.lst.psbsd <- list()
    
    for(i in 1:nrow(rf.ro.a1)){
      
      a.pnt <- cbind(rf.ro.a1[i,5],rf.ro.a1[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      ro.rf.lst.psbsd[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    ro.a1.mat.rep.psbsd <- matrix(0, nrow=nrow(ro.rf.lst.psbsd[[1]]), ncol = (length(ro.rf.lst.psbsd)+1))
    
    for (i in 1:length(ro.rf.lst.psbsd)){
      
      a <- as.matrix(ro.rf.lst.psbsd[[i]]$value, ncol=1)
      ro.a1.mat.rep.psbsd[,i] <- a
      
      if(i == length(ro.rf.lst.psbsd)) {
        
        b <- as.matrix(ro.rf.lst.psbsd[[i]]$datez.char, ncol=1)
        ro.a1.mat.rep.psbsd[,i+1] <- b
        
      }
      
      
    }
    
    
    
    ro.a1.df.rep.psbsd <- as.data.frame(ro.a1.mat.rep.psbsd)
    ro.a1.df.rep.psbsd$it <- "psbsd"
    
    
    
    
  }
  
  
  
  
  tiflist.all.ndvi <- tiflist.all.ndvi.ps2sd
  datez.u.all <- datez.u.all.ps2sd
  
  
  # Makes df: "ro.a1.df.rep.ps2sd", keep repetitions in observations
  
  # rf.ro.a1 
  
  {
    
    
    ro.rf.lst.ps2sd <- list()
    
    for(i in 1:nrow(rf.ro.a1)){
      
      a.pnt <- cbind(rf.ro.a1[i,5],rf.ro.a1[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      ro.rf.lst.ps2sd[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    ro.a1.mat.rep.ps2sd <- matrix(0, nrow=nrow(ro.rf.lst.ps2sd[[1]]), ncol = (length(ro.rf.lst.ps2sd)+1))
    
    for (i in 1:length(ro.rf.lst.ps2sd)){
      
      a <- as.matrix(ro.rf.lst.ps2sd[[i]]$value, ncol=1)
      ro.a1.mat.rep.ps2sd[,i] <- a
      
      if(i == length(ro.rf.lst.ps2sd)) {
        
        b <- as.matrix(ro.rf.lst.ps2sd[[i]]$datez.char, ncol=1)
        ro.a1.mat.rep.ps2sd[,i+1] <- b
        
      }
      
      
    }
    
    
    
    ro.a1.df.rep.ps2sd <- as.data.frame(ro.a1.mat.rep.ps2sd)
    ro.a1.df.rep.ps2sd$it <- "ps2sd"
    
    
    
    
  }
  
  
  
  
  
}


#  FOR B3, Willow

{
  
  tiflist.all.ndvi <- tiflist.all.ndvi.ps2
  datez.u.all <- datez.u.all.ps2
  
  # Makes df: "wl.b3.df.rep.ps2", keep repetitions in observations
  
  # rf.wl.b3 
  
  {
    
    
    wl.rf.lst.ps2 <- list()
    
    for(i in 1:nrow(rf.wl.b3)){
      
      a.pnt <- cbind(rf.wl.b3[i,5],rf.wl.b3[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      wl.rf.lst.ps2[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    wl.b3.mat.rep.ps2 <- matrix(0, nrow=nrow(wl.rf.lst.ps2[[1]]), ncol = (length(wl.rf.lst.ps2)+1))
    
    for (i in 1:length(wl.rf.lst.ps2)){
      
      a <- as.matrix(wl.rf.lst.ps2[[i]]$value, ncol=1)
      wl.b3.mat.rep.ps2[,i] <- a
      
      if(i == length(wl.rf.lst.ps2)) {
        
        b <- as.matrix(wl.rf.lst.ps2[[i]]$datez.char, ncol=1)
        wl.b3.mat.rep.ps2[,i+1] <- b
        
      }
      
      
    }
    
    
    
    wl.b3.df.rep.ps2 <- as.data.frame(wl.b3.mat.rep.ps2)
    wl.b3.df.rep.ps2$it <- "ps2"
    
    
    
    
  }
  
  
  
  
  tiflist.all.ndvi <- tiflist.all.ndvi.psbsd
  datez.u.all <- datez.u.all.psbsd
  
  
  # Makes df: "wl.b3.df.rep.psbsd", keep repetitions in observations
  
  # rf.wl.b3 
  
  {
    
    
    wl.rf.lst.psbsd <- list()
    
    for(i in 1:nrow(rf.wl.b3)){
      
      a.pnt <- cbind(rf.wl.b3[i,5],rf.wl.b3[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      wl.rf.lst.psbsd[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    wl.b3.mat.rep.psbsd <- matrix(0, nrow=nrow(wl.rf.lst.psbsd[[1]]), ncol = (length(wl.rf.lst.psbsd)+1))
    
    for (i in 1:length(wl.rf.lst.psbsd)){
      
      a <- as.matrix(wl.rf.lst.psbsd[[i]]$value, ncol=1)
      wl.b3.mat.rep.psbsd[,i] <- a
      
      if(i == length(wl.rf.lst.psbsd)) {
        
        b <- as.matrix(wl.rf.lst.psbsd[[i]]$datez.char, ncol=1)
        wl.b3.mat.rep.psbsd[,i+1] <- b
        
      }
      
      
    }
    
    
    
    wl.b3.df.rep.psbsd <- as.data.frame(wl.b3.mat.rep.psbsd)
    wl.b3.df.rep.psbsd$it <- "psbsd"
    
    
    
    
  }
  
  
  
  
  tiflist.all.ndvi <- tiflist.all.ndvi.ps2sd
  datez.u.all <- datez.u.all.ps2sd
  
  
  # Makes df: "wl.b3.df.rep.ps2sd", keep repetitions in observations
  
  # rf.wl.b3 
  
  {
    
    
    wl.rf.lst.ps2sd <- list()
    
    for(i in 1:nrow(rf.wl.b3)){
      
      a.pnt <- cbind(rf.wl.b3[i,5],rf.wl.b3[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      wl.rf.lst.ps2sd[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    wl.b3.mat.rep.ps2sd <- matrix(0, nrow=nrow(wl.rf.lst.ps2sd[[1]]), ncol = (length(wl.rf.lst.ps2sd)+1))
    
    for (i in 1:length(wl.rf.lst.ps2sd)){
      
      a <- as.matrix(wl.rf.lst.ps2sd[[i]]$value, ncol=1)
      wl.b3.mat.rep.ps2sd[,i] <- a
      
      if(i == length(wl.rf.lst.ps2sd)) {
        
        b <- as.matrix(wl.rf.lst.ps2sd[[i]]$datez.char, ncol=1)
        wl.b3.mat.rep.ps2sd[,i+1] <- b
        
      }
      
      
    }
    
    
    
    wl.b3.df.rep.ps2sd <- as.data.frame(wl.b3.mat.rep.ps2sd)
    wl.b3.df.rep.ps2sd$it <- "ps2sd"
    
    
    
    
  }
  
  
  
  
  
}


#  FOR A3, Prairie Grass, closest to wetland style grass

{
  
  tiflist.all.ndvi <- tiflist.all.ndvi.ps2
  datez.u.all <- datez.u.all.ps2
  
  # Makes df: "pga3.a3.df.rep.ps2", keep repetitions in observations
  
  # rf.pg.a3 
  
  {
    
    
    pga3.rf.lst.ps2 <- list()
    
    for(i in 1:nrow(rf.pg.a3)){
      
      a.pnt <- cbind(rf.pg.a3[i,5],rf.pg.a3[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      pga3.rf.lst.ps2[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    pga3.a3.mat.rep.ps2 <- matrix(0, nrow=nrow(pga3.rf.lst.ps2[[1]]), ncol = (length(pga3.rf.lst.ps2)+1))
    
    for (i in 1:length(pga3.rf.lst.ps2)){
      
      a <- as.matrix(pga3.rf.lst.ps2[[i]]$value, ncol=1)
      pga3.a3.mat.rep.ps2[,i] <- a
      
      if(i == length(pga3.rf.lst.ps2)) {
        
        b <- as.matrix(pga3.rf.lst.ps2[[i]]$datez.char, ncol=1)
        pga3.a3.mat.rep.ps2[,i+1] <- b
        
      }
      
      
    }
    
    
    
    pga3.a3.df.rep.ps2 <- as.data.frame(pga3.a3.mat.rep.ps2)
    pga3.a3.df.rep.ps2$it <- "ps2"
    
    
    
    
  }
  
  
  
  
  tiflist.all.ndvi <- tiflist.all.ndvi.psbsd
  datez.u.all <- datez.u.all.psbsd
  
  
  # Makes df: "pga3.a3.df.rep.psbsd", keep repetitions in observations
  
  # rf.pg.a3 
  
  {
    
    
    pga3.rf.lst.psbsd <- list()
    
    for(i in 1:nrow(rf.pg.a3)){
      
      a.pnt <- cbind(rf.pg.a3[i,5],rf.pg.a3[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      pga3.rf.lst.psbsd[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    pga3.a3.mat.rep.psbsd <- matrix(0, nrow=nrow(pga3.rf.lst.psbsd[[1]]), ncol = (length(pga3.rf.lst.psbsd)+1))
    
    for (i in 1:length(pga3.rf.lst.psbsd)){
      
      a <- as.matrix(pga3.rf.lst.psbsd[[i]]$value, ncol=1)
      pga3.a3.mat.rep.psbsd[,i] <- a
      
      if(i == length(pga3.rf.lst.psbsd)) {
        
        b <- as.matrix(pga3.rf.lst.psbsd[[i]]$datez.char, ncol=1)
        pga3.a3.mat.rep.psbsd[,i+1] <- b
        
      }
      
      
    }
    
    
    
    pga3.a3.df.rep.psbsd <- as.data.frame(pga3.a3.mat.rep.psbsd)
    pga3.a3.df.rep.psbsd$it <- "psbsd"
    
    
    
    
  }
  
  
  
  
  tiflist.all.ndvi <- tiflist.all.ndvi.ps2sd
  datez.u.all <- datez.u.all.ps2sd
  
  
  # Makes df: "pga3.a3.df.rep.ps2sd", keep repetitions in observations
  
  # rf.pg.a3 
  
  {
    
    
    pga3.rf.lst.ps2sd <- list()
    
    for(i in 1:nrow(rf.pg.a3)){
      
      a.pnt <- cbind(rf.pg.a3[i,5],rf.pg.a3[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      pga3.rf.lst.ps2sd[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    pga3.a3.mat.rep.ps2sd <- matrix(0, nrow=nrow(pga3.rf.lst.ps2sd[[1]]), ncol = (length(pga3.rf.lst.ps2sd)+1))
    
    for (i in 1:length(pga3.rf.lst.ps2sd)){
      
      a <- as.matrix(pga3.rf.lst.ps2sd[[i]]$value, ncol=1)
      pga3.a3.mat.rep.ps2sd[,i] <- a
      
      if(i == length(pga3.rf.lst.ps2sd)) {
        
        b <- as.matrix(pga3.rf.lst.ps2sd[[i]]$datez.char, ncol=1)
        pga3.a3.mat.rep.ps2sd[,i+1] <- b
        
      }
      
      
    }
    
    
    
    pga3.a3.df.rep.ps2sd <- as.data.frame(pga3.a3.mat.rep.ps2sd)
    pga3.a3.df.rep.ps2sd$it <- "ps2sd"
    
    
    
    
  }
  
  
  
  
  
}


#  FOR C4, Prairie Grass, more like mixed forbs

{
  
  tiflist.all.ndvi <- tiflist.all.ndvi.ps2
  datez.u.all <- datez.u.all.ps2
  
  # Makes df: "pgc4.c4.df.rep.ps2", keep repetitions in observations
  
  # rf.pg.c4 
  
  {
    
    
    pgc4.rf.lst.ps2 <- list()
    
    for(i in 1:nrow(rf.pg.c4)){
      
      a.pnt <- cbind(rf.pg.c4[i,5],rf.pg.c4[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      pgc4.rf.lst.ps2[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    pgc4.c4.mat.rep.ps2 <- matrix(0, nrow=nrow(pgc4.rf.lst.ps2[[1]]), ncol = (length(pgc4.rf.lst.ps2)+1))
    
    for (i in 1:length(pgc4.rf.lst.ps2)){
      
      a <- as.matrix(pgc4.rf.lst.ps2[[i]]$value, ncol=1)
      pgc4.c4.mat.rep.ps2[,i] <- a
      
      if(i == length(pgc4.rf.lst.ps2)) {
        
        b <- as.matrix(pgc4.rf.lst.ps2[[i]]$datez.char, ncol=1)
        pgc4.c4.mat.rep.ps2[,i+1] <- b
        
      }
      
      
    }
    
    
    
    pgc4.c4.df.rep.ps2 <- as.data.frame(pgc4.c4.mat.rep.ps2)
    pgc4.c4.df.rep.ps2$it <- "ps2"
    
    
    
    
  }
  
  
  
  
  tiflist.all.ndvi <- tiflist.all.ndvi.psbsd
  datez.u.all <- datez.u.all.psbsd
  
  
  # Makes df: "pgc4.c4.df.rep.psbsd", keep repetitions in observations
  
  # rf.pg.c4 
  
  {
    
    
    pgc4.rf.lst.psbsd <- list()
    
    for(i in 1:nrow(rf.pg.c4)){
      
      a.pnt <- cbind(rf.pg.c4[i,5],rf.pg.c4[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      pgc4.rf.lst.psbsd[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    pgc4.c4.mat.rep.psbsd <- matrix(0, nrow=nrow(pgc4.rf.lst.psbsd[[1]]), ncol = (length(pgc4.rf.lst.psbsd)+1))
    
    for (i in 1:length(pgc4.rf.lst.psbsd)){
      
      a <- as.matrix(pgc4.rf.lst.psbsd[[i]]$value, ncol=1)
      pgc4.c4.mat.rep.psbsd[,i] <- a
      
      if(i == length(pgc4.rf.lst.psbsd)) {
        
        b <- as.matrix(pgc4.rf.lst.psbsd[[i]]$datez.char, ncol=1)
        pgc4.c4.mat.rep.psbsd[,i+1] <- b
        
      }
      
      
    }
    
    
    
    pgc4.c4.df.rep.psbsd <- as.data.frame(pgc4.c4.mat.rep.psbsd)
    pgc4.c4.df.rep.psbsd$it <- "psbsd"
    
    
    
    
  }
  
  
  
  
  tiflist.all.ndvi <- tiflist.all.ndvi.ps2sd
  datez.u.all <- datez.u.all.ps2sd
  
  
  # Makes df: "pgc4.c4.df.rep.ps2sd", keep repetitions in observations
  
  # rf.pg.c4 
  
  {
    
    
    pgc4.rf.lst.ps2sd <- list()
    
    for(i in 1:nrow(rf.pg.c4)){
      
      a.pnt <- cbind(rf.pg.c4[i,5],rf.pg.c4[i,6]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      pgc4.rf.lst.ps2sd[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    pgc4.c4.mat.rep.ps2sd <- matrix(0, nrow=nrow(pgc4.rf.lst.ps2sd[[1]]), ncol = (length(pgc4.rf.lst.ps2sd)+1))
    
    for (i in 1:length(pgc4.rf.lst.ps2sd)){
      
      a <- as.matrix(pgc4.rf.lst.ps2sd[[i]]$value, ncol=1)
      pgc4.c4.mat.rep.ps2sd[,i] <- a
      
      if(i == length(pgc4.rf.lst.ps2sd)) {
        
        b <- as.matrix(pgc4.rf.lst.ps2sd[[i]]$datez.char, ncol=1)
        pgc4.c4.mat.rep.ps2sd[,i+1] <- b
        
      }
      
      
    }
    
    
    
    pgc4.c4.df.rep.ps2sd <- as.data.frame(pgc4.c4.mat.rep.ps2sd)
    pgc4.c4.df.rep.ps2sd$it <- "ps2sd"
    
    
    
    
  }
  
  
  
  
  
}


#  FOR B1, Fallow Land, more like dead/dying mixed forbs

{
  
  tiflist.all.ndvi <- tiflist.all.ndvi.ps2
  datez.u.all <- datez.u.all.ps2
  
  # Makes df: "fl.b1.df.rep.ps2", keep repetitions in observations
  
  # rf.fl.b1 
  
  {
    
    
    fl.rf.lst.ps2 <- list()
    
    for(i in 1:nrow(rf.fl.b1)){
      
      a.pnt <- cbind(rf.fl.b1[i,1],rf.fl.b1[i,2]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      fl.rf.lst.ps2[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    fl.b1.mat.rep.ps2 <- matrix(0, nrow=nrow(fl.rf.lst.ps2[[1]]), ncol = (length(fl.rf.lst.ps2)+1))
    
    for (i in 1:length(fl.rf.lst.ps2)){
      
      a <- as.matrix(fl.rf.lst.ps2[[i]]$value, ncol=1)
      fl.b1.mat.rep.ps2[,i] <- a
      
      if(i == length(fl.rf.lst.ps2)) {
        
        b <- as.matrix(fl.rf.lst.ps2[[i]]$datez.char, ncol=1)
        fl.b1.mat.rep.ps2[,i+1] <- b
        
      }
      
      
    }
    
    
    
    fl.b1.df.rep.ps2 <- as.data.frame(fl.b1.mat.rep.ps2)
    fl.b1.df.rep.ps2$it <- "ps2"
    
    
    
    
  }
  
  
  
  tiflist.all.ndvi <- tiflist.all.ndvi.psbsd
  datez.u.all <- datez.u.all.psbsd
  
  
  # Makes df: "fl.b1.df.rep.psbsd", keep repetitions in observations
  
  # rf.fl.b1 
  
  {
    
    
    fl.rf.lst.psbsd <- list()
    
    for(i in 1:nrow(rf.fl.b1)){
      
      a.pnt <- cbind(rf.fl.b1[i,1],rf.fl.b1[i,2]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      fl.rf.lst.psbsd[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    fl.b1.mat.rep.psbsd <- matrix(0, nrow=nrow(fl.rf.lst.psbsd[[1]]), ncol = (length(fl.rf.lst.psbsd)+1))
    
    for (i in 1:length(fl.rf.lst.psbsd)){
      
      a <- as.matrix(fl.rf.lst.psbsd[[i]]$value, ncol=1)
      fl.b1.mat.rep.psbsd[,i] <- a
      
      if(i == length(fl.rf.lst.psbsd)) {
        
        b <- as.matrix(fl.rf.lst.psbsd[[i]]$datez.char, ncol=1)
        fl.b1.mat.rep.psbsd[,i+1] <- b
        
      }
      
      
    }
    
    
    
    fl.b1.df.rep.psbsd <- as.data.frame(fl.b1.mat.rep.psbsd)
    fl.b1.df.rep.psbsd$it <- "psbsd"
    
    
    
    
  }
  
  

  
  tiflist.all.ndvi <- tiflist.all.ndvi.ps2sd
  datez.u.all <- datez.u.all.ps2sd
  
  
  # Makes df: "fl.b1.df.rep.ps2sd", keep repetitions in observations
  
  # rf.fl.b1 
  
  {
    
    
    fl.rf.lst.ps2sd <- list()
    
    for(i in 1:nrow(rf.fl.b1)){
      
      a.pnt <- cbind(rf.fl.b1[i,1],rf.fl.b1[i,2]) |> vect(crs="+proj=longlat")
      
      # point sample the same crs as the raster
      
      a.pnt.prj <- crds(project(a.pnt, new.crs))
      
      res <- 6 # This resolution grabs 4 pixels with the sample point in center (i.e., a "4 square")
      
      #Create a single cell around that point and expand (if needed):
      b <- rep(a.pnt.prj, each=2) + c(-res, res) / 2 
      grid1 <- rast(ext(b), crs=crs(new.crs), ncol=1, nrow=1)
      
      grid2 <- as.polygons(grid1)
      
      # check to see if point and polygon make sense
      # plot(grid2)
      # points(a.pnt.prj)
      
      
      c1.lst <- list()
      
      for(j in 1:length(tiflist.all.ndvi)){
        
        aa <- rast(tiflist.all.ndvi[j])
        bb <- try(crop(x=aa, y=grid2), silent = TRUE)
        
        if("try-error" %in% class(bb)){
          bb <- rast(matrix(NA,nrow=res/3,ncol=res/3),crs=crs(new.crs)) # edit on 9/18/25
        }
        
        c1.lst[[j]] <- bb
        
      }
      
      
      
      c1.d.all <- matrix(0, nrow = ((res/3) * (res/3)),ncol=length(c1.lst))
      for(k in 1:length(c1.lst)){
        skip_to_next <- FALSE
        
        tryCatch({
          
          aaa <- as.matrix(c1.lst[[k]])
          bbb <- c(aaa)
          # bbb[bbb < 0] <- NA # edit on 9/18/25
          
          c1.d.all[,k] <- bbb},
          error = function(e) { skip_to_next <- TRUE})
        
        if(skip_to_next) { next }  
        
      }
      
      
      c1.df <- as.data.frame(gather(as.data.frame(c1.d.all)))
      planet.datez <- cbind.data.frame(datez.u.all)
      planet.datez$datez <- as.Date(planet.datez$datez.u.all, "%Y%m%d")
      c1.df$datez.char <- rep(datez.u.all, each = nrow(c1.d.all))
      c1.df$datez <- as.Date(c1.df$datez.char, "%Y%m%d") 
      
      c1.df.ndvi.6 <- c1.df
      #change 0's and NaN from planet data into NA's
      c1.df.ndvi.6$value[c1.df.ndvi.6$value == 0] <- NA
      c1.df.ndvi.6$value[is.nan(c1.df.ndvi.6$value)] <- NA 
      
      c1.df.ndvi.6$monthz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y%m%d"),"%m")
      c1.df.ndvi.6$yearz <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y")
      c1.df.ndvi.6$yr_mo <- format(as.Date(c1.df.ndvi.6$datez, format="%Y-%m-%d"),"%Y-%m")
      
      
      fl.rf.lst.ps2sd[[i]] <- c1.df.ndvi.6
      
      
      
      
      
    }
    
    
    
    fl.b1.mat.rep.ps2sd <- matrix(0, nrow=nrow(fl.rf.lst.ps2sd[[1]]), ncol = (length(fl.rf.lst.ps2sd)+1))
    
    for (i in 1:length(fl.rf.lst.ps2sd)){
      
      a <- as.matrix(fl.rf.lst.ps2sd[[i]]$value, ncol=1)
      fl.b1.mat.rep.ps2sd[,i] <- a
      
      if(i == length(fl.rf.lst.ps2sd)) {
        
        b <- as.matrix(fl.rf.lst.ps2sd[[i]]$datez.char, ncol=1)
        fl.b1.mat.rep.ps2sd[,i+1] <- b
        
      }
      
      
    }
    
    
    
    fl.b1.df.rep.ps2sd <- as.data.frame(fl.b1.mat.rep.ps2sd)
    fl.b1.df.rep.ps2sd$it <- "ps2sd"
    
    
    
    
  }
  
  
  
  
  
}







# Plot to be sure that NDVI looks good  

{

check1 <- rbind.data.frame(pf.c1.df.rep.ps2,pf.c1.df.rep.psbsd,pf.c1.df.rep.ps2sd)

check1$datez <- as.Date(check1$V34, format = "%Y%m%d" ) 
check1$it <- as.factor(check1$it)
  
check1_long <- check1 %>%
  pivot_longer(
    cols = -c(it,datez,V34),
    names_to = "locs_pnts",
    values_to = "ndvi"
  )
  
check1_long$ndvi <- as.numeric(check1_long$ndvi)  

ggplot(check1_long) + geom_point(aes(x=datez,y=ndvi,color=it, shape = locs_pnts)) + 
  xlab("Date (d)") + ylab("NDVI") + labs(color="Instrument \nType") + scale_color_viridis_d()


}



### Things I need to export...

### All my reponse variables...

pf_c1_df_rep_all <- rbind.data.frame(pf.c1.df.rep.ps2,pf.c1.df.rep.psbsd,pf.c1.df.rep.ps2sd)

tm_a2_df_rep_all <- rbind.data.frame(tm.a2.df.rep.ps2,tm.a2.df.rep.psbsd,tm.a2.df.rep.ps2sd)

ro_a1_df_rep_all <- rbind.data.frame(ro.a1.df.rep.ps2,ro.a1.df.rep.psbsd,ro.a1.df.rep.ps2sd)

wl_b3_df_rep_all <- rbind.data.frame(wl.b3.df.rep.ps2,wl.b3.df.rep.psbsd,wl.b3.df.rep.ps2sd)

pga3_a3_df_rep_all <- rbind.data.frame(pga3.a3.df.rep.ps2,pga3.a3.df.rep.psbsd,pga3.a3.df.rep.ps2sd)

pgc4_c4_df_rep_all <- rbind.data.frame(pgc4.c4.df.rep.ps2,pgc4.c4.df.rep.psbsd,pgc4.c4.df.rep.ps2sd)

fl_b1_df_rep_all <- rbind.data.frame(fl.b1.df.rep.ps2,fl.b1.df.rep.psbsd,fl.b1.df.rep.ps2sd)

### All my covariates...

# tranA_phd_avg, tranB_phd_avg, tranC_phd_avg,
# tranA_postdoc_avg, tranB_postdoc_avg, tranC_postdoc_avg,
# q.data, met.data

  
  
  # EXPORT covariates and y dataset seperately
  
  # EXPORT EXPORT EXPORT!!!
  
  write.csv(pf_c1_df_rep_all, 
            "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/pf_c1_df_rep_all.csv")
  write.csv(tm_a2_df_rep_all, 
            "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/tm_a2_df_rep_all.csv")
  write.csv(ro_a1_df_rep_all, 
            "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/ro_a1_df_rep_all.csv") 
  write.csv(wl_b3_df_rep_all, 
            "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/wl_b3_df_rep_all.csv") 
  write.csv(pga3_a3_df_rep_all, 
            "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/pga3_a3_df_rep_all.csv") 
  write.csv(pgc4_c4_df_rep_all, 
            "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/pgc4_c4_df_rep_all.csv") 
  write.csv(fl_b1_df_rep_all, 
            "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/fl_b1_df_rep_all.csv") 
  
  
  # tranA_phd_avg, tranB_phd_avg, tranC_phd_avg,
  # tranA_postdoc_avg, tranB_postdoc_avg, tranC_postdoc_avg,
  # q.data, met.data
  
  write.csv(tranA_phd_avg, 
            "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/tranA_phd_avg.csv")
  write.csv(tranB_phd_avg, 
            "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/tranB_phd_avg.csv")
  write.csv(tranC_phd_avg, 
            "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/tranC_phd_avg.csv")

  write.csv(tranA_postdoc_avg, 
            "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/tranA_postdoc_avg.csv")
  write.csv(tranB_postdoc_avg, 
            "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/tranB_postdoc_avg.csv")
  write.csv(tranC_postdoc_avg, 
            "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/tranC_postdoc_avg.csv")


  write.csv(q.data, 
            "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/q_data.csv")
  write.csv(met.data, 
            "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/met_data.csv")





