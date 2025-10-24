###
### 10172025
### Use imputeTS (if necessary) and Rssa for SSA and (i.e., separate signal and noise)
### and then run RF on all the responses and covars
###


# For response-data manipulation
library(imputeTS)
library(Rssa)

# Split-apply-combine and plot packages
library(tidyverse)

#plotting
library(gridExtra)

# RF Packages for regression
library(vip) # for variable importance plots
library(rsample) # for resampling procedures
library(ranger) # a c++ implementation of random forest

###
###
### Here are the Y's, pfc1 only (standardized), no repetitions, daily means
###
###




pf_c1_df_rep_all <- read.csv( 
  "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/pf_c1_df_rep_all.csv",header = TRUE)

pg_c4_df_rep_all <- read.csv( 
  "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/pgc4_c4_df_rep_all.csv",header = TRUE)

wl_b3_df_rep_all <- read.csv( 
  "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/wl_b3_df_rep_all.csv",header = TRUE)

tm_a2_df_rep_all <- read.csv( 
  "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/tm_a2_df_rep_all.csv",header = TRUE)

ro_a1_df_rep_all <- read.csv( 
  "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/ro_a1_df_rep_all.csv",header = TRUE)

pg_a3_df_rep_all <- read.csv( 
  "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/pga3_a3_df_rep_all.csv",header = TRUE)

fl_b1_df_rep_all <- read.csv( 
  "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/fl_b1_df_rep_all.csv",header = TRUE)


all_raw_y <- list(pf_c1_df_rep_all, pg_c4_df_rep_all, wl_b3_df_rep_all, tm_a2_df_rep_all,
                  ro_a1_df_rep_all,  pg_a3_df_rep_all, fl_b1_df_rep_all)

veg.name <- c("Populus deltoides", "Mixed Forbs", "Narrowleaf Willow", "Tamarix ramnossima",
              "Russian Olive", "Prairie Grass", "Fallow")

dy_std_y <- list()

for(i in 1:length(all_raw_y)){
  
  pf_c1_df_rep_all <- all_raw_y[[i]]
  pf_c1_df_rep_all$datez.char <- as.character(pf_c1_df_rep_all[,c(colnames(pf_c1_df_rep_all)[ncol(pf_c1_df_rep_all) - 1])])
  pf_c1_df_rep_all$datez <- as.Date(pf_c1_df_rep_all$datez.char, format = "%Y%m%d") 
  pf_c1_df_rep_all$it <- as.factor(pf_c1_df_rep_all$it)
  
  colstokeep <- c("X","it","datez","datez.char",colnames(pf_c1_df_rep_all)[ncol(pf_c1_df_rep_all) - 3])
  
  pf_c1_df_rep_all_long <- pf_c1_df_rep_all %>%
    pivot_longer(
      cols = -c("X","it","datez","datez.char",colnames(pf_c1_df_rep_all)[ncol(pf_c1_df_rep_all) - 3]),
      names_to = "locs_pnts",
      values_to = "ndvi"
    )
  
  pf_c1_df_rep_all_long$ndvi <- as.numeric(pf_c1_df_rep_all_long$ndvi)  
  
  pfc1 <- pf_c1_df_rep_all_long[,-c(1)]
  
  ### Take mean of NDVI by date && instrument type:
  ### This reduces the 4 repetitions in the observation down to the average obs per day,
  ### per instrument type
  
  pfc1 <- pfc1 %>%
    group_by(it, datez) %>%
    summarize(
      ndvi_avg = mean(ndvi, na.rm = TRUE)
    )
  
  pfc1 <- pfc1 %>% mutate_all(~ifelse(is.nan(.), NA, .))
  
  # A question remains here, do I standardize before or after I take the daily means for ps2 and psbsd overlaps?
  
  pfc1_ps2 <- pfc1 %>% filter(it == "ps2")
  pfc1_psbsd <- pfc1 %>% filter(it == "psbsd")
  pfc1_ps2sd <- pfc1 %>% filter(it == "ps2sd")
  
  pfc1_ps2$ndvi_std <- (pfc1_ps2$ndvi_avg - mean(pfc1_ps2$ndvi_avg,na.rm=TRUE)) / sd(pfc1_ps2$ndvi_avg, na.rm=TRUE)
  pfc1_psbsd$ndvi_std <- (pfc1_psbsd$ndvi_avg - mean(pfc1_psbsd$ndvi_avg,na.rm=TRUE)) / sd(pfc1_psbsd$ndvi_avg, na.rm=TRUE)
  pfc1_ps2sd$ndvi_std <- (pfc1_ps2sd$ndvi_avg - mean(pfc1_ps2sd$ndvi_avg,na.rm=TRUE)) / sd(pfc1_ps2sd$ndvi_avg, na.rm=TRUE)
  
  
  pfc1 <- rbind.data.frame(pfc1_ps2, pfc1_psbsd)
  
  pfc1 <- pfc1[,-c(3)]
  
  pfc1 <- pfc1 %>%
    group_by(datez) %>%
    summarize(
      ndvi_std_avg = mean(ndvi_std, na.rm = TRUE)
    )
  
  # now combine and only keep the unique dates from pfc1_ps2sd that pfc1 doesnt have...
  
  pfc1_ps2sd <- pfc1_ps2sd[,-c(1,3)]
  colnames(pfc1_ps2sd) <- c("datez","ndvi_std_avg")
  dates_to_remove <- c(pfc1$datez)
  pfc1_ps2sd <- pfc1_ps2sd %>%
    filter(!datez %in% dates_to_remove)
  
  pfc1 <- rbind.data.frame(pfc1, pfc1_ps2sd)
  pfc1$ndvi_std_avg[is.nan(pfc1$ndvi_std_avg)] <- NA
  
  #any(duplicated(pfc1$datez)) #, FALSE, ensure I have only unique days left
  
  colnames(pfc1) <- c("datez.int","NDVI_std")
  
  #plot(y=pfc1$NDVI_std,x=pfc1$datez)
  
  pfc1$datez <- as.Date(pfc1$datez.int, origin = "1970-01-01")
  
  # Drop data that comes before 2017-01-01
  
  pfc1 <- pfc1[pfc1$datez >= "2017-01-01", ]
  
  #plot(y=pfc1$NDVI_std,x=pfc1$datez)
  
  dy_std_y[[i]] <- pfc1
  
  
}


### Weekly Data Second...

wk_std_y <- list()

for(i in 1:length(dy_std_y)){
  
  a <- dy_std_y[[i]]
  
  a$wk_year <- strftime(a$datez, format = "%Y-%V")
  
  a_wk_avgs <- a %>%
    group_by(wk_year) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
  
  a_wk_avgs$datez <- as.Date(a_wk_avgs$datez.int, origin = "1970-01-01")
  
  a_wk_avgs$NDVI_std[is.nan(a_wk_avgs$NDVI_std)] <- NA
  
  wk_std_y[[i]] <- a_wk_avgs
  
}


### Monthly Data Third...

mo_std_y <- list()

for(i in 1:length(dy_std_y)){
  
  a <- dy_std_y[[i]]
  
  a$mo_year <- strftime(a$datez, format = "%Y-%m")
  
  a_mo_avgs <- a %>%
    group_by(mo_year) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
  
  a_mo_avgs$datez <-  as.Date(paste0(a_mo_avgs$mo_year, "-01"), format = "%Y-%m-%d")
  
  a_mo_avgs$NDVI_std[is.nan(a_mo_avgs$NDVI_std)] <- NA
  
  mo_std_y[[i]] <- a_mo_avgs
  
  
}




# Make daily, weekly, and monthly ssa lists

dy.ssa.lst <- list()

for(i in 1:length(dy_std_y)){
  
  ### Make complete date vectors...can include missing obs for given date value
  
  # Define the start and end dates
  start_date <- min(dy_std_y[[i]][["datez"]])
  end_date <- max(dy_std_y[[i]][["datez"]])
  
  # Create a date vector including all days between the two dates
  date_vector <- as.data.frame(seq(from = start_date, to = end_date, by = "days"))
  colnames(date_vector) <- c("datez")
  
  dy.ssa <- full_join(x=date_vector, y=dy_std_y[[i]], by="datez")
  
  
  
  if(sum(is.na(dy.ssa$NDVI_std)) > 0) {
    
    test.dy.ssa <- ssa(dy.ssa$NDVI_std, L=1095, force.decompose = FALSE)
    
    # Use igapfill() to impute the missing values
    dy.filled <- igapfill(test.dy.ssa, groups = list(1:4)) 
    
    dy.filled2 <- imputeTS::na_kalman(dy.ssa$NDVI_std, model = "auto.arima", smooth = TRUE )
    
    # plot(dy.filled) +
    # lines(dy.ssa$NDVI_std, col="red")
    
    dy.final <- ssa(dy.filled, L = 1095, kind = c("toeplitz-ssa"))
    
    dy.final2 <- ssa(dy.filled2, L = 1095, kind = c("toeplitz-ssa"))
    
    # Plot the eigenvectors (elementary series) to visually identify the annual component
    # This shows the first two eigenvectors capture the periodicity
    
    # plot(dy.final, type = "vectors", idx = 1:24) # Adjust idx based on your needs
    # plot(dy.final2, type = "vectors", idx = 1:24) # Adjust idx based on your needs
    
    # Reconstruct the annual cycle (e.g., using eigentriples 1 and 2)
    annual_component <- reconstruct(dy.final, groups = list(annual = 1:4))
    
    annual_component2 <- reconstruct(dy.final2, groups = list(annual = 1:4))
    
    # Subtract the annual component from the original series
    deseasonalized_ndvi <- dy.ssa$NDVI_std - annual_component$annual
    
    deseasonalized_ndvi2 <- dy.ssa$NDVI_std - annual_component2$annual
    
    dy.ssa$deseas_std_NDVI <- deseasonalized_ndvi
    
    dy.ssa$dekalm_std_NDVI <- deseasonalized_ndvi2
    
    # acf(dy.ssa$deseas_std_NDVI, lag.max = 730, na.action = na.pass)
    # acf(dy.ssa$dekalm_std_NDVI, lag.max = 730, na.action = na.pass)
    
    
    # dy.ssa$mo_year <- strftime(dy.ssa$datez, format = "%b")
    # dy.ssa$Month <- as.factor(dy.ssa$mo_year)
    # ggplot(dy.ssa) + geom_point(aes(x=datez, y=deseas_std_NDVI,color=Month)) +
    # scale_colour_viridis_d()
    
    
    
  }else{
    
    test.dy.ssa <- ssa(dy.ssa$NDVI_std, L=1095, , kind = c("toeplitz-ssa"))
    
    dy.final <- test.dy.ssa
    
    # Plot the eigenvectors (elementary series) to visually identify the annual component
    # This shows the first two eigenvectors capture the periodicity
    # plot(dy.final, type = "vectors", idx = 1:24) # Adjust idx based on your needs
    
    # Reconstruct the annual cycle (e.g., using eigentriples 1 and 2)
    annual_component <- reconstruct(dy.final, groups = list(annual = 1:4))
    
    # Subtract the annual component from the original series
    deseasonalized_ndvi <- dy.ssa$NDVI_std - annual_component$annual
    dy.ssa$deseas_std_NDVI <- deseasonalized_ndvi
    
  }
  
  
  dy.ssa.lst[[i]] <- dy.ssa
  
  
}


wk.ssa.lst <- list()

for(i in 1:length(wk_std_y)){
  
  # Define the start and end dates
  start_date <- min(wk_std_y[[i]][["datez"]])
  end_date <- max(wk_std_y[[i]][["datez"]])
  
  # Create a date vector including all days between the two dates
  date_vector <- as.data.frame(seq(from = start_date, to = end_date, by = "weeks"))
  colnames(date_vector) <- c("datez")
  date_vector$wk_year <- strftime(date_vector$datez, format = "%Y-%V")
  
  wk.ssa <- full_join(x=date_vector, y=wk_std_y[[i]], by="wk_year")
  wk.ssa <- wk.ssa[, -c(3,5)]
  colnames(wk.ssa) <- c("datez","wk_year","NDVI_std")
  
  wk.ssa <- wk.ssa[!is.na(wk.ssa$datez),]
  
  
  if(sum(is.na(wk.ssa$NDVI_std)) > 0) {
    
    test.wk.ssa <- ssa(wk.ssa$NDVI_std, L=156, force.decompose = FALSE)
    
    # Use igapfill() to impute the missing values
    wk.filled <- igapfill(test.wk.ssa, groups = list(1:4)) 
    
    wk.filled2 <- imputeTS::na_kalman(wk.ssa$NDVI_std, model = "auto.arima", smooth = TRUE )
    
    # plot(wk.filled) +
    # lines(wk.ssa$NDVI_std, col="red")
    
    wk.final <- ssa(wk.filled, L = 156, kind = c("toeplitz-ssa"))
    
    wk.final2 <- ssa(wk.filled2, L = 156, kind = c("toeplitz-ssa"))
    
    # Plot the eigenvectors (elementary series) to visually identify the annual component
    # This shows the first two eigenvectors capture the periodicity
    #plot(wk.final, type = "vectors", idx = 1:24) # Adjust idx based on your needs
    
    # Reconstruct the annual cycle (e.g., using eigentriples 1 and 2)
    # 
    annual_component <- reconstruct(wk.final, groups = list(annual = 1:4))
    
    annual_component2 <- reconstruct(wk.final2, groups = list(annual = 1:4))
    
    # Subtract the annual component from the original series
    deseasonalized_ndvi <- wk.ssa$NDVI_std - annual_component$annual
    
    deseasonalized_ndvi2 <- wk.ssa$NDVI_std - annual_component2$annual
    
    wk.ssa$deseas_std_NDVI <- deseasonalized_ndvi
    
    wk.ssa$dekalm_std_NDVI <- deseasonalized_ndvi2
    
    # acf(wk.ssa$deseas_std_NDVI, lag.max = 156, na.action = na.pass)
    # wk.ssa$mo_year <- strftime(wk.ssa$datez_all, format = "%b")
    # wk.ssa$Month <- as.factor(wk.ssa$mo_year)
    # ggplot(wk.ssa) + geom_point(aes(x=datez_all, y=deseas_std_NDVI,color=Month)) +
    #   scale_colour_viridis_d()
    
    
    
  }else{
    
    test.wk.ssa <- ssa(wk.ssa$NDVI_std, L=156, , kind = c("toeplitz-ssa"))
    
    wk.final <- test.wk.ssa
    
    # Plot the eigenvectors (elementary series) to visually identify the annual component
    # This shows the first two eigenvectors capture the periodicity
    # plot(wk.final, type = "vectors", idx = 1:24) # Adjust idx based on your needs
    
    # Reconstruct the annual cycle (e.g., using eigentriples 1 and 2)
    annual_component <- reconstruct(wk.final, groups = list(annual = 1:4))
    
    # Subtract the annual component from the original series
    deseasonalized_ndvi <- wk.ssa$NDVI_std - annual_component$annual
    wk.ssa$deseas_std_NDVI <- deseasonalized_ndvi
    
  }
  
  
  wk.ssa.lst[[i]] <- wk.ssa
  
  
}


mo.ssa.lst <- list()

for(i in 1:length(mo_std_y)){

# Define the start and end dates
start_date <- min(mo_std_y[[i]][["datez"]])
end_date <- max(mo_std_y[[i]][["datez"]])

# Create a date vector including all days between the two dates
date_vector <- as.data.frame(seq(from = start_date, to = end_date, by = "months"))
colnames(date_vector) <- c("datez")

mo.ssa <- full_join(x=date_vector, y=mo_std_y[[i]], by="datez")
mo.ssa <- mo.ssa[, -c(3)]
colnames(mo.ssa) <- c("datez","mo_year","NDVI_std")

mo.ssa <- mo.ssa[!is.na(mo.ssa$datez),]


if(sum(is.na(mo.ssa$NDVI_std)) > 0) {
  
  test.mo.ssa <- ssa(mo.ssa$NDVI_std, L=36, force.decompose = FALSE)
  
  # Use igapfill() to impute the missing values
  mo.filled <- igapfill(test.mo.ssa, groups = list(1:3)) 
  
  mo.filled2 <- imputeTS::na_kalman(mo.ssa$NDVI_std, model = "auto.arima", smooth = TRUE )
  
  mo.final <- ssa(mo.filled, L = 36, kind = c("toeplitz-ssa"))
  
  mo.final2 <- ssa(mo.filled2, L = 36, kind = c("toeplitz-ssa"))
  
  # Plot the eigenvectors (elementary series) to visually identify the annual component
  # This shows the first two eigenvectors capture the periodicity
  # plot(mo.final, type = "vectors", idx = 1:24) # Adjust idx based on your needs
  
  # Reconstruct the annual cycle (e.g., using eigentriples 1 and 2)
  annual_component <- reconstruct(mo.final, groups = list(annual = 1:2))
  
  annual_component2 <- reconstruct(mo.final2, groups = list(annual = 1:2))
  
  # Subtract the annual component from the original series
  deseasonalized_ndvi <- mo.ssa$NDVI_std - annual_component$annual
  
  deseasonalized_ndvi2 <- mo.ssa$NDVI_std - annual_component2$annual
  
  mo.ssa$deseas_std_NDVI <- deseasonalized_ndvi
  
  mo.ssa$dekalm_std_NDVI <- deseasonalized_ndvi2
  
  
}else{
  
  test.mo.ssa <- ssa(mo.ssa$NDVI_std, L=36, , kind = c("toeplitz-ssa"))
  
  mo.final <- test.mo.ssa
  
  # Plot the eigenvectors (elementary series) to visually identify the annual component
  # This shows the first two eigenvectors capture the periodicity
  # plot(mo.final, type = "vectors", idx = 1:24) # Adjust idx based on your needs
  
  # Reconstruct the annual cycle (e.g., using eigentriples 1 and 2)
  annual_component <- reconstruct(mo.final, groups = list(annual = 1:2))
  
  # Subtract the annual component from the original series
  deseasonalized_ndvi <- mo.ssa$NDVI_std - annual_component$annual
  mo.ssa$deseas_std_NDVI <- deseasonalized_ndvi
  
}

mo.ssa.lst[[i]] <- mo.ssa


}




###
###
### Here are the X's
###
###


# 2016-01-01 to 2024-12-31, 3288 days

met_data <- read.csv( 
  "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/met_data.csv",header = TRUE)
colnames(met_data)[colnames(met_data) == "datez"] <- "datez.char"
met_data$datez <- as.Date(met_data$datez.char, format = "%Y-%m-%d")


# combine the X's into a single dataframe, and plot
# Only Keep Temp, Precip, ASCE ETr, and datez.char

met_data_reduced <- met_data[,c(4,5,7,9,11,13,14,16,17,20,22,38)] 

# plot(met_data_reduced$Avg.Temp_deg.C)
# plot(met_data_reduced$RH.Max_.)
# plot(met_data_reduced$Vapor.Pressure_kPa)
# plot(met_data_reduced$Liquid.Precip_mm) # remove negative values, erroneous
# plot(met_data_reduced$Wind.Run_km)
# plot(met_data_reduced$Gust.Speed_m.s)
# plot(met_data_reduced$Solar.Rad_W.m.2) # remove negative values, erroneous
# plot(met_data_reduced$ASCE.ETr_mm)

met_data_reduced$Liquid.Precip_mm <- replace(met_data_reduced$Liquid.Precip_mm, met_data_reduced$Liquid.Precip_mm < 0, 0)
met_data_reduced$Solar.Rad_W.m.2 <- replace(met_data_reduced$Solar.Rad_W.m.2, met_data_reduced$Solar.Rad_W.m.2 < 0, NA)







q_data <- read.csv( 
  "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/q_data.csv",header = TRUE)
colnames(q_data)[colnames(q_data) == "datez"] <- "datez.char"
q_data$datez <- as.Date(q_data$datez.char, format = "%Y-%m-%d")


q_data$q_cms <- q_data$DISCHRG.Value * 0.0283168
q_data$gage_ht_m <- q_data$GAGE_HT2.Value * 0.3048 ### Keep GAGE HEIGHT #2!!!!!
q_data$q_temp_C <- (q_data$WATTEMP.Value - 32) * (5/9)

# this dataframe keeps specific conductance (salinity measure) in the dataframe
#q_data_reduced <- q_data[,c(24,40,42:44)]


q_data_reduced <- q_data[,c(24,40,42:44)]











tranA_phd_avg <- read.csv( 
  "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/tranA_phd_avg.csv",header = TRUE)
tranA_postdoc_avg <- read.csv( 
  "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/tranA_postdoc_avg.csv",header = TRUE)

tranB_phd_avg <- read.csv( 
  "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/tranB_phd_avg.csv",header = TRUE)
tranB_postdoc_avg <- read.csv( 
  "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/tranB_postdoc_avg.csv",header = TRUE)

tranC_phd_avg <- read.csv( 
  "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/tranC_phd_avg.csv",header = TRUE)
tranC_postdoc_avg <- read.csv( 
  "C:/Users/19139/Desktop/NSF_Hypoth2_repeat/Results_ArkansasR_09212025/tranC_postdoc_avg.csv",header = TRUE)


date_df <- as.data.frame(seq(as.Date("2016-01-01"), as.Date("2024-12-31"), by = "days"), ncol=1)
colnames(date_df) <- c("datez")


# make list of gw df's in same order as veg-y lists


# pf_c1_df_rep_all, pg_c4_df_rep_all, wl_b3_df_rep_all, tm_a2_df_rep_all,
# ro_a1_df_rep_all,  pg_a3_df_rep_all, fl_b1_df_rep_all

{


c1_phd <- tranC_phd_avg[,c(2,5,7)]
c1_phd$C1_Temp_C <- (c1_phd$C1_Temp_F - 32) * (5/9)
c1_phd <- c1_phd[,c(1,3:4)] # only keep dtgw, gw_temp for now

c1_post <- tranC_postdoc_avg[,c(2,5,7)] # only keep dtgw, gw_temp for now

c1 <- rbind.data.frame(c1_phd,c1_post)
colnames(c1)[colnames(c1) == "datez"] <- "datez.int" 
c1$datez <- as.Date(c1$datez.int, origin = "1970-01-01")

# combine with date vector to create 3288 rows, with missing values
c1_full <- full_join(c1, date_df, by = "datez") %>%
  complete(datez = seq(min(.$datez), max(.$datez), by = "1 day"))

# I only need data from 2016-01-01 to 2024-12-31

c1_full <- c1_full[c1_full$datez < "2025-01-01", ]




c4_phd <- tranC_phd_avg[,c(2,14,16)]
c4_phd$C4_Temp_C <- (c4_phd$C4_Temp_F - 32) * (5/9)
c4_phd <- c4_phd[,c(1,3:4)] # only keep dtgw, gw_temp for now

c4_post <- tranC_postdoc_avg[,c(2,14,16)] # only keep dtgw, gw_temp for now

c4 <- rbind.data.frame(c4_phd,c4_post)
colnames(c4)[colnames(c4) == "datez"] <- "datez.int" 
c4$datez <- as.Date(c4$datez.int, origin = "1970-01-01")

# combine with date vector to create 3288 rows, with missing values
c4_full <- full_join(c4, date_df, by = "datez") %>%
  complete(datez = seq(min(.$datez), max(.$datez), by = "1 day"))

# I only need data from 2016-01-01 to 2024-12-31

c4_full <- c4_full[c4_full$datez < "2025-01-01", ]




b3_phd <- tranB_phd_avg[,c(2,11,13)]
b3_phd$B3_Temp_C <- (b3_phd$B3_Temp_F - 32) * (5/9)
b3_phd <- b3_phd[,c(1,3:4)] # only keep dtgw, gw_temp for now

b3_post <- tranB_postdoc_avg[,c(2,11,13)] # only keep dtgw, gw_temp for now

b3 <- rbind.data.frame(b3_phd,b3_post)
colnames(b3)[colnames(b3) == "datez"] <- "datez.int" 
b3$datez <- as.Date(b3$datez.int, origin = "1970-01-01")

# combine with date vector to create 3288 rows, with missing values
b3_full <- full_join(b3, date_df, by = "datez") %>%
  complete(datez = seq(min(.$datez), max(.$datez), by = "1 day"))

# I only need data from 2016-01-01 to 2024-12-31

b3_full <- b3_full[b3_full$datez < "2025-01-01", ]








a2_phd <- tranA_phd_avg[,c(2,8,10)]
a2_phd$A2_Temp_C <- (a2_phd$A2_Temp_F - 32) * (5/9)
a2_phd <- a2_phd[,c(1,3:4)] # only keep dtgw, gw_temp for now

a2_post <- tranA_postdoc_avg[,c(2,8,10)] # only keep dtgw, gw_temp for now

a2 <- rbind.data.frame(a2_phd,a2_post)
colnames(a2)[colnames(a2) == "datez"] <- "datez.int" 
a2$datez <- as.Date(a2$datez.int, origin = "1970-01-01")

# combine with date vector to create 3288 rows, with missing values
a2_full <- full_join(a2, date_df, by = "datez") %>%
  complete(datez = seq(min(.$datez), max(.$datez), by = "1 day"))

# I only need data from 2016-01-01 to 2024-12-31

a2_full <- a2_full[a2_full$datez < "2025-01-01", ]





a1_phd <- tranA_phd_avg[,c(2,5,7)]
a1_phd$A1_Temp_C <- (a1_phd$A1_Temp_F - 32) * (5/9)
a1_phd <- a1_phd[,c(1,3:4)] # only keep dtgw, gw_temp for now

a1_post <- tranA_postdoc_avg[,c(2,5,7)] # only keep dtgw, gw_temp for now

a1 <- rbind.data.frame(a1_phd,a1_post)
colnames(a1)[colnames(a1) == "datez"] <- "datez.int" 
a1$datez <- as.Date(a1$datez.int, origin = "1970-01-01")

# combine with date vector to create 3288 rows, with missing values
a1_full <- full_join(a1, date_df, by = "datez") %>%
  complete(datez = seq(min(.$datez), max(.$datez), by = "1 day"))

# I only need data from 2016-01-01 to 2024-12-31

a1_full <- a1_full[a1_full$datez < "2025-01-01", ]






a3_phd <- tranA_phd_avg[,c(2,11,13)]
a3_phd$A3_Temp_C <- (a3_phd$A3_Temp_F - 32) * (5/9)
a3_phd <- a3_phd[,c(1,3:4)] # only keep dtgw, gw_temp for now

a3_post <- tranA_postdoc_avg[,c(2,11,13)] # only keep dtgw, gw_temp for now

a3 <- rbind.data.frame(a3_phd,a3_post)
colnames(a3)[colnames(a3) == "datez"] <- "datez.int" 
a3$datez <- as.Date(a3$datez.int, origin = "1970-01-01")

# combine with date vector to create 3288 rows, with missing values
a3_full <- full_join(a3, date_df, by = "datez") %>%
  complete(datez = seq(min(.$datez), max(.$datez), by = "1 day"))

# I only need data from 2016-01-01 to 2024-12-31

a3_full <- a3_full[a3_full$datez < "2025-01-01", ]





b1_phd <- tranB_phd_avg[,c(2,5,7)]
b1_phd$B1_Temp_C <- (b1_phd$B1_Temp_F - 32) * (5/9)
b1_phd <- b1_phd[,c(1,3:4)] # only keep dtgw, gw_temp for now

b1_post <- tranB_postdoc_avg[,c(2,5,7)] # only keep dtgw, gw_temp for now

b1 <- rbind.data.frame(b1_phd,b1_post)
colnames(b1)[colnames(b1) == "datez"] <- "datez.int" 
b1$datez <- as.Date(b1$datez.int, origin = "1970-01-01")

# combine with date vector to create 3288 rows, with missing values
b1_full <- full_join(b1, date_df, by = "datez") %>%
  complete(datez = seq(min(.$datez), max(.$datez), by = "1 day"))

# I only need data from 2016-01-01 to 2024-12-31

b1_full <- b1_full[b1_full$datez < "2025-01-01", ]



}


gw_lst <- list(c1_full, c4_full, b3_full, a2_full, a1_full, a3_full, b1_full)



met_data_reduced$datez <- as.Date(met_data_reduced$datez.char, format = "%Y-%m-%d")
q_data_reduced$datez <- as.Date(q_data_reduced$datez.char, format = "%Y-%m-%d")

q_met_all <- full_join(met_data_reduced, q_data_reduced, by = "datez") %>%
  complete(datez = seq(min(.$datez), max(.$datez), by = "1 day"))


# Full join: For any date within this range that is missing from the joined data frame
# "complete" inserts an NA row


### Daily X's List


dy_std_x <- list()

for(i in 1:length(gw_lst)){
  
  # combine with date vector to create 3288 rows, with missing values
  x_full <- full_join(q_met_all, gw_lst[[i]], by = "datez") %>%
    complete(datez = seq(min(.$datez), max(.$datez), by = "1 day"))
  
  x_full_datez <- x_full$datez
  
  # Add new covariate: SW_temp minus GW_temp
  
  x_full$swT_m_gwT <- x_full[,18] - x_full[,21] 
  
  x_full <- x_full[,-c(1,13,15,19)]
  
  x_full_std <- scale(x_full)
  
  x_full_std <- as.data.frame(x_full_std)
  
  x_full_std$datez <- x_full_datez
  
  
  # Drop data that comes before 2017-01-01
  
  x_full_std <- x_full_std[x_full_std$datez >= "2017-01-01", ]
  
  
  dy_std_x[[i]] <- x_full_std
  
  
}




### Weekly X's List


wk_std_x <- list()

for(i in 1:length(dy_std_x)){
  
  x_wk <- dy_std_x[[i]]
  
  x_wk$wk_year <- strftime(x_wk$datez, format = "%Y-%V")
  
  x_wk_avgs <- x_wk %>%
    group_by(wk_year) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
  
  year_week_string <- x_wk_avgs$wk_year
  year <- as.numeric(substr(year_week_string, 1, 4))
  week_number <- as.numeric(substr(year_week_string, 6, 7))
  first_day_of_year <- ymd(paste0(year, "-01-01"))
  
  x_wk_avgs$datez <- floor_date(first_day_of_year + weeks(week_number - 1),
                                unit = "week",
                                week_start = 1)
  
  # Replace NaN with NA in numeric columns using dplyr
  x_wk_avgs <- x_wk_avgs %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))
  
  
  wk_std_x[[i]] <- x_wk_avgs
  
}




### Monthly X's List


mo_std_x <- list()

for(i in 1:length(dy_std_x)){
  
  x_mo <- dy_std_x[[i]]
  
  x_mo$mo_year <- strftime(x_mo$datez, format = "%Y-%m")
  
  x_mo_avgs <- x_mo %>%
    group_by(mo_year) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
  
  x_mo_avgs$datez <- as.Date(paste0(x_mo_avgs$mo_year, "-01"), format = "%Y-%m-%d")
  
  # Replace NaN with NA in numeric columns using dplyr
  x_mo_avgs <- x_mo_avgs %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))
  
  mo_std_x[[i]] <- x_mo_avgs
  
}






# Right here is where I could implement the code to do away with variables I dont need

# Or I could wait, make an exhaustive amount of lags per each then filter each one down






# Start with the monthly covariates first

mo_covs_x <- list()

datez_mo_cov_x <- list()

for(i in 1:length(mo_std_x)){
  
  x_mo_avgs <- mo_std_x[[i]]
  
  # Identify numeric columns
  numeric_cols <- names(x_mo_avgs)[sapply(x_mo_avgs, is.numeric)]
  
  # NAme lag t = 0 for Xs, 
  
  x_mo_avgs_l0 <- x_mo_avgs %>%
    mutate(across(all_of(numeric_cols), .names = "{.col}_l0"))
  
  x_mo_avgs_l0 <- x_mo_avgs_l0[,c(1,20:38)]
  
  #x_mo_avgs$datez <- x_mo_avgs$datez
  
  
  # Identify numeric columns
  numeric_cols <- names(x_mo_avgs)[sapply(x_mo_avgs, is.numeric)]
  
  
  # Make lag t = 1 for Xs, 
  
  x_mo_lag1 <- x_mo_avgs %>%
    mutate(across(all_of(numeric_cols), ~lag(., n = 1), .names = "{.col}_l1"))
  x_mo_lag1$datez <- x_mo_avgs$datez
  
  # Make average of lagged 1-2-3-4 month and take average (this tries to capture intra-seasonal)

  x_mo_lag1234avg <- x_mo_avgs %>%
    mutate(across(all_of(numeric_cols), ~ {
      lag_1 <- lag(., n = 1)
      lag_2 <- lag(., n = 2)
      lag_3 <- lag(., n = 3)
      lag_4 <- lag(., n = 4)
      (lag_1 + lag_2 + lag_3 + lag_4) / 4
    }, .names = "{.col}_l_1_4avg"))

  x_mo_lag1234avg$datez <- x_mo_avgs$datez

  # Make average of lagged 5-6-7 month and take average (this tries to capture intra-seasonal)
  
  x_mo_lag4567avg <- x_mo_avgs %>%
    mutate(across(all_of(numeric_cols), ~ {
      lag_4 <- lag(., n = 4)
      lag_5 <- lag(., n = 5)
      lag_6 <- lag(., n = 6)
      lag_7 <- lag(., n = 7)
      (lag_4 + lag_5 + lag_6 + lag_7) / 4
    }, .names = "{.col}_l_4_7avg"))
  
  x_mo_lag4567avg$datez <- x_mo_avgs$datez
  
  
  x_mo_avgs <- right_join(x=x_mo_avgs, y=x_mo_lag1[,c(20:38)], by = "datez")
  x_mo_avgs <- right_join(x=x_mo_avgs, y=x_mo_lag1234avg[,c(20:38)], by = "datez")
  x_mo_avgs <- right_join(x=x_mo_avgs, y=x_mo_lag4567avg[,c(20:38)], by = "datez")
  
  
  # Replace with _l0 column names
  
  x_mo_avgs <- x_mo_avgs[,-c(2:19)]
  
  x_mo_avgs <- cbind.data.frame(x_mo_avgs, x_mo_avgs_l0[,c(3:20)])
  
  
  
  
  # Join X and y datasets, remove missing NA's in y dataset using Right Join instead of Full Join
  
  y_x_mo_avgs <- right_join(x=x_mo_avgs, y=mo.ssa.lst[[i]], by = "datez")
  
  # Make column for Months
  
  y_x_mo_avgs$Month_l0 <- as.factor(format(y_x_mo_avgs$datez, "%b"))
  
  y_x_mo_avgs <- y_x_mo_avgs %>%
    mutate(
      month = month(datez),
      season_l0 = case_when(
        month %in% c(12, 1, 2) ~ "Winter",
        month %in% c(3, 4, 5) ~ "Spring",
        month %in% c(6, 7, 8) ~ "Summer",
        TRUE ~ "Fall" # Covers months 9, 10, 11
      )
    )
  
  y_x_mo_avgs$season_l0  <- as.factor(y_x_mo_avgs$season_l0)
  
  # Remove all winter time observations (based on my own data's evidence of dormant state)
  
  y_x_mo_avgs <- y_x_mo_avgs %>%
    filter(!Month_l0 %in% c("Dec", "Jan", "Feb")) %>%
    mutate(Month_l0 = droplevels(Month_l0)) # Drop unused levels
  
  # Replace NaN with NA in numeric columns using dplyr
  
  y_x_mo_avgs <- y_x_mo_avgs %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))
  
  #Drop all rows if NA value in NDVI column
  
  y_x_mo_avgs <- y_x_mo_avgs[!is.na(y_x_mo_avgs$deseas_std_NDVI), ]
  
  # # Store Datez Column
  # 
   datez_mo_cov_x[[i]] <- y_x_mo_avgs[,c("datez")]
  # 
  # #Drop Dates columns
  # 
  # y_x_mo_avgs_cln <- y_x_mo_avgs[,-c(1,11)]
  
  
  mo_covs_x[[i]] <- y_x_mo_avgs
  
}


# Start with the weekly covariates second

wk_covs_x <- list()

datez_wk_cov_x <- list()

for(i in 1:length(wk_std_x)){
  
  x_wk_avgs <- wk_std_x[[i]]
  
  
  # Identify numeric columns
  numeric_cols <- names(x_wk_avgs)[sapply(x_wk_avgs, is.numeric)]
  
  # NAme lag t = 0 for Xs, 
  
  x_wk_avgs_l0 <- x_wk_avgs %>%
    mutate(across(all_of(numeric_cols), .names = "{.col}_l0"))
  
  x_wk_avgs_l0 <- x_wk_avgs_l0[,c(1,20:38)]
  
  
  
  
  # Identify numeric columns
  numeric_cols <- names(x_wk_avgs)[sapply(x_wk_avgs, is.numeric)]
  
  
  
  # Make lag t = 1 for Xs, 
  
  x_wk_lag1 <- x_wk_avgs %>%
    mutate(across(all_of(numeric_cols), ~lag(., n = 1), .names = "{.col}_l1"))
  x_wk_lag1$datez <- x_wk_avgs$datez
  
  # Make average of lagged 1_4 weeks and take average (this tries to capture intra-monthly
  x_wk_lag1_4avg <- x_wk_avgs %>%
    mutate(across(all_of(numeric_cols), ~ {
      lag_1 <- lag(., n = 1)
      lag_2 <- lag(., n = 2)
      lag_3 <- lag(., n = 3)
      lag_4 <- lag(., n = 4)
      (lag_1 + lag_2 + lag_3 + lag_4) / 4
    }, .names = "{.col}_l_1_4avg"))
  
  x_wk_lag1_4avg$datez <- x_wk_avgs$datez
  
  # Make average of lagged 5-8 weeks (this tries to capture intra-smonthly)
  
  x_wk_lag5_8avg <- x_wk_avgs %>%
    mutate(across(all_of(numeric_cols), ~ {
      lag_5 <- lag(., n = 5)
      lag_6 <- lag(., n = 6)
      lag_7 <- lag(., n = 7)
      lag_8 <- lag(., n = 8)
      (lag_5 + lag_6 + lag_7 + lag_8) / 4
    }, .names = "{.col}_l_5_8avg"))
  
  x_wk_lag5_8avg$datez <- x_wk_avgs$datez
  
  
  x_wk_avgs  <- right_join(x=x_wk_avgs, y=x_wk_lag1[,c(1,20:38)], by = "wk_year")
  x_wk_avgs  <- right_join(x=x_wk_avgs, y=x_wk_lag1_4avg[,c(1,20:38)], by = "wk_year")
  x_wk_avgs  <- right_join(x=x_wk_avgs, y=x_wk_lag5_8avg[,c(1,20:38)], by = "wk_year")
  

  
  # Replace with _l0 column names
  
  x_wk_avgs <- x_wk_avgs[,-c(2:19)]
  
  x_wk_avgs <- cbind.data.frame(x_wk_avgs, x_wk_avgs_l0[,c(3:20)])
  
  
  
  y_x_wk_avgs <- right_join(x=x_wk_avgs, y=wk.ssa.lst[[i]], by = "wk_year")
  
  
  
  year_week_string <- y_x_wk_avgs$wk_year
  
  year <- as.numeric(substr(year_week_string, 1, 4))
  week_number <- as.numeric(substr(year_week_string, 6, 7))
  first_day_of_year <- ymd(paste0(year, "-01-01"))
  
  y_x_wk_avgs$datez <- floor_date(first_day_of_year + weeks(week_number - 1), 
                                  unit = "week", 
                                  week_start = 1)
  
  y_x_wk_avgs$Month_l0 <- as.factor(format(y_x_wk_avgs$datez , "%b"))
  
  
  y_x_wk_avgs <- y_x_wk_avgs %>%
    mutate(
      month = month(datez),
      season_l0 = case_when(
        month %in% c(12, 1, 2) ~ "Winter",
        month %in% c(3, 4, 5) ~ "Spring",
        month %in% c(6, 7, 8) ~ "Summer",
        TRUE ~ "Fall" # Covers months 9, 10, 11
      )
    )
  
  y_x_wk_avgs$season_l0  <- as.factor(y_x_wk_avgs$season_l0)
  
  
  
  y_x_wk_avgs <- y_x_wk_avgs %>%
    filter(!Month_l0 %in% c("Dec", "Jan", "Feb")) %>%
    mutate(Month_l0 = droplevels(Month_l0)) # Drop unused levels
  
  
  y_x_wk_avgs <- y_x_wk_avgs %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))
  
  y_x_wk_avgs <- y_x_wk_avgs[!is.na(y_x_wk_avgs$deseas_std_NDVI), ]
  
  #Keep/store the dates column
  datez_wk_cov_x[[i]] <- y_x_wk_avgs[,c("datez","wk_year")]
  
  # y_x_wk_avgs_cln <- y_x_wk_avgs[,-c(1,11,12,22,32,34)]
  
  wk_covs_x[[i]] <- y_x_wk_avgs
  
}


# Start with the Daily covariates third

dy_covs_x <- list()

datez_dy_cov_x <- list()

for(i in 1:length(dy_std_x)){
  
  x_dy <- dy_std_x[[i]]
  
  # Identify numeric columns
  numeric_cols <- names(x_dy)[sapply(x_dy, is.numeric)]
  
  # NAme lag t = 0 for Xs, 
  
  x_dy_l0 <- x_dy %>%
    mutate(across(all_of(numeric_cols), .names = "{.col}_l0"))
  
  x_dy_l0 <- x_dy_l0[,c(19:37)]
  
  
  
  # Identify numeric columns
  numeric_cols <- names(x_dy)[sapply(x_dy, is.numeric)]
  
  # Make lag t = 1 for Xs, 
  
  x_dy_lag1 <- x_dy %>%
    mutate(across(all_of(numeric_cols), ~lag(., n = 1), .names = "{.col}_l1"))
  x_dy_lag1$datez <- x_dy$datez
  
  
  # Compute lags 3 days and their average for each numeric column (intra-weekly)
  # This kinda gets at management scenario controlling weekly flows...
  x_dy_lag1_4avg <- x_dy %>%
    mutate(across(all_of(numeric_cols), ~ {
      lag_1 <- lag(., n = 1)
      lag_2 <- lag(., n = 2)
      lag_3 <- lag(., n = 3)
      lag_4 <- lag(., n = 4)
      (lag_1 + lag_2 + lag_3 + lag_4) / 4
    }, .names = "{.col}_l_1_4avg"))
  
  x_dy_lag1_4avg$datez <- x_dy$datez
  
  # Compute lags 7 days and their average for each numeric column (intra-monthly)
  # This kinda gets at management scenario controlling weekly flows...
  x_dy_lag4_7avg <- x_dy %>%
    mutate(across(all_of(numeric_cols), ~ {
      lag_4 <- lag(., n = 4)
      lag_5 <- lag(., n = 5)
      lag_6 <- lag(., n = 6)
      lag_7 <- lag(., n = 7)
      (lag_4 + lag_5 + lag_6 + lag_7) / 4
    }, .names = "{.col}_l_4_7avg"))
  
  x_dy_lag4_7avg$datez <- x_dy$datez
  
  
  x_dy_avgs <- right_join(x=x_dy, y=x_dy_lag1[,c(19:37)], by = "datez")
  x_dy_avgs <- right_join(x=x_dy_avgs, y=x_dy_lag1_4avg[,c(19:37)], by = "datez")
  x_dy_avgs <- right_join(x=x_dy_avgs, y=x_dy_lag4_7avg[,c(19:37)], by = "datez")
  
  
  
  # Replace with _l0 column names
  
  x_dy_avgs <- x_dy_avgs[,-c(1:18)]
  
  x_dy_avgs <- cbind.data.frame(x_dy_avgs, x_dy_l0[,c(2:19)])
  
  
  
  y_x_dy_avgs <- right_join(x=x_dy_avgs, y=dy.ssa.lst[[i]], by = "datez")
  
  # Add Month col
  y_x_dy_avgs$Month_l0 <- as.factor(format(y_x_dy_avgs$datez, "%b"))
  
  
  y_x_dy_avgs <- y_x_dy_avgs %>%
    mutate(
      month = month(datez),
      season_l0 = case_when(
        month %in% c(12, 1, 2) ~ "Winter",
        month %in% c(3, 4, 5) ~ "Spring",
        month %in% c(6, 7, 8) ~ "Summer",
        TRUE ~ "Fall" # Covers months 9, 10, 11
      )
    )
  
  y_x_dy_avgs$season_l0  <- as.factor(y_x_dy_avgs$season_l0)
  
  
  
  
  # Drop all rows with winter time months in them...
  y_x_dy_avgs <- y_x_dy_avgs %>%
    filter(!Month_l0 %in% c("Dec", "Jan", "Feb")) %>%
    mutate(Month_l0 = droplevels(Month_l0)) # Drop unused levels
  
  
  y_x_dy_avgs <- y_x_dy_avgs %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))
  
  
  # get rid of rows that have NA's in the response variable
  y_x_dy_avgs <- y_x_dy_avgs[!is.na(y_x_dy_avgs$deseas_std_NDVI), ]
  
  #Keep/store the dates column
  datez_dy_cov_x[[i]] <- y_x_dy_avgs[,c("datez")]
  
  # #drop the dates column
  # y_x_dy_avgs_cln <- y_x_dy_avgs[,-c(10)]
  
  
  dy_covs_x[[i]] <- y_x_dy_avgs
  
}








### Implement code to go inside each "mo_covs_x" and get variables I need

# All the vars for MONTHLY...

# "Avg.Temp_deg.C"     "Max.Temp_deg.C"     "Min.Temp_deg.C"     "RH.Max_."          
# "RH.Min_."           "Vapor.Pressure_kPa" "Liquid.Precip_mm"   "Wind.Run_km"       
# "Gust.Speed_m.s"     "Solar.Rad_W.m.2"    "ASCE.ETr_mm"        "COND.Value"        
# "q_cms"              "gage_ht_m"          "q_temp_C"           "deseas_std_NDVI"
# "Month_l0"

# letters for other covariates in data

# "_DTGW_m"
# "_Temp_C"
# "swT_m_gwT"

# letters for lags in monthly data

# _l0
# _l1
# _l_1_4avg
# _l_4_7avg


# So, which variables do I want to use (MONTHLY):

namez_covars <- c("Avg.Temp_deg.C","Liquid.Precip_mm","ASCE.ETr_mm",
                  "_DTGW_m","_Temp_C","swT_m_gwT","Month_l0","season_l0")

#namez_lagz <- c("_l0","_l1","_l234avg","_l567avg")

select_y <- c("deseas_std_NDVI")

namez_lagz <- c("_l0","_l1")

y_x_mo_fin <- list()

for(i in 1:length(mo_covs_x)){
  
  df1 <- mo_covs_x[[i]]
  
  y_df <- df1 %>%
    select(matches(select_y))
  
  # Select columns based on variables I wanted, but all the lags included
  # So, just pick the micromet, discharge gage and dtgw stuff at t=0
  
  df2 <- df1 %>%
    select(matches(namez_covars))
  
  
  df3 <- df2 %>%
    select(contains(namez_lagz))


  df4 <- cbind.data.frame(y_df, df3)
  
  
  y_x_mo_fin[[i]] <- df4

 
  
}






# So, which variables do I want to use (WEEKLY):

namez_covars <- c("Avg.Temp_deg.C","Liquid.Precip_mm","ASCE.ETr_mm",
                  "_DTGW_m","_Temp_C","swT_m_gwT","Month_l0","season_l0")

#namez_lagz <- c("_l0","_l1","_l_1_4avg")

#select_y <- c("deseas_std_NDVI","dekalm_std_NDVI")

select_y <- c("deseas_std_NDVI")

namez_lagz <- c("_l0","_l1")

y_x_wk_fin <- list()

for(i in 1:length(wk_covs_x)){
  
  df1 <- wk_covs_x[[i]]
  
  y_df <- df1 %>%
    select(matches(select_y))
  
  # Select columns based on variables I wanted, but all the lags included
  # So, just pick the micromet, discharge gage and dtgw stuff at t=0
  
  df2 <- df1 %>%
    select(matches(namez_covars))
  
  
  df3 <- df2 %>%
    select(contains(namez_lagz))
  
  
  df4 <- cbind.data.frame(y_df, df3)
  
  
  y_x_wk_fin[[i]] <- df4
  
  
  
}





# So, which variables do I want to use (DAILY):

namez_covars <- c("Avg.Temp_deg.C","Liquid.Precip_mm","ASCE.ETr_mm",
                  "_DTGW_m","_Temp_C","swT_m_gwT","Month_l0","season_l0")

#namez_lagz <- c("_l0","_l1","_l_1_4avg","_l_4_7avg")

#select_y <- c("deseas_std_NDVI","dekalm_std_NDVI")

select_y <- c("deseas_std_NDVI")

namez_lagz <- c("_l0","_l1")

y_x_dy_fin <- list()


for(i in 1:length(dy_covs_x)){
  
  df1 <- dy_covs_x[[i]]
  
  y_df <- df1 %>%
    select(matches(select_y))
  
  # Select columns based on variables I wanted, but all the lags included
  # So, just pick the micromet, discharge gage and dtgw stuff at t=0
  
  df2 <- df1 %>%
    select(matches(namez_covars))
  
  
  df3 <- df2 %>%
    select(contains(namez_lagz))
  
  
  df4 <- cbind.data.frame(y_df, df3)
  
  
  y_x_dy_fin[[i]] <- df4
  
  
  
}



















### Test plot to make sure things looks normal after combining into lists

# test_plot <- as.data.frame(dy_covs_x[[1]])
# 
# ggplot(data=test_plot) + geom_point(aes(x=c(1:1047),y=deseas_std_NDVI))
# ggplot(data=test_plot) + geom_point(aes(x=c(1:1047),y=Avg.Temp_deg.C))
# ggplot(data=test_plot) + geom_point(aes(x=c(1:1047),y=C1_DTGW_m))

###
### Monthly RF process...
###

# Stratified sampling with the rsample package
seedz <- set.seed(82)
# lag is useful if you have lagged predictors

obsVspred.lst.mo <- list()
impurz.lst.mo <- list()
permz.lst.mo <- list()
import.lst.mo <- list()

for(i in 1:length(y_x_mo_fin)){
  
  split_mo <- initial_time_split(y_x_mo_fin[[i]], prop = 4/5, lag = 1)
  #split_mo <- initial_split(data = y_x_mo_fin[[i]], prop = 0.75)
  
  my_train  <- training(split_mo)
  my_test   <- testing(split_mo)
  

    
    # number of features
    n_features <- length(setdiff(names(my_train), "deseas_std_NDVI" ))
    
    # train a default random forest model
    # THIS IS JUST A SIMPLE MODEL TO MAKE COMPARISON TOO....
    dummy1 <- ranger(
      deseas_std_NDVI ~ ., 
      data = my_train,
      mtry = floor(n_features / 3),
      respect.unordered.factors = 'order', # this is recommended for categorical vars
      na.action = "na.omit",
      seed = seedz
    )
    
    # get OOB RMSE
    default_rmse <- sqrt(dummy1$prediction.error)
    
    
    
    
    
    
    # Make hyperperameter grid search to find optimal model
    
    # minimum node size is often use as 5 for regression problems
    # mtry should be centered on number-of-features divided by 3
    # sample fraction should be anywhere from 70% to 80%
    
    # create hyperparameter grid
    hyper_grid <- expand.grid(
      mtry = floor(seq(from=2,to=n_features,length.out=5)), # If many relevant predictors, lower m-try works better
      min.node.size = c(1, 3, 5, 10), # 5 for most regression problems
      replace = c(TRUE, FALSE),                               
      sample.fraction = c(0.7,0.75,0.80),                       
      rmse = NA                                               
    )
    
    #nrow(hyper_grid)
    
    
    
    # execute full cartesian grid search
    
    for(j in seq_len(nrow(hyper_grid))) {
      # fit model for ith hyperparameter combination
      fit <- ranger(
        formula         = deseas_std_NDVI ~ ., 
        data            = my_train, 
        num.trees       = n_features * 10,
        mtry            = hyper_grid$mtry[j],
        min.node.size   = hyper_grid$min.node.size[j],
        replace         = hyper_grid$replace[j],
        sample.fraction = hyper_grid$sample.fraction[j],
        verbose         = TRUE,
        seed            = seedz,
        respect.unordered.factors = 'order',
        na.action = "na.omit",
        write.forest = TRUE
      )
      # export OOB error 
      hyper_grid$rmse[j] <- sqrt(fit$prediction.error)
    }
    
    
    
    
    
    
    # assess top 10 models
    top10 <- hyper_grid %>%
      arrange(rmse) %>%
      mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
      head(10)
    
    
    
    
    # For the ranger package...
    # Once you identify the optimal parameters from the grid search:
    # 1) Rerun the model with those hyperparameter values
    # 2) AND crank up the number of trees to create more stable estimates of feature importance
    
    
    # re-run model with impurity-based variable importance
    rf_impurity <- ranger(
      formula = deseas_std_NDVI ~ ., 
      data = my_train, 
      num.trees = 2000, # can be exhaustive with this pick, now that best params are known
      mtry = top10$mtry[1],
      min.node.size = top10$min.node.size[1],
      sample.fraction = top10$sample.fraction[1],
      replace = top10$replace[1],
      importance = "impurity",
      respect.unordered.factors = "order",
      na.action = "na.omit",
      verbose = FALSE,
      seed  = seedz
    )
    
    # re-run model with permutation-based variable importance
    rf_permutation <- ranger(
      formula = deseas_std_NDVI ~ ., 
      data = my_train, 
      num.trees = 2000,
      mtry = top10$mtry[1],
      min.node.size = top10$min.node.size[1],
      sample.fraction = top10$sample.fraction[1],
      replace = top10$replace[1],
      importance = "permutation",
      respect.unordered.factors = "order",
      na.action = "na.omit",
      verbose = FALSE,
      seed  = seedz
    )
    
    
    
    # plot of the top #n_features most important variables based on impurity and permutation
    
    
    impurz <- vip::vip(rf_impurity, num_features = 10, bar = FALSE,geom = "point") + labs(x="Predictor", y="Importance (Impurity)",subtitle="(a)")
    permz <- vip::vip(rf_permutation, num_features = 10, bar = FALSE, geom = "point") + labs(x="Predictor", y="Importance (Permutation)",subtitle="(a)")

    # Assuming 'model' is your trained tree-based model (e.g., from ranger, xgboost, rpart)
    importance_scores <- vi(rf_impurity, method = "model")
    sorted_importance <- importance_scores[order(importance_scores$Importance, decreasing = TRUE), ]
    top_10_features <- sorted_importance$Variable[1:10]
    
    
    
    fit_pred <- predict(fit, data = my_test)
    fit_pred_vals <- as.data.frame(fit_pred$predictions)
    fit_pred_vals$Observed <- my_test$deseas_std_NDVI
    fit_pred_vals$Month <- my_test$Month_l0
    colnames(fit_pred_vals) <- c("Predicted","Observed","Month")
    
    fit_pred_vals$PFG <- veg.name[i]
    
    import.lst.mo[[i]] <- top_10_features
    
    impurz.lst.mo[[i]] <- impurz
    permz.lst.mo[[i]] <- permz
    
    obsVspred.lst.mo[[i]] <- fit_pred_vals
  
  
} 

obsVspred.df.mo <- do.call(rbind, obsVspred.lst.mo)

obsVspred.df.mo$PFG <-as.factor(obsVspred.df.mo$PFG)


r_squared_values <- obsVspred.df.mo %>%
  group_by(PFG) %>%
  do({
    model <- lm(Predicted ~ Observed, data = .)
    data.frame(r_squared = summary(model)$r.squared,
               intercept = coef(model)[1],
               slope = coef(model)[2])
  })


p.mo <- ggplot(obsVspred.df.mo, aes(x = Observed, y = Predicted, color = PFG)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(data = r_squared_values,
            aes(x = max(obsVspred.df.mo$Observed) * 0.8,
                y = max(obsVspred.df.mo$Predicted) * 0.15 - (as.numeric(PFG) - 1) * max(obsVspred.df.mo$Predicted) * 0.25, # Adjust y for each group
                label = paste0("~R^{2} == ", round(r_squared, 2)),color=PFG),
            inherit.aes = FALSE, parse = TRUE) + # Prevent inheriting color from main plot
  labs(title = "(a)",
       x = "Observed",
       y = "Predicted") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_viridis_d() + theme(legend.position = "top")









# What does it look like to extract the top 10 vars and run it again...
# This would involve examining the importance/impurity graph as we decreas no. of vars

  
select_y <- c("deseas_std_NDVI")

y_x_mo_fin2 <- list()

for(i in 1:length(mo_covs_x)){
  
  df1 <- mo_covs_x[[i]]
  
  y_df <- df1 %>%
    select(matches(select_y))
  
  # Select columns based on variables I wanted, but all the lags included
  # So, just pick the micromet, discharge gage and dtgw stuff at t=0
  
  namez_exact <- c(import.lst.mo[[i]])
  
  df2 <- df1 %>%
    select(matches(namez_exact))
  
  
  df3 <- cbind.data.frame(y_df, df2)
  
  
  y_x_mo_fin2[[i]] <- df3
  
  
  
}




obsVspred.lst.mo2 <- list()
impurz.lst.mo2 <- list()
permz.lst.mo2 <- list()
import.lst.mo2 <- list()

for(i in 1:length(y_x_mo_fin2)){
  
  split_mo <- initial_time_split(y_x_mo_fin2[[i]], prop = 4/5, lag = 1)
  #split_mo <- initial_split(data = y_x_mo_fin2[[i]], prop = 0.75)
  
  my_train  <- training(split_mo)
  my_test   <- testing(split_mo)
  
  
  
  # number of features
  n_features <- length(setdiff(names(my_train), "deseas_std_NDVI" ))
  
  # train a default random forest model
  # THIS IS JUST A SIMPLE MODEL TO MAKE COMPARISON TOO....
  dummy1 <- ranger(
    deseas_std_NDVI ~ ., 
    data = my_train,
    mtry = floor(n_features / 3),
    respect.unordered.factors = 'order', # this is recommended for categorical vars
    na.action = "na.omit",
    seed = seedz
  )
  
  # get OOB RMSE
  default_rmse <- sqrt(dummy1$prediction.error)
  
  
  
  
  
  
  # Make hyperperameter grid search to find optimal model
  
  # minimum node size is often use as 5 for regression problems
  # mtry should be centered on number-of-features divided by 3
  # sample fraction should be anywhere from 70% to 80%
  
  # create hyperparameter grid
  hyper_grid <- expand.grid(
    mtry = floor(seq(from=2,to=n_features,length.out=5)), # If many relevant predictors, lower m-try works better
    min.node.size = c(1, 3, 5, 10), # 5 for most regression problems
    replace = c(TRUE, FALSE),                               
    sample.fraction = c(0.7,0.75,0.80),                       
    rmse = NA                                               
  )
  
  #nrow(hyper_grid)
  
  
  
  # execute full cartesian grid search
  
  for(j in seq_len(nrow(hyper_grid))) {
    # fit model for ith hyperparameter combination
    fit <- ranger(
      formula         = deseas_std_NDVI ~ ., 
      data            = my_train, 
      num.trees       = n_features * 10,
      mtry            = hyper_grid$mtry[j],
      min.node.size   = hyper_grid$min.node.size[j],
      replace         = hyper_grid$replace[j],
      sample.fraction = hyper_grid$sample.fraction[j],
      verbose         = TRUE,
      seed            = seedz,
      respect.unordered.factors = 'order',
      na.action = "na.omit",
      write.forest = TRUE
    )
    # export OOB error 
    hyper_grid$rmse[j] <- sqrt(fit$prediction.error)
  }
  
  
  
  
  
  
  # assess top 10 models
  top10 <- hyper_grid %>%
    arrange(rmse) %>%
    mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
    head(10)
  
  
  
  
  # For the ranger package...
  # Once you identify the optimal parameters from the grid search:
  # 1) Rerun the model with those hyperparameter values
  # 2) AND crank up the number of trees to create more stable estimates of feature importance
  
  
  # re-run model with impurity-based variable importance
  rf_impurity <- ranger(
    formula = deseas_std_NDVI ~ ., 
    data = my_train, 
    num.trees = 2000, # can be exhaustive with this pick, now that best params are known
    mtry = top10$mtry[1],
    min.node.size = top10$min.node.size[1],
    sample.fraction = top10$sample.fraction[1],
    replace = top10$replace[1],
    importance = "impurity",
    respect.unordered.factors = "order",
    na.action = "na.omit",
    verbose = FALSE,
    seed  = seedz
  )
  
  # re-run model with permutation-based variable importance
  rf_permutation <- ranger(
    formula = deseas_std_NDVI ~ ., 
    data = my_train, 
    num.trees = 2000,
    mtry = top10$mtry[1],
    min.node.size = top10$min.node.size[1],
    sample.fraction = top10$sample.fraction[1],
    replace = top10$replace[1],
    importance = "permutation",
    respect.unordered.factors = "order",
    na.action = "na.omit",
    verbose = FALSE,
    seed  = seedz
  )
  
  
  
  # plot of the top #n_features most important variables based on impurity and permutation
  
  
  impurz <- vip::vip(rf_impurity, num_features = 10, bar = FALSE,geom = "point") + labs(x="Predictor", y="Importance (Impurity)",subtitle="(b)")
  permz <- vip::vip(rf_permutation, num_features = 10, bar = FALSE, geom = "point") + labs(x="Predictor", y="Importance (Permutation)",subtitle="(b)")
  
  # Assuming 'model' is your trained tree-based model (e.g., from ranger, xgboost, rpart)
  importance_scores <- vi(rf_impurity, method = "model")
  sorted_importance <- importance_scores[order(importance_scores$Importance, decreasing = TRUE), ]
  top_10_features <- sorted_importance$Variable[1:10]
  
  
  
  fit_pred <- predict(fit, data = my_test)
  fit_pred_vals <- as.data.frame(fit_pred$predictions)
  fit_pred_vals$Observed <- my_test$deseas_std_NDVI
  # fit_pred_vals$Month <- my_test$Month_l0
  # colnames(fit_pred_vals) <- c("Predicted","Observed","Month")
  
  colnames(fit_pred_vals) <- c("Predicted","Observed")
  
  fit_pred_vals$PFG <- veg.name[i]
  
  import.lst.mo2[[i]] <- top_10_features
  
  impurz.lst.mo2[[i]] <- impurz
  permz.lst.mo2[[i]] <- permz
  
  obsVspred.lst.mo2[[i]] <- fit_pred_vals
  
  
} 


obsVspred.df.mo2 <- do.call(rbind, obsVspred.lst.mo2)

obsVspred.df.mo2$PFG <-as.factor(obsVspred.df.mo2$PFG)


r_squared_values <- obsVspred.df.mo2 %>%
  group_by(PFG) %>%
  do({
    model <- lm(Predicted ~ Observed, data = .)
    data.frame(r_squared = summary(model)$r.squared,
               intercept = coef(model)[1],
               slope = coef(model)[2])
  })


p.mo2 <- ggplot(obsVspred.df.mo2, aes(x = Observed, y = Predicted, color = PFG)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(data = r_squared_values,
            aes(x = max(obsVspred.df.mo$Observed) * 0.8,
                y = max(obsVspred.df.mo$Predicted) * 0.15 - (as.numeric(PFG) - 1) * max(obsVspred.df.mo$Predicted) * 0.25, # Adjust y for each group
                label = paste0("~R^{2} == ", round(r_squared, 2)),color=PFG),
            inherit.aes = FALSE, parse = TRUE) + # Prevent inheriting color from main plot
  labs(title = "(b)",
       x = "Observed",
       y = "Predicted") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_viridis_d() + theme(legend.position = "top")



gridExtra::grid.arrange(p.mo, p.mo2,ncol = 2)


gridExtra::grid.arrange(impurz.lst.mo[[1]], impurz.lst.mo2[[1]],ncol = 2)














###
### Weekly RF process...
###

# Stratified sampling with the rsample package
seedz <- set.seed(82)
# lag is useful if you have lagged predictors

obsVspred.lst.wk <- list()
impurz.lst.wk <- list()
permz.lst.wk <- list()
import.lst.wk <- list()

for(i in 1:length(y_x_wk_fin)){
  
  split_mo <- initial_time_split(y_x_wk_fin[[i]], prop = 4/5, lag = 1)
  #split_mo <- initial_split(data = y_x_wk_fin[[i]], prop = 0.75)
  
  my_train  <- training(split_mo)
  my_test   <- testing(split_mo)
  
  
  
  # number of features
  n_features <- length(setdiff(names(my_train), "deseas_std_NDVI" ))
  
  # train a default random forest model
  # THIS IS JUST A SIMPLE MODEL TO MAKE COMPARISON TOO....
  dummy1 <- ranger(
    deseas_std_NDVI ~ ., 
    data = my_train,
    mtry = floor(n_features / 3),
    respect.unordered.factors = 'order', # this is recommended for categorical vars
    na.action = "na.omit",
    seed = seedz
  )
  
  # get OOB RMSE
  default_rmse <- sqrt(dummy1$prediction.error)
  
  
  
  
  
  
  # Make hyperperameter grid search to find optimal model
  
  # minimum node size is often use as 5 for regression problems
  # mtry should be centered on number-of-features divided by 3
  # sample fraction should be anywhere from 70% to 80%
  
  # create hyperparameter grid
  hyper_grid <- expand.grid(
    mtry = floor(seq(from=2,to=n_features,length.out=5)), # If many relevant predictors, lower m-try works better
    min.node.size = c(1, 3, 5, 10), # 5 for most regression problems
    replace = c(TRUE, FALSE),                               
    sample.fraction = c(0.7,0.75,0.80),                       
    rmse = NA                                               
  )
  
  #nrow(hyper_grid)
  
  
  
  # execute full cartesian grid search
  
  for(j in seq_len(nrow(hyper_grid))) {
    # fit model for ith hyperparameter combination
    fit <- ranger(
      formula         = deseas_std_NDVI ~ ., 
      data            = my_train, 
      num.trees       = n_features * 10,
      mtry            = hyper_grid$mtry[j],
      min.node.size   = hyper_grid$min.node.size[j],
      replace         = hyper_grid$replace[j],
      sample.fraction = hyper_grid$sample.fraction[j],
      verbose         = TRUE,
      seed            = seedz,
      respect.unordered.factors = 'order',
      na.action = "na.omit",
      write.forest = TRUE
    )
    # export OOB error 
    hyper_grid$rmse[j] <- sqrt(fit$prediction.error)
  }
  
  
  
  
  
  
  # assess top 10 models
  top10 <- hyper_grid %>%
    arrange(rmse) %>%
    mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
    head(10)
  
  
  
  
  # For the ranger package...
  # Once you identify the optimal parameters from the grid search:
  # 1) Rerun the model with those hyperparameter values
  # 2) AND crank up the number of trees to create more stable estimates of feature importance
  
  
  # re-run model with impurity-based variable importance
  rf_impurity <- ranger(
    formula = deseas_std_NDVI ~ ., 
    data = my_train, 
    num.trees = 2000, # can be exhaustive with this pick, now that best params are known
    mtry = top10$mtry[1],
    min.node.size = top10$min.node.size[1],
    sample.fraction = top10$sample.fraction[1],
    replace = top10$replace[1],
    importance = "impurity",
    respect.unordered.factors = "order",
    na.action = "na.omit",
    verbose = FALSE,
    seed  = seedz
  )
  
  # re-run model with permutation-based variable importance
  rf_permutation <- ranger(
    formula = deseas_std_NDVI ~ ., 
    data = my_train, 
    num.trees = 2000,
    mtry = top10$mtry[1],
    min.node.size = top10$min.node.size[1],
    sample.fraction = top10$sample.fraction[1],
    replace = top10$replace[1],
    importance = "permutation",
    respect.unordered.factors = "order",
    na.action = "na.omit",
    verbose = FALSE,
    seed  = seedz
  )
  
  
  
  # plot of the top #n_features most important variables based on impurity and permutation
  
  
  impurz <- vip::vip(rf_impurity, num_features = 10, bar = FALSE,geom = "point") + labs(x="Predictor", y="Importance (Impurity)",subtitle="(a)")
  permz <- vip::vip(rf_permutation, num_features = 10, bar = FALSE, geom = "point") + labs(x="Predictor", y="Importance (Permutation)",subtitle="(a)")
  
  
  # Assuming 'model' is your trained tree-based model (e.g., from ranger, xgboost, rpart)
  importance_scores <- vi(rf_impurity, method = "model")
  sorted_importance <- importance_scores[order(importance_scores$Importance, decreasing = TRUE), ]
  top_10_features <- sorted_importance$Variable[1:10]
  
  
  
  fit_pred <- predict(fit, data = my_test)
  fit_pred_vals <- as.data.frame(fit_pred$predictions)
  fit_pred_vals$Observed <- my_test$deseas_std_NDVI
  fit_pred_vals$Month <- my_test$Month_l0
  colnames(fit_pred_vals) <- c("Predicted","Observed","Month")
  
  fit_pred_vals$PFG <- veg.name[i]
  
  import.lst.wk[[i]] <- top_10_features
  
  impurz.lst.wk[[i]] <- impurz
  permz.lst.wk[[i]] <- permz
  
  obsVspred.lst.wk[[i]] <- fit_pred_vals
  
  
} 


obsVspred.df.wk <- do.call(rbind, obsVspred.lst.wk)

obsVspred.df.wk$PFG <-as.factor(obsVspred.df.wk$PFG)


r_squared_values <- obsVspred.df.wk %>%
  group_by(PFG) %>%
  do({
    model <- lm(Predicted ~ Observed, data = .)
    data.frame(r_squared = summary(model)$r.squared,
               intercept = coef(model)[1],
               slope = coef(model)[2])
  })


p.wk <- ggplot(obsVspred.df.wk, aes(x = Observed, y = Predicted, color = PFG)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(data = r_squared_values,
            aes(x = max(obsVspred.df.wk$Observed) * 0.8,
                y = max(obsVspred.df.wk$Predicted) * 0.15 - (as.numeric(PFG) - 1) * max(obsVspred.df.wk$Predicted) * 0.25, # Adjust y for each group
                label = paste0("~R^{2} == ", round(r_squared, 2)),color=PFG),
            inherit.aes = FALSE, parse = TRUE) + # Prevent inheriting color from main plot
  labs(title = "(a)",
       x = "Observed",
       y = "Predicted") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_viridis_d() + theme(legend.position = "top")







namez_covars <- c("Avg.Temp_deg.C","Liquid.Precip_mm","ASCE.ETr_mm",
                  "_DTGW_m","_Temp_C","swT_m_gwT","Month_l0","season_l0")

select_y <- c("dekalm_std_NDVI")

namez_lagz <- c("_l0","_l1")

y_x_wk_fin2 <- list()

for(i in 1:length(wk_covs_x)){
  
  df1 <- wk_covs_x[[i]]
  
  y_df <- df1 %>%
    select(matches(select_y))
  
  # Select columns based on variables I wanted, but all the lags included
  # So, just pick the micromet, discharge gage and dtgw stuff at t=0
  
  df2 <- df1 %>%
    select(matches(namez_covars))
  
  
  df3 <- df2 %>%
    select(contains(namez_lagz))
  
  
  df4 <- cbind.data.frame(y_df, df3)
  
  
  y_x_wk_fin2[[i]] <- df4
  
  
  
}



obsVspred.lst.wk2 <- list()
impurz.lst.wk2 <- list()
permz.lst.wk2 <- list()
import.lst.wk2 <- list()

for(i in 1:length(y_x_wk_fin2)){
  
  split_mo <- initial_time_split(y_x_wk_fin2[[i]], prop = 4/5, lag = 1)
  #split_mo <- initial_split(data = y_x_wk_fin2[[i]], prop = 0.75)
  
  my_train  <- training(split_mo)
  my_test   <- testing(split_mo)
  
  
  
  # number of features
  n_features <- length(setdiff(names(my_train), "dekalm_std_NDVI" ))
  
  # train a default random forest model
  # THIS IS JUST A SIMPLE MODEL TO MAKE COMPARISON TOO....
  dummy1 <- ranger(
    dekalm_std_NDVI ~ ., 
    data = my_train,
    mtry = floor(n_features / 3),
    respect.unordered.factors = 'order', # this is recommended for categorical vars
    na.action = "na.omit",
    seed = seedz
  )
  
  # get OOB RMSE
  default_rmse <- sqrt(dummy1$prediction.error)
  
  
  
  
  
  
  # Make hyperperameter grid search to find optimal model
  
  # minimum node size is often use as 5 for regression problems
  # mtry should be centered on number-of-features divided by 3
  # sample fraction should be anywhere from 70% to 80%
  
  # create hyperparameter grid
  hyper_grid <- expand.grid(
    mtry = floor(seq(from=2,to=n_features,length.out=5)), # If many relevant predictors, lower m-try works better
    min.node.size = c(1, 3, 5, 10), # 5 for most regression problems
    replace = c(TRUE, FALSE),                               
    sample.fraction = c(0.7,0.75,0.80),                       
    rmse = NA                                               
  )
  
  #nrow(hyper_grid)
  
  
  
  # execute full cartesian grid search
  
  for(j in seq_len(nrow(hyper_grid))) {
    # fit model for ith hyperparameter combination
    fit <- ranger(
      formula         = dekalm_std_NDVI ~ ., 
      data            = my_train, 
      num.trees       = n_features * 10,
      mtry            = hyper_grid$mtry[j],
      min.node.size   = hyper_grid$min.node.size[j],
      replace         = hyper_grid$replace[j],
      sample.fraction = hyper_grid$sample.fraction[j],
      verbose         = TRUE,
      seed            = seedz,
      respect.unordered.factors = 'order',
      na.action = "na.omit",
      write.forest = TRUE
    )
    # export OOB error 
    hyper_grid$rmse[j] <- sqrt(fit$prediction.error)
  }
  
  
  
  
  
  
  # assess top 10 models
  top10 <- hyper_grid %>%
    arrange(rmse) %>%
    mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
    head(10)
  
  
  
  
  # For the ranger package...
  # Once you identify the optimal parameters from the grid search:
  # 1) Rerun the model with those hyperparameter values
  # 2) AND crank up the number of trees to create more stable estimates of feature importance
  
  
  # re-run model with impurity-based variable importance
  rf_impurity <- ranger(
    formula = dekalm_std_NDVI ~ ., 
    data = my_train, 
    num.trees = 2000, # can be exhaustive with this pick, now that best params are known
    mtry = top10$mtry[1],
    min.node.size = top10$min.node.size[1],
    sample.fraction = top10$sample.fraction[1],
    replace = top10$replace[1],
    importance = "impurity",
    respect.unordered.factors = "order",
    na.action = "na.omit",
    verbose = FALSE,
    seed  = seedz
  )
  
  # re-run model with permutation-based variable importance
  rf_permutation <- ranger(
    formula = dekalm_std_NDVI ~ ., 
    data = my_train, 
    num.trees = 2000,
    mtry = top10$mtry[1],
    min.node.size = top10$min.node.size[1],
    sample.fraction = top10$sample.fraction[1],
    replace = top10$replace[1],
    importance = "permutation",
    respect.unordered.factors = "order",
    na.action = "na.omit",
    verbose = FALSE,
    seed  = seedz
  )
  
  
  
  # plot of the top #n_features most important variables based on impurity and permutation
  
  
  impurz <- vip::vip(rf_impurity, num_features = 10, bar = FALSE,geom = "point") + labs(x="Predictor", y="Importance (Impurity)",subtitle="(b)")
  permz <- vip::vip(rf_permutation, num_features = 10, bar = FALSE, geom = "point") + labs(x="Predictor", y="Importance (Permutation)",subtitle="(b)")
  
  # Assuming 'model' is your trained tree-based model (e.g., from ranger, xgboost, rpart)
  importance_scores <- vi(rf_impurity, method = "model")
  sorted_importance <- importance_scores[order(importance_scores$Importance, decreasing = TRUE), ]
  top_10_features <- sorted_importance$Variable[1:10]
  
  
  
  fit_pred <- predict(fit, data = my_test)
  fit_pred_vals <- as.data.frame(fit_pred$predictions)
  fit_pred_vals$Observed <- my_test$dekalm_std_NDVI
  # fit_pred_vals$Month <- my_test$Month_l0
  # colnames(fit_pred_vals) <- c("Predicted","Observed","Month")
  
  colnames(fit_pred_vals) <- c("Predicted","Observed")
  
  fit_pred_vals$PFG <- veg.name[i]
  
  import.lst.wk2[[i]] <- top_10_features
  
  impurz.lst.wk2[[i]] <- impurz
  permz.lst.wk2[[i]] <- permz
  
  obsVspred.lst.wk2[[i]] <- fit_pred_vals
  
  
} 


obsVspred.df.wk2 <- do.call(rbind, obsVspred.lst.wk2)

obsVspred.df.wk2$PFG <-as.factor(obsVspred.df.wk2$PFG)


r_squared_values <- obsVspred.df.wk2 %>%
  group_by(PFG) %>%
  do({
    model <- lm(Predicted ~ Observed, data = .)
    data.frame(r_squared = summary(model)$r.squared,
               intercept = coef(model)[1],
               slope = coef(model)[2])
  })


p.wk2 <- ggplot(obsVspred.df.wk2, aes(x = Observed, y = Predicted, color = PFG)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(data = r_squared_values,
            aes(x = max(obsVspred.df.wk$Observed) * 0.8,
                y = max(obsVspred.df.wk$Predicted) * 0.15 - (as.numeric(PFG) - 1) * max(obsVspred.df.wk$Predicted) * 0.25, # Adjust y for each group
                label = paste0("~R^{2} == ", round(r_squared, 2)),color=PFG),
            inherit.aes = FALSE, parse = TRUE) + # Prevent inheriting color from main plot
  labs(title = "(b)",
       x = "Observed",
       y = "Predicted") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_viridis_d() + theme(legend.position = "top")



gridExtra::grid.arrange(p.wk, p.wk2, ncol = 2)


gridExtra::grid.arrange(impurz.lst.wk[[1]], impurz.lst.wk2[[1]],ncol = 2)




















###
### Daily RF process...
###

# Stratified sampling with the rsample package
seedz <- set.seed(82)
# lag is useful if you have lagged predictors

obsVspred.lst.dy <- list()
impurz.lst.dy <- list()
permz.lst.dy <- list()
import.lst.dy <- list()

for(i in 1:length(y_x_dy_fin)){
  
  split_mo <- initial_time_split(y_x_dy_fin[[i]], prop = 4/5, lag = 1)
  #split_mo <- initial_split(data = y_x_dy_fin[[i]], prop = 0.75)
  
  my_train  <- training(split_mo)
  my_test   <- testing(split_mo)
  
  
  
  # number of features
  n_features <- length(setdiff(names(my_train), "deseas_std_NDVI" ))
  
  # train a default random forest model
  # THIS IS JUST A SIMPLE MODEL TO MAKE COMPARISON TOO....
  dummy1 <- ranger(
    deseas_std_NDVI ~ ., 
    data = my_train,
    mtry = floor(n_features / 3),
    respect.unordered.factors = 'order', # this is recommended for categorical vars
    na.action = "na.omit",
    seed = seedz
  )
  
  # get OOB RMSE
  default_rmse <- sqrt(dummy1$prediction.error)
  
  
  
  
  
  
  # Make hyperperameter grid search to find optimal model
  
  # minimum node size is often use as 5 for regression problems
  # mtry should be centered on number-of-features divided by 3
  # sample fraction should be anywhere from 70% to 80%
  
  # create hyperparameter grid
  hyper_grid <- expand.grid(
    mtry = floor(seq(from=2,to=n_features,length.out=5)), # If many relevant predictors, lower m-try works better
    min.node.size = c(1, 3, 5, 10), # 5 for most regression problems
    replace = c(TRUE, FALSE),                               
    sample.fraction = c(0.7,0.75,0.80),                       
    rmse = NA                                               
  )
  
  #nrow(hyper_grid)
  
  
  
  # execute full cartesian grid search
  
  for(j in seq_len(nrow(hyper_grid))) {
    # fit model for ith hyperparameter combination
    fit <- ranger(
      formula         = deseas_std_NDVI ~ ., 
      data            = my_train, 
      num.trees       = n_features * 10,
      mtry            = hyper_grid$mtry[j],
      min.node.size   = hyper_grid$min.node.size[j],
      replace         = hyper_grid$replace[j],
      sample.fraction = hyper_grid$sample.fraction[j],
      verbose         = TRUE,
      seed            = seedz,
      respect.unordered.factors = 'order',
      na.action = "na.omit",
      write.forest = TRUE
    )
    # export OOB error 
    hyper_grid$rmse[j] <- sqrt(fit$prediction.error)
  }
  
  
  
  
  
  
  # assess top 10 models
  top10 <- hyper_grid %>%
    arrange(rmse) %>%
    mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
    head(10)
  
  
  
  
  # For the ranger package...
  # Once you identify the optimal parameters from the grid search:
  # 1) Rerun the model with those hyperparameter values
  # 2) AND crank up the number of trees to create more stable estimates of feature importance
  
  
  # re-run model with impurity-based variable importance
  rf_impurity <- ranger(
    formula = deseas_std_NDVI ~ ., 
    data = my_train, 
    num.trees = 2000, # can be exhaustive with this pick, now that best params are known
    mtry = top10$mtry[1],
    min.node.size = top10$min.node.size[1],
    sample.fraction = top10$sample.fraction[1],
    replace = top10$replace[1],
    importance = "impurity",
    respect.unordered.factors = "order",
    na.action = "na.omit",
    verbose = FALSE,
    seed  = seedz
  )
  
  # re-run model with permutation-based variable importance
  rf_permutation <- ranger(
    formula = deseas_std_NDVI ~ ., 
    data = my_train, 
    num.trees = 2000,
    mtry = top10$mtry[1],
    min.node.size = top10$min.node.size[1],
    sample.fraction = top10$sample.fraction[1],
    replace = top10$replace[1],
    importance = "permutation",
    respect.unordered.factors = "order",
    na.action = "na.omit",
    verbose = FALSE,
    seed  = seedz
  )
  
  
  
  # plot of the top #n_features most important variables based on impurity and permutation
  
  
  impurz <- vip::vip(rf_impurity, num_features = 10, bar = FALSE,geom = "point") + labs(x="Predictor", y="Importance (Impurity)",subtitle="(a)")
  permz <- vip::vip(rf_permutation, num_features = 10, bar = FALSE, geom = "point") + labs(x="Predictor", y="Importance (Permutation)",subtitle="(a)")
  
  
  # Assuming 'model' is your trained tree-based model (e.g., from ranger, xgboost, rpart)
  importance_scores <- vi(rf_impurity, method = "model")
  sorted_importance <- importance_scores[order(importance_scores$Importance, decreasing = TRUE), ]
  top_10_features <- sorted_importance$Variable[1:10]
  
  
  
  fit_pred <- predict(fit, data = my_test)
  fit_pred_vals <- as.data.frame(fit_pred$predictions)
  fit_pred_vals$Observed <- my_test$deseas_std_NDVI
  fit_pred_vals$Month <- my_test$Month_l0
  colnames(fit_pred_vals) <- c("Predicted","Observed","Month")
  
  fit_pred_vals$PFG <- veg.name[i]
  
  import.lst.dy[[i]] <- top_10_features
  
  impurz.lst.dy[[i]] <- impurz
  permz.lst.dy[[i]] <- permz
  
  obsVspred.lst.dy[[i]] <- fit_pred_vals
  
  
} 


obsVspred.df.dy <- do.call(rbind, obsVspred.lst.dy)

obsVspred.df.dy$PFG <-as.factor(obsVspred.df.dy$PFG)


r_squared_values <- obsVspred.df.dy %>%
  group_by(PFG) %>%
  do({
    model <- lm(Predicted ~ Observed, data = .)
    data.frame(r_squared = summary(model)$r.squared,
               intercept = coef(model)[1],
               slope = coef(model)[2])
  })


p.dy <- ggplot(obsVspred.df.dy, aes(x = Observed, y = Predicted, color = PFG)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(data = r_squared_values,
            aes(x = max(obsVspred.df.dy$Observed) * 0.8,
                y = max(obsVspred.df.dy$Predicted) * 0.15 - (as.numeric(PFG) - 1) * max(obsVspred.df.dy$Predicted) * 0.25, # Adjust y for each group
                label = paste0("~R^{2} == ", round(r_squared, 2)),color=PFG),
            inherit.aes = FALSE, parse = TRUE) + # Prevent inheriting color from main plot
  labs(title = "(a)",
       x = "Observed",
       y = "Predicted") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_viridis_d() + theme(legend.position = "top")










# So, which variables do I want to use (DAILY):

namez_covars <- c("Avg.Temp_deg.C","Liquid.Precip_mm","ASCE.ETr_mm",
                  "_DTGW_m","_Temp_C","swT_m_gwT","Month_l0","season_l0")

#namez_lagz <- c("_l0","_l1","_l_1_4avg","_l_4_7avg")

#select_y <- c("deseas_std_NDVI","dekalm_std_NDVI")

select_y <- c("dekalm_std_NDVI")

namez_lagz <- c("_l0","_l1")

y_x_dy_fin2 <- list()


for(i in 1:length(dy_covs_x)){
  
  df1 <- dy_covs_x[[i]]
  
  y_df <- df1 %>%
    select(matches(select_y))
  
  # Select columns based on variables I wanted, but all the lags included
  # So, just pick the micromet, discharge gage and dtgw stuff at t=0
  
  df2 <- df1 %>%
    select(matches(namez_covars))
  
  
  df3 <- df2 %>%
    select(contains(namez_lagz))
  
  
  df4 <- cbind.data.frame(y_df, df3)
  
  
  y_x_dy_fin2[[i]] <- df4
  
  
  
}





obsVspred.lst.dy2 <- list()
impurz.lst.dy2 <- list()
permz.lst.dy2 <- list()
import.lst.dy2 <- list()

for(i in 1:length(y_x_dy_fin2)){
  
  split_mo <- initial_time_split(y_x_dy_fin2[[i]], prop = 4/5, lag = 1)
  #split_mo <- initial_split(data = y_x_dy_fin2[[i]], prop = 0.75)
  
  my_train  <- training(split_mo)
  my_test   <- testing(split_mo)
  
  
  
  # number of features
  n_features <- length(setdiff(names(my_train), "dekalm_std_NDVI" ))
  
  # train a default random forest model
  # THIS IS JUST A SIMPLE MODEL TO MAKE COMPARISON TOO....
  dummy1 <- ranger(
    dekalm_std_NDVI ~ ., 
    data = my_train,
    mtry = floor(n_features / 3),
    respect.unordered.factors = 'order', # this is recommended for categorical vars
    na.action = "na.omit",
    seed = seedz
  )
  
  # get OOB RMSE
  default_rmse <- sqrt(dummy1$prediction.error)
  
  
  
  
  
  
  # Make hyperperameter grid search to find optimal model
  
  # minimum node size is often use as 5 for regression problems
  # mtry should be centered on number-of-features divided by 3
  # sample fraction should be anywhere from 70% to 80%
  
  # create hyperparameter grid
  hyper_grid <- expand.grid(
    mtry = floor(seq(from=2,to=n_features,length.out=5)), # If many relevant predictors, lower m-try works better
    min.node.size = c(1, 3, 5, 10), # 5 for most regression problems
    replace = c(TRUE, FALSE),                               
    sample.fraction = c(0.7,0.75,0.80),                       
    rmse = NA                                               
  )
  
  #nrow(hyper_grid)
  
  
  
  # execute full cartesian grid search
  
  for(j in seq_len(nrow(hyper_grid))) {
    # fit model for ith hyperparameter combination
    fit <- ranger(
      formula         = dekalm_std_NDVI ~ ., 
      data            = my_train, 
      num.trees       = n_features * 10,
      mtry            = hyper_grid$mtry[j],
      min.node.size   = hyper_grid$min.node.size[j],
      replace         = hyper_grid$replace[j],
      sample.fraction = hyper_grid$sample.fraction[j],
      verbose         = TRUE,
      seed            = seedz,
      respect.unordered.factors = 'order',
      na.action = "na.omit",
      write.forest = TRUE
    )
    # export OOB error 
    hyper_grid$rmse[j] <- sqrt(fit$prediction.error)
  }
  
  
  
  
  
  
  # assess top 10 models
  top10 <- hyper_grid %>%
    arrange(rmse) %>%
    mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
    head(10)
  
  
  
  
  # For the ranger package...
  # Once you identify the optimal parameters from the grid search:
  # 1) Rerun the model with those hyperparameter values
  # 2) AND crank up the number of trees to create more stable estimates of feature importance
  
  
  # re-run model with impurity-based variable importance
  rf_impurity <- ranger(
    formula = dekalm_std_NDVI ~ ., 
    data = my_train, 
    num.trees = 2000, # can be exhaustive with this pick, now that best params are known
    mtry = top10$mtry[1],
    min.node.size = top10$min.node.size[1],
    sample.fraction = top10$sample.fraction[1],
    replace = top10$replace[1],
    importance = "impurity",
    respect.unordered.factors = "order",
    na.action = "na.omit",
    verbose = FALSE,
    seed  = seedz
  )
  
  # re-run model with permutation-based variable importance
  rf_permutation <- ranger(
    formula = dekalm_std_NDVI ~ ., 
    data = my_train, 
    num.trees = 2000,
    mtry = top10$mtry[1],
    min.node.size = top10$min.node.size[1],
    sample.fraction = top10$sample.fraction[1],
    replace = top10$replace[1],
    importance = "permutation",
    respect.unordered.factors = "order",
    na.action = "na.omit",
    verbose = FALSE,
    seed  = seedz
  )
  
  
  
  # plot of the top #n_features most important variables based on impurity and permutation
  
  
  impurz <- vip::vip(rf_impurity, num_features = 10, bar = FALSE,geom = "point") + labs(x="Predictor", y="Importance (Impurity)",subtitle="(b)")
  permz <- vip::vip(rf_permutation, num_features = 10, bar = FALSE, geom = "point") + labs(x="Predictor", y="Importance (Permutation)",subtitle="(b)")
  
  # Assuming 'model' is your trained tree-based model (e.g., from ranger, xgboost, rpart)
  importance_scores <- vi(rf_impurity, method = "model")
  sorted_importance <- importance_scores[order(importance_scores$Importance, decreasing = TRUE), ]
  top_10_features <- sorted_importance$Variable[1:10]
  
  
  
  fit_pred <- predict(fit, data = my_test)
  fit_pred_vals <- as.data.frame(fit_pred$predictions)
  fit_pred_vals$Observed <- my_test$dekalm_std_NDVI
  # fit_pred_vals$Month <- my_test$Month_l0
  # colnames(fit_pred_vals) <- c("Predicted","Observed","Month")
  
  colnames(fit_pred_vals) <- c("Predicted","Observed")
  
  fit_pred_vals$PFG <- veg.name[i]
  
  import.lst.dy2[[i]] <- top_10_features
  
  impurz.lst.dy2[[i]] <- impurz
  permz.lst.dy2[[i]] <- permz
  
  obsVspred.lst.dy2[[i]] <- fit_pred_vals
  
  
} 


obsVspred.df.dy2 <- do.call(rbind, obsVspred.lst.dy2)

obsVspred.df.dy2$PFG <-as.factor(obsVspred.df.dy2$PFG)


r_squared_values <- obsVspred.df.dy2 %>%
  group_by(PFG) %>%
  do({
    model <- lm(Predicted ~ Observed, data = .)
    data.frame(r_squared = summary(model)$r.squared,
               intercept = coef(model)[1],
               slope = coef(model)[2])
  })


p.dy2 <- ggplot(obsVspred.df.dy2, aes(x = Observed, y = Predicted, color = PFG)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(data = r_squared_values,
            aes(x = max(obsVspred.df.dy$Observed) * 0.8,
                y = max(obsVspred.df.dy$Predicted) * 0.15 - (as.numeric(PFG) - 1) * max(obsVspred.df.dy$Predicted) * 0.25, # Adjust y for each group
                label = paste0("~R^{2} == ", round(r_squared, 2)),color=PFG),
            inherit.aes = FALSE, parse = TRUE) + # Prevent inheriting color from main plot
  labs(title = "(b)",
       x = "Observed",
       y = "Predicted") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_color_viridis_d() + theme(legend.position = "top")



gridExtra::grid.arrange(p.dy, p.dy2, ncol = 2)


gridExtra::grid.arrange(impurz.lst.dy[[1]], impurz.lst.dy2[[1]],ncol = 2)












