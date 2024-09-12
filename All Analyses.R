# R code for the analysis of sex-ratio reaction norms and their correlations with nest temperature proxies

# Load the R packages used in this analysis

library(tidyverse)
library(phytools)
library(embryogrowth)
library(raster)
library(sf)
library(brms)
library(curl)

####################################################################################################

# Fit sex-ratio reaction norms to obtain Tpiv and TRT estimates

# Begin by reading in the ROSIE database

data <- read.csv("https://raw.githubusercontent.com/calebkrueger/ROSIE/main/Database.csv")

# Filter to include data from constant incubation experiments without chemical manipulations
# Exclude species with GSD, multiple instances of the same data, incomplete data, and studies using unreliable sexing techniques

data <- data %>% filter(SDM == "TSD",
                        Incubation_Setup == "Constant",
                        !Sexing_Method %in% c("Morphometrics", "Morphometrics, macroscopy"),
                        Data_Elsewhere %in% c(0,NA),
                        !Source == "Janzen, F. J., Unpublished data.", #Included in Krueger & Janzen 2023 but missing Data_Elsewhere == 1 tag
                        is.na(Chemical_Treatment),
                        !is.na(Mean_Temp),
                        !is.na(Males),
                        !is.na(Females))

# Remove irrelevant variables

data <- data %>% dplyr::select(Family,
                               Species,
                               Latitude,
                               Longitude,
                               Latitude_Range,
                               Longitude_Range,
                               Mean_Temp,
                               Males,
                               Females,
                               Source)

# Create new variable of species + coordinates to indicate unique "populations"

data <- data %>% mutate(sp.coords = as.character(interaction(Species, Latitude, Longitude, sep = ", ", drop = T)))

# Pelomedusa subrufa has no location information, creating NAs, so manually assign a sp.coords value

data[is.na(data$sp.coords),]$sp.coords <- "Pelomedusa_subrufa, NA, NA"
data$sp.coords <- as.factor(data$sp.coords)

# Identify and remove species with no mixed-sex data - cannot fit reaction norms

data %>% 
  mutate(spt = as.character(interaction(sp.coords, Mean_Temp, sep = ", ", drop = T))) %>%
  group_by(spt) %>%
  summarize(Males = sum(Males), Females = sum(Females), sp.coords = unique(sp.coords)) %>%
  group_by(spt) %>%
  summarize(Mix = length(Males) - length(Males[Males == 0]) - length(Females[Females == 0]), sp.coords = unique(sp.coords)) %>%
  filter(Mix > 0) -> keepers

data %>% filter(sp.coords %in% keepers$sp.coords) -> data

# Remove some populations with too few data to fit reaction norms

data %>% filter(!sp.coords == "Chrysemys_picta, 47.5, -95.21666667") -> data
data %>% filter(!sp.coords == "Chelydra_serpentina, 44.4305135, -79.226305") -> data
data %>% filter(!sp.coords == "Chelydra_serpentina, 46.443014, -93.3672395") -> data

# Exclude some data on multiple species

data <- data[!data$Species == "Graptemys_pseudogeographica/ouachitensis",]

# Remove some irrelevant data from lower temps so models can fit

data %>% filter(!(Species == "Pelomedusa_subrufa" & Mean_Temp < 26)) -> data
data %>% filter(!(Species == "Pelusios_castaneus" & Mean_Temp < 26)) -> data
data %>% filter(!(sp.coords == "Chelydra_serpentina, 40.7114365, -74.4717875" & Mean_Temp < 26)) -> data

# Remove erroneous data showing females at low temperatures in Chrysemys picta

data %>% filter(!(Species == "Chrysemys_picta" & Mean_Temp == 20)) -> data

# Match up data from Fort Ripley, MN with different location accuracies to avoid replicating population

data[data$sp.coords == "Chelydra_serpentina, 46.202449, -94.4147375", c("Latitude_Range","Longitude_Range")] <- c(rep(0.255304, 8), rep(0.163537, 8))

# Read in Thomson et al. (2021) phylogeny

tree <- read.tree("Thomson.nwk")

# Add missing species based on divergence estimates from timetree.org

tree <- bind.tip(tree,
                 tip.label = "Podocnemis_lewyana",
                 where = which(tree$tip.label == "Podocnemis_unifilis"),
                 position = 15.62)

tree <- ladderize(tree)

# Drop data for species not found in phylogeny

data <- data[data$Species %in% tree$tip.label,]

# Create dataframe to hold pivotal temperature and TRT values and credible intervals

data %>%
  dplyr::select(sp.coords,
                Family,
                Species,
                Latitude,
                Longitude,
                Latitude_Range,
                Longitude_Range) %>%
  distinct() %>%
  add_column(Tpiv = NA,
             Tpiv_SD = NA,
             TRT = NA,
             TRT_SD = NA,
             logTRT_SD = NA,
             TRT_Coverage = NA,
             prop = NA,
             resolution = NA,
             Source = NA) -> df

# Identify populations that need to be fit as a pattern Ia vs pattern II reaction norm
# NOTE: This does NOT mean these species are classified as this pattern, but rather that the available data best represent one pattern vs the other

forII <- c("Chelydra_serpentina, 43.537978, -82.349297",
           "Chelydra_serpentina, 45.591395, -78.5189905",
           "Chelydra_serpentina, 39.0422895, -86.2459755",
           "Chelydra_serpentina, 41.930585, -90.1112055",
           "Chelydra_serpentina, 44.065051, -79.49139",
           "Chelydra_serpentina, 47.1935815, -95.221055",
           "Chelydra_serpentina, 46.0272355, -86.621559",
           "Chelydra_serpentina, 32.6390315, -92.5755055",
           "Chelydra_serpentina, 29.9251295, -85.0121755",
           "Chelydra_serpentina, 26.4088845, -80.8418345",
           "Clemmys_guttata, 36.864652, -79.0698825",
           "Heosemys_grandis, 9.929448, 102.2701795",
           "Kinosternon_baurii, 26.4088845, -80.8418345",
           "Kinosternon_baurii, 29.9251295, -85.0121755",
           "Kinosternon_creaseri, 19.6773545, -88.741656",
           "Kinosternon_flavescens, 41.675582, -102.277737",
           "Kinosternon_leucostomum, 8.813487, -84.9906775",
           "Kinosternon_scorpioides, 19.019158, -93.0359735",
           "Kinosternon_stejnegeri, 30.150685, -111.4038485",
           "Kinosternon_subrubrum, 32.2377005, -94.1030405",
           "Kinosternon_subrubrum, 32.9534275, -82.9752375",
           "Macrochelys_temminckii, 29.9251295, -85.0121755",
           "Mauremys_sinensis, 22.899816, 113.6410245",
           "Melanochelys_trijuga, 18.267259, 86.746364",
           "Pelomedusa_subrufa, NA, NA",
           "Pelusios_castaneus, 1.397743, 2.266259",
           "Podocnemis_expansa, -9.9714295, -50.094228",
           "Rhinoclemmys_areolata, 18.36754, -91.8317355",
           "Sternotherus_carinatus, 32.28193, -93.2729625",
           "Sternotherus_minor, 29.9251295, -85.0121755",
           "Sternotherus_odoratus, 34.371102, -86.243556",
           "Sternotherus_odoratus, 39.0422895, -86.2459755",
           "Terrapene_carolina, 39.0422895, -86.2459755",
           "Terrapene_carolina, 27.7713485, -82.172115")

pops.II <- data %>% filter(sp.coords %in% forII)
pops.Ia <- data %>% filter(!sp.coords %in% forII)

# Fit reaction norms

par(mfrow = c(2,3))
set.seed(123)
for(i in unique(pops.II$sp.coords)){
  log<-tsd(pops.II[pops.II$sp.coords == i,],
           males = pops.II[pops.II$sp.coords == i, "Males"],
           females = pops.II[pops.II$sp.coords == i, "Females"], 
           temperatures = pops.II[pops.II$sp.coords == i, "Mean_Temp"],
           males.freq = F,
           equation = "logistic",
           parameters.initial = c(P_low = min(pops.II[pops.II$sp.coords == i, "Mean_Temp"]) + 2,
                                  S_low = 0.9,
                                  P_high = max(pops.II[pops.II$sp.coords == i, "Mean_Temp"]) - 2,
                                  S_high = -0.9),
           replicate.CI = 0)
  log.p <- tsd_MHmcmc_p(log, accept = T)
  log.p$Min[1] <- 0
  log.p$Max[1] <- 100
  log.p$Min[2] <- 0
  log.p$Max[2] <- 100
  log.p$Min[3] <- 0
  log.p$Max[3] <- 100
  log.p$Prior2[3] <- 10
  log.p$Prior1[4] <- 0
  log.p$Prior2[4] <- 10
  log.p$Min[4] <- -100
  log.p$Max[4] <- 0
  log.result <- tsd_MHmcmc(result = log, parametersMCMC = log.p, n.iter = 1e4, adaptive = T)
  CI <- P_TRT(x = log, resultmcmc = log.result, replicate.CI = 1e4)
  df[df$sp.coords == i, "Tpiv"] <- as.numeric(CI$P_TRT_quantiles[2,8])
  df[df$sp.coords == i, "Tpiv_SD"] <- as.numeric(sd(CI$P_TRT[,8]))
  df[df$sp.coords == i, "TRT"] <- as.numeric(CI$P_TRT_quantiles[2,7])
  df[df$sp.coords == i, "TRT_SD"] <- as.numeric(sd(CI$P_TRT[,7]))
  df[df$sp.coords == i, "logTRT_SD"] <- as.numeric(sd(log(CI$P_TRT[,7])))
  pops.II[pops.II$sp.coords == i,] %>%
    filter(Mean_Temp <= CI$P_TRT_quantiles[2,5]) %>%
    summarize(max(Mean_Temp)) %>%
    as.numeric() -> min
  if(min == -Inf) {min <- CI$P_TRT_quantiles[2,5]} else {}
  pops.II[pops.II$sp.coords == i,] %>%
    filter(Mean_Temp >= CI$P_TRT_quantiles[2,6]) %>%
    summarize(min(Mean_Temp)) %>%
    as.numeric() -> max
  if(max == Inf) {max <- CI$P_TRT_quantiles[2,6]} else {}
  pops.II[pops.II$sp.coords == i,] %>%
    filter(Mean_Temp <= max) %>%
    filter(Mean_Temp >= min) -> sub
  if(max(sub$Mean_Temp) > CI$P_TRT_quantiles[2,6] & min(sub$Mean_Temp) < CI$P_TRT_quantiles[2,5]) {
    df[df$sp.coords == i, "prop"] <- 1
  } else {
    if(max(sub$Mean_Temp) < CI$P_TRT_quantiles[2,6] & min(sub$Mean_Temp) > CI$P_TRT_quantiles[2,5]) {
      df[df$sp.coords == i, "prop"] <- (max(sub$Mean_Temp) - min(sub$Mean_Temp)) / (CI$P_TRT_quantiles[2,7])
    } else {
      if(max(sub$Mean_Temp) > CI$P_TRT_quantiles[2,6] & min(sub$Mean_Temp) > CI$P_TRT_quantiles[2,5]) {
        df[df$sp.coords == i, "prop"] <- (CI$P_TRT_quantiles[2,6] - min(sub$Mean_Temp)) / (CI$P_TRT_quantiles[2,7])
      } else {
        if(max(sub$Mean_Temp) < CI$P_TRT_quantiles[2,6] & min(sub$Mean_Temp) < CI$P_TRT_quantiles[2,5]) {
          df[df$sp.coords == i, "prop"] <- (max(sub$Mean_Temp) - CI$P_TRT_quantiles[2,5]) / (CI$P_TRT_quantiles[2,7])
        } else { df[df$sp.coords == i, "prop"] <- 0 }
      }
    }
  }
  pops.II[pops.II$sp.coords == i,] %>%
    filter(Mean_Temp <= max) %>%
    filter(Mean_Temp >= min) %>%
    summarize((length(unique(Mean_Temp)))/(max(Mean_Temp) - min(Mean_Temp))) %>%
    as.numeric() -> df[df$sp.coords == i, "resolution"]
  df[df$sp.coords == i, "TRT_Coverage"] <- df[df$sp.coords == i, "resolution"] * df[df$sp.coords == i, "prop"]
  pops.II[pops.II$sp.coords == i,] %>%
    arrange(Source) %>%
    summarize(Source = unique(Source)) -> sources
  paste(unlist(sources), collapse = " ") -> df[df$sp.coords == i, "Source"]
  # Code for visualizing the reaction norm fit
  #plot(log, resultmcmc = log.result, main = i, xlim = c(18,35), show.PTRT = F, show.observations = T)
  #abline(h = 0.50, lty = "dotted")
  #abline(v = as.numeric(CI$P_TRT_quantiles[1,5]), lty = "dashed", col = "red")
  #abline(v = as.numeric(CI$P_TRT_quantiles[2,5]), lwd = 2, col = "red")
  #abline(v = as.numeric(CI$P_TRT_quantiles[3,5]), lty = "dashed", col = "red")
  #abline(v = as.numeric(CI$P_TRT_quantiles[1,6]), lty = "dashed", col = "blue")
  #abline(v = as.numeric(CI$P_TRT_quantiles[2,6]), lwd = 2, col = "blue")
  #abline(v = as.numeric(CI$P_TRT_quantiles[3,6]), lty = "dashed", col = "blue")
  #abline(v = as.numeric(CI$P_TRT_quantiles[1,8]), lty = "dashed", col = "black")
  #abline(v = as.numeric(CI$P_TRT_quantiles[2,8]), lwd = 2, col = "black")
  #abline(v = as.numeric(CI$P_TRT_quantiles[3,8]), lty = "dashed", col = "black")
  #plot(log.result, parameters = "P_high", xlim = c(0, 50), legend = F)
  #plot(log.result, parameters = "S_high", xlim = c(-25, 0), legend = F)
}

set.seed(123)
for(i in unique(pops.Ia$sp.coords)){
  log<-tsd(pops.Ia[pops.Ia$sp.coords == i,],
           males = pops.Ia[pops.Ia$sp.coords == i, "Males"],
           females = pops.Ia[pops.Ia$sp.coords == i, "Females"], 
           temperatures = pops.Ia[pops.Ia$sp.coords == i, "Mean_Temp"],
           males.freq = F,
           equation = "logistic",
           replicate.CI = 0)
  log.p <- tsd_MHmcmc_p(log, accept = T)
  log.p$Min[1] <- 0
  log.p$Max[1] <- 100
  log.p$Min[2] <- -100
  log.p$Max[2] <- 0
  log.p$Prior2[1] <- 10
  log.p$Prior1[2] <- 0
  log.p$Prior2[2] <- 10
  log.result <- tsd_MHmcmc(result = log, parametersMCMC = log.p, n.iter = 1e4, adaptive = T)
  CI <- P_TRT(x = log, resultmcmc = log.result, replicate.CI = 1e4)
  df[df$sp.coords == i, "Tpiv"] <- as.numeric(CI$P_TRT_quantiles[2,4])
  df[df$sp.coords == i, "Tpiv_SD"] <- as.numeric(sd(CI$P_TRT[,4]))
  df[df$sp.coords == i, "TRT"] <- as.numeric(CI$P_TRT_quantiles[2,3])
  df[df$sp.coords == i, "TRT_SD"] <- as.numeric(sd(CI$P_TRT[,3]))
  df[df$sp.coords == i, "logTRT_SD"] <- as.numeric(sd(log(CI$P_TRT[,3])))
  pops.Ia[pops.Ia$sp.coords == i,] %>%
    filter(Mean_Temp <= CI$P_TRT_quantiles[2,1]) %>%
    summarize(max(Mean_Temp)) %>%
    as.numeric() -> min
  if(min == -Inf) {min <- CI$P_TRT_quantiles[2,1]} else {}
  pops.Ia[pops.Ia$sp.coords == i,] %>%
    filter(Mean_Temp >= CI$P_TRT_quantiles[2,2]) %>%
    summarize(min(Mean_Temp)) %>%
    as.numeric() -> max
  if(max == Inf) {max <- CI$P_TRT_quantiles[2,2]} else {}
  pops.Ia[pops.Ia$sp.coords == i,] %>%
    filter(Mean_Temp <= max) %>%
    filter(Mean_Temp >= min) -> sub
  if(max(sub$Mean_Temp) > CI$P_TRT_quantiles[2,2] & min(sub$Mean_Temp) < CI$P_TRT_quantiles[2,1]) {
    df[df$sp.coords == i, "prop"] <- 1
  } else {
    if(max(sub$Mean_Temp) < CI$P_TRT_quantiles[2,2] & min(sub$Mean_Temp) > CI$P_TRT_quantiles[2,1]) {
      df[df$sp.coords == i, "prop"] <- (max(sub$Mean_Temp) - min(sub$Mean_Temp)) / (CI$P_TRT_quantiles[2,3])
    } else {
      if(max(sub$Mean_Temp) > CI$P_TRT_quantiles[2,2] & min(sub$Mean_Temp) > CI$P_TRT_quantiles[2,1]) {
        df[df$sp.coords == i, "prop"] <- (CI$P_TRT_quantiles[2,2] - min(sub$Mean_Temp)) / (CI$P_TRT_quantiles[2,3])
      } else {
        if(max(sub$Mean_Temp) < CI$P_TRT_quantiles[2,2] & min(sub$Mean_Temp) < CI$P_TRT_quantiles[2,1]) {
          df[df$sp.coords == i, "prop"] <- (max(sub$Mean_Temp) - CI$P_TRT_quantiles[2,1]) / (CI$P_TRT_quantiles[2,3])
        } else { df[df$sp.coords == i, "prop"] <- 0 }
      }
    }
  }
  pops.Ia[pops.Ia$sp.coords == i,] %>%
    filter(Mean_Temp <= max) %>%
    filter(Mean_Temp >= min) %>%
    summarize((length(unique(Mean_Temp)))/(max(Mean_Temp) - min(Mean_Temp))) %>%
    as.numeric() -> df[df$sp.coords == i, "resolution"]
  df[df$sp.coords == i, "TRT_Coverage"] <- df[df$sp.coords == i, "resolution"] * df[df$sp.coords == i, "prop"]
  pops.Ia[pops.Ia$sp.coords == i,] %>%
    arrange(Source) %>%
    summarize(Source = unique(Source)) -> sources
  paste(unlist(sources), collapse = " ") -> df[df$sp.coords == i, "Source"]
  # Code for visualizing the reaction norm fit
  #plot(log, resultmcmc = log.result, main = i, xlim = c(18,35), show.PTRT = F, show.observations = T)
  #abline(h = 0.95, lty = "dotted")
  #bline(h = 0.05, lty = "dotted")
  #abline(v = as.numeric(CI$P_TRT_quantiles[1,1]), lty = "dashed", col = "red")
  #abline(v = as.numeric(CI$P_TRT_quantiles[2,1]), lwd = 2, col = "red")
  #abline(v = as.numeric(CI$P_TRT_quantiles[3,1]), lty = "dashed", col = "red")
  #abline(v = as.numeric(CI$P_TRT_quantiles[1,2]), lty = "dashed", col = "blue")
  #abline(v = as.numeric(CI$P_TRT_quantiles[2,2]), lwd = 2, col = "blue")
  #abline(v = as.numeric(CI$P_TRT_quantiles[3,2]), lty = "dashed", col = "blue")
  #abline(v = as.numeric(CI$P_TRT_quantiles[1,4]), lty = "dashed", col = "black")
  #abline(v = as.numeric(CI$P_TRT_quantiles[2,4]), lwd = 2, col = "black")
  #abline(v = as.numeric(CI$P_TRT_quantiles[3,4]), lty = "dashed", col = "black")
  #plot(log.result, parameters = "P", xlim = c(0, 50), legend = F)
  #plot(log.result, parameters = "S", xlim = c(-25, 0), legend = F)
}

# Add data on mean annual temperatures from CHELSA for each population

# First, download and aggregate the CHELSA temperature raster

curl_download("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio1_1981-2010_V.2.1.tif",
              "./CHELSA_MAT.tif")
tavg <- raster("CHELSA_MAT.tif")
tavg <- aggregate(tavg, 5)
tavg <- (tavg * 0.1) - 273.15

# Extract mean annual temperature at each population's coordinates

# Define bounding box around each population

df$Latitude_Max <- df$Latitude + (df$Latitude_Range/2)
df$Latitude_Min <- df$Latitude - (df$Latitude_Range/2)
df$Longitude_Max <- df$Longitude + (df$Longitude_Range/2)
df$Longitude_Min <- df$Longitude - (df$Longitude_Range/2)

df %>%
  dplyr::select(sp.coords, Latitude_Max, Latitude_Min, Longitude_Max, Longitude_Min) %>%
  dplyr::filter(!is.na(Latitude_Max)) -> lat.long

coords <- st_as_sf(lat.long,
                   coords = c("Longitude_Max", "Latitude_Max"),
                   crs = 4326)[,c(1,4)]
coords <- rbind(coords, st_as_sf(lat.long,
                                 coords = c("Longitude_Min", "Latitude_Max"),
                                 crs = 4326)[,c(1,4)])
coords <- rbind(coords, st_as_sf(lat.long,
                                 coords = c("Longitude_Max", "Latitude_Min"),
                                 crs = 4326)[,c(1,4)])
coords <- rbind(coords, st_as_sf(lat.long,
                                 coords = c("Longitude_Min", "Latitude_Min"),
                                 crs = 4326)[,c(1,4)])

coords %>%
  group_by(sp.coords) %>%
  summarise() %>%
  st_cast("POLYGON") -> polys

polys <- st_transform(polys, crs(tavg))

# Extract mean and SD of mean annual temperature values within bounding box

polys$Temps <- as.vector(unlist(raster::extract(tavg, polys, fun = mean, na.rm = TRUE), use.names = FALSE))
polys$Temps_SD <- as.vector(unlist(raster::extract(tavg, polys, fun = sd, na.rm = TRUE), use.names = FALSE))

polys %>%
  base::as.data.frame() %>%
  dplyr::select(sp.coords, Temps, Temps_SD) %>%
  left_join(x = df, by = "sp.coords") -> df

df[is.na(df$Temps_SD),"Temps_SD"] <- 0

# Log-transform select variables for analysis

df$logTRT <- log(df$TRT)
df$logCoverage <- log(df$TRT_Coverage)

# Create two data frames - one of populations with well-defined Tpiv and another with well-defined TRT

# List of species to drop from Tpiv analyses

tpiv <- c("Kinosternon_baurii, 29.9251295, -85.0121755",
          "Kinosternon_creaseri, 19.6773545, -88.741656",
          "Kinosternon_leucostomum, 8.813487, -84.9906775",
          "Terrapene_carolina, 27.7713485, -82.172115",
          "Astrochelys_radiata, -24.4213925, 45.2127545",
          "Caretta_caretta, -23.7442075, 113.5620435",
          "Chelonia_mydas, 5.993292, -54.8125825",
          "Chelonia_mydas, 20.5694525, -157.6830405",
          "Chelonoidis_niger, -1.377946, -89.6822405",
          "Chelydra_serpentina, 46.202449, -94.4147375",
          "Chelydra_serpentina, 45.575434, -78.683943",
          "Chrysemys_picta, 46.0272355, -86.621559",
          "Chrysemys_picta, 37.223365, -80.484596",
          "Geochelone_elegans, 16.7827945, 76.7525175",
          "Geochelone_platynota, 22.00165, 95.308987",
          "Malaclemys_terrapin, 39.162667, -74.6824625",
          "Manouria_impressa, 15.525993, 101.089441",
          "Peltocephalus_dumerilianus, -0.406805556, -63.44736111",
          "Podocnemis_expansa, -1.008516, -71.64148",
          "Podocnemis_lewyana, 7.8470115, -74.437783",
          "Podocnemis_unifilis, -1.4369195, -70.710111",
          "Podocnemis_unifilis, -11.955611, -53.509643",
          "Podocnemis_vogli, 4.903391, -70.4720135",
          "Pseudemys_texana, 31.0514145, -98.7093235",
          "Stigmochelys_pardalis, -11.916498, 30.377801",
          "Terrapene_carolina, 38.905738, -77.033739",
          "Testudo_graeca, 32.656925, 35.481161",
          "Trachemys_scripta, 30.3259625, -97.7720965",
          "Trachemys_scripta, 30.5040795, -90.4592465")

# List of species to drop for TRT analyses

trt <- c("Chelydra_serpentina, 44.065051, -79.49139",
         "Chelydra_serpentina, 29.9251295, -85.0121755",
         "Chelydra_serpentina, 26.4088845, -80.8418345",
         "Clemmys_guttata, 36.864652, -79.0698825",
         "Heosemys_grandis, 9.929448, 102.2701795",
         "Kinosternon_baurii, 26.4088845, -80.8418345",
         "Kinosternon_baurii, 29.9251295, -85.0121755",
         "Kinosternon_creaseri, 19.6773545, -88.741656",
         "Kinosternon_flavescens, 41.675582, -102.277737",
         "Kinosternon_leucostomum, 8.813487, -84.9906775",
         "Kinosternon_subrubrum, 32.2377005, -94.1030405",
         "Kinosternon_subrubrum, 32.9534275, -82.9752375",
         "Macrochelys_temminckii, 29.9251295, -85.0121755",
         "Melanochelys_trijuga, 18.267259, 86.746364",
         "Pelomedusa_subrufa, NA, NA",
         "Pelusios_castaneus, 1.397743, 2.266259",
         "Rhinoclemmys_areolata, 18.36754, -91.8317355",
         "Sternotherus_minor, 29.9251295, -85.0121755",
         "Sternotherus_odoratus, 34.371102, -86.243556",
         "Terrapene_carolina, 39.0422895, -86.2459755",
         "Terrapene_carolina, 27.7713485, -82.172115",
         "Astrochelys_radiata, -24.4213925, 45.2127545",
         "Caretta_caretta, -23.7442075, 113.5620435",
         "Chelonia_mydas, 5.993292, -54.8125825",
         "Chelonia_mydas, 20.5694525, -157.6830405",
         "Chelonoidis_niger, -1.377946, -89.6822405",
         "Chelydra_serpentina, 46.202449, -94.4147375",
         "Chelydra_serpentina, 45.575434, -78.683943",
         "Chrysemys_picta, 46.0272355, -86.621559",
         "Chrysemys_picta, 37.223365, -80.484596",
         "Geochelone_elegans, 16.7827945, 76.7525175",
         "Geochelone_platynota, 22.00165, 95.308987",
         "Kinosternon_leucostomum, 17.4911865, -93.04148575", #Removed for being a pattern II species without a defined TRT
         "Malaclemys_terrapin, 39.162667, -74.6824625",
         "Malayemys_macrocephala, 14.3532445, 100.447429",
         "Manouria_impressa, 15.525993, 101.089441",
         "Peltocephalus_dumerilianus, -0.406805556, -63.44736111",
         "Podocnemis_expansa, -1.008516, -71.64148",
         "Podocnemis_lewyana, 7.8470115, -74.437783",
         "Podocnemis_unifilis, -1.4369195, -70.710111",
         "Podocnemis_unifilis, -11.955611, -53.509643",
         "Podocnemis_vogli, 4.903391, -70.4720135",
         "Pseudemys_texana, 31.0514145, -98.7093235",
         "Stigmochelys_pardalis, -11.916498, 30.377801",
         "Terrapene_carolina, 38.905738, -77.033739",
         "Testudo_graeca, 32.656925, 35.481161",
         "Trachemys_scripta, 30.3259625, -97.7720965",
         "Trachemys_scripta, 30.5040795, -90.4592465")

df.tpiv <- df %>% filter(!sp.coords %in% tpiv)
df.trt <- df %>% filter(!sp.coords %in% trt)

tree.tpiv <- drop.tip(tree, setdiff(tree$tip.label, df.tpiv$Species))
tree.trt <- drop.tip(tree, setdiff(tree$tip.label, df.trt$Species))
tree.tpiv <- ladderize(tree.tpiv)
tree.trt <- ladderize(tree.trt)

# Calculate species mean and mean-centered values for coverage

df.trt$smlogCoverage <- sapply(split(df.trt$logCoverage, df.trt$Species), mean)[df.trt$Species]
df.trt$sclogCoverage <- df.trt$logCoverage - df.trt$smlogCoverage

df.tpiv$Source <- as.numeric(as.factor(df.tpiv$Source))
df.trt$Source <- as.numeric(as.factor(df.trt$Source))

rm(list = c("CI","coords","data","keepers","lat.long","log","log.p","log.result","polys","pops.Ia","pops.II","sources","sub","tavg","forII","max","min"))

####################################################################################################

# First, assess relationship between TRT coverage and TRT estimates
# Also calculate phylogenetic heritability of TRT

df.trt %>%
  mutate(Species2 = Species) -> mcmc.data

vcv <- vcv.phylo(tree.trt, corr = T)

bf1 <- brmsformula(logTRT | se(logTRT_SD, sigma = TRUE) ~ smlogCoverage + sclogCoverage + (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))

prior <- c(prior(normal(0, 5), "b"),
           prior(normal(0, 5), "Intercept"),
           prior(student_t(3, 0, 2.5), "sd"),
           prior(student_t(3, 0, 2.5), "sigma"))

trt <- brm(bf1,
           data = mcmc.data,
           data2 = list(vcv = vcv),
           prior = prior,
           chains = 4,
           cores = 4,
           iter = 5000,
           thin = 10,
           warmup = 1000,
           control = list(adapt_delta = 0.99),
           seed = 123)

mcmc_plot(trt, type = "trace")
mcmc_plot(trt, type = "acf")
mcmc_plot(trt, type = "rhat")
mcmc_plot(trt, type = "neff")

rstan::check_hmc_diagnostics(trt$fit)

pp_check(trt, type = "intervals")
pp_check(trt, ndraws = 25)

summary(trt)

# Calculate the phylogenetic heritability of logTRT

trt %>%
  as_tibble() %>%
  dplyr::select(sigb = sd_Species__Intercept,
                sigb2 = sd_Species2__Intercept,
                sigb3 = sd_Source__Intercept,
                sige = sigma) %>%
  mutate(h2 = sigb^2 / (sigb^2 + sigb2^2 + sigb3^2 + sige^2)) %>%
  pull(h2) %>%
  quantile(probs = c(0.025,0.5,0.975))

####################################################################################################

# Calculate phylogenetic heritability of Tpiv

df.tpiv %>%
  mutate(Species2 = Species) -> mcmc.data

vcv <- vcv.phylo(tree.tpiv, corr = T)

bf1 <- brmsformula(Tpiv | se(Tpiv_SD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))

prior <- c(prior(normal(30, 5), "Intercept"),
           prior(student_t(3, 0, 2.5), "sd"),
           prior(student_t(3, 0, 2.5), "sigma"))

tpiv <- brm(bf1,
            data = mcmc.data,
            data2 = list(vcv = vcv),
            prior = prior,
            chains = 4,
            cores = 4,
            iter = 5000,
            thin = 10,
            warmup = 1000,
            control = list(adapt_delta = 0.99),
            seed = 123)

mcmc_plot(tpiv, type = "trace")
mcmc_plot(tpiv, type = "acf")
mcmc_plot(tpiv, type = "rhat")
mcmc_plot(tpiv, type = "neff")

rstan::check_hmc_diagnostics(tpiv$fit)

pp_check(tpiv, type = "intervals")
pp_check(tpiv, ndraws = 25)

summary(tpiv)

tpiv %>%
  as_tibble() %>%
  dplyr::select(sigb = sd_Species__Intercept,
                sigb2 = sd_Species2__Intercept,
                sigb3 = sd_Source__Intercept,
                sige = sigma) %>%
  mutate(h2 = sigb^2 / (sigb^2 + sigb2^2 + sigb3^2 + sige^2)) %>%
  pull(h2) %>%
  quantile(probs = c(0.025,0.5,0.975))

####################################################################################################

# Fit model to check for correlations between Tpiv and logTRT

df.trt %>%
  mutate(Species2 = Species) -> mcmc.data

vcv <- vcv.phylo(tree.trt, corr = T)

bf1 <- brmsformula(logTRT | se(logTRT_SD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))
bf2 <- brmsformula(logCoverage ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))
bf3 <- brmsformula(Tpiv | se(Tpiv_SD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))

prior <- c(prior(normal(0, 5), "Intercept", resp = "logCoverage"),
           prior(normal(0, 5), "Intercept", resp = "logTRT"),
           prior(normal(30, 5), "Intercept", resp = "Tpiv"),
           prior(student_t(3, 0, 2.5), "sd", resp = "logCoverage"),
           prior(student_t(3, 0, 2.5), "sd", resp = "logTRT"),
           prior(student_t(3, 0, 2.5), "sd", resp = "Tpiv"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "logCoverage"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "logTRT"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "Tpiv"))

tpiv.trt <- brm(bf1 + bf2 + bf3 + set_rescor(TRUE),
                data = mcmc.data,
                data2 = list(vcv = vcv),
                prior = prior,
                chains = 4,
                cores = 4,
                iter = 5000,
                thin = 10,
                warmup = 1000,
                control = list(adapt_delta = 0.99),
                seed = 123)

mcmc_plot(tpiv.trt, type = "trace")
mcmc_plot(tpiv.trt, type = "acf")
mcmc_plot(tpiv.trt, type = "rhat")
mcmc_plot(tpiv.trt, type = "neff")

rstan::check_hmc_diagnostics(tpiv.trt$fit)

pp_check(tpiv.trt, type = "intervals", resp = "Tpiv")
pp_check(tpiv.trt, ndraws = 25, resp = "Tpiv")
pp_check(tpiv.trt, type = "intervals", resp = "logTRT")
pp_check(tpiv.trt, ndraws = 25, resp = "logTRT")
pp_check(tpiv.trt, type = "intervals", resp = "logCoverage")
pp_check(tpiv.trt, ndraws = 25, resp = "logCoverage")

summary(tpiv.trt)

as.data.frame(tpiv.trt) -> tpiv.trt.df

tpiv.trt.df$pc.TRT.Coverage <- (tpiv.trt.df$rescor__logTRT__logCoverage - (tpiv.trt.df$rescor__logCoverage__Tpiv*tpiv.trt.df$rescor__logTRT__Tpiv)) / sqrt((1-tpiv.trt.df$rescor__logCoverage__Tpiv^2)*(1-tpiv.trt.df$rescor__logTRT__Tpiv^2))
tpiv.trt.df$pc.TRT.Tpiv <- (tpiv.trt.df$rescor__logTRT__Tpiv - (tpiv.trt.df$rescor__logTRT__logCoverage*tpiv.trt.df$rescor__logCoverage__Tpiv)) / sqrt((1-tpiv.trt.df$rescor__logCoverage__Tpiv^2)*(1-tpiv.trt.df$rescor__logTRT__logCoverage^2))
tpiv.trt.df$pc.Tpiv.Coverage <- (tpiv.trt.df$rescor__logCoverage__Tpiv - (tpiv.trt.df$rescor__logTRT__logCoverage*tpiv.trt.df$rescor__logTRT__Tpiv)) / sqrt((1-tpiv.trt.df$rescor__logTRT__logCoverage^2)*(1-tpiv.trt.df$rescor__logTRT__Tpiv^2))
tpiv.trt.df$ppc.TRT.Coverage <- (tpiv.trt.df$cor_Species__logTRT_Intercept__logCoverage_Intercept - (tpiv.trt.df$cor_Species__logTRT_Intercept__Tpiv_Intercept*tpiv.trt.df$cor_Species__logCoverage_Intercept__Tpiv_Intercept)) / sqrt((1-tpiv.trt.df$cor_Species__logTRT_Intercept__Tpiv_Intercept^2)*(1-tpiv.trt.df$cor_Species__logCoverage_Intercept__Tpiv_Intercept^2))
tpiv.trt.df$ppc.TRT.Tpiv <- (tpiv.trt.df$cor_Species__logTRT_Intercept__Tpiv_Intercept - (tpiv.trt.df$cor_Species__logTRT_Intercept__logCoverage_Intercept*tpiv.trt.df$cor_Species__logCoverage_Intercept__Tpiv_Intercept)) / sqrt((1-tpiv.trt.df$cor_Species__logTRT_Intercept__logCoverage_Intercept^2)*(1-tpiv.trt.df$cor_Species__logCoverage_Intercept__Tpiv_Intercept^2))
tpiv.trt.df$ppc.Tpiv.Coverage <- (tpiv.trt.df$cor_Species__logCoverage_Intercept__Tpiv_Intercept - (tpiv.trt.df$cor_Species__logTRT_Intercept__logCoverage_Intercept*tpiv.trt.df$cor_Species__logTRT_Intercept__Tpiv_Intercept)) / sqrt((1-tpiv.trt.df$cor_Species__logTRT_Intercept__logCoverage_Intercept^2)*(1-tpiv.trt.df$cor_Species__logTRT_Intercept__Tpiv_Intercept^2))

apply(tpiv.trt.df[,630:635], 2, quantile, probs = c(0.025,0.5,0.975), na.rm = T)

####################################################################################################

# Fit model to check for correlations between logTRT and mean annual temps

df.trt %>%
  filter(!is.na(Temps)) %>%
  filter(Latitude_Range < 2) %>%
  mutate(Species2 = Species) -> mcmc.data

tree.mcmc <- drop.tip(tree.trt, setdiff(tree.trt$tip.label, mcmc.data$Species))
vcv <- vcv.phylo(tree.mcmc, corr = T)

bf1 <- brmsformula(logTRT | se(logTRT_SD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))
bf2 <- brmsformula(logCoverage ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))
bf3 <- brmsformula(Temps | se(Temps_SD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))

prior <- c(prior(normal(0, 5), "Intercept", resp = "logCoverage"),
           prior(normal(0, 5), "Intercept", resp = "logTRT"),
           prior(normal(20, 5), "Intercept", resp = "Temps"),
           prior(student_t(3, 0, 2.5), "sd", resp = "logCoverage"),
           prior(student_t(3, 0, 2.5), "sd", resp = "logTRT"),
           prior(student_t(3, 0, 2.5), "sd", resp = "Temps"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "logCoverage"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "logTRT"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "Temps"))

trt.temps <- brm(bf1 + bf2 + bf3 + set_rescor(TRUE),
                 data = mcmc.data,
                 data2 = list(vcv = vcv),
                 prior = prior,
                 chains = 4,
                 cores = 4,
                 iter = 5000,
                 thin = 10,
                 warmup = 1000,
                 control = list(adapt_delta = 0.99),
                 seed = 123)

mcmc_plot(trt.temps, type = "trace")
mcmc_plot(trt.temps, type = "acf")
mcmc_plot(trt.temps, type = "rhat")
mcmc_plot(trt.temps, type = "neff")

rstan::check_hmc_diagnostics(trt.temps$fit)

pp_check(trt.temps, type = "intervals", resp = "logTRT")
pp_check(trt.temps, ndraws = 25, resp = "logTRT")
pp_check(trt.temps, type = "intervals", resp = "Temps")
pp_check(trt.temps, ndraws = 25, resp = "Temps")
pp_check(trt.temps, type = "intervals", resp = "logCoverage")
pp_check(trt.temps, ndraws = 25, resp = "logCoverage")

summary(trt.temps)

as.data.frame(trt.temps) -> trt.temps.df

trt.temps.df$pc.TRT.Coverage <- (trt.temps.df$rescor__logTRT__logCoverage - (trt.temps.df$rescor__logCoverage__Temps*trt.temps.df$rescor__logTRT__Temps)) / sqrt((1-trt.temps.df$rescor__logCoverage__Temps^2)*(1-trt.temps.df$rescor__logTRT__Temps^2))
trt.temps.df$pc.TRT.Temps <- (trt.temps.df$rescor__logTRT__Temps - (trt.temps.df$rescor__logTRT__logCoverage*trt.temps.df$rescor__logCoverage__Temps)) / sqrt((1-trt.temps.df$rescor__logCoverage__Temps^2)*(1-trt.temps.df$rescor__logTRT__logCoverage^2))
trt.temps.df$pc.Temps.Coverage <- (trt.temps.df$rescor__logCoverage__Temps - (trt.temps.df$rescor__logTRT__logCoverage*trt.temps.df$rescor__logTRT__Temps)) / sqrt((1-trt.temps.df$rescor__logTRT__logCoverage^2)*(1-trt.temps.df$rescor__logTRT__Temps^2))
trt.temps.df$ppc.TRT.Coverage <- (trt.temps.df$cor_Species__logTRT_Intercept__logCoverage_Intercept - (trt.temps.df$cor_Species__logTRT_Intercept__Temps_Intercept*trt.temps.df$cor_Species__logCoverage_Intercept__Temps_Intercept)) / sqrt((1-trt.temps.df$cor_Species__logTRT_Intercept__Temps_Intercept^2)*(1-trt.temps.df$cor_Species__logCoverage_Intercept__Temps_Intercept^2))
trt.temps.df$ppc.TRT.Temps <- (trt.temps.df$cor_Species__logTRT_Intercept__Temps_Intercept - (trt.temps.df$cor_Species__logTRT_Intercept__logCoverage_Intercept*trt.temps.df$cor_Species__logCoverage_Intercept__Temps_Intercept)) / sqrt((1-trt.temps.df$cor_Species__logTRT_Intercept__logCoverage_Intercept^2)*(1-trt.temps.df$cor_Species__logCoverage_Intercept__Temps_Intercept^2))
trt.temps.df$ppc.Temps.Coverage <- (trt.temps.df$cor_Species__logCoverage_Intercept__Temps_Intercept - (trt.temps.df$cor_Species__logTRT_Intercept__logCoverage_Intercept*trt.temps.df$cor_Species__logTRT_Intercept__Temps_Intercept)) / sqrt((1-trt.temps.df$cor_Species__logTRT_Intercept__logCoverage_Intercept^2)*(1-trt.temps.df$cor_Species__logTRT_Intercept__Temps_Intercept^2))

apply(trt.temps.df[,429:434], 2, quantile, probs = c(0.025,0.5,0.975), na.rm = T)

####################################################################################################

# Fit model to check for correlations between Tpiv and mean annual temps

df.tpiv %>%
  filter(!is.na(Temps)) %>%
  filter(Latitude_Range < 2) %>%
  mutate(Species2 = Species) -> mcmc.data

tree.mcmc <- drop.tip(tree.tpiv, setdiff(tree.tpiv$tip.label, mcmc.data$Species))
vcv <- vcv.phylo(tree.mcmc, corr = T)

bf1 <- brmsformula(Tpiv | se(Tpiv_SD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))
bf2 <- brmsformula(Temps | se(Temps_SD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))

prior <- c(prior(normal(30, 5), "Intercept", resp = "Tpiv"),
           prior(normal(20, 5), "Intercept", resp = "Temps"),
           prior(student_t(3, 0, 2.5), "sd", resp = "Tpiv"),
           prior(student_t(3, 0, 2.5), "sd", resp = "Temps"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "Tpiv"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "Temps"))

tpiv.temps <- brm(bf1 + bf2 + set_rescor(TRUE),
                  data = mcmc.data,
                  data2 = list(vcv = vcv),
                  prior = prior,
                  chains = 4,
                  cores = 4,
                  iter = 5000,
                  thin = 10,
                  warmup = 1000,
                  control = list(adapt_delta = 0.99),
                  seed = 123)

mcmc_plot(tpiv.temps, type = "trace")
mcmc_plot(tpiv.temps, type = "acf")
mcmc_plot(tpiv.temps, type = "rhat")
mcmc_plot(tpiv.temps, type = "neff")

rstan::check_hmc_diagnostics(tpiv.temps$fit)

pp_check(tpiv.temps, type = "intervals", resp = "Tpiv")
pp_check(tpiv.temps, ndraws = 25, resp = "Tpiv")
pp_check(tpiv.temps, type = "intervals", resp = "Temps")
pp_check(tpiv.temps, ndraws = 25, resp = "Temps")

summary(tpiv.temps)

####################################################################################################

# Plot data

png(filename = "Tpiv MAT Correlation.png", res = 300, width = 2900, height = 2100)

df.tpiv %>%
  filter(!is.na(Temps)) %>%
  filter(Latitude_Range < 2) %>%
  ggplot(aes(x = Temps, y = Tpiv, color = Family, group = Species)) +
  geom_point(size = 3.5) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              se = F,
              show.legend = FALSE,
              linewidth = 1.25) +
  theme_classic() +
  labs(x = expression(Mean~Annual~Temperature~"("*degree*"C)"),
       y = expression(T[piv]~"("*degree*"C)")) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        legend.position = c(0.13,0.8)) +
  scale_x_continuous(breaks = c(5,10,15,20,25,30)) +
  scale_color_brewer(palette = "Paired")

dev.off()

png(filename = "TRT Coverage Correlation.png", res = 300, width = 2900, height = 2100)

df.trt %>%
  ggplot(aes(x = logCoverage, y = logTRT, color = Temps)) +
  geom_point(size = 3.5) +
  theme_classic() +
  labs(x = expression("log(Coverage)"),
       y = expression("log(TRT)"),
       color = expression("MAT ("*degree*"C)")) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        legend.position = c(0.1,0.15)) +
  scale_color_gradientn(colors = c("#1F78B4","#A6CEE3","#FB9A99","#E31A1C")) 

dev.off()

## END OF MAIN ANALYSES

####################################################################################################

### SUPPLEMENTAL ANALYSES ###
## Examine influence of temperature data on correlations

# Read in temp during warmest quarter from CHELSA
# and mean annual soil temps (0-5cm) from Soil Temp Database

curl_download("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio10_1981-2010_V.2.1.tif",
              "./CHELSA_QTEMPS.tif")
curl_download("https://zenodo.org/records/7134169/files/SBIO1_0_5cm_Annual_Mean_Temperature.tif",
              "./SOILTEMPS.tif")
curl_download("https://zenodo.org/records/7134169/files/PCA_int_ext_0_5cm.tif",
              "./SOILTEMPSMASK.tif")

tq <- raster("CHELSA_QTEMPS.tif")
stavg <- raster("SOILTEMPS.tif")
mask <- raster("SOILTEMPSMASK.tif")

tq <- aggregate(tq, 5)
stavg <- aggregate(stavg, 5)
mask <- aggregate(mask, 5)

mask[mask < 0.9] <- NA

stavg <- mask(stavg, mask)

tq <- (tq * 0.1) - 273.15

# Extract temperature at each population's coordinates

df$Latitude_Max <- df$Latitude + (df$Latitude_Range/2)
df$Latitude_Min <- df$Latitude - (df$Latitude_Range/2)
df$Longitude_Max <- df$Longitude + (df$Longitude_Range/2)
df$Longitude_Min <- df$Longitude - (df$Longitude_Range/2)

df %>%
  dplyr::select(sp.coords, Latitude_Max, Latitude_Min, Longitude_Max, Longitude_Min) %>%
  dplyr::filter(!is.na(Latitude_Max)) -> lat.long

coords <- st_as_sf(lat.long,
                   coords = c("Longitude_Max", "Latitude_Max"),
                   crs = 4326)[,c(1,4)]
coords <- rbind(coords, st_as_sf(lat.long,
                                 coords = c("Longitude_Min", "Latitude_Max"),
                                 crs = 4326)[,c(1,4)])
coords <- rbind(coords, st_as_sf(lat.long,
                                 coords = c("Longitude_Max", "Latitude_Min"),
                                 crs = 4326)[,c(1,4)])
coords <- rbind(coords, st_as_sf(lat.long,
                                 coords = c("Longitude_Min", "Latitude_Min"),
                                 crs = 4326)[,c(1,4)])

coords %>%
  group_by(sp.coords) %>%
  summarise() %>%
  st_cast("POLYGON") -> polys

polys <- st_transform(polys, crs(tq))

polys$qTemps <- as.vector(unlist(raster::extract(tq, polys, fun = mean, na.rm = TRUE), use.names = FALSE))
polys$qTemps_SD <- as.vector(unlist(raster::extract(tq, polys, fun = sd, na.rm = TRUE), use.names = FALSE))
polys$SoilTemps <- as.vector(unlist(raster::extract(stavg, polys, fun = mean, na.rm = TRUE), use.names = FALSE))
polys$SoilTemps_SD <- as.vector(unlist(raster::extract(stavg, polys, fun = sd, na.rm = TRUE), use.names = FALSE))

polys %>%
  base::as.data.frame() %>%
  dplyr::select(sp.coords, qTemps, qTemps_SD, SoilTemps, SoilTemps_SD) %>%
  left_join(x = df, by = "sp.coords") -> df

df[is.na(df$qTemps_SD),"qTemps_SD"] <- 0
df[is.na(df$SoilTemps_SD),"SoilTemps_SD"] <- 0

# List of species to drop from Tpiv analyses

tpiv <- c("Kinosternon_baurii, 29.9251295, -85.0121755",
          "Kinosternon_creaseri, 19.6773545, -88.741656",
          "Kinosternon_leucostomum, 8.813487, -84.9906775",
          "Terrapene_carolina, 27.7713485, -82.172115",
          "Astrochelys_radiata, -24.4213925, 45.2127545",
          "Caretta_caretta, -23.7442075, 113.5620435",
          "Chelonia_mydas, 5.993292, -54.8125825",
          "Chelonia_mydas, 20.5694525, -157.6830405",
          "Chelonoidis_niger, -1.377946, -89.6822405",
          "Chelydra_serpentina, 46.202449, -94.4147375",
          "Chelydra_serpentina, 45.575434, -78.683943",
          "Chrysemys_picta, 46.0272355, -86.621559",
          "Chrysemys_picta, 37.223365, -80.484596",
          "Geochelone_elegans, 16.7827945, 76.7525175",
          "Geochelone_platynota, 22.00165, 95.308987",
          "Malaclemys_terrapin, 39.162667, -74.6824625",
          "Manouria_impressa, 15.525993, 101.089441",
          "Peltocephalus_dumerilianus, -0.406805556, -63.44736111",
          "Podocnemis_expansa, -1.008516, -71.64148",
          "Podocnemis_lewyana, 7.8470115, -74.437783",
          "Podocnemis_unifilis, -1.4369195, -70.710111",
          "Podocnemis_unifilis, -11.955611, -53.509643",
          "Podocnemis_vogli, 4.903391, -70.4720135",
          "Pseudemys_texana, 31.0514145, -98.7093235",
          "Stigmochelys_pardalis, -11.916498, 30.377801",
          "Terrapene_carolina, 38.905738, -77.033739",
          "Testudo_graeca, 32.656925, 35.481161",
          "Trachemys_scripta, 30.3259625, -97.7720965",
          "Trachemys_scripta, 30.5040795, -90.4592465")

# List of species to drop for TRT analyses

trt <- c("Chelydra_serpentina, 44.065051, -79.49139",
         "Chelydra_serpentina, 29.9251295, -85.0121755",
         "Chelydra_serpentina, 26.4088845, -80.8418345",
         "Clemmys_guttata, 36.864652, -79.0698825",
         "Heosemys_grandis, 9.929448, 102.2701795",
         "Kinosternon_baurii, 26.4088845, -80.8418345",
         "Kinosternon_baurii, 29.9251295, -85.0121755",
         "Kinosternon_creaseri, 19.6773545, -88.741656",
         "Kinosternon_flavescens, 41.675582, -102.277737",
         "Kinosternon_leucostomum, 8.813487, -84.9906775",
         "Kinosternon_subrubrum, 32.2377005, -94.1030405",
         "Kinosternon_subrubrum, 32.9534275, -82.9752375",
         "Macrochelys_temminckii, 29.9251295, -85.0121755",
         "Melanochelys_trijuga, 18.267259, 86.746364",
         "Pelomedusa_subrufa, NA, NA",
         "Pelusios_castaneus, 1.397743, 2.266259",
         "Rhinoclemmys_areolata, 18.36754, -91.8317355",
         "Sternotherus_minor, 29.9251295, -85.0121755",
         "Sternotherus_odoratus, 34.371102, -86.243556",
         "Terrapene_carolina, 39.0422895, -86.2459755",
         "Terrapene_carolina, 27.7713485, -82.172115",
         "Astrochelys_radiata, -24.4213925, 45.2127545",
         "Caretta_caretta, -23.7442075, 113.5620435",
         "Chelonia_mydas, 5.993292, -54.8125825",
         "Chelonia_mydas, 20.5694525, -157.6830405",
         "Chelonoidis_niger, -1.377946, -89.6822405",
         "Chelydra_serpentina, 46.202449, -94.4147375",
         "Chelydra_serpentina, 45.575434, -78.683943",
         "Chrysemys_picta, 46.0272355, -86.621559",
         "Chrysemys_picta, 37.223365, -80.484596",
         "Geochelone_elegans, 16.7827945, 76.7525175",
         "Geochelone_platynota, 22.00165, 95.308987",
         "Kinosternon_leucostomum, 17.4911865, -93.04148575", #Removed for being a pattern II species without a defined TRT
         "Malaclemys_terrapin, 39.162667, -74.6824625",
         "Malayemys_macrocephala, 14.3532445, 100.447429",
         "Manouria_impressa, 15.525993, 101.089441",
         "Peltocephalus_dumerilianus, -0.406805556, -63.44736111",
         "Podocnemis_expansa, -1.008516, -71.64148",
         "Podocnemis_lewyana, 7.8470115, -74.437783",
         "Podocnemis_unifilis, -1.4369195, -70.710111",
         "Podocnemis_unifilis, -11.955611, -53.509643",
         "Podocnemis_vogli, 4.903391, -70.4720135",
         "Pseudemys_texana, 31.0514145, -98.7093235",
         "Stigmochelys_pardalis, -11.916498, 30.377801",
         "Terrapene_carolina, 38.905738, -77.033739",
         "Testudo_graeca, 32.656925, 35.481161",
         "Trachemys_scripta, 30.3259625, -97.7720965",
         "Trachemys_scripta, 30.5040795, -90.4592465")

df.tpiv <- df %>% filter(!sp.coords %in% tpiv)
df.trt <- df %>% filter(!sp.coords %in% trt)

tree.tpiv <- drop.tip(tree, setdiff(tree$tip.label, df.tpiv$Species))
tree.trt <- drop.tip(tree, setdiff(tree$tip.label, df.trt$Species))
tree.tpiv <- ladderize(tree.tpiv)
tree.trt <- ladderize(tree.trt)

# Calculate species mean and mean-centered values for coverage

df.trt$smlogCoverage <- sapply(split(df.trt$logCoverage, df.trt$Species), mean)[df.trt$Species]
df.trt$sclogCoverage <- df.trt$logCoverage - df.trt$smlogCoverage

df.tpiv$Source <- as.numeric(as.factor(df.tpiv$Source))
df.trt$Source <- as.numeric(as.factor(df.trt$Source))

####################################################################################################

# Repeat TRT analysis using temperatures during warmest quarter of the year

df.trt %>%
  filter(!is.na(qTemps)) %>%
  filter(Latitude_Range < 2) %>%
  mutate(Species2 = Species) -> mcmc.data

tree.mcmc <- drop.tip(tree.trt, setdiff(tree.trt$tip.label, mcmc.data$Species))
vcv <- vcv.phylo(tree.mcmc, corr = T)

bf1 <- brmsformula(logTRT | se(logTRT_SD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))
bf2 <- brmsformula(logCoverage ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))
bf3 <- brmsformula(qTemps | se(qTemps_SD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))

prior <- c(prior(normal(0, 5), "Intercept", resp = "logCoverage"),
           prior(normal(0, 5), "Intercept", resp = "logTRT"),
           prior(normal(20, 5), "Intercept", resp = "qTemps"),
           prior(student_t(3, 0, 2.5), "sd", resp = "logCoverage"),
           prior(student_t(3, 0, 2.5), "sd", resp = "logTRT"),
           prior(student_t(3, 0, 2.5), "sd", resp = "qTemps"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "logCoverage"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "logTRT"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "qTemps"))

trt.qtemps <- brm(bf1 + bf2 + bf3 + set_rescor(TRUE),
                  data = mcmc.data,
                  data2 = list(vcv = vcv),
                  prior = prior,
                  chains = 4,
                  cores = 4,
                  iter = 5000,
                  thin = 10,
                  warmup = 1000,
                  control = list(adapt_delta = 0.99),
                  seed = 123)

mcmc_plot(trt.qtemps, type = "trace")
mcmc_plot(trt.qtemps, type = "acf")
mcmc_plot(trt.qtemps, type = "rhat")
mcmc_plot(trt.qtemps, type = "neff")

rstan::check_hmc_diagnostics(trt.qtemps$fit)

pp_check(trt.qtemps, type = "intervals", resp = "logTRT")
pp_check(trt.qtemps, ndraws = 25, resp = "logTRT")
pp_check(trt.qtemps, type = "intervals", resp = "qTemps")
pp_check(trt.qtemps, ndraws = 25, resp = "qTemps")
pp_check(trt.qtemps, type = "intervals", resp = "logCoverage")
pp_check(trt.qtemps, ndraws = 25, resp = "logCoverage")

summary(trt.qtemps)

as.data.frame(trt.qtemps) -> trt.qtemps.df

trt.qtemps.df$pc.TRT.Coverage <- (trt.qtemps.df$rescor__logTRT__logCoverage - (trt.qtemps.df$rescor__logCoverage__qTemps*trt.qtemps.df$rescor__logTRT__qTemps)) / sqrt((1-trt.qtemps.df$rescor__logCoverage__qTemps^2)*(1-trt.qtemps.df$rescor__logTRT__qTemps^2))
trt.qtemps.df$pc.TRT.qTemps <- (trt.qtemps.df$rescor__logTRT__qTemps - (trt.qtemps.df$rescor__logTRT__logCoverage*trt.qtemps.df$rescor__logCoverage__qTemps)) / sqrt((1-trt.qtemps.df$rescor__logCoverage__qTemps^2)*(1-trt.qtemps.df$rescor__logTRT__logCoverage^2))
trt.qtemps.df$pc.qTemps.Coverage <- (trt.qtemps.df$rescor__logCoverage__qTemps - (trt.qtemps.df$rescor__logTRT__logCoverage*trt.qtemps.df$rescor__logTRT__qTemps)) / sqrt((1-trt.qtemps.df$rescor__logTRT__logCoverage^2)*(1-trt.qtemps.df$rescor__logTRT__qTemps^2))
trt.qtemps.df$ppc.TRT.Coverage <- (trt.qtemps.df$cor_Species__logTRT_Intercept__logCoverage_Intercept - (trt.qtemps.df$cor_Species__logTRT_Intercept__qTemps_Intercept*trt.qtemps.df$cor_Species__logCoverage_Intercept__qTemps_Intercept)) / sqrt((1-trt.qtemps.df$cor_Species__logTRT_Intercept__qTemps_Intercept^2)*(1-trt.qtemps.df$cor_Species__logCoverage_Intercept__qTemps_Intercept^2))
trt.qtemps.df$ppc.TRT.qTemps <- (trt.qtemps.df$cor_Species__logTRT_Intercept__qTemps_Intercept - (trt.qtemps.df$cor_Species__logTRT_Intercept__logCoverage_Intercept*trt.qtemps.df$cor_Species__logCoverage_Intercept__qTemps_Intercept)) / sqrt((1-trt.qtemps.df$cor_Species__logTRT_Intercept__logCoverage_Intercept^2)*(1-trt.qtemps.df$cor_Species__logCoverage_Intercept__qTemps_Intercept^2))
trt.qtemps.df$ppc.qTemps.Coverage <- (trt.qtemps.df$cor_Species__logCoverage_Intercept__qTemps_Intercept - (trt.qtemps.df$cor_Species__logTRT_Intercept__logCoverage_Intercept*trt.qtemps.df$cor_Species__logTRT_Intercept__qTemps_Intercept)) / sqrt((1-trt.qtemps.df$cor_Species__logTRT_Intercept__logCoverage_Intercept^2)*(1-trt.qtemps.df$cor_Species__logTRT_Intercept__qTemps_Intercept^2))

apply(trt.qtemps.df[,429:434], 2, quantile, probs = c(0.025,0.5,0.975), na.rm = T)

####################################################################################################

# Repeat TRT analysis using soil temperatures

df.trt %>%
  filter(!is.na(SoilTemps)) %>%
  filter(Latitude_Range < 2) %>%
  mutate(Species2 = Species) -> mcmc.data

tree.mcmc <- drop.tip(tree.trt, setdiff(tree.trt$tip.label, mcmc.data$Species))
vcv <- vcv.phylo(tree.mcmc, corr = T)

bf1 <- brmsformula(logTRT | se(logTRT_SD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))
bf2 <- brmsformula(logCoverage ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))
bf3 <- brmsformula(SoilTemps | se(SoilTemps_SD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))

prior <- c(prior(normal(0, 5), "Intercept", resp = "logCoverage"),
           prior(normal(0, 5), "Intercept", resp = "logTRT"),
           prior(normal(20, 5), "Intercept", resp = "SoilTemps"),
           prior(student_t(3, 0, 2.5), "sd", resp = "logCoverage"),
           prior(student_t(3, 0, 2.5), "sd", resp = "logTRT"),
           prior(student_t(3, 0, 2.5), "sd", resp = "SoilTemps"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "logCoverage"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "logTRT"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "SoilTemps"))

trt.soiltemps <- brm(bf1 + bf2 + bf3 + set_rescor(TRUE),
                     data = mcmc.data,
                     data2 = list(vcv = vcv),
                     prior = prior,
                     chains = 4,
                     cores = 4,
                     iter = 5000,
                     thin = 10,
                     warmup = 1000,
                     control = list(adapt_delta = 0.99),
                     seed = 123)

mcmc_plot(trt.soiltemps, type = "trace")
mcmc_plot(trt.soiltemps, type = "acf")
mcmc_plot(trt.soiltemps, type = "rhat")
mcmc_plot(trt.soiltemps, type = "neff")

rstan::check_hmc_diagnostics(trt.soiltemps$fit)

pp_check(trt.soiltemps, type = "intervals", resp = "logTRT")
pp_check(trt.soiltemps, ndraws = 25, resp = "logTRT")
pp_check(trt.soiltemps, type = "intervals", resp = "SoilTemps")
pp_check(trt.soiltemps, ndraws = 25, resp = "SoilTemps")
pp_check(trt.soiltemps, type = "intervals", resp = "logCoverage")
pp_check(trt.soiltemps, ndraws = 25, resp = "logCoverage")

summary(trt.soiltemps)

as.data.frame(trt.soiltemps) -> trt.soiltemps.df

trt.soiltemps.df$pc.TRT.Coverage <- (trt.soiltemps.df$rescor__logTRT__logCoverage - (trt.soiltemps.df$rescor__logCoverage__SoilTemps*trt.soiltemps.df$rescor__logTRT__SoilTemps)) / sqrt((1-trt.soiltemps.df$rescor__logCoverage__SoilTemps^2)*(1-trt.soiltemps.df$rescor__logTRT__SoilTemps^2))
trt.soiltemps.df$pc.TRT.SoilTemps <- (trt.soiltemps.df$rescor__logTRT__SoilTemps - (trt.soiltemps.df$rescor__logTRT__logCoverage*trt.soiltemps.df$rescor__logCoverage__SoilTemps)) / sqrt((1-trt.soiltemps.df$rescor__logCoverage__SoilTemps^2)*(1-trt.soiltemps.df$rescor__logTRT__logCoverage^2))
trt.soiltemps.df$pc.SoilTemps.Coverage <- (trt.soiltemps.df$rescor__logCoverage__SoilTemps - (trt.soiltemps.df$rescor__logTRT__logCoverage*trt.soiltemps.df$rescor__logTRT__SoilTemps)) / sqrt((1-trt.soiltemps.df$rescor__logTRT__logCoverage^2)*(1-trt.soiltemps.df$rescor__logTRT__SoilTemps^2))
trt.soiltemps.df$ppc.TRT.Coverage <- (trt.soiltemps.df$cor_Species__logTRT_Intercept__logCoverage_Intercept - (trt.soiltemps.df$cor_Species__logTRT_Intercept__SoilTemps_Intercept*trt.soiltemps.df$cor_Species__logCoverage_Intercept__SoilTemps_Intercept)) / sqrt((1-trt.soiltemps.df$cor_Species__logTRT_Intercept__SoilTemps_Intercept^2)*(1-trt.soiltemps.df$cor_Species__logCoverage_Intercept__SoilTemps_Intercept^2))
trt.soiltemps.df$ppc.TRT.SoilTemps <- (trt.soiltemps.df$cor_Species__logTRT_Intercept__SoilTemps_Intercept - (trt.soiltemps.df$cor_Species__logTRT_Intercept__logCoverage_Intercept*trt.soiltemps.df$cor_Species__logCoverage_Intercept__SoilTemps_Intercept)) / sqrt((1-trt.soiltemps.df$cor_Species__logTRT_Intercept__logCoverage_Intercept^2)*(1-trt.soiltemps.df$cor_Species__logCoverage_Intercept__SoilTemps_Intercept^2))
trt.soiltemps.df$ppc.SoilTemps.Coverage <- (trt.soiltemps.df$cor_Species__logCoverage_Intercept__SoilTemps_Intercept - (trt.soiltemps.df$cor_Species__logTRT_Intercept__logCoverage_Intercept*trt.soiltemps.df$cor_Species__logTRT_Intercept__SoilTemps_Intercept)) / sqrt((1-trt.soiltemps.df$cor_Species__logTRT_Intercept__logCoverage_Intercept^2)*(1-trt.soiltemps.df$cor_Species__logTRT_Intercept__SoilTemps_Intercept^2))

apply(trt.soiltemps.df[,399:404], 2, quantile, probs = c(0.025,0.5,0.975), na.rm = T)

####################################################################################################

# Repeat Tpiv analysis using temperatures during warmest quarter of the year

df.tpiv %>%
  filter(!is.na(qTemps)) %>%
  filter(Latitude_Range < 2) %>%
  mutate(Species2 = Species) -> mcmc.data

tree.mcmc <- drop.tip(tree.tpiv, setdiff(tree.tpiv$tip.label, mcmc.data$Species))
vcv <- vcv.phylo(tree.mcmc, corr = T)

bf1 <- brmsformula(Tpiv | se(Tpiv_SD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))
bf2 <- brmsformula(qTemps | se(qTemps_SD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))

prior <- c(prior(normal(30, 5), "Intercept", resp = "Tpiv"),
           prior(normal(25, 5), "Intercept", resp = "qTemps"),
           prior(student_t(3, 0, 2.5), "sd", resp = "Tpiv"),
           prior(student_t(3, 0, 2.5), "sd", resp = "qTemps"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "Tpiv"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "qTemps"))

tpiv.qtemps <- brm(bf1 + bf2 + set_rescor(TRUE),
                   data = mcmc.data,
                   data2 = list(vcv = vcv),
                   prior = prior,
                   chains = 4,
                   cores = 4,
                   iter = 5000,
                   thin = 10,
                   warmup = 1000,
                   control = list(adapt_delta = 0.99),
                   seed = 123)

mcmc_plot(tpiv.qtemps, type = "trace")
mcmc_plot(tpiv.qtemps, type = "acf")
mcmc_plot(tpiv.qtemps, type = "rhat")
mcmc_plot(tpiv.qtemps, type = "neff")

rstan::check_hmc_diagnostics(tpiv.qtemps$fit)

pp_check(tpiv.qtemps, type = "intervals", resp = "Tpiv")
pp_check(tpiv.qtemps, ndraws = 25, resp = "Tpiv")
pp_check(tpiv.qtemps, type = "intervals", resp = "qTemps")
pp_check(tpiv.qtemps, ndraws = 25, resp = "qTemps")

summary(tpiv.qtemps)

####################################################################################################

# Repeat Tpiv analysis using soil temperatures

df.tpiv %>%
  filter(!is.na(SoilTemps)) %>%
  filter(Latitude_Range < 2) %>%
  mutate(Species2 = Species) -> mcmc.data

tree.mcmc <- drop.tip(tree.tpiv, setdiff(tree.tpiv$tip.label, mcmc.data$Species))
vcv <- vcv.phylo(tree.mcmc, corr = T)

bf1 <- brmsformula(Tpiv | se(Tpiv_SD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))
bf2 <- brmsformula(SoilTemps | se(SoilTemps_SD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))

prior <- c(prior(normal(30, 5), "Intercept", resp = "Tpiv"),
           prior(normal(25, 5), "Intercept", resp = "SoilTemps"),
           prior(student_t(3, 0, 2.5), "sd", resp = "Tpiv"),
           prior(student_t(3, 0, 2.5), "sd", resp = "SoilTemps"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "Tpiv"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "SoilTemps"))

tpiv.soiltemps <- brm(bf1 + bf2 + set_rescor(TRUE),
                      data = mcmc.data,
                      data2 = list(vcv = vcv),
                      prior = prior,
                      chains = 4,
                      cores = 4,
                      iter = 5000,
                      thin = 10,
                      warmup = 1000,
                      control = list(adapt_delta = 0.99),
                      seed = 123)

mcmc_plot(tpiv.soiltemps, type = "trace")
mcmc_plot(tpiv.soiltemps, type = "acf")
mcmc_plot(tpiv.soiltemps, type = "rhat")
mcmc_plot(tpiv.soiltemps, type = "neff")

rstan::check_hmc_diagnostics(tpiv.soiltemps$fit)

pp_check(tpiv.soiltemps, type = "intervals", resp = "Tpiv")
pp_check(tpiv.soiltemps, ndraws = 25, resp = "Tpiv")
pp_check(tpiv.soiltemps, type = "intervals", resp = "SoilTemps")
pp_check(tpiv.soiltemps, ndraws = 25, resp = "SoilTemps")

summary(tpiv.soiltemps)

