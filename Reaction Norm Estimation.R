# R code for fitting sex-ratio reaction norms to turtle "populations"

# Load packages

library(tidyverse) # v2.0.0
library(phytools) # v2.3-0
library(embryogrowth) # v8.3-4

####################################################################################################

## Fit sex-ratio reaction norms ##

# Begin by reading in the ROSIE database

data <- read.csv("https://raw.githubusercontent.com/calebkrueger/ROSIE/main/Database.csv")

# Filter to include data from constant incubation experiments without chemical manipulations
# Exclude species with GSD, multiple instances of the same data, incomplete data, and studies using unreliable sexing techniques

data <- data %>% filter(SDM == "TSD",
                        Incubation_Setup == "Constant",
                        !Sexing_Method %in% c("Morphometrics", "Morphometrics, macroscopy"),
                        Data_Elsewhere %in% c(0,NA),
                        !Source == "Janzen, F. J., Unpublished data.", # Included in Krueger & Janzen 2023 but missing Data_Elsewhere == 1 tag
                        is.na(Chemical_Treatment),
                        !is.na(Mean_Temp),
                        !is.na(Males),
                        !is.na(Females))

# Retain relevant variables

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
             Tpiv_MAD = NA,
             TRT = NA,
             TRT_MAD = NA,
             logTRT_MAD = NA,
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
  log.p$Prior2[3] <- 10
  log.p$Prior1[4] <- 0
  log.p$Prior2[4] <- 10
  log.p$Min[1] <- 0
  log.p$Max[1] <- 100
  log.p$Min[2] <- 0
  log.p$Max[2] <- 100
  log.p$Min[3] <- 0
  log.p$Max[3] <- 100
  log.p$Min[4] <- -100
  log.p$Max[4] <- 0
  log.result <- tsd_MHmcmc(result = log, parametersMCMC = log.p, n.iter = 1e4, adaptive = T)
  CI <- P_TRT(x = log, resultmcmc = log.result, replicate.CI = 1e4)
  df[df$sp.coords == i, "Tpiv"] <- as.numeric(median(CI$P_TRT[,8]))
  df[df$sp.coords == i, "Tpiv_MAD"] <- as.numeric(mad(CI$P_TRT[,8]))
  df[df$sp.coords == i, "TRT"] <- as.numeric(median(CI$P_TRT[,7]))
  df[df$sp.coords == i, "TRT_MAD"] <- as.numeric(mad(CI$P_TRT[,7]))
  df[df$sp.coords == i, "logTRT_MAD"] <- as.numeric(mad(log(CI$P_TRT[,7])))
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
  plot(log, resultmcmc = log.result, main = i, xlim = c(18,35), show.PTRT = F, show.observations = T)
  abline(h = 0.50, lty = "dotted")
  abline(h = 0.95, lty = "dotted")
  plot(log.result, parameters = "P_high", xlim = c(0, 50), legend = F, breaks = seq(0, 100, l = 100))
  plot(log.result, parameters = "S_high", xlim = c(-25, 0), legend = F, breaks = seq(-100, 0, l = 200))
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
  log.p$Prior2[1] <- 10
  log.p$Prior1[2] <- 0
  log.p$Prior2[2] <- 10
  log.p$Min[1] <- 0
  log.p$Max[1] <- 100
  log.p$Min[2] <- -100
  log.p$Max[2] <- 0
  log.result <- tsd_MHmcmc(result = log, parametersMCMC = log.p, n.iter = 1e4, adaptive = T)
  CI <- P_TRT(x = log, resultmcmc = log.result, replicate.CI = 1e4)
  df[df$sp.coords == i, "Tpiv"] <- as.numeric(median(CI$P_TRT[,4]))
  df[df$sp.coords == i, "Tpiv_MAD"] <- as.numeric(mad(CI$P_TRT[,4]))
  df[df$sp.coords == i, "TRT"] <- as.numeric(median(CI$P_TRT[,3]))
  df[df$sp.coords == i, "TRT_MAD"] <- as.numeric(mad(CI$P_TRT[,3]))
  df[df$sp.coords == i, "logTRT_MAD"] <- as.numeric(mad(log(CI$P_TRT[,3])))
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
  plot(log, resultmcmc = log.result, main = i, xlim = c(18,35), show.PTRT = F, show.observations = T)
  plot(log.result, parameters = "P", xlim = c(0, 50), legend = F, breaks = seq(0, 100, l = 100))
  plot(log.result, parameters = "S", xlim = c(-25, 0), legend = F, breaks = seq(-100, 0, l = 200))
}

# Log-transform select variables for analysis

df$logTRT <- log(df$TRT)
df$logCoverage <- log(df$TRT_Coverage)

# Write to file

write.csv(df, 'rn_out.csv', row.names = F)
