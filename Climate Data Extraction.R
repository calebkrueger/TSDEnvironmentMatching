# R code for working with climate rasters to obtain temperature variables

# Load packages

library(tidyverse) # v2.0.0
library(phytools) # v2.3-0
library(sf) #1.0-16
library(terra) #1.7-78

terraOptions(memmax = 8)

####################################################################################################

# Read in output of reaction norm estimation R script

df <- read.csv('rn_out.csv')

# Read in Thomson et al. (2021) phylogeny

tree <- read.tree("Thomson.nwk")

# Add missing species based on divergence estimates from timetree.org

tree <- bind.tip(tree,
                 tip.label = "Podocnemis_lewyana",
                 where = which(tree$tip.label == "Podocnemis_unifilis"),
                 position = 15.62)

tree <- ladderize(tree)

####################################################################################################

## Add temperature data ##

# Files downloaded beforehand
# See 'CHELSA_Monthly_*.txt' files for download urls

# Begin with mean temperature during the estimated thermosensitive period (TSP) months

# Set wd to folder with files from "CHELSA_Monthly_Temps_1981-2010.txt"

rastlist <- list.files(pattern = '.tif',
                       all.files = T,
                       full.names = F)

tavg <- lapply(rastlist, rast)

# Read in the mask

lsmask <- rast("landseamask.sdat")

# Next, read in data on nesting seasons and incubation durations by population
# Data are on github and in supplementary materials from publisher

nests <- read.csv("nesting.csv", header = T)
nests <- nests[!nests$tsp_months=="",]

# Extract mean temperature during incubation months at each population's coordinates

# Define bounding box around each population

df$Latitude_Max <- df$Latitude + (df$Latitude_Range/2)
df$Latitude_Min <- df$Latitude - (df$Latitude_Range/2)
df$Longitude_Max <- df$Longitude + (df$Longitude_Range/2)
df$Longitude_Min <- df$Longitude - (df$Longitude_Range/2)

df %>%
  dplyr::select(sp.coords, Latitude_Max, Latitude_Min, Longitude_Max, Longitude_Min) %>%
  dplyr::filter(!is.na(Latitude_Max)) -> lat.long

coords <- st_as_sf(lat.long,
                   coords = c("Longitude_Min", "Latitude_Max"),
                   crs = 4326)[,c(1,4)]
coords <- rbind(coords, st_as_sf(lat.long,
                                 coords = c("Longitude_Max", "Latitude_Max"),
                                 crs = 4326)[,c(1,4)])
coords <- rbind(coords, st_as_sf(lat.long,
                                 coords = c("Longitude_Max", "Latitude_Min"),
                                 crs = 4326)[,c(1,4)])
coords <- rbind(coords, st_as_sf(lat.long,
                                 coords = c("Longitude_Min", "Latitude_Min"),
                                 crs = 4326)[,c(1,4)])
coords <- rbind(coords, st_as_sf(lat.long,
                                 coords = c("Longitude_Min", "Latitude_Max"),
                                 crs = 4326)[,c(1,4)])

coords %>%
  group_by(sp.coords) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") -> polys

polys <- st_transform(polys, crs(tavg[[1]]))

polys$Tavg <- NA
polys$Tavg_MAD <- NA

# Select relevant monthly rasters by population and calculate median and MAD of temperatures

temp <- list()
temp.mask <- NA
iter <- 1
out <- NA

for(i in nests$sp.coords){
  print(i)
  months <- unlist(strsplit(nests[nests$sp.coords==i,"tsp_months"], ","))
  for(j in months){
    temp.mask <- crop(lsmask, polys[polys$sp.coords==i, "geometry"])
    temp[[iter]] <- mask(crop(tavg[[as.numeric(j)]], polys[polys$sp.coords==i, "geometry"]), temp.mask)
    iter <- iter + 1
  }
  out <- app(sds(temp), mean)
  polys[polys$sp.coords==i, "Tavg"] <- global(out, median, na.rm = T)[[1]]
  polys[polys$sp.coords==i, "Tavg_MAD"] <- global(out, mad, na.rm = T)[[1]]
  temp <- list()
  iter <- 1
  # For populations on minor islands not represented in masked rasters, use ocean surface temps:
  if(is.na(polys[polys$sp.coords==i,]$Tavg)==T){
    for(j in months){
      temp[[iter]] <- crop(tavg[[as.numeric(j)]], polys[polys$sp.coords==i,"geometry"])
      iter <- iter + 1
    }
    out <- app(sds(temp), mean)
    polys[polys$sp.coords==i, "Tavg"] <- global(out, median, na.rm = T)[[1]]
    polys[polys$sp.coords==i, "Tavg_MAD"] <- global(out, mad, na.rm = T)[[1]]
    temp <- list()
    iter <- 1
  } else {}
}

polys %>%
  base::as.data.frame() %>%
  dplyr::select(sp.coords, Tavg, Tavg_MAD) %>%
  left_join(x = df, by = "sp.coords") -> df

# Next, intra-annual variation in temperature

# Set wd to folder with files from "CHELSA_Monthly_Max_Min_1981-2010.txt"

rastlist <- list.files(pattern = 'CHELSA_tasmax',
                       all.files = T,
                       full.names = F)

tmax <- lapply(rastlist, rast)

rastlist <- list.files(pattern = 'CHELSA_tasmin',
                       all.files = T,
                       full.names = F)

tmin <- lapply(rastlist, rast)

polys$wyt <- NA
polys$wyt_MAD <- NA

# Select relevant monthly rasters by population and warmest average daily high and coldest mean daily low
# Intra-annual variation will be represented as the max minus the min

maxrasts <- list()
minrasts <- list()
temp.mask <- NA
iter <- 1
hi <- NA
lo <- NA
out <- NA

for(i in nests$sp.coords){
  print(i)
  months <- unlist(strsplit(nests[nests$sp.coords==i,"tsp_months"], ","))
  dat <- data.frame(months, maxtemp = rep(NA, length(months)), mintemp = rep(NA, length(months)))
  for(j in months){
    temp.mask <- crop(lsmask, polys[polys$sp.coords==i, "geometry"])
    maxrasts[[iter]] <- mask(crop(tmax[[as.numeric(j)]], polys[polys$sp.coords==i,"geometry"]), temp.mask)
    minrasts[[iter]] <- mask(crop(tmin[[as.numeric(j)]], polys[polys$sp.coords==i,"geometry"]), temp.mask)
    dat[dat$months==j,"maxtemp"] <- global(maxrasts[[iter]], mean, na.rm = T)[[1]]
    dat[dat$months==j,"mintemp"] <- global(minrasts[[iter]], mean, na.rm = T)[[1]]
    iter <- iter + 1
  }
  if(is.na(dat$maxtemp)[1]==T){
    for(j in months){
      maxrasts[[iter]] <- crop(tmax[[as.numeric(j)]], polys[polys$sp.coords==i,"geometry"])
      minrasts[[iter]] <- crop(tmin[[as.numeric(j)]], polys[polys$sp.coords==i,"geometry"])
      dat[dat$months==j,"maxtemp"] <- global(maxrasts[[iter]], mean, na.rm = T)[[1]]
      dat[dat$months==j,"mintemp"] <- global(minrasts[[iter]], mean, na.rm = T)[[1]]
      iter <- iter + 1
    }
    hi <- crop(tmax[[as.numeric(dat[which.max(dat$maxtemp),]$months)]], polys[polys$sp.coords==i,"geometry"])
    lo <- crop(tmin[[as.numeric(dat[which.min(dat$mintemp),]$months)]], polys[polys$sp.coords==i,"geometry"])
    out <- hi - lo
    polys[polys$sp.coords==i, "wyt"] <- global(out, median, na.rm = T)[[1]]
    polys[polys$sp.coords==i, "wyt_MAD"] <- global(out, mad, na.rm = T)[[1]]
    maxrasts <- list()
    minrasts <- list()
    iter <- 1
  } else {
    hi <- mask(crop(tmax[[as.numeric(dat[which.max(dat$maxtemp),]$months)]], polys[polys$sp.coords==i,"geometry"]), temp.mask)
    lo <- mask(crop(tmin[[as.numeric(dat[which.min(dat$mintemp),]$months)]], polys[polys$sp.coords==i,"geometry"]), temp.mask)
    out <- hi - lo
    polys[polys$sp.coords==i, "wyt"] <- global(out, median, na.rm = T)[[1]]
    polys[polys$sp.coords==i, "wyt_MAD"] <- global(out, mad, na.rm = T)[[1]]
    maxrasts <- list()
    minrasts <- list()
    iter <- 1
  }
}

polys %>%
  base::as.data.frame() %>%
  dplyr::select(sp.coords, wyt, wyt_MAD) %>%
  left_join(x = df, by = "sp.coords") -> df

# Last, interannual variation in temperature

# Calculated as standard deviation of mean temperature during TSP months

# Set wd to folder with files from "CHELSA_Monthly_Temps_by_Year_1981-2010.txt"

rastlist <- list.files(pattern = 'CHELSA_tas_',
                       all.files = T,
                       full.names = F)

temps <- lapply(rastlist, rast)

for(i in 1:length(temps)){
  names(temps)[[i]] <- names(temps[[i]])
}

# Calculate SD for TSP months per population

polys$byt <- NA
polys$byt_MAD <- NA

temp <- list()
yr <- list()
yr.temp <- list()
temp.mask <- NA
iter <- 1
out <- NA

for(i in nests$sp.coords){
  months <- str_pad(unlist(strsplit(nests[nests$sp.coords==i,"tsp_months"], ",")), 2, pad = "0")
  rast.names <- paste("CHELSA_tas_", months, sep = "")
  temp <- Filter(function(x) isTRUE(substr(names(x), 1, 13) %in% rast.names), temps)
  temp.mask <- crop(lsmask, polys[polys$sp.coords==i, "geometry"])
  temp <- lapply(temp, crop, polys[polys$sp.coords==i, "geometry"])
  temp <- lapply(temp, mask, temp.mask)
  temp <- lapply(temp, function(x) {x * 0.1 - 273.15})
  for(j in 1981:2010){
    yr[[iter]] <- app(sds(temp[grep(j, names(temp))]), mean)
    iter <- iter + 1
  }
  out <- app(sds(yr), sd)
  polys[polys$sp.coords==i, "byt"] <- global(out, median, na.rm = T)[[1]]
  polys[polys$sp.coords==i, "byt_MAD"] <- global(out, mad, na.rm = T)[[1]]
  iter <- 1
  if(is.na(polys[polys$sp.coords==i,]$byt)==T){
    temp <- Filter(function(x) isTRUE(substr(names(x), 1, 13) %in% rast.names), temps)
    temp <- lapply(temp, crop, polys[polys$sp.coords==i, "geometry"])
    temp <- lapply(temp, function(x) {x * 0.1 - 273.15})
    for(j in 1981:2010){
      yr[[iter]] <- app(sds(temp[grep(j, names(temp))]), mean)
      iter <- iter + 1
    }
    out <- app(sds(yr), sd)
    polys[polys$sp.coords==i, "byt"] <- global(out, median, na.rm = T)[[1]]
    polys[polys$sp.coords==i, "byt_MAD"] <- global(out, mad, na.rm = T)[[1]]
    iter <- 1
  } else {}
}

polys %>%
  base::as.data.frame() %>%
  dplyr::select(sp.coords, byt, byt_MAD) %>%
  left_join(x = df, by = "sp.coords") -> df

# Write to file

write.csv(df, "rasters_out.csv", row.names = F)
