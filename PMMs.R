# R code for fitting phylogenetic mixed models to reaction norm and climate data

# Load packages

library(tidyverse) # v2.0.0
library(phytools) # v2.3-0
library(brms) #2.21.0

####################################################################################################

# Read in output of Climate Data Extraction R script

df <- read.csv('rasters_out.csv')

# Read in Thomson et al. (2021) phylogeny

tree <- read.tree("Thomson.nwk")

# Add missing species based on divergence estimates from timetree.org

tree <- bind.tip(tree,
                 tip.label = "Podocnemis_lewyana",
                 where = which(tree$tip.label == "Podocnemis_unifilis"),
                 position = 15.62)

tree <- ladderize(tree)

####################################################################################################

# Create two data frames - one of populations with well-defined Tpiv and another with well-defined TRT

# List of species to drop from Tpiv analyses

tpiv <- c("Kinosternon_baurii, 29.9251295, -85.0121755",
          "Kinosternon_creaseri, 19.6773545, -88.741656",
          "Kinosternon_leucostomum, 8.813487, -84.9906775",
          "Terrapene_carolina, 27.7713485, -82.172115",
          "Astrochelys_radiata, -24.4213925, 45.2127545",
          "Caretta_caretta, -23.7442075, 113.5620435",
          "Chelonia_mydas, 20.5694525, -157.6830405",
          "Chelydra_serpentina, 46.202449, -94.4147375",
          "Chelydra_serpentina, 45.575434, -78.683943",
          "Chrysemys_picta, 46.0272355, -86.621559",
          "Chrysemys_picta, 37.223365, -80.484596",
          "Geochelone_elegans, 16.7827945, 76.7525175",
          "Geochelone_platynota, 22.00165, 95.308987",
          "Kinosternon_flavescens, 40.635338, -99.870176", # Potentially pattern II with one temp in lower TRT and one in upper TRT
          "Malaclemys_terrapin, 39.162667, -74.6824625",
          "Manouria_impressa, 15.525993, 101.089441",
          "Peltocephalus_dumerilianus, -0.406805556, -63.44736111",
          "Podocnemis_expansa, -1.008516, -71.64148",
          "Podocnemis_unifilis, -1.4369195, -70.710111",
          "Podocnemis_vogli, 4.903391, -70.4720135",
          "Stigmochelys_pardalis, -11.916498, 30.377801",
          "Terrapene_carolina, 38.905738, -77.033739",
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
         "Kinosternon_flavescens, 40.635338, -99.870176", # Potentially pattern II with one temp in lower TRT and one in upper TRT
         "Kinosternon_leucostomum, 17.4911865, -93.04148575", # Pattern II species without a defined TRT
         "Macrochelys_temminckii, 34.7522945, -90.963746", # Pattern II species without a defined TRT
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

# First, assess relationship between TRT coverage and TRT estimates
# Also calculate phylogenetic heritability of TRT

df.trt %>%
  mutate(Species2 = Species) -> mcmc.data

vcv <- vcv.phylo(tree.trt, corr = T)

bf1 <- brmsformula(logTRT | se(logTRT_MAD, sigma = TRUE) ~ smlogCoverage + sclogCoverage + (1|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))

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

# Calculate proportion of total variance explained by fixed effects

x <- VarCorr(trt, summary = F)
y.pred <- posterior_epred(trt, re.form = NA)

r2 <- rep(NA, nrow(y.pred))

for(i in 1:length(r2)){
  r2[i] <- var(y.pred[i,]) / (var(y.pred[i,]) + x$residual__$sd[i]^2 + x$Source$sd[i]^2 + x$Species$sd[i]^2 + x$Species2$sd[i]^2)
}

quantile(r2, probs = c(0.025, 0.5, 0.975))

####################################################################################################

# Calculate phylogenetic heritability of Tpiv

df.tpiv %>%
  mutate(Species2 = Species) -> mcmc.data

vcv <- vcv.phylo(tree.tpiv, corr = T)

bf1 <- brmsformula(Tpiv | se(Tpiv_MAD, sigma = TRUE) ~ (1|gr(Species, cov = vcv)) + (1|Species2) + (1|Source))

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

bf1 <- brmsformula(logTRT | se(logTRT_MAD, sigma = TRUE) ~ smlogCoverage + sclogCoverage + (1|b|gr(Species, cov = vcv)) + (1|c|Species2) + (1|Source))
bf2 <- brmsformula(Tpiv | se(Tpiv_MAD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|c|Species2) + (1|Source))

prior <- c(prior(normal(0, 5), "Intercept", resp = "logTRT"),
           prior(normal(30, 5), "Intercept", resp = "Tpiv"),
           prior(student_t(3, 0, 2.5), "sd", resp = "logTRT"),
           prior(student_t(3, 0, 2.5), "sd", resp = "Tpiv"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "logTRT"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "Tpiv"))

tpiv.trt <- brm(bf1 + bf2 + set_rescor(TRUE),
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

summary(tpiv.trt)

# Calculate correlations

cor.mat <- VarCorr(tpiv.trt, summary = F)
mat <- NA

rhos <- data.frame(phy = rep(NA, 1600),
                   nonphy = rep(NA, 1600),
                   sp = rep(NA, 1600),
                   pop = rep(NA, 1600),
                   ind = rep(NA, 1600),
                   tot = rep(NA, 1600))

for(i in 1:1600){
  mat <- solve(cor.mat$Species$cov[i,,])
  rhos[i,]$phy <- (-mat[1,2] / sqrt(mat[1,1] * mat[2,2]))
  mat <- solve(cor.mat$Species2$cov[i,,])
  rhos[i,]$nonphy <- (-mat[1,2] / sqrt(mat[1,1] * mat[2,2]))
  mat <- solve(cor.mat$residual__$cov[i,,])
  rhos[i,]$pop <- (-mat[1,2] / sqrt(mat[1,1] * mat[2,2]))
  mat <- solve(cor.mat$Species$cov[i,,] + cor.mat$Species2$cov[i,,])
  rhos[i,]$sp <- (-mat[1,2] / sqrt(mat[1,1] * mat[2,2]))
  mat <- solve(cor.mat$Species2$cov[i,,] + cor.mat$residual__$cov[i,,])
  rhos[i,]$ind <- (-mat[1,2] / sqrt(mat[1,1] * mat[2,2]))
  mat <- solve(cor.mat$Species$cov[i,,] + cor.mat$Species2$cov[i,,] + cor.mat$residual__$cov[i,,])
  rhos[i,]$tot <- (-mat[1,2] / sqrt(mat[1,1] * mat[2,2]))
}

rhos %>%
  apply(2, quantile, probs = c(0.025, 0.5, 0.975))

# Calculate proportion of posterior that is negative

apply(rhos, 2, function(x){sum(x<0)/length(x)})

####################################################################################################

# Fit model to check for correlations between Tpiv, TRT, and temp variables

df.trt %>%
  filter(!is.na(Tavg)) %>%
  mutate(Species2 = Species) -> mcmc.data

tree.mcmc <- drop.tip(tree.trt, setdiff(tree.trt$tip.label, mcmc.data$Species))
vcv <- vcv.phylo(tree.mcmc, corr = T)

bf1 <- brmsformula(Tpiv | se(Tpiv_MAD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|c|Species2) + (1|Source))
bf2 <- brmsformula(logTRT | se(logTRT_MAD, sigma = TRUE) ~ smlogCoverage + sclogCoverage + (1|b|gr(Species, cov = vcv)) + (1|c|Species2) + (1|Source))
bf3 <- brmsformula(Tavg | se(Tavg_MAD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|c|Species2) + (1|Source))
bf4 <- brmsformula(wyt | se(wyt_MAD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|c|Species2) + (1|Source))
bf5 <- brmsformula(byt | se(byt_MAD, sigma = TRUE) ~ (1|b|gr(Species, cov = vcv)) + (1|c|Species2) + (1|Source))

prior <- c(prior(normal(30, 5), "Intercept", resp = "Tpiv"),
           prior(normal(0, 5), "Intercept", resp = "logTRT"),
           prior(normal(25, 5), "Intercept", resp = "Tavg"),
           prior(normal(10, 5), "Intercept", resp = "wyt"),
           prior(normal(0, 5), "Intercept", resp = "byt"),
           prior(student_t(3, 0, 2.5), "sd", resp = "Tpiv"),
           prior(student_t(3, 0, 2.5), "sd", resp = "logTRT"),
           prior(student_t(3, 0, 2.5), "sd", resp = "Tavg"),
           prior(student_t(3, 0, 2.5), "sd", resp = "wyt"),
           prior(student_t(3, 0, 2.5), "sd", resp = "byt"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "Tpiv"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "logTRT"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "Tavg"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "wyt"),
           prior(student_t(3, 0, 2.5), "sigma", resp = "byt"))

tpiv.trt.temps <- brm(bf1 + bf2 + bf3 + bf4 + bf5 + set_rescor(TRUE),
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

mcmc_plot(tpiv.trt.temps, type = "trace")
mcmc_plot(tpiv.trt.temps, type = "acf")
mcmc_plot(tpiv.trt.temps, type = "rhat")
mcmc_plot(tpiv.trt.temps, type = "neff")

rstan::check_hmc_diagnostics(tpiv.trt.temps$fit)

pp_check(tpiv.trt.temps, type = "intervals", resp = "Tpiv")
pp_check(tpiv.trt.temps, ndraws = 25, resp = "Tpiv")
pp_check(tpiv.trt.temps, type = "intervals", resp = "logTRT")
pp_check(tpiv.trt.temps, ndraws = 25, resp = "logTRT")
pp_check(tpiv.trt.temps, type = "intervals", resp = "Tavg")
pp_check(tpiv.trt.temps, ndraws = 25, resp = "Tavg")
pp_check(tpiv.trt.temps, type = "intervals", resp = "wyt")
pp_check(tpiv.trt.temps, ndraws = 25, resp = "wyt")
pp_check(tpiv.trt.temps, type = "intervals", resp = "byt")
pp_check(tpiv.trt.temps, ndraws = 25, resp = "byt")

summary(tpiv.trt.temps)

# Partial correlation between Tpiv, TRT, and temps:

cor.mat <- VarCorr(tpiv.trt.temps, summary = F)
mat <- NA

rhos.tpiv.tavg <- data.frame(phy = rep(NA, 1600),
                             nonphy = rep(NA, 1600),
                             pop = rep(NA, 1600),
                             sp = rep(NA, 1600),
                             ind = rep(NA, 1600),
                             tot = rep(NA, 1600))
rhos.tpiv.wyt <- data.frame(phy = rep(NA, 1600),
                            nonphy = rep(NA, 1600),
                            pop = rep(NA, 1600),
                            sp = rep(NA, 1600),
                            ind = rep(NA, 1600),
                            tot = rep(NA, 1600))
rhos.tpiv.byt <- data.frame(phy = rep(NA, 1600),
                            nonphy = rep(NA, 1600),
                            pop = rep(NA, 1600),
                            sp = rep(NA, 1600),
                            ind = rep(NA, 1600),
                            tot = rep(NA, 1600))
rhos.trt.tavg <- data.frame(phy = rep(NA, 1600),
                            nonphy = rep(NA, 1600),
                            pop = rep(NA, 1600),
                            sp = rep(NA, 1600),
                            ind = rep(NA, 1600),
                            tot = rep(NA, 1600))
rhos.trt.wyt <- data.frame(phy = rep(NA, 1600),
                           nonphy = rep(NA, 1600),
                           pop = rep(NA, 1600),
                           sp = rep(NA, 1600),
                           ind = rep(NA, 1600),
                           tot = rep(NA, 1600))
rhos.trt.byt <- data.frame(phy = rep(NA, 1600),
                           nonphy = rep(NA, 1600),
                           pop = rep(NA, 1600),
                           sp = rep(NA, 1600),
                           ind = rep(NA, 1600),
                           tot = rep(NA, 1600))

for(i in 1:1600){
  mat <- solve(cor.mat$Species$cov[i,,])
  rhos.tpiv.tavg[i,]$phy <- (-mat[1,3] / sqrt(mat[1,1] * mat[3,3]))
  rhos.tpiv.wyt[i,]$phy <- (-mat[1,4] / sqrt(mat[1,1] * mat[4,4]))
  rhos.tpiv.byt[i,]$phy <- (-mat[1,5] / sqrt(mat[1,1] * mat[5,5]))
  rhos.trt.tavg[i,]$phy <- (-mat[2,3] / sqrt(mat[2,2] * mat[3,3]))
  rhos.trt.wyt[i,]$phy <- (-mat[2,4] / sqrt(mat[2,2] * mat[4,4]))
  rhos.trt.byt[i,]$phy <- (-mat[2,5] / sqrt(mat[2,2] * mat[5,5]))
  mat <- solve(cor.mat$Species2$cov[i,,])
  rhos.tpiv.tavg[i,]$nonphy <- (-mat[1,3] / sqrt(mat[1,1] * mat[3,3]))
  rhos.tpiv.wyt[i,]$nonphy <- (-mat[1,4] / sqrt(mat[1,1] * mat[4,4]))
  rhos.tpiv.byt[i,]$nonphy <- (-mat[1,5] / sqrt(mat[1,1] * mat[5,5]))
  rhos.trt.tavg[i,]$nonphy <- (-mat[2,3] / sqrt(mat[2,2] * mat[3,3]))
  rhos.trt.wyt[i,]$nonphy <- (-mat[2,4] / sqrt(mat[2,2] * mat[4,4]))
  rhos.trt.byt[i,]$nonphy <- (-mat[2,5] / sqrt(mat[2,2] * mat[5,5]))
  mat <- solve(cor.mat$residual__$cov[i,,])
  rhos.tpiv.tavg[i,]$pop <- (-mat[1,3] / sqrt(mat[1,1] * mat[3,3]))
  rhos.tpiv.wyt[i,]$pop <- (-mat[1,4] / sqrt(mat[1,1] * mat[4,4]))
  rhos.tpiv.byt[i,]$pop <- (-mat[1,5] / sqrt(mat[1,1] * mat[5,5]))
  rhos.trt.tavg[i,]$pop <- (-mat[2,3] / sqrt(mat[2,2] * mat[3,3]))
  rhos.trt.wyt[i,]$pop <- (-mat[2,4] / sqrt(mat[2,2] * mat[4,4]))
  rhos.trt.byt[i,]$pop <- (-mat[2,5] / sqrt(mat[2,2] * mat[5,5]))
  mat <- solve(cor.mat$Species$cov[i,,] + cor.mat$Species2$cov[i,,])
  rhos.tpiv.tavg[i,]$sp <- (-mat[1,3] / sqrt(mat[1,1] * mat[3,3]))
  rhos.tpiv.wyt[i,]$sp <- (-mat[1,4] / sqrt(mat[1,1] * mat[4,4]))
  rhos.tpiv.byt[i,]$sp <- (-mat[1,5] / sqrt(mat[1,1] * mat[5,5]))
  rhos.trt.tavg[i,]$sp <- (-mat[2,3] / sqrt(mat[2,2] * mat[3,3]))
  rhos.trt.wyt[i,]$sp <- (-mat[2,4] / sqrt(mat[2,2] * mat[4,4]))
  rhos.trt.byt[i,]$sp <- (-mat[2,5] / sqrt(mat[2,2] * mat[5,5]))
  mat <- solve(cor.mat$Species2$cov[i,,] + cor.mat$residual__$cov[i,,])
  rhos.tpiv.tavg[i,]$ind <- (-mat[1,3] / sqrt(mat[1,1] * mat[3,3]))
  rhos.tpiv.wyt[i,]$ind <- (-mat[1,4] / sqrt(mat[1,1] * mat[4,4]))
  rhos.tpiv.byt[i,]$ind <- (-mat[1,5] / sqrt(mat[1,1] * mat[5,5]))
  rhos.trt.tavg[i,]$ind <- (-mat[2,3] / sqrt(mat[2,2] * mat[3,3]))
  rhos.trt.wyt[i,]$ind <- (-mat[2,4] / sqrt(mat[2,2] * mat[4,4]))
  rhos.trt.byt[i,]$ind <- (-mat[2,5] / sqrt(mat[2,2] * mat[5,5]))
  mat <- solve(cor.mat$Species$cov[i,,] + cor.mat$Species2$cov[i,,] + cor.mat$residual__$cov[i,,])
  rhos.tpiv.tavg[i,]$tot <- (-mat[1,3] / sqrt(mat[1,1] * mat[3,3]))
  rhos.tpiv.wyt[i,]$tot <- (-mat[1,4] / sqrt(mat[1,1] * mat[4,4]))
  rhos.tpiv.byt[i,]$tot <- (-mat[1,5] / sqrt(mat[1,1] * mat[5,5]))
  rhos.trt.tavg[i,]$tot <- (-mat[2,3] / sqrt(mat[2,2] * mat[3,3]))
  rhos.trt.wyt[i,]$tot <- (-mat[2,4] / sqrt(mat[2,2] * mat[4,4]))
  rhos.trt.byt[i,]$tot <- (-mat[2,5] / sqrt(mat[2,2] * mat[5,5]))
}

rhos.tpiv.tavg %>%
  apply(2, quantile, probs = c(0.025, 0.5, 0.975))
rhos.tpiv.wyt %>%
  apply(2, quantile, probs = c(0.025, 0.5, 0.975))
rhos.tpiv.byt %>%
  apply(2, quantile, probs = c(0.025, 0.5, 0.975))
rhos.trt.tavg %>%
  apply(2, quantile, probs = c(0.025, 0.5, 0.975))
rhos.trt.wyt %>%
  apply(2, quantile, probs = c(0.025, 0.5, 0.975))
rhos.trt.byt %>%
  apply(2, quantile, probs = c(0.025, 0.5, 0.975))

# Calculate proportion of posterior that is negative

apply(rhos.tpiv.tavg, 2, function(x){sum(x<0)/length(x)})
apply(rhos.tpiv.wyt, 2, function(x){sum(x<0)/length(x)})
apply(rhos.tpiv.byt, 2, function(x){sum(x<0)/length(x)})
apply(rhos.trt.tavg, 2, function(x){sum(x<0)/length(x)})
apply(rhos.trt.wyt, 2, function(x){sum(x<0)/length(x)})
apply(rhos.trt.byt, 2, function(x){sum(x<0)/length(x)})
