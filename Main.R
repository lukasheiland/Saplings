packages <- c("dplyr", "tidyr", "magrittr", "raster", "sf", "sdmTMB", "parallel", "ggplot2", "cowplot")
sapply(packages, FUN = library, character.only = T)

# Sourcing -------------------------------------
source("Functions.R")
datapath <- "Data/"


# Load and prepare data ---------------------------------------------------
B <- readRDS(file.path(datapath, "DE_BWI_big_abund.rds"))
S <- readRDS(file.path(datapath, "DE_BWI_small_abund.rds"))
E <- readRDS(file.path(datapath, "DE_BWI_Env_sf.rds"))
Stack <- readRDS(file.path(datapath, "Predictorstack.rds")) ## variables in stack: names(Stack) %>% dput()
# st_crs(Stack) == st_crs(E)
Newdata <- as.data.frame(Stack, xy = T) %>%
  rename(X = x, Y = y)

D <- combineData(B, S, E)
rm(B, S, E)

#### Selection of species and predictors -------
taxa <- D$tax %>% unique() %>% as.character()
species <- taxa[c(13, 15, 32, 25)]
predictor <- c("phCaCl_esdacc", "sand_esdact", "cn_esdacc", "tPeriodic2010_mh", "precPeriodic2010_mh",
               "ba_ha_plot")
predictor_rf <- c("ba_ha_species")
predictor_all <- c(predictor, predictor_rf)

#### Subset -------
data_species <- lapply(species, subsetData, D = D, predictor = predictor_all, n_plots = 10000)
names(data_species) <- species



# Fit and predict ---------------------------------------------------
fits <- mclapply(data_species, fitModel, predictor_fixed = predictor, predictor_randomfield = predictor_rf, mc.cores = getOption("mc.cores", 4L))
predictions <- mclapply(fits, predictFit, Newdata = Newdata, mc.cores = getOption("mc.cores", 4L))


# Plot and save ---------------------------------------------------
# vars <- c("est", "est_non_rf", "est_rf", "omega_s", "zeta_s_ba_ha_species")

plots <- lapply(predictions, plotMap, var = est)
plotgrid <- cowplot::plot_grid(plotlist = plots, ncol = 2)
ggsave("Maps.pdf", plotgrid, width = 20, units = "cm")

plotgrid_spatial_random_intercept <- lapply(predictions, plotMap, var = omega_s) %>% cowplot::plot_grid(plotlist = .)
plotgrid_spatial_random_intercept
# plotgrid_spatial_ba_ha_species <- lapply(predictions, plotMap, var = zeta_s_ba_ha_species) %>% cowplot::plot_grid(plotlist = .)


