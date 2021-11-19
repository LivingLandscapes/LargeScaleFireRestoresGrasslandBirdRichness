##############################################################################
#### Large-scale fire management restores grassland bird richness in #########
#### private lands. Roberts et al. in Ecological Solutions & Evidence. #######
##############################################################################

##### NOTE: Because fire history (i.e., fire perimeters) are proprietary data,
#           the exact locations of fires cannot be shared. This means some  
#           figures and tests described in the manuscript cannot be reproduced 
#           here. Those interested in fire history data can contact the 
#           authors to request access.

#=============================================================================
## Preparations

# # Clean environment?
# rm(list=ls())

# List of packages necessary to run this script:
require(librarian)
shelf(tidyverse, data.table, cowplot, gridExtra, grid, here, sp, rgdal,
      raster, rgeos, mapproj, maptools, vroom, vegan, mgcv, MuMIn, gratia)

# Set paths
DataPath <- "LoessCanyons_BBS_Data"
ResultsPath <- "LoessCanyons_BBS_Results"
FigurePath <- "LoessCanyons_BBS_Figures"

# Source functions
source(here("mydraw.gam.R"))

# Set random seed
set.seed(90101)

#=============================================================================
## Loading data

# Load Loess Canyons ecoregion boundary as a SpatialPolygonsDataFrame
loessCanyons_boundary <- readOGR(dsn = here(DataPath,
                                            "LoessCanyons_BUL_Boundary"),
                                 layer = "Loess_canyons")

# Make tree cover raster images into a data.table for plotting predictions
treeTifFolder <- here(DataPath, 
                      "LoessCanyons_BBS_TreeCover_Images")
treeImages <- rbindlist(lapply(list.files(treeTifFolder),
                               function(X) {
                                 rast <- raster(here(treeTifFolder, X))
                                 rast <- projectRaster(rast,
                                                       crs = CRS("+proj=utm +zone=14 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"),
                                                       res = c(400, 400))
                                 temp <- as.data.frame(rast, xy = TRUE)
                                 setnames(temp, names(temp), c("easting", "northing", "mean"))
                                 temp$Year <- str_sub(X, 24, 27)
                                 return(na.omit(temp))
                               }))
treeImages[ , "Year_Num" := ifelse(Year == 2010, 1, 
                                   ifelse(Year == 2013, 3, 7))]


# Load bird survey + tree cover + (summarized) fire history data 
# for hierarchical generalized additive model
dat_grass <- fread(here(DataPath, "LoessCanyonsBBS_DataRaw.csv"))
##### Data column descriptions: 
# - Route: name of breeding bird survey route. 
# - Stop: ID for each stop along survey route+
# - Year: year of survey
# - Rich_Grass: grassland bird species richness 
# - Route_factor: name of breeding bird survey route. 
# - Year_Num: survey years re-numbered from 1 (2010) to 7 (2016)
# - Burned: if this stop was burned on or after the given year, 1, else 0
# - Year_Burned: the year that a given stop was burned
# - easting/northing: easting/northing coordinates in UTM 
# - count: number of 30x30m pixels within 400m of each stop. Number of pixels
#          vary because some are masked because they were not 'rangeland' 
#          pixels. See Methods/Data Collection/Tree Cover for details.
# - mean: mean percent tree cover across all 30x30m pixels within 400m of
#         each stop
# - stdDev: standard deviation of percent tree cover across all 30x30m pixels 
#           within 400m of each stop
# - TSF: years-since-fire

# Make sure that Route_factor column is a factor
dat_grass[ , "Route_factor" := as.factor(Route_factor)]

#=============================================================================
## Analysis
#=============================================================================

# Generalized additive model
gam_fit <- 
  gam(Rich_Grass ~ 
        te(Year_Num, easting, northing) + s(mean) + 
        s(Route_factor, bs = "re"),
      family = poisson,
      data = dat_grass,
      method = "REML")
# summary(gam_fit)
# appraise(gam_fit)
# draw(gam_fit)

# Plot the route factor fit
gam_plots <- mydraw.gam(gam_fit)
gam_TreeEffect_plot <- 
  gam_plots[[1]] +
  scale_x_continuous(breaks = seq(0, 50, 10)) +
  xlab("Percent Tree Cover (mean)") +
  ylab("Effect on Grassland\nBird Richness") +
  theme_bw() +
  theme(plot.title = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
# ggsave(plot = gam_TreeEffect_plot,
#        filename = here(FigurePath, "LoessCanyons_BBS_FigS2.jpeg"),
#        dpi = 600,
#        width = 6)

# Create a data.table for predicting grassland bird richness onto a map
nd <- 
  treeImages %>%
  mutate(Route_factor = "CallahanJeffrey")

# Predict grassland bird richness
predFit <- predict(gam_fit, 
                   newdata = nd,
                   type = "response",
                   se.fit = TRUE)
predFit_dt <- data.table(fit = predFit$fit,
                         CI = predFit$se.fit * 1.645,
                         nd)

# Get info on differences in richness by pixel
predFit_dt$ID <- rep(1:as.numeric(predFit_dt[ , .N, by = "Year"][1,2]), 3)

# Color ramps
treeBreaks <- seq(0, 
                  max(predFit_dt$mean) + 1, 
                  1)
treeRamp <- colorRampPalette(c("gold",
                               "goldenrod",
                               "goldenrod2", 
                               "goldenrod4",
                               "brown",
                               "brown4",
                               "black"))(length(treeBreaks) - 1)

predBreaks <- seq(0, 
                  max(predFit_dt$fit) + 0.05, 
                  0.05)
predRamp <- colorRampPalette(c("black", 
                               "darkred",
                               "red",
                               "orange", 
                               "palegoldenrod", 
                               "greenyellow",
                               "green",
                               "green3"))(length(predBreaks) - 1)

# Plot the maps
treeMap <-
  ggplot() +
  coord_equal() +
  geom_tile(data = predFit_dt,
            mapping = aes(x = easting, y = northing, fill = mean)) +
  scale_fill_gradientn(colors = treeRamp,
                       limits = c(0, 94),
                       name = "Mean Percent\nTree Cover") +
  facet_wrap(~ Year) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.box.margin = unit(rep(0, 4), "cm"),
        legend.title = element_text(size = 11),
        strip.text = element_text(size = 12)) +
  ylab("Northing") +
  xlab("Easting")
# treeMap
predMap <- 
  ggplot() +
  coord_equal() +
  geom_tile(data = predFit_dt,
            mapping = aes(x = easting, y = northing, fill = fit)) +
  scale_fill_viridis_c(
    name = "Grassland\nBird Richness",
    breaks = seq(1, 7, 1)) +
  facet_wrap(~ Year) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.box.margin = unit(rep(0, 4), "cm"),
        legend.title = element_text(size = 11),
        strip.text = element_text(size = 12)) +
  ylab("Northing") +
  xlab("Easting")
# ggsave(here(FigurePath, "PredictedMap_2.jpeg"),
#        dpi = 600,
#        width = 11, height = 5)


#========================================================================
## Plot grassland bird richness + tree cover change maps

# Calculate the change in richness from 2010 - 2016
PixelDiff <- predFit_dt[Year %in% c(2010, 2016) , list(diff = diff(fit)), by = "ID"]
# range(PixelDiff$diff)
# hist(PixelDiff$diff)

# How much increasing?
increase_df <- 
  data.frame(MagnitudeOfIncrease = c("< 0", "> 0", "> 0.50", "> 1", "> 2"),
             PercentOfCanyonsIncreasing = c(nrow(PixelDiff[diff < 0])/nrow(PixelDiff),
                                            nrow(PixelDiff[diff > 0])/nrow(PixelDiff),
                                            nrow(PixelDiff[diff > 0.5])/nrow(PixelDiff),
                                            nrow(PixelDiff[diff > 1])/nrow(PixelDiff),
                                            nrow(PixelDiff[diff > 2])/nrow(PixelDiff))) %>%
  mutate(HectaresIncreasing = PercentOfCanyonsIncreasing * nrow(PixelDiff) * 400^2 * 0.0001,
         PercentOfCanyonsIncreasing = round(PercentOfCanyonsIncreasing, 2) * 100)
# increase_df
# write.csv(increase_df, 
#           file = here(ResultsPath, "IncreaseRichnessTable.csv"),
#           row.names = FALSE)

# Grassland bird richness change for plotting
change_dt <- copy(PixelDiff)
change_dt <- left_join(change_dt, predFit_dt[Year == 2016])

# Grassland bird richness change map
changeMap <- 
  ggplot() +
  coord_equal() +
  geom_tile(data = change_dt,
            mapping = aes(x = easting, y = northing, fill = diff)) +
  scale_fill_gradient2(
    name = "Change in\ngrassland\nbird richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 10),
        legend.box.margin = unit(rep(0, 4), "cm"),
        legend.title = element_text(size = 9)) +
  ylab("Northing") +
  xlab("Easting")
# changeMap
# ggsave(plot = changeMap,
#        filename = here(FigurePath, "ChangeMap.jpeg"),
#        dpi = 600)

# Difference in tree cover from 2010 - 2016
TreeDiff <- copy(predFit_dt)[Year %in% c(2010, 2016) , list(diff = diff(mean)), by = "ID"]
treeChange_dt <- left_join(TreeDiff, 
                           predFit_dt[Year == 2010,
                                      c("mean", "ID",
                                        "easting", "northing")])
treeChange_dt <- left_join(treeChange_dt, 
                           predFit_dt[Year == 2016,
                                      c("ID")])
treeChange_dt[ , "Year_Burned" := as.factor(Year_Burned)]

# 
treeChangeMap <- 
  ggplot() +
  coord_equal() +
  geom_tile(data = treeChange_dt,
            mapping = aes(x = easting, y = northing, fill = diff)) +
  scale_fill_gradient2(name = "Change in\npercent mean\ntree cover") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 10),
        legend.box.margin = unit(rep(0, 4), "cm"),
        legend.title = element_text(size = 9)) +
  ylab("Northing") +
  xlab("Easting")
# treeChangeMap
# ggsave(plot = treeChangeMap,
#        filename = here(FigurePath, "TreeChangeMap.jpeg"),
#        dpi = 600)

# Combine the bird/tree change plots
changePlots <- 
  plot_grid(treeChangeMap, changeMap, 
            ncol = 1,
            align = "hv",
            labels = c("A", "B"),
            hjust = -2)
# ggsave(plot = changePlots,
#        filename = here(FigurePath, "LoessCanyons_BBS_Fig3.jpeg"),
#        dpi = 600,
#        height = 8,
#        width = 6)

#=============================================================================
## Response of richness to change in tree cover + years-since-fire

setnames(treeChange_dt, "diff", "TreeDiff")
compareDT <- left_join(change_dt, 
                       treeChange_dt[,c("ID", "TreeDiff", 
                                        "easting", "northing")])

##### NOTE: Because fire history (i.e., fire perimeters) are proprietary data,
#           the exact locations of fires cannot be shared.
# #
# tst_cor_sample <- 
#   compareDT %>% 
#   sample_n(nrow(compareDT) * 0.5)
# 
# gam_changeFit <- gam(diff ~ te(TSF, TreeDiff),
#                      method = "REML",
#                      data = tst_cor_sample)
# gam_changePlot <- 
#   mydraw.gam(gam_changeFit)[[1]] + 
#   # theme_minimal() +
#   theme(plot.title = element_blank(),
#         axis.title = element_text(size = 12),
#         axis.text = element_text(size = 11)) +
#   scale_fill_viridis_b() +
#   scale_x_continuous(breaks = c(-1, 0, 1, 3, 5, 10, 14)) +
#   xlab("Years-since-fire") +
#   ylab("Change in\npercent mean tree cover")
# ggsave(plot = gam_changePlot,
#        filename = here(FigurePath, "LoessCanyons_BBS_FigS4.jpeg"),
#        dpi = 600,
#        height = 4,
#        width = 7)
# summary(gam_changeFit)
# appraise(gam_changeFit)
# draw(gam_changeFit)

# Pearson's Correlation Coefficient between Mean % tree cover in 2016 and 
# TreeDiff was only 0.59.
cor(compareDT$mean, compareDT$TreeDiff)

#========================================================================
# Plot a supplemental figue: confidence map

confidenceMap <- 
  ggplot() +
  coord_equal() +
  geom_tile(data = predFit_dt,
            mapping = aes(x = easting, y = northing, 
                          fill = CI)) +
  scale_fill_gradient(
    name = "Predicted error\nin grassland\nbird richness\n(90% confidence)") +
  facet_wrap(~ Year) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 8),
        legend.box.margin = unit(rep(0, 4), "cm"),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 12)) +
  ylab("Northing") +
  xlab("Easting")
# confidenceMap
# ggsave(plot = confidenceMap,
#        filename = here(FigurePath, "LoessCanyons_BBS_FigS3.jpeg"),
#        dpi = 600,
#        width = 8)

