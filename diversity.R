## Crown Mountain Diversity Modelling ###

# Includes both richness and abundance

{
library(tidyverse)
library(wildrtrax)

library(lme4)
library(glmmTMB)
library(unmarked)
library(terra)
library(MuMIn)                # model selection
library(effects)            
library(emmeans)              # for estimating marginal means
}

# Prep Diversity Data --------------------------------------------------------------

dat <- read_csv("dat/BBS_combined.csv") |> 
  mutate(Year = factor(Year))

# Summarize data for richness and abundance per visit (a survey)   
# This needs to include songbirds and woodpeckers. 

div_dat <- dat |>
  filter(Order %in% c("Passeriformes", "Piciformes")) |> 
  group_by(location, datetime, Year, SurveyID, DUR, Latitude, Longitude) |> 
  summarize(n_species = n_distinct(`Common Name`),
            n_detections = sum(count)) |> 
  ungroup()

# Add Covariates ---------------------------------------------------------

# Survey-level covariates, like time of day and day of the year are needed
# to estimate richness and abundance. There is a function in wildrtrax that
# can be used get these.

# Its an internal function and needs to be called with ":::"
Sys.setenv(WT_USERNAME = 'Jeffm', WT_PASSWORD = 'Engineer592')
wt_auth()

# These fields are required for WildRTrax
div_dat <- div_dat |> 
  mutate(recording_date_time = datetime,
         task_duration = DUR*60) |> 
  rename(latitude = Latitude,
         longitude = Longitude)

vars <- wildrtrax:::.make_x(div_dat)

div_dat <- bind_cols(div_dat, vars) |> 
  dplyr::select("location", "datetime",  "Year", "SurveyID", "DUR", "latitude",
         "longitude", "TSSR", "JDAY", "LCC2",
         "LCC4", "TREE", "n_species", "n_detections")  

# Now add elevation
div_dat <- st_as_sf(div_dat, 
                       coords = c("longitude", "latitude"), crs = 4326,
                    remove = FALSE)
dem <- rast("dat/demz11.tiff")
elev <- terra::extract(dem, div_dat)
div_dat$elevation <- elev$demz11

# Add Canopy Height
canopy <- rast("C:/Users/jeff.matheson/OneDrive - Tetra Tech, Inc/Documents/GIS_Projects/Crown_Mountain/EO/ETH_GlobalCanopyHeight_10m_2020_N48W117_Map.tif")
# Create 50 m buffer around each point
buffers <- st_buffer(div_dat, dist = 50)
# Calculate mean canopy height within each buffer
library(exactextractr)
mean_canopy <- exactextractr::exact_extract(canopy, buffers, 'mean')
#Add mean canopy height to dataframe
div_dat$canopy_height <- mean_canopy

st_write(div_dat, "dat/div_dat_sp.gpkg", append = FALSE)

write_csv(div_dat, "dat/div_dat.csv")


# Naive --------------------------------------------------------------------

div_dat <- read_csv("dat/div_dat.csv")

table(div_dat$Year, div_dat$LCC2)

richness_estimates <- div_dat |> 
  group_by(location, Year) |> 
  summarize(mean_raw = mean(n_species))

# Global mean survey richness
mean(richness_estimates$mean_raw)

# GLMM ---------------------------------------------------------------------

library(glmmTMB)

# Richness
div_dat$DUR <- as.numeric(div_dat$DUR)
div_dat$elevation_scaled <- as.numeric(scale(div_dat$elevation))


mod1 <- glmmTMB(n_species ~ Year + elevation_scaled + canopy_height + (1|DUR) + 
                  (1|location) + (1|JDAY) + (1|TSSR), 
              data = div_dat, 
              family = (poisson("log")),
              na.action = na.fail)
mod1

d <- dredge(mod1)
d

ae <- allEffects(mod1)
plot(ae, residuals="TRUE")
ae

# Try various models to determine what random effects are worth including.
mod2 <- glmmTMB(n_species ~ elevation_scaled + canopy_height + (1|DUR) + 
                  (1|location) + (1|TSSR), 
                data = div_dat, 
                family = (poisson("log")),
                na.action = na.fail)

mod3 <- glmmTMB(n_species ~ elevation_scaled + canopy_height + (1|DUR) + 
                  (1|location) + (1|JDAY), 
                data = div_dat, 
                family = (poisson("log")),
                na.action = na.fail)

mod4 <- glmmTMB(n_species ~ elevation_scaled + canopy_height + 
                  (1|location) + (1|TSSR), 
                data = div_dat, 
                family = (poisson("log")),
                na.action = na.fail)

mod5 <- glmmTMB(n_species ~ elevation_scaled + canopy_height + 
                  (1|location), 
                data = div_dat, 
                family = (poisson("log")),
                na.action = na.fail)

mod6 <- glmmTMB(n_species ~ (1|Year) + elevation_scaled + canopy_height +
                  (1|location), 
                data = div_dat, 
                family = (poisson("log")),
                na.action = na.fail)
dredge(mod6)


anova(mod1, mod2, mod3, mod4, mod5, mod6)

mod_top <- mod5

# Predict richness for each location surveyed. 

richness_glmm <- div_dat |> 
  mutate(rich_est = predict(mod_top, newdata = div_dat, 
                                 type = "response",
                                 re.form = NA)) |> 
  group_by(location) |> 
    summarize(mean_glmm = mean(rich_est))
  
richness_estimates <- richness_estimates |> 
  left_join(richness_glmm, by = "location")

mean(richness_estimates$mean_glmm)


# N-Mixture ----------------------------------------------------------------

library(unmarked)
library(tidyr)
library(dplyr)
library(tidyverse)

# Create visit order within site-year
div_dat <- div_dat %>%
  group_by(location, Year) %>%
  mutate(visit = row_number()) %>%
  ungroup()
table(div_dat$visit)

# Pivot to wide format for counts
y_wide <- div_dat %>%
  dplyr::select(location, Year, visit, n_species) %>%
  pivot_wider(names_from = visit, values_from = n_species)

# Site covariates (take first row per site-year)
site_covs <- div_dat %>%
  group_by(location, Year) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(elev_scaled = scale(elevation),
         mon_year = Year - 2014) |> 
  dplyr::select(location, latitude, longitude, Year, mon_year, elevation, elev_scaled, canopy_height)

# Observation covariates
site_order <- y_wide$location
y <- y_wide %>%
  dplyr::select(-location) %>%
  as.matrix()

# JDAY
JDAY_mat <- div_dat %>%
  dplyr::select(location, visit, JDAY) %>%
  pivot_wider(names_from = visit, values_from = JDAY, values_fill = NA) %>%
  slice(match(site_order, location)) %>%
  dplyr::select(-location) %>%
  as.matrix()

# DUR
DUR_mat <- div_dat %>%
  dplyr::select(location, visit, DUR) %>%
  pivot_wider(names_from = visit, values_from = DUR, values_fill = NA) %>%
  slice(match(site_order, location)) %>%
  dplyr::select(-location) %>%
  as.matrix()

# TSSR
TSSR_mat <- div_dat %>%
  dplyr::select(location, visit, TSSR) %>%
  pivot_wider(names_from = visit, values_from = TSSR, values_fill = NA) %>%
  slice(match(site_order, location)) %>%
  dplyr::select(-location) %>%
  as.matrix()

# Combine into a list
obs_covs_list <- list(JDAY = JDAY_mat,
                      DUR  = DUR_mat,
                      TSSR = TSSR_mat)

# Build unmarked frame
umf <- unmarkedFramePCount(y = as.matrix(y_wide[,-c(1:2)]),
                           siteCovs = site_covs[,-c(1:1)],
                           obsCovs = obs_covs_list)


fm1 <- pcount(~ JDAY + TSSR + DUR ~ elev_scaled + canopy_height, data = umf, K = 60)
summary(fm1)

fm2 <- pcount(~ JDAY + DUR ~ elev_scaled + canopy_height, data = umf, K = 60)

fm3 <- pcount(~ DUR ~ elev_scaled + canopy_height, data = umf, K = 60)

fm4 <- pcount(~ JDAY ~ elev_scaled + canopy_height, data = umf, K = 60)
summary(fm4)

fm5 <- pcount(~1 ~ elev_scaled + canopy_height, data = umf, K = 60)
summary(fm5)

fm6 <- pcount(~1 ~ mon_year + elev_scaled + canopy_height, data = umf, K = 60)
summary(fm6)

fm7 <- pcount(~ JDAY + TSSR + DUR ~ mon_year + elev_scaled + canopy_height, data = umf, K = 60)
summary(fm7)

fm8 <- pcount(~ JDAY + DUR ~ mon_year + elev_scaled + canopy_height, data = umf, K = 60)
summary(fm8)

fm9 <- pcount(~ JDAY + TSSR ~ mon_year + elev_scaled + canopy_height, data = umf, K = 60)
summary(fm9)

fm10 <- pcount(~ TSSR ~ mon_year + elev_scaled + canopy_height, data = umf, K = 60)
summary(fm10)

# Look for best supported model
dredge(fm7)

AIC(fm1, fm2, fm3, fm4, fm5, fm6, fm7, fm8, fm9, fm10)

fm_top <- fm10

# Function returning three fit-statistics.
fitstats <- function(fm) {
  y    <- getY(fm)
  mu   <- fitted(fm)           # expected counts per survey (same shape as y)
  res  <- y - mu
  SSE  <- sum(res^2, na.rm=TRUE)
  FT   <- sum((sqrt(y) - sqrt(mu))^2, na.rm=TRUE)
  CHI  <- sum((res^2) / pmax(mu, 1e-6), na.rm=TRUE)  # Pearson χ²
  c(SSE = SSE, FT = FT, CHI = CHI)
}

pb <- parboot(fm_top, fitstats, nsim = 200, report = 1)
# c-hat from Pearson:
c_hat <-c_hat pb@t0["CHI"] / mean(pb@t.star[,"CHI"], na.rm=TRUE)
c_hat

# Model fit is ok
# Quick check to see if nbinom might be better. No. 
fm11 <- pcount(~ TSSR ~ mon_year + elev_scaled + canopy_height, 
               data = umf, K = 60, mixture  = "NB")
summary(fm11)
AIC(fm10, fm11)


# Predict for each location -------------------------------------------------

# First create location dataframe
loc_preds <- div_dat %>%
  group_by(location) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(location, latitude, longitude, Year, elevation, canopy_height)

loc_preds$elev_scaled <- as.numeric(scale(loc_preds$elevation))
loc_preds$mon_year <- as.numeric(as.character(loc_preds$Year)) - 2014

# Now predict
richness_nmixture <- predict(fm_top, type = "state", 
                     newdata = loc_preds, backTransform = T)

# Add to richness dataframe
richness_estimates$mean_nmixture <- richness_nmixture$Predicted
  
colMeans(richness_estimates[, 2:4], na.rm = TRUE)

richness_estimates %>%
  group_by(Year) |> 
  summarise(across(c("mean_raw", "mean_glmm", "mean_nmixture"), ~ mean(.x, na.rm = TRUE)))

plot(richness_estimates$mean_raw, richness_estimates$mean_glmm)
plot(richness_estimates$mean_glmm, richness_estimates$mean_nmixture)

cor_mat <- cor(richness_estimates[3:5],
               use = "pairwise.complete.obs",
               method = "spearman")
cor_mat

richness_estimates <- richness_estimates |> 
  left_join(loc_preds, by = c("location", "Year"))

richness_estimates <- st_as_sf(richness_estimates,
                              coords = c("longitude", "latitude"), crs = 4326)

write_csv(richness_estimates, "dat/richness_estimates.csv")
st_write(richness_estimates, "dat/richness_estimates.gpkg")


# Map ----------------------------------------------------------------------

library(sf)
library(leaflet)
library(RColorBrewer)
library(scales)

richness_estimates <- st_read("dat/richness_estimates.gpkg")

fp_path <- "C:/Users/jeff.matheson/OneDrive - Tetra Tech, Inc/Documents/GIS_Projects/Crown_Mountain/footprint/Project_Footprint.shp"

fp_path <- "C:/Users/jeff.matheson/OneDrive - Tetra Tech, Inc/Documents/R_Projects/cm_birds/dat/Project_Footprint/Project_Footprint.shp"

fp <- st_read(fp_path) |> st_transform(crs = 4326) |> st_zm()

plot(fp)


# Define the palette using RColorBrewer
pal <- colorNumeric(
  palette = (brewer.pal(9, "YlOrRd")),  # blue → green → red
  domain = richness_estimates$mean_nmixture
)

leaflet(richness_estimates) |>
  addTiles() |>
  addPolygons(
    data = fp,
    fillOpacity = 0.5
  ) |>
  addCircleMarkers(
    radius = ~scales::rescale(mean_nmixture^2, c(2, 20)),  # size scaling
#    radius = ~mean_nmixture,
    color = ~pal(mean_nmixture),                   # apply color palette
    stroke = FALSE,
    fillOpacity = 0.5
  ) |>
  addLegend(
    "bottomright",
    pal = pal,
    values = ~mean_nmixture,
    title = "Predicted"
  )





leaflet() |>
  addTiles() |> 
  addPolygons(data = fp)

