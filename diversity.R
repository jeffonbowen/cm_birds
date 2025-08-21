## Crown Mountain Diversity Modelling ###

# Includes both richness and abundance

library(tidyverse)
library(wildrtrax)
library(unmarked)

library(MuMIn)                # model selection
library(effects)            
library(emmeans)              # for estimating marginal means
library(hms)

# Prep Data --------------------------------------------------------------

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
  select("location", "datetime",  "Year", "SurveyID", "DUR", "latitude",
         "longitude", "TSSR", "JDAY", "LCC2",
         "LCC4", "TREE", "n_species", "n_detections")  

write_csv(div_dat, "dat/div_dat.csv")

div_dat_sp <- st_as_sf(div_dat, 
                       coords = c("longitude", "latitude"), crs = 4326)

st_write(div_dat_sp, "dat/div_dat_sp.gpkg")

# Naive --------------------------------------------------------------------

table(div_dat$Year, div_dat$LCC2)




# GLMM ---------------------------------------------------------------------






# N-Mixture ----------------------------------------------------------------



library(elevatr)
library(sf)
library(raster)

# Define AOI
bbox <- st_bbox(c(xmin = -114.9, ymin = 49.6, xmax = -114.7, ymax = 49.9), crs = st_crs(4326))
aoi <- st_as_sfc(bbox)
aoi <- st_sf(geometry = aoi)

# Plot AOI
plot(aoi, main = "Area of Interest")

# Get DEM (zoom level 11)
dem <- get_elev_raster(aoi, z = 11, clip = "bbox")

# Plot DEM
plot(dem, main = "Digital Elevation Model")


