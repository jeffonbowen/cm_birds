### Process Data ###


# Setup -------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(sf)
library(mapview)


# Load species name reference list
species_ref <- read_excel("dat/SpeciesList_ref.xlsx", sheet = 1) |> 
  select(`Common Name` = `English Name`, `Bird Code`, Order, `Breeding Bird`) |> 
  dplyr::filter(!is.na(`Breeding Bird`))


# Load the dillon data and assemble ----------------------------------------

dat_survey <- read_excel("dat/dillon/A2. Survey data.xlsx", sheet = "reformat") |> 
  rename(location = "Point Count Name")

# Load the BBS data
dat1 <- read_excel("dat/dillon/A3. BBS Data, July 4-July 11 2019.xlsx", sheet = "wide") |> 
  pivot_longer(!(1:4), 
               names_to = "Common Name", 
               values_to = "count") |> 
  filter(!is.na(count))

dat2 <- read_excel("dat/dillon/A3. BBS Data, June 29-July 3 2014.xlsx", sheet = "wide") |> 
  pivot_longer(!(1:4), 
               names_to = "Common Name", 
               values_to = "count") |> 
  filter(!is.na(count))

dat3 <- read_excel("dat/dillon/A3. BBS Data, June 30-July 7 2017.xlsx", sheet = "wide") |> 
  pivot_longer(!(1:4), 
               names_to = "Common Name", 
               values_to = "count") |> 
  filter(!is.na(count))

dat4 <- read_excel("dat/dillon/A3. BBS Data, June 30-July 8 2018.xlsx", sheet = "wide") |> 
  pivot_longer(!(1:4), 
               names_to = "Common Name", 
               values_to = "count") |> 
  filter(!is.na(count))

dat5 <- read_excel("dat/dillon/A3. BBS Data, June 5-10 2014.xlsx", sheet = "wide") |> 
  pivot_longer(!(1:4), 
               names_to = "Common Name", 
               values_to = "count") |> 
  filter(!is.na(count))

dat6 <- read_excel("dat/dillon/A3. BBS Data, June 5-June 11 2017.xlsx", sheet = "wide") |> 
  pivot_longer(!(1:4), 
               names_to = "Common Name", 
               values_to = "count") |> 
  filter(!is.na(count))

dat_sp <- bind_rows(dat1, dat2, dat3, dat4, dat5, dat6) |> 
  rename(location = "Survey Point Name") 

dat <- dat_sp |> 
  left_join(dat_survey, by = "location") |> 
  mutate(
    # Extract time part from Time_Start
    time_only = format(as.POSIXct(`Time Start`), "%H:%M:%S"),
    # Combine correct Date with extracted time
    datetime = as.POSIXct(paste(Date, time_only), format = "%Y-%m-%d %H:%M:%S"),
    Year = year(datetime),
    DUR = 7,
    location = paste0(Year, "-", location),
    SurveyID = paste0(location, "_", datetime)) |> 
  filter(!is.na(Easting)) |> 
  st_as_sf(coords = c("Easting", "Northing"), crs = 32611) |>
  st_transform(crs = 4326) |> 
  select(location, datetime, Year, SurveyID, DUR, `Common Name`, count)
  
dat$Longitude = st_coordinates(dat)[, 1]
dat$Latitude = st_coordinates(dat)[, 2]

mapview(dat)  

# Fix a few species names and remove unidentified species

dat <- dat |> 
  filter(!grepl("^Unid", `Common Name`)) |> 
  mutate(`Common Name` = case_when(
    `Common Name` == "Clark's nuthatch" ~ "Clark's Nutcracker",
    `Common Name` == "Cassin's finch" ~ "Cassin's Finch",
    `Common Name` == "Gray Jay" ~ "Canada Jay",
    `Common Name` == "Harry Woodpecker" ~ "Hairy Woodpecker",
    `Common Name` == "Northern Flicker" ~ "Northern Flicker",
    `Common Name` == "MacGillvery's Warbler" ~ "MacGillivray's Warbler",
    `Common Name` == "Pacific-slope Flycatcher" ~ "Western Flycatcher",
    `Common Name` == "Red-winged BlackBird" ~ "Red-winged Blackbird",
    `Common Name` == "Tenesse Warbler" ~ "Tennessee Warbler",
    `Common Name` == "Vesper's Sparrow" ~ "Vesper Sparrow",
    `Common Name` == "Violet-Green Swallow" ~ "Violet-green Swallow",
    `Common Name` == "Western-wood Pewee" ~ "Western Wood-Pewee",
    `Common Name` == "Yellow-Rumped Warbler" ~ "Yellow-rumped Warbler",
    TRUE ~ `Common Name`
  ))

dat <- dat |> left_join(species_ref, by = "Common Name") |> 
  select(location, datetime, Year, SurveyID, DUR, Latitude, Longitude, 
         `Common Name`, `Bird Code`, count, Order)

qa <- dat |> 
  group_by(`Bird Code`, `Common Name`, Order) |> 
  summarize(n = n(), .groups = "drop") |> 
  arrange(`Bird Code`) |> 
  filter(is.na(`Bird Code`))

write_csv(dat, "dat/dillon/dillon_bird_data.csv")

dat_dillon <- read_csv("dat/dillon/dillon_bird_data.csv") |> 
  select(!geometry)


# TT Data -----------------------------------------------------------------

surveys <- read_excel("dat/tt/NWP Data Songbird Data Summer 2024 (RAW DATA).xlsx",
                      sheet = 1) |> 
  mutate(
    time_only = format(as.POSIXct(Time), "%H:%M:%S"),
    datetime = as.POSIXct(paste(Date, time_only), format = "%Y-%m-%d %H:%M:%S"),
    Year = 2024,
    DUR = 10)|> 
  select(location = `Sample Station Label`,
         datetime, Year, SurveyID, Latitude, Longitude, DUR)

obs <- read_excel("dat/tt/NWP Data Songbird Data Summer 2024 (RAW DATA).xlsx",
                  sheet = 2) |> 
  select(SurveyID, `Common Name`, `Total Count`) |>
  group_by(SurveyID, `Common Name`) |>
  summarize(count = sum(`Total Count`)) |> 
  filter(!is.na(`Common Name`))

dat_tt <- obs |> left_join(surveys, by = "SurveyID") |> 
  filter(!is.na(location))

dat_tt <- dat_tt |> left_join(species_ref, by = "Common Name") |> 
  mutate(
    `Common Name` = case_when(
      `Common Name` == "Gray Jay" ~ "Canada Jay",
      TRUE ~ `Common Name`),
    `Bird Code` = case_when(
    `Common Name` == "NONE" ~ "NONE",
    TRUE ~ `Bird Code`),
    location = paste0(Year, "-", location)) |> 
  select(location, datetime, Year, SurveyID, DUR, Latitude, Longitude, 
       `Common Name`, `Bird Code`, count, Order)



# Combine -----------------------------------------------------------------

dat_combined <- bind_rows(dat_dillon, dat_tt)

write_csv(dat_combined, "dat/BBS_combined.csv")


# View --------------------------------------------------------------------

dat_combined <- read_csv("dat/BBS_combined.csv")

dat_sp <- dat_combined |> 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

mapview(dat_sp, zcol = "Year")
