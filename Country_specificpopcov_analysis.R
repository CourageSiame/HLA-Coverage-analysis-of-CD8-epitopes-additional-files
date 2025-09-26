# Install required packages if not already installed
# You can run these lines once, then comment them out
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("maps") # For world map data
# install.packages("tidyverse")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(maps)
library(tidyverse)
library(ggrepel)

# Step 1: Read the CSV file
data <- read.csv("epitope_data_countries.csv")

# Define the group of epitopes
selected_epitopes <- c("CSPep5", "CSPep13", "TRAPep1", "CSPep11")
#selected_epitopes <- "CSPep5"
# Filter data and convert to percentage
transformeddata <- data %>%
  filter(Epitopes %in% selected_epitopes) %>%
  mutate(across(c(2:79), function(x) x * 100))

# Step 2: Reshape the data to long format
long_data <- transformeddata %>%
  pivot_longer(
    cols = -Epitopes,
    names_to = "Country",
    values_to = "Coverage"
  ) %>%
  rename(Epitope = Epitopes)

# Step 3: Standardize country names and filter
# The key change is here: we first clean the country names by replacing periods with spaces.
long_data <- long_data %>%
  mutate(Country = str_replace_all(Country, "\\.", " ")) %>%
  mutate(
    Standard_Country = case_when(
      Country %in% c("England", "Scotland", "Wales", "Ireland Northern") ~ "UK",
      Country == "Ireland South" ~ "Ireland",
      Country == "Ivory Coast" ~ "Ivory Coast",
      Country == "Macedonia" ~ "North Macedonia",
      Country == "United States" ~ "USA",
      Country == "Guinea Bissau" ~ "Guinea-Bissau",
      TRUE ~ Country # This now works correctly for all other countries with spaces.
    )
  ) %>%
  filter(!is.na(Standard_Country) & !is.na(Coverage))

# Step 4: Load world map data
world <- map_data("world")

# Step 5: Join the data with the map
map_data_joined <- world %>%
  left_join(long_data, by = c("region" = "Standard_Country"), relationship = "many-to-many") %>%
  filter(!is.na(Coverage))

# Step 6: Get centroids for country labels
centroids <- map_data_joined %>%
  group_by(Epitope, region) %>%
  summarise(
    long = mean(long),
    lat = mean(lat),
    .groups = 'drop'
  ) %>%
  rename(Standard_Country = region)

# Step 7: Order the facets
map_data_joined$Epitope <- factor(map_data_joined$Epitope, levels = sort(unique(map_data_joined$Epitope)))
centroids$Epitope <- factor(centroids$Epitope, levels = sort(unique(centroids$Epitope)))

# Step 8: Create the plot
Plot_TPEPS <- ggplot(map_data_joined, aes(x = long, y = lat)) +
  geom_polygon(aes(fill = Coverage, group = group), color = "black", size = 0.1) +
  geom_text_repel(data = centroids, aes(x = long, y = lat, label = Standard_Country),
                  size = 2, color = "black", min.segment.length = 0.5,
                  max.overlaps = Inf) +
  coord_fixed(1.3) +
  facet_wrap(~ Epitope, ncol = 2, dir = "h") +
  scale_fill_gradient(low = "white", high = "darkred", na.value = "grey90", limits = c(0, 100)) +
  theme_minimal() +
  labs(
    title = "",
    fill = "%PC",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12)
  )

print(Plot_TPEPS)

# Step 9: Save the plot
ggsave("hla_coverage_worldmap2.jpeg", width = 16, height = 12, dpi = 300)


#####Visualizing all epitopes coverage############

#Step 1: Transform the data into percentages
# converting the values to percentage
scalar <- 100  
datatransformed <- data %>% mutate(across(c(2:79), function(x) x * scalar))


# Step 2: Reshape the data to long format
# The original data is wide (epitopes as rows, countries as columns).
# We use pivot_longer to make it long: one row per epitope-country combination.
# 'names_to = "Country"' takes the country names from column headers.
# 'values_to = "Coverage"' stores the coverage values.
longData <- datatransformed %>%
  pivot_longer(
    cols = -Epitopes, # All columns except the first (Epitopes)
    names_to = "Country",
    values_to = "Coverage"
  ) %>%
  rename(Epitope = Epitopes) # Rename the first column to 'Epitope' for clarity


# Step 3: Standardize country names for matching with map data and create short forms
# We use mutate with case_when to create a new column 'Standard_Country' with corrected names
# and a 'Short_Name' column for the short forms to be used in the legend.
# Entries with NA coverage will be filtered out.
longData <- longData %>%
  mutate(Country = str_replace_all(Country, "\\.", " ")) %>%
  mutate(
    Standard_Country = case_when(
      Country %in% c("England", "Scotland", "Wales", "Ireland Northern") ~ "UK",
      Country == "Ireland South" ~ "Ireland",
      Country == "Ivory Coast" ~ "Ivory Coast",
      Country == "Macedonia" ~ "North Macedonia",
      Country == "United States" ~ "USA",
      Country == "Guinea Bissau" ~ "Guinea-Bissau",
      TRUE ~ Country # This now works correctly for all other countries with spaces.
    )
  ) %>%
  filter(!is.na(Standard_Country) & !is.na(Coverage))

# Step 4: Load world map data
# map_data("world") provides a data frame with longitude, latitude, and region (country) for plotting polygons.
World <- map_data("world")

# Step 5: Join the data with the map
# We left-join the map data with our long_data on 'region' == 'Standard_Country'.
# This adds coverage values to the map polygons where countries match.
mapData_joined <- World %>%
  left_join(longData, by = c("region" = "Standard_Country"), relationship = "many-to-many") %>%
  # Filter out map regions that have no corresponding data (i.e., countries with NA coverage)
  filter(!is.na(Coverage))

# Step 6: Create the plot

Plot_all <- ggplot(mapData_joined, aes(x = long, y = lat, group = group, fill = Coverage)) +
  geom_polygon(color = "black", size = 0.1) + # Draw country borders
  coord_fixed(1.3) + # Fix aspect ratio for world map
  facet_wrap(~ Epitope, ncol = 4) + # Facet by epitope, 4 columns for layout
  scale_fill_gradient(low = "white", high = "darkred", na.value = "grey90") + # Color scale; grey for missing
  theme_minimal() + # Clean theme
  labs(
    title = "",
    fill = "%PC",
    x = NULL,
    y = NULL
  ) + # Labels; remove x/y axis titles
  theme(
    axis.text = element_blank(), # Remove axis text/ticks
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )  # Remove grid lines
# Add the legend as a text annotation
print(Plot_all)

# Step 7: Save the plot (optional)
# This saves the plot to a PNG file. Adjust width/height as needed.
ggsave("countrySpecificCov_all_epitopes.png", plot = Plot_all, width = 16, height = 12, dpi = 300)

