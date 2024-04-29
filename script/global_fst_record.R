rm(list=ls());gc()
graphics.off()
library(ggplot2)
library(dplyr)
library(reshape2)
v1.1 = read.csv("../data/GlobalAlienSpeciesFirstRecordDatabase_v1.1.csv")
v1.2 = read.csv("../data/GlobalAlienSpeciesFirstRecordDatabase_v1.2.csv")
v2 = read.csv("../data/GlobalAlienSpeciesFirstRecordDatabase_v2.csv")

v3.1 = read.csv("../data/GlobalAlienSpeciesFirstRecordDatabase_v3.1_freedata.csv")


### All species ###
# Assuming all dataframes have a 'FirstRecord' column

# Select only the necessary columns (FirstRecord) and add a version column
nsp_v1.1 <- v1.1 %>% select(FirstRecord) %>% mutate(version = "v1.1_20200228")
nsp_v1.2 <- v1.2 %>% select(FirstRecord) %>% mutate(version = "v1.2_20200228")
nsp_v2 <- v2 %>% select(FirstRecord) %>% mutate(version = "v2_20210322")
nsp_v3.1 <- v3.1 %>% select(FirstRecord) %>% mutate(version = "v3.1_20231025")
# Combine the dataframes
combined_df <- rbind(nsp_v1.1, nsp_v1.2, nsp_v2, nsp_v3.1)

# Filter for years after 1600
combined_df <- combined_df %>% 
  filter(FirstRecord >= 1700)

# Count the number of species per year for each source
species_count_df <- combined_df %>%
  group_by(version, FirstRecord) %>%
  summarise(Count = n(), .groups = 'drop')

# Calculate cumulative count
species_cumulative_df <- species_count_df %>%
  group_by(version) %>%
  mutate(CumulativeCount = cumsum(Count))

# Plot using ggplot2
ggplot(species_count_df, aes(x = FirstRecord, y = Count, color = version)) +
  geom_line() +
  labs(x = "Year", y = "Number of Species", 
       title = "Number of Species First Seen After 1700 by Version") +
  theme_bw() +
  scale_color_brewer(palette = "Set1")  # Use a color palette for clarity
ggsave("../graph/Number of Species First Seen After 1700 by version.png")
# Plot using ggplot2
ggplot(species_cumulative_df, aes(x = FirstRecord, y = CumulativeCount, color = version)) +
  geom_line() +
  labs(x = "Year", y = "Cumulative Number of Species", 
       title = "Cumulative Number of Species Over Time by Version") +
  theme_bw() +
  scale_color_brewer(palette = "Set1") 
ggsave("../graph/Cumulative Number of Species Over Time by Version.png") 


###Now let's examine each lifeform###
#for version 1, let's use v1.2
lifeforms_v1_2 <- unique(v1.2$LifeForm)
lifeforms_v2 <- unique(v2$LifeForm)
lifeforms_v3_1 <- unique(v3.1$LifeForm)

# Find the common LifeForms
common_lifeforms <- Reduce(intersect, list(lifeforms_v1_2, lifeforms_v2, lifeforms_v3_1))

# View the common LifeForms
print(common_lifeforms)

# Function to prepare data for a given LifeForm
prepare_data_for_lifeform <- function(lifeform, dfs, df_versions) {
  # Prepare each dataframe
  prepared_dfs <- lapply(seq_along(dfs), function(i) {
    df <- dfs[[i]] %>%
      filter(LifeForm == lifeform, FirstRecord > 1600) %>%
      group_by(FirstRecord) %>%
      summarise(Count = n(), .groups = 'drop') %>%
      mutate(CumulativeCount = cumsum(Count),
             version = df_versions[i])
    return(df)
  })
  
  # Combine the dataframes
  combined_df <- do.call(rbind, prepared_dfs)
  return(combined_df)
}

# List of dataframes and their versions
dfs <- list(v1.2, v2, v3.1)
df_versions <- c("v1.2", "v2", "v3.1")

# Loop through each common LifeForm and plot
for(lifeform in common_lifeforms) {
  lifeform_data <- prepare_data_for_lifeform(lifeform, dfs, df_versions)
  title = paste("Species Count for", lifeform)
  # Define the maximum for scaling the secondary axis
  max_count <- max(lifeform_data$Count, na.rm = TRUE)
  max_cumulative <- max(lifeform_data$CumulativeCount, na.rm = TRUE)
  scale_factor <- max_cumulative / max_count
  
  # Plotting
  gg <- ggplot(lifeform_data, aes(x = FirstRecord)) +
    geom_point(aes(y = Count, color = version), size = 3) +
    geom_line(aes(y = CumulativeCount / scale_factor, color = version), size = 1.2) +
    scale_y_continuous(name = "Number of Species", 
                       sec.axis = sec_axis(~ . * scale_factor, name = "Cumulative Count")) +
    labs(title = title,
         x = "Year") +
    theme_bw() +
    scale_color_brewer(palette = "Set1")
  ggsave(paste0("../graph/lifeform/", title, ".png"))
  print(gg)
}

### Let's examine each regions###
region_v1_2 <- unique(v1.2$Region)
region_v2 <- unique(v2$Region)
region_v3_1 <- unique(v3.1$Region)

# Find the common regions
common_regions <- Reduce(intersect, list(region_v1_2, region_v2, region_v3_1))

prepare_data_for_region <- function(region, dfs, df_versions) {
  # Prepare each dataframe
  prepared_dfs <- lapply(seq_along(dfs), function(i) {
    df <- dfs[[i]] %>%
      filter(Region == region, FirstRecord > 1600) %>%
      group_by(FirstRecord) %>%
      summarise(Count = n(), .groups = 'drop') %>%
      mutate(CumulativeCount = cumsum(Count),
             version = df_versions[i])
    return(df)
  })
  
  # Combine the dataframes
  combined_df <- do.call(rbind, prepared_dfs)
  return(combined_df)
}

for(region in common_regions) {
  region_data <- prepare_data_for_region(region, dfs, df_versions)
  title = paste("Species Count for", region)
  # Define the maximum for scaling the secondary axis
  max_count <- max(region_data$Count, na.rm = TRUE)
  max_cumulative <- max(region_data$CumulativeCount, na.rm = TRUE)
  scale_factor <- max_cumulative / max_count
  
  # Plotting
  gg <- ggplot(region_data, aes(x = FirstRecord)) +
    geom_point(aes(y = Count, color = version), size = 3) +
    geom_line(aes(y = CumulativeCount / scale_factor, color = version), size = 1.2) +
    scale_y_continuous(name = "Number of Species", 
                       sec.axis = sec_axis(~ . * scale_factor, name = "Cumulative Count")) +
    labs(title = title,
         x = "Year") +
    theme_bw() +
    scale_color_brewer(palette = "Set1")
  ggsave(paste0("../graph/region/", title, ".png"))
  print(gg)
}

#regional species
regional_sp_2016 = readRDS("../data/regional_sp_2016.rds")

rownames(regional_sp_2016)[rownames(regional_sp_2016) == "Australia/New Zealand"] <- "Australia and New Zealand"

regional_sp_df <- as.data.frame(regional_sp_2016, stringsAsFactors = FALSE)
regional_sp_df$Region <- rownames(regional_sp_df)

# Now, use the melt function from reshape2 to convert this to a long format
regional_sp_long <- melt(regional_sp_df, id.vars = "Region", variable.name = "Species", value.name = "Year")

# Remove NA values if there are any
regional_sp_long <- na.omit(regional_sp_long)

# Convert Year to a numeric value if it's not already
regional_sp_long$Year <- as.numeric(as.character(regional_sp_long$Year))
regional_sp_long <- regional_sp_long[regional_sp_long$Year >= 1970,]
# Aggregate the data by year for each region
yearly_counts_2016 <- aggregate(cbind(Count = Year) ~ Region + Year, data = regional_sp_long, FUN = length)

# Plot the trends for each region
for(region_2016 in unique(yearly_counts_2016$Region)) {
  # Subset the data for the current region
  region_data_2016 <- subset(yearly_counts_2016, Region == region)
  title <- paste("Yearly Trend of Species Occurrences in", region)
  # Create the plot
  p <- ggplot(region_data, aes(x = Year, y = Count)) +
    geom_line() +
    labs(title = title,
         x = "Year",
         y = "Number of Species Occurrences") +
    theme_bw()
  
  # Define the title for use in file name
  
  
  # Save the plot as a PNG file
  ggsave(filename = paste0("../graph/ecoregion_2016/", title, ".png"), plot = p, width = 10, height = 6)
}

regional_sp_2023 = readRDS("../data/regional_sp_2023.rds")
regional_sp_df <- as.data.frame(regional_sp_2023, stringsAsFactors = FALSE)
regional_sp_df$Region <- rownames(regional_sp_df)

# Now, use the melt function from reshape2 to convert this to a long format
regional_sp_long <- melt(regional_sp_df, id.vars = "Region", variable.name = "Species", value.name = "Year")

# Remove NA values if there are any
regional_sp_long <- na.omit(regional_sp_long)

# Convert Year to a numeric value if it's not already
regional_sp_long$Year <- as.numeric(as.character(regional_sp_long$Year))
regional_sp_long <- regional_sp_long[regional_sp_long$Year >= 1970,]
# Aggregate the data by year for each region
yearly_counts_2023 <- aggregate(cbind(Count = Year) ~ Region + Year, data = regional_sp_long, FUN = length)

# Plot the trends for each region
for(region in unique(yearly_counts_2023$Region)) {
  # Subset the data for the current region
  region_data_2023 <- subset(yearly_counts_2023, Region == region)
  title <- paste("Yearly Trend of Species Occurrences in", region)
  # Create the plot
  p <- ggplot(region_data_2023, aes(x = Year, y = Count)) +
    geom_line() +
    labs(title = title,
         x = "Year",
         y = "Number of Species Occurrences") +
    theme_bw()
  
  # Define the title for use in file name
  
  
  # Save the plot as a PNG file
  ggsave(filename = paste0("../graph/ecoregion_2023/", title, ".png"), plot = p, width = 10, height = 6)
}

yearly_counts_2016$DatasetYear <- '2016'
yearly_counts_2023$DatasetYear <- '2023'


# Combine the datasets
combined_data <- rbind(yearly_counts_2016, yearly_counts_2023)

# Plot the trends for each region
unique_regions <- union(unique(yearly_counts_2016$Region), unique(yearly_counts_2023$Region))
for(region in unique_regions) {
  # Subset the data for the current region from combined data
  region_data_combined <- subset(combined_data, Region == region)
  
  # Create the plot with two lines
  p <- ggplot(region_data_combined, aes(x = Year, y = Count, color = DatasetYear)) +
    geom_line() +
    labs(title = paste("Yearly Trend of Species Occurrences in", region),
         x = "Year",
         y = "Number of Species Occurrences") +
    theme_bw() +
    scale_color_manual(values = c("2016" = "blue", "2023" = "red"))
  
  # Save the plot as a PNG file
  file_name <- paste0("../graph/ecoregion_combined/", gsub(" ", "_", gsub("/", "_", region)), ".png")
  ggsave(filename = file_name, plot = p, width = 10, height = 6)
}

combined_data$Year <- as.numeric(combined_data$Year)

# Plot the trends for all regions combined using facet_wrap
p <- ggplot(combined_data, aes(x = Year, y = Count, color = DatasetYear)) +
  geom_line() +
  facet_wrap(~Region, scales = "free_y") + # Use free_y if the counts vary significantly between regions
  geom_vline(xintercept = c(1996, 2002), linetype = "dashed", color = "black") + # Add vertical lines for the years 1996 and 2002
  labs(title = "Yearly Trend of Species Occurrences by Region",
       x = "Year",
       y = "Number of Species Occurrences") +
  theme_bw() +
  scale_color_manual(values = c("2016" = "blue", "2023" = "red"))
ggsave(filename = "../graph/ecoregion_combined/Yearly Trend of Species Occurrences by Region.png", plot = p, width = 10, height = 6)



