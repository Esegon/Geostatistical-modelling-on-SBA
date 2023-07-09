library(sf)
library(ggplot2)
# Read the shapefile
main_datas <- st_read("C:/Users/Esegon/Desktop/Geostatistic/merged_data/main_dataset.shp")

# Rename the columns
col_names <- c("education", "head_of_family", "wealth_index", "residence",
               "woman_last_birth_age", "work_status", "husband_education_level", "birth_order")

col_names_old <- c("hghs___", "s______", "wlth_nd", "residnc", "Ag_____", "wrk_stt", "hsbnd__", "brth_br")

for (i in seq_along(col_names_old)) {
  col_idx <- which(names(main_datas) == col_names_old[i])
  if (length(col_idx) > 0) {
    names(main_datas)[col_idx] <- col_names[i]
  }
}

# Check for NA values in the specified columns
na_check <- rowSums(is.na(main_datas[, col_names]))

# Create the SBA column
main_datas$SBA <- ifelse(na_check == 0, 1, 2)

# Save the modified shapefile
# st_write(main_datas, "C:/Users/Esegon/Desktop/Geostatistic/merged_data/main_dataset_with_SBAa1.shp")

#main_dat <- st_read("C:/Users/Esegon/Desktop/Geostatistic/merged_data/main_dataset_with_SBAa1.shp")

#View(main_datas)
#plot(main_datas)


replace_na_with_zero <- function(data) {
  # Check for NA values in the specified columns and replace with 0
  columns <- c("education", "head_of_family", "wealth_index", "residence", "woman_last_birth_age",
               "work_status", "husband_education_level", "birth_order","trvl_tm")
  
  for (col in columns) {
    data[[col]][is.na(data[[col]])] <- 0  # Replace NA values with 0
  }
  
  return(data)
}

# Usage:
updated_data <- replace_na_with_zero(main_datas)
View(updated_data)
# --------------------------------------------------------------------------------------
# check for NA Values 
# Assuming your data is stored in a dataframe called "updated_data"
columns <- c("education", "head_of_family", "wealth_index", "residence",
             "woman_last_birth_age", "work_status", "husband_education_level", "birth_order", "trvl_tm")

# Check for NA values in the specified columns
na_counts <- colSums(is.na(updated_data[columns]))

# Print the NA counts for each column
for (i in seq_along(na_counts)) {
  cat("Column", names(na_counts)[i], "has", na_counts[i], "NA values.\n")
}
library(sf)

file_path <- "C:/Users/Esegon/Desktop/Geostatistic/night/ntl_dataf.shp"
# Read the shapefile
modi <- st_read(file_path)
# Add a common ID column to updated_data
updated_data$ID <- 1:nrow(updated_data)

# Add a common ID column to modi
modi$ID <- 1:nrow(modi)

# Perform the merge based on common ID attribute
merged_dataset <- cbind(updated_data, modi[match(updated_data$ID, modi$ID), ])

# View the merged dataset
View(merged_dataset)

# Alternatively, you can print the structure of the merged dataset
str(merged_dataset)
library(dplyr)

# Rename the column in the merged_dataset and assign it to a new variable
final_merge <- merged_dataset %>% 
  rename(nighttime = geometry.1)

View(final_merge)

# ----------------------------------------------------------------------------------------

# NTL Spatial Binarization 

library(sf)

# Define the bin size
bin_size <- 25

# Create a new column 'NTL' in the dataset
final_merge$NTL <- 0

# Iterate through the nighttime column
for (i in 1:length(final_merge$nighttime)) {
  # Get the spatial object from the nighttime data point
  sfc_point <- final_merge$nighttime[[i]]
  
  # Check if sfc_point is not NULL
  if (!is.null(sfc_point)) {
    # Extract the x and y coordinates from the sfc_POINT object
    coords <- st_coordinates(sfc_point)
    x <- coords[1]
    y <- coords[2]
    
    # Check for missing values or NaN in the coordinates
    if (!is.na(x) && !is.na(y) && !is.nan(x) && !is.nan(y)) {
      # Determine the bin indices for the spatial coordinates
      bin_x <- floor(x / bin_size)
      bin_y <- floor(y / bin_size)
      
      # Assign the binary value based on the presence within the bin
      if (bin_x >= 0 && bin_y >= 0) {
        final_merge$NTL[i] <- 1
      }
    }
  }
}

# Display the updated dataset with the 'NTL' column
View(final_merge)

# -----------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------
library(spdep)
library(spatialreg)

# Step 1: Subset the "final_merge" data using the desired columns and limit it to the desired number of rows
columns <- c("education", "head_of_family", "wealth_index", "residence", "woman_last_birth_age", "work_status", "husband_education_level", "birth_order", "trvl_tm", "NTL", "SBA")
subset_data <- final_merge[1:1800, columns]

# Step 2: Adjust neighbor search parameters
k <- 20  # Adjust the value of k for the number of nearest neighbors

# Create the neighbor list using knn2nb with the adjusted k value
nb <- knn2nb(knearneigh(subset_data$geometry, k = k))

# Step 3: Create spatial weights matrix
w <- nb2listw(nb, style = "W")

# Step 4: Perform spatial autoregressive modeling
sar_model <- lagsarlm(SBA ~ education + head_of_family + wealth_index + residence + woman_last_birth_age + work_status + husband_education_level + birth_order + trvl_tm + NTL, data = subset_data, listw = w)

# Step 5: Inspect the results of the spatial autoregressive model
summary(sar_model)

# Step 6: Calculate the accuracy metrics
predicted <- residuals(sar_model) + fitted(sar_model)
actual <- subset_data$SBA
# Calculate the total sum of squares (TSS)
tss <- sum((subset_data$SBA - mean(subset_data$SBA))^2)

# Calculate the residual sum of squares (RSS)
rss <- sum((subset_data$SBA - predicted)^2)

# Calculate the R-squared value
r_squared <- 1 - (rss / tss)

rmse <- sqrt(mean((predicted - actual)^2))
mae <- mean(abs(predicted - actual))
coefficients <- coef(sar_model)


# Print the R-squared value
cat("R-squared (R²):", r_squared, "\n")
cat("Root Mean Squared Error (RMSE):", rmse, "\n")
cat("Mean Absolute Error (MAE):", mae, "\n")
cat("Coefficient Matrix:\n")
print(coefficients)

# Calculate AIC and BIC
aic <- AIC(sar_model)
bic <- BIC(sar_model)

# Print AIC and BIC
cat("AIC:", aic, "\n")
cat("BIC:", bic, "\n")

# ------------------------------------------------------------------------------------------------------------------

# 

# Generate predictions using the SAR model
predicted <- predict(sar_model, newdata = subset_data, type = "response", listw = w)

# Print the predictions
print(predicted)

# -------------------------------------------------------------------------------------------------------------------

# map of predicted SBA rates
library(ggplot2)
library(sf)

# Read the county administrative boundary shapefile
county_shapefile <- st_read("C:/Users/Esegon/Desktop/Geostatistic/kenyan-counties/County.shp")

# Convert the subset_data to an sf object if it's not already in sf format
if (!inherits(subset_data, "sf")) {
  subset_data_sf <- st_as_sf(subset_data, coords = c("longitude", "latitude"), crs = st_crs(4326))
} else {
  subset_data_sf <- subset_data
}

# Add the predicted values as a new column in the sf object
subset_data_sf$predicted <- predicted

# Plot the map with the county boundaries and predicted values
ggplot() +
  geom_sf(data = county_shapefile, fill = "transparent", color = "black") +
  geom_sf(data = subset_data_sf, aes(fill = predicted, color = predicted)) +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(fill = "Predicted SBA Rate", color = "Predicted SBA Rate") +
  theme_minimal()

# -------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(sf)

# Load the county administrative boundary shapefile
counties <- st_read("C:/Users/Esegon/Desktop/Geostatistic/kenyan-counties/County.shp")

# Convert the subset_data to an sf object if it's not already in sf format
if (!inherits(subset_data, "sf")) {
  subset_data_sf <- st_as_sf(subset_data, coords = c("longitude", "latitude"), crs = st_crs(4326))
} else {
  subset_data_sf <- subset_data
}

# Add the predicted values as a new column in the sf object
subset_data_sf$predicted <- predicted

# Perform spatial join to associate predicted values with county boundaries
joined_data <- st_join(counties, subset_data_sf)

# Plot the heat map
ggplot() +
  geom_sf(data = joined_data, aes(fill = predicted), color = "black") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(fill = "Predicted SBA Rate") +
  theme_minimal() +
  ggtitle("Heatmap: Predicted SBA Rate") +
  theme(plot.title = element_text(hjust = 0.5))

# -----------------------------------------------------------------------------------------------
library(ggplot2)
library(gridExtra)
library(cowplot)

# Create bar plots for each variable with different colors
plot_education <- ggplot(subset_data, aes(x = education, y = predicted, fill = education)) +
  geom_bar(stat = "identity") +
  labs(x = "Education", y = "SBA Rate")

plot_head_of_family <- ggplot(subset_data, aes(x = head_of_family, y = predicted, fill = head_of_family)) +
  geom_bar(stat = "identity") +
  labs(x = "Head of Family", y = "SBA Rate")

plot_wealth_index <- ggplot(subset_data, aes(x = wealth_index, y = predicted, fill = wealth_index)) +
  geom_bar(stat = "identity") +
  labs(x = "Wealth Index", y = "SBA Rate")

plot_residence <- ggplot(subset_data, aes(x = residence, y = predicted, fill = residence)) +
  geom_bar(stat = "identity") +
  labs(x = "Residence", y = "SBA Rate")

plot_woman_last_birth_age <- ggplot(subset_data, aes(x = woman_last_birth_age, y = predicted, fill = woman_last_birth_age)) +
  geom_bar(stat = "identity") +
  labs(x = "Woman Last Birth Age", y = "SBA Rate")

plot_work_status <- ggplot(subset_data, aes(x = work_status, y = predicted, fill = work_status)) +
  geom_bar(stat = "identity") +
  labs(x = "Work Status", y = "SBA Rate")

plot_husband_education_level <- ggplot(subset_data, aes(x = husband_education_level, y = predicted, fill = husband_education_level)) +
  geom_bar(stat = "identity") +
  labs(x = "Husband Education Level", y = "SBA Rate")

plot_birth_order <- ggplot(subset_data, aes(x = birth_order, y = predicted, fill = birth_order)) +
  geom_bar(stat = "identity") +
  labs(x = "Birth Order", y = "SBA Rate")

plot_trvl_tm <- ggplot(subset_data, aes(x = trvl_tm, y = predicted, fill = trvl_tm)) +
  geom_bar(stat = "identity") +
  labs(x = "Travel Time", y = "SBA Rate")

plot_NTL <- ggplot(subset_data, aes(x = NTL, y = predicted, fill = NTL)) +
  geom_bar(stat = "identity") +
  labs(x = "NTL", y = "SBA Rate")

# Arrange the plots in a grid
grid_arrange <- plot_grid(plot_education, plot_head_of_family, plot_wealth_index, plot_residence,
                          plot_woman_last_birth_age, plot_work_status, plot_husband_education_level,
                          plot_birth_order, plot_trvl_tm, plot_NTL,
                          ncol = 2, labels = c("Education", "Head of Family", "Wealth Index", "Residence",
                                               "Woman Last Birth Age", "Work Status", "Husband Education Level",
                                               "Birth Order", "Travel Time", "NTL"),
                          label_size = 8)

# Set the title and center it
plot_title <- ggdraw() +
  draw_label("Predicted SBA Rate by Variable",
             fontface = "bold", size = 14, hjust = 0.5)

# Combine the title and grid of plots
final_plot <- plot_grid(plot_title, grid_arrange, nrow = 2, rel_heights = c(0.1, 0.9))

# Display the final plot
final_plot

# ---------------------------------------------------------------------------------
library(ggplot2)
library(sf)

# Load the county administrative boundary shapefile
counties <- st_read("C:/Users/Esegon/Desktop/Geostatistic/kenyan-counties/County.shp")

# Convert the subset_data to an sf object if it's not already in sf format
if (!inherits(subset_data, "sf")) {
  subset_data_sf <- st_as_sf(subset_data, coords = c("longitude", "latitude"), crs = st_crs(4326))
} else {
  subset_data_sf <- subset_data
}

# Add the predicted values as a new column in the sf object
subset_data_sf$predicted <- (10 - predicted)

# Perform spatial join to associate predicted values with county boundaries
joined_data <- st_join(counties, subset_data_sf)

# Plot the heat map
ggplot() +
  geom_sf(data = joined_data, aes(fill = predicted), color = NA) +
  scale_fill_gradientn(
    colours = c("blue", "cyan", "green", "yellow", "red"),
    name = "Predicted SBA Rate (%)",
    labels = scales::percent_format()
  ) +
  labs(title = "Spatial Distribution of Predicted SBA Rate (Heat Map)",
       caption = "Source: Your Data Source") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) 
# ------------------------------------------------------------------------------------

library(ggplot2)
library(sf)
library(scales)

# Load the county administrative boundary shapefile
counties <- st_read("C:/Users/Esegon/Desktop/Geostatistic/kenyan-counties/County.shp")

# Convert the subset_data to an sf object if it's not already in sf format
if (!inherits(subset_data, "sf")) {
  subset_data_sf <- st_as_sf(subset_data, coords = c("longitude", "latitude"), crs = st_crs(4326))
} else {
  subset_data_sf <- subset_data
}

# Add the predicted values as a new column in the sf object
subset_data_sf$predicted <- predicted

# Perform spatial join to associate predicted values with county boundaries
joined_data <- st_join(counties, subset_data_sf)

# Plot the point-heat map
ggplot() +
  geom_sf(data = joined_data, aes(fill = predicted), color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = median(joined_data$predicted),
                       labels = percent_format(scale = 100)) +
  labs(title = "Predicted SBA Rate",
       fill = "Predicted SBA Rate (%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
# -------------------------------------------------------------------------------------------------

library(ggplot2)
library(sf)

# Load the county administrative boundary shapefile
counties <- st_read("C:/Users/Esegon/Desktop/Geostatistic/kenyan-counties/County.shp")

# Convert the subset_data to an sf object if it's not already in sf format
if (!inherits(subset_data, "sf")) {
  subset_data_sf <- st_as_sf(subset_data, coords = c("geometry"), crs = st_crs(4326))
} else {
  subset_data_sf <- subset_data
}

# Add the predicted values as a new column in the sf object
subset_data_sf$predicted <- predicted

# Convert predicted values to percentages
subset_data_sf$predicted_percentage <- subset_data_sf$predicted * 100

# Create a heatmap plot
heatmap_plot <- ggplot() +
  geom_sf(data = counties, fill = "white", color = "black") +
  geom_sf(data = subset_data_sf, aes(fill = predicted_percentage), size = 3) +
  scale_fill_gradient(low = "blue", high = "red", labels = scales::percent_format()) +
  labs(fill = "Predicted SBA Rate (%)", title = "Heatmap for Predicted SBA Rates") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Display the heatmap plot
print(heatmap_plot)

# -----------------------------------------------------------------------------------------
library(ggplot2)
library(sf)
library(gdistance)

# Read the county administrative boundary shapefile
county_shapefile <- st_read("C:/Users/Esegon/Desktop/Geostatistic/kenyan-counties/County.shp")

# Convert the subset_data to an sf object if it's not already in sf format
if (!inherits(subset_data, "sf")) {
  subset_data_sf <- st_as_sf(subset_data, coords = c("longitude", "latitude"), crs = st_crs(4326))
} else {
  subset_data_sf <- subset_data
}

# Add the predicted values as a new column in the sf object
subset_data_sf$predicted <- predicted

# Calculate the distance matrix
dist_matrix <- st_distance(subset_data_sf)

# Get the indices of 5 nearest neighbors for each point
k <- 5
nearest_neighbors <- apply(dist_matrix, 1, function(row) order(row)[2:(k+1)])

# Create a subset of the data with only the 5 nearest neighbors
subset_data_neighbors <- subset_data_sf[unlist(nearest_neighbors), ]

# Plot the map with the county boundaries and predicted values
ggplot() +
  geom_sf(data = county_shapefile, fill = "transparent", color = "black") +
  geom_sf(data = subset_data_sf, aes(fill = predicted, color = predicted)) +
  geom_sf(data = subset_data_neighbors, aes(fill = predicted), color = "black", alpha = 0.6) +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_color_gradient(low = "blue", high = "red") +
  labs(fill = "Predicted SBA Rate", color = "Predicted SBA Rate") +
  theme_minimal()

# ---------------------------------------------------------------------------------

library(ggplot2)
library(sf)

# Load the county administrative boundary shapefile
counties <- st_read("C:/Users/Esegon/Desktop/Geostatistic/kenyan-counties/County.shp")

# Convert the subset_data to an sf object if it's not already in sf format
if (!inherits(subset_data, "sf")) {
  subset_data_sf <- st_as_sf(subset_data, coords = c("longitude", "latitude"), crs = st_crs(4326))
} else {
  subset_data_sf <- subset_data
}

# Add the predicted values as a new column in the sf object
subset_data_sf$predicted <- predicted  # Subtract predicted values from 100 to reverse the scaling
#View(subset_data_sf$predicted)
# Perform spatial join to associate predicted values with county boundaries
joined_data <- st_join(counties, subset_data_sf)

# Plot the heat map
ggplot() +
  geom_sf(data = joined_data, aes(fill = predicted), color = NA) +
  scale_fill_gradientn(
    colours = c("red", "yellow", "green", "cyan", "blue"),
    name = "Predicted SBA Rate (%)",
    labels = scales::percent_format(),
    limits = (c(0, 2))/2
  ) +
  labs(title = "Spatial Distribution of Predicted SBA Rate (Heat Map)",
       caption = "Source: Your Data Source") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
# ---------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(sf)

# Load the county administrative boundary shapefile
counties <- st_read("C:/Users/Esegon/Desktop/Geostatistic/kenyan-counties/County.shp")

# Convert the subset_data to an sf object if it's not already in sf format
if (!inherits(subset_data, "sf")) {
  subset_data_sf <- st_as_sf(subset_data, coords = c("longitude", "latitude"), crs = st_crs(4326))
} else {
  subset_data_sf <- subset_data
}

# Add the predicted values as a new column in the sf object
subset_data_sf$predicted <- (10 - predicted)

# Perform spatial join to associate predicted values with county boundaries
joined_data <- st_join(counties, subset_data_sf)

# Specify the original column names
new_columns <- c("education", "head_of_family", "wealth_index", "residence",
                      "woman_last_birth_age", "work_status", "husband_education_level", "birth_order")

# Specify the new column names
#new_columns <- c("educated", "head_of_family", "rich", "urban_rural",
                 #"woman_last_birth_age", "work_status", "husband_education_level", "birth_order")

# Create a separate heat map for each column and store it in the heatmaps list
heatmaps <- list()

for (i in seq_along(original_columns)) {
  column <- original_columns[i]
  new_column <- new_columns[i]
  
  heat_map <- ggplot() +
    geom_sf(data = joined_data, aes(fill = .data[[new_column]]), color = NA) +
    scale_fill_gradientn(
      colours = c("blue", "cyan", "green", "yellow", "red"),
      name = paste("Predicted", new_column),
      labels = scales::percent_format()
    ) +
    labs(title = paste("Spatial Distribution of Predicted", new_column, "(Heat Map)"),
         caption = "Source: Your Data Source") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  heatmaps[[new_column]] <- heat_map
  
  # Save the heat map as a shapefile
  #output_file <- paste0("C:/Users/Esegon/Desktop/Geostatistic/practicals/", new_column, "_heatmap.shp")
  #st_write(joined_data, output_file)
}

# Save the heat maps as images
for (i in seq_along(original_columns)) {
  column <- original_columns[i]
  new_column <- new_columns[i]
  
  output_file <- paste0("C:/Users/Esegon/Desktop/Geostatistic/practicals/", new_column, "_heatmap.png")
  ggsave(output_file, heatmaps[[new_column]], width = 8, height = 6)
}

# ----------------------------------------------------------------------------------------
 


