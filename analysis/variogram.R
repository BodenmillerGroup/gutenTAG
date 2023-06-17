library(spatial)
library(gstat)
library(spatstat)
library(sp)
library(ggplot2)

# 1. Construct Variogram function ####
constructVariogram <- function(df = Final_intensity_matrix_targeted, coords = Location_pixels){

  # Include make.names to remove problematic characters, such as spaces, +, etc.
  names(df) <- make.names(names(df))
  names <- colnames(df)

  # Create list of variograms for each marker
  variograms <- list()
  for (i in 1:ncol(df)){

    marker <- data.frame(df[i])
    colnames(marker) <- names[i]

    # Assign coordinates
    marker$x <- coords$x
    marker$y <- coords$y
    coordinates(marker) <- ~x+y

    # Construct a formula from the current column name
    formula_string <- paste0(names[i], "~1")
    formula_object <- as.formula(formula_string)

    # Calculate variogram for the marker
    marker_variogram <- variogram(formula_object, marker)

    # Store the result in the list
    variograms[[names[i]]] <- marker_variogram

  }

  return(variograms)

}



# Apply function to Final Intensity Matrix
vgrams <- constructVariogram(df = Final_intensity_matrix_targeted, coords = Location_pixels)

names <- colnames(Final_intensity_matrix_targeted)

# Specify an initial model to which you fit the sample variogram
initial_model <- vgm(psill = 100, model = "Exp", range = 25, nugget = 10)

# Fit a variogram model to each variogram in the list
fitted_models <- lapply(vgrams, function(vgram) fit.variogram(vgram, initial_model))


# 2. Compute predicted semivariances for each fitted model ####

# Initialise list
predicted_values <- list()
for (k in 1:length(fitted_models)){

  # Compute predicted semivariance values of fitted model at each distance
  predicted_vals <- variogramLine(object = fitted_models[[k]], maxdist = max(vgrams[[k]]$dist))

  # Add each set of predicted values to list of predicted values
  predicted_values[[names[k]]] <- predicted_vals

}




# 3. Create dataframe for plotting all variograms as a facet plot ####

# First for the empirical variogram

# Add index name to the gamma column in each dataframe in the list
for(i in seq_along(vgrams)) {
  names(vgrams[[i]])[3] <- paste0("gamma_", i)
}

# Remove unnecessary columns
short_vgrams <-  list()
for (i in 1:length(vgrams)){

  # obtain distance and gamma column for marker
  short_vgram <- vgrams[[i]][c(2,3)]
  # add to list
  short_vgrams[[names[i]]] <- short_vgram

}

# Merge to dataframe by dist
vgram_df <- Reduce(function(df1, df2) merge(df1, df2, by = "dist"), short_vgrams)
distances <- vgram_df["dist"]
vgram_df <- vgram_df[-1]
colnames(vgram_df) <- names
vgram_df <- bind_cols(distances, vgram_df)

# Then for the fitted values

# Add a name to the gamma column in each dataframe in the list
for(i in seq_along(predicted_values)) {
  names(predicted_values[[i]])[2] <- paste0("gamma_", i)
}

# Combine all dataframes in the list into a single dataframe
gamma_df <- Reduce(function(df1, df2) merge(df1, df2, by = "dist"), predicted_values)

distances <- gamma_df["dist"]
gamma_df <- gamma_df[-1]
colnames(gamma_df) <- names
gamma_df <- bind_cols(distances, gamma_df)



#ggplot(data = vgram_df) +
#  geom_point(aes(x = dist, y = PD1), color = "red") +
#  geom_line(data = gamma_df, aes(x = dist, y = PD1), color = "darkblue") +
#  ggtitle(paste("Variogram")) +
#  ylab("Semivariance") +
#  xlab("Distance")


# Reshape data to long format for ggplot2
vgram_long <- tidyr::pivot_longer(vgram_df, cols = -dist, names_to = "feature", values_to = "gamma")
gamma_long <- tidyr::pivot_longer(gamma_df, cols = -dist, names_to = "feature", values_to = "gamma_pred")

# Combine into a single dataframe
plot_data <- dplyr::full_join(vgram_long, gamma_long, by = c("dist", "feature"))

# 4. Plot single sample variogram and fitted model.  ####

# Need to specify which marker you want to plot by indexing
ggplot(data = vgrams[[29]]) +
  geom_point(aes(x = dist, y = gamma_29), color = "red") +
  geom_line(data = predicted_values[[29]], aes(x = dist, y = gamma_29), color = "darkblue") +
  ggtitle(paste("Variogram of", names(vgrams[3]))) +
  ylab("Semivariance") +
  xlab("Distance")


# 5. Create facet plot for each variogram ####
gg_variogram <- ggplot(data = plot_data) +
  geom_point(aes(x = dist, y = gamma), color = "red") +
  geom_line(aes(x = dist, y = gamma_pred), color = "darkblue") +
  facet_wrap(~ feature, scales = "free_y") +
  ggtitle("Variograms") +
  ylab("Semivariance") +
  xlab("Distance")


# 6. Save ####
variogram_dir = paste(sample_path,"/", "variograms",sep = "")
# if directory does not exist, create it
if( !dir.exists(variogram_dir)) {
  dir.create(variogram_dir)
}
print(variogram_dir)

# Save facet plot
ggsave(filename = "variograms.png",
       plot = gg_variogram,
       device = "png",
       path = variogram_dir)

# save table
write.table(vgram_df, file = paste(variogram_dir,"/", "variogram_df.txt", sep = ""), sep="\t", quote = FALSE, row.names = FALSE)


