# Step 1: Install and load the necessary libraries ####
{
library(umap)
library(ggplot2)
}

# Step 2: Load the MALDI-imaging data into a data frame ####
data <- Final_intensity_matrix_targeted

# Step 3: Preprocess the data if necessary ####

# Step 4: UMAP ####
umap_representation <- umap(data)
umap_coordinates <- umap_representation$layout

# Step 5: Visualize the UMAP results ####
umap_df <- data.frame(UMAP1 = umap_coordinates[,1],
                      UMAP2 = umap_coordinates[,2])

# Step 6. Create Intensity-normalised dataframe with UMAP coordinates and all marker intensities  ####

# normalise intensity values between 0 and 1
norm_intensity <- data.frame(apply(Final_intensity_matrix_targeted, MARGIN = 2, function(x) data.frame(x/max(x))))
colnames(norm_intensity) <- colnames(Final_intensity_matrix_targeted)

# bind umap data to normalised intensity data
norm_umap_intensity_df <- data.frame(do.call(cbind, c(umap_df, norm_intensity)))

# Step 7. Plot ####

#### a. Plot individually ####

# plot individual UMAP with intensity overlay
ggplot(norm_umap_intensity_df) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = VIM), alpha = 0.6) +
  scale_color_viridis_c() +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "UMAP on MALDI-imaging data",
       x = "UMAP Dimension 1",
       y = "UMAP Dimension 2",
       color = "Intensity")


#### b. Make wide dataframe long ####

# original dataframe: rows as pixels, features (UMAPs and channels) as columns
# new dataframe: UMAP1 and 2 as first two columns, channel as column 3, intensity as channel 4.
# structure: for each UMAP coordinate (pixel), there is a row for every channel and a corresponding intensity
long_norm <- tidyr::pivot_longer(norm_umap_intensity_df, cols = -c(UMAP1, UMAP2), names_to = "channel", values_to = "intensity")

#### c. Plot facet umap plot ####
facet_intensity_umap <- ggplot(long_norm, aes(x = UMAP1, y = UMAP2, color = intensity)) +
  geom_point(size = 0.25) +
  scale_color_viridis_c() +
  facet_wrap(~ channel) + # Create a separate plot for each channel
  #theme_minimal() +
  labs(title = "UMAP plots colored by intensity for each channel",
       x = "UMAP1", y = "UMAP2", color = "Intensity")


# Step 8. Save ####
umap_dir = paste(sample_path,"/", "umap",sep = "")
# if directory does not exist, create it
if( !dir.exists(umap_dir)) {
  dir.create(umap_dir)
}
print(umap_dir)

# Save facet plot of UMAP intensities
ggsave(filename = "all_umap.png",
       plot = facet_intensity_umap,
       device = "png",
       path = umap_dir)

# Save UMAP coordinates
write.table(umap_df, file = paste(sample_path,"/", "umap_coordinates.txt", sep = ""), sep="\t", quote = FALSE, row.names = FALSE)






