# MALDI-imaging integration with Shiny ####

# Load packages
library(cytomapper)
library(imager)

# 2. Create CytoImageList object ####

# Convert flattened dataframe column to intensity matrix
curateMatrix <- function(dataframe, channel, Location_data){
  Matrix_image = matrix(0, ncol = max(Location_data$y), nrow=max(Location_data$x))
  x = dataframe[,channel]
  Matrix_image[as.matrix(Location_data)] = x
  return(Matrix_image)
}

# Create list of matrices from intensity dataframe columns
matrix_list <- lapply(colnames(Final_intensity_matrix_targeted), curateMatrix, Location_data = Location_data, dataframe = Final_intensity_matrix_targeted)
# Create multidimensional array
my_array <- do.call(abind, list(matrix_list, along = 3))
# Create multichannel image from array
my_image <- Image(my_array)
# Create CytoImageList object from image
my_image <- CytoImageList(my_image)
# Assign channel names
channelNames(my_image) <- colnames(Final_intensity_matrix_targeted)

# 3. Plot pixels ####

# Plot channels

plotPixels(my_image,
           colour_by = c("CD68", "PanCK"))


plotPixels(my_image, colour_by = c("PR.AB", "CK5", "Cytokeratin7"),
           colour = list(PR.AB = c("black", "red"), CK5  = c("black", "cyan"), Cytokeratin7 = c("black", "green")))


