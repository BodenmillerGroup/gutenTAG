## Cytomapper testing

# Install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("cytomapper")

# Load packages
library(cytomapper)
library(imager)

# Test 1: initial test ####

# List of matrices for each channel of Final intensity matrix targetted
test <- lapply(Final_intensity_matrix_targeted, function(x) matrix(x, nrow = range(Location_data$x)[2], ncol = range(Location_data$y)[2], byrow = TRUE))
# Convert list to 3 dimensional array
test1 <- do.call(abind, list(test, along = 3))

# Convert array to image
test2 <- apply(test1, MARGIN = c(1,2), Image)

# Create CytoImageList object
CytoImageList(image1 = test2)
CytoImageList(test2)



# Test 2: Create list of one channel images and make CytoImageList ####


# List of matrices
matrix_list <- lapply(Final_intensity_matrix_targeted, function(x) matrix(x, nrow = range(Location_data$x)[2], ncol = range(Location_data$y)[2], byrow = TRUE))
# List of Images
image_list <- lapply(matrix_list, Image)
# Create CytoImageList
kidneyImage <- CytoImageList(image_list)





# Test 3:  Create a 33 channel image ####

# adjust location data
Location_data$x = Location_data$x - min(Location_data$x) + 1
Location_data$y = Location_data$y - min(Location_data$y) + 1

# List of matrices for each channel of Final intensity matrix targetted
matrix_list <- lapply(Final_intensity_matrix_targeted, function(x) matrix(x, nrow = max(Location_data$x), ncol = max(Location_data$y), byrow = TRUE))
# Convert list to 3 dimensional array
my_array <- do.call(abind, list(matrix_list, along = 3))

# Convert array to n channel image
my_image <- Image(my_array)
# Create CytoImageList from Image
kidney_image <- CytoImageList(my_image)


# Plot pixels
imageData(kidney_image)[,,1]
plotPixels(kidney_image, colour_by = "HLA-ABC")



# Test 4: use old method of making matrix
# New method

curateMatrix <- function(dataframe, channel, Location_data){

  Matrix_image = matrix(0, ncol = max(Location_data$y), nrow=max(Location_data$x))
  x = dataframe[,channel]
  #x_max = quantile(x, probs = quantile_lim)
  #x[x>x_max] = x_max
  Matrix_image[as.matrix(Location_data)] = x

  return(Matrix_image)
}



matrix_list <- lapply(Final_intensity_matrix_targeted, curateMatrix)
my_array <- do.call(abind, list(matrix_list, along = 3))
my_image <- Image(my_array)
kidney_image <- CytoImageList(my_image)

plotPixels(kidney_image, colour_by = c("Caveolin-1", "HLA-ABC"))



# ChatGPT

curateMatrix_GPT <- function(channel, Location_data, dataframe) {
  Matrix_image = matrix(0, ncol = max(Location_data$y), nrow = max(Location_data$x))
  x = dataframe[,channel]
  Matrix_image[cbind(Location_data$x, Location_data$y)] = x
  return(Matrix_image)
}


curateMatrix <- function(dataframe, channel, Location_data){
  Matrix_image = matrix(0, ncol = max(Location_data$y), nrow=max(Location_data$x))
  x = dataframe[,channel]
  Matrix_image[as.matrix(Location_data)] = x
  return(Matrix_image)
}


matrix_list <- lapply(colnames(Final_intensity_matrix_targeted), curateMatrix, Location_data = Location_data, dataframe = Final_intensity_matrix_targeted)
my_array <- do.call(abind, list(matrix_list, along = 3))
my_image <- Image(my_array)
kidney_image <- CytoImageList(my_image)
channelNames(kidney_image) <- colnames(Final_intensity_matrix_targeted)

plotPixels(kidney_image,
           colour_by = c("Collagen-1A1", "HLA-ABC", "Actin-αSM-(SMA)"))

