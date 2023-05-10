#


# 1). Loading packages and functions ####
{cat("Loading of the libraries....  ")
  suppressMessages(library(Cardinal))
  suppressMessages(library(fftw))
  suppressMessages(library(N2R))
  suppressMessages(library(Matrix))
  suppressMessages(library(igraph))
  suppressMessages(library(splus2R))
  suppressMessages(library(EBImage))
  suppressMessages(library(imager))
  suppressMessages(library(FactoMineR))
  suppressMessages(library(hexbin))
  suppressMessages(library(matrixStats))
  suppressMessages(library(mclust))
  suppressMessages(library(dplyr))
  suppressMessages(library(matter))
  cat("... done ! \n")
}

# Load in functions from helper script
source("/mnt/msi_volume/Rscripts/maldi-processing/code/helper_functions.R")

# 2. Read in data
Path_to_imzml_file <- "/mnt/msi_volume/experiments/080523_dilution/13.imzML"
region13 <- readMSIData(Path_to_imzml_file)
region <- "region13"

# 3. Create custom peakAnnotation df ####
Name_vector <- c("uncharged_peptide_plus1", "charged_peptide_plus1", "uncharged_peptide_plus2", "charged_peptide_plus2")
Mass_vector <- as.numeric(c(919.584, 961.628, 460.298, 481.335))

peakAnnotation <- data.frame(cbind(Name_vector, Mass_vector))
rownames(peakAnnotation) <- NULL
peakAnnotation$Mass_vector <- as.numeric(peakAnnotation$Mass_vector)
peakAnnotation = dplyr::arrange(peakAnnotation, peakAnnotation$Mass_vector)

typeof(peakAnnotation$Mass_vector)

# 4. Pre-processing ####

peakPre <- region13 %>%
  Cardinal::normalize(method = "tic") %>%
  smoothSignal(method = "gaussian", plot=FALSE) %>%
  reduceBaseline(method="locmin") %>%
  process(BPPARAM = MulticoreParam(workers = 4))


# 5. Processing to align ####

snr = 2
window_width = 50
refList =  peakAnnotation$Mass_vector
Location_pixels = as.data.frame(pData(peakPre))[,c("x","y")]

peakProc <- peakPre %>%
  peakPick(method="mad", SNR = snr,window = window_width) %>%
  peakAlign(ref = refList, tolerance = 0.2, units = "mz", method="diff") %>%
  process()

# retrieve intensity matter matrix
mat_mat <- peakProc@imageData$data$intensity
# convert to matrix then to dataframe
my_matrix <- t(as.matrix(mat_mat))
my_df <- data.frame(my_matrix)
colnames(my_df) <- peakAnnotation$Name_vector


# 6. Crop location coordinates  ####

x_values <- Location_pixels$x
y_values <- Location_pixels$y

# crop at desired x value position. Edit this value for different crops
x_values <- x_values[x_values > 160]

# crop y_values accordingly
y_values <- y_values[1:length(x_values)]

# remake Location_pixels_crop dataframe
Location_pixels_crop <- data.frame(matrix(NA, nrow = length(x_values), ncol = 0))

Location_pixels_crop$x <- x_values
Location_pixels_crop$y <- y_values
Location_pixels_crop <- cbind(x_values, y_values)
Location_pixels_crop <- data.frame(Location_pixels_crop)


# 7. Replace zeros with NA for entire DF ####
df_na <- my_df
for (i in 1:ncol(my_df)){

  temp_channel <- my_df[,i]
  temp_channel[temp_channel == 0] <- NA
  df_na[,i] <- temp_channel

}


# test plot of uncropped dataframe

Plot_channel_drop = function(image, channel_number=4, quantile_lim = 0.99) {

  Matrix_image = matrix(NA,ncol = max(Location_pixels$y), nrow=max(Location_pixels$x))

  # intensity vector
  x = image[,channel_number]
  # get intensity value for 99th percentile most intense pixels
  x_max = quantile(x, probs = quantile_lim)
  # set all pixel values greater than x_max to that of x_max
  x[x>x_max] = x_max
  x[x == 0] <- NA
  Matrix_image[as.matrix(Location_pixels)] = x
  # convert matrix to image
  Matrix_image = as.cimg(Matrix_image-min(Matrix_image, na.rm = TRUE))

  # add colour channels to the image
  Matrix_image = add.color(Matrix_image,simple = TRUE)

  # set red and blue channels to zero to get only green
  R(Matrix_image) <- 0
  B(Matrix_image) <- 0
  #G(Matrix_image) <- 0

  plot((Matrix_image))

}

Plot_channel(my_df, 1)


# 8. Crop intensity dataframe according to cropped pixel coordinates ####

# initialise dataframe
my_cropped_df <- data.frame(matrix(data = 0, nrow = length(x_values), ncol = 1))
# crop dataframe
for (k in 1:ncol(my_df)){

  # uses df_na, the dataframe with zeros replaced by NAs so the mean isn't weighted
  temp_channel <- my_df[,k]
  temp_channel <- data.frame(temp_channel[1:length(x_values)])
  my_cropped_df[,k] <- temp_channel

}
# set column names
colnames(my_cropped_df) <- colnames(my_df)


# 9. Calculate mean for each column ####
mean_vector <- c()
for (j in 1:ncol(my_df)){

  temp_df <- data.frame(my_df[,j])
  temp_mean <- colMeans((temp_df), na.rm = TRUE)
  mean_vector <- data.frame(rbind(mean_vector, temp_mean))

}
rownames(mean_vector) <- colnames(my_df)
colnames(mean_vector) <- region
mean_df <- data.frame(mean_vector)

# Convert to DF, add column, fill column with mean values for crop2

mean_df["region14_crop2"] <- NA

for (l in 1:ncol(my_cropped_df)){

  temp_df <- data.frame(my_cropped_df[,l])
  mean_df[,2][l] <- colMeans((temp_df), na.rm = TRUE)

}

region14_mean_df<- mean_df

# 10. Plot image ####

Plot_channel_cropped = function(image, channel_number=4, quantile_lim = 0.99) {

  Matrix_image = matrix(NA,ncol = max(Location_pixels_crop$y), nrow=max(Location_pixels_crop$x))

  # intensity vector
  x = image[,channel_number]
  # get intensity value for 99th percentile most intense pixels
  x_max = quantile(x, probs = quantile_lim)
  # set all pixel values greater than x_max to that of x_max
  x[x>x_max] = x_max
  #x[x == 0] <- NA
  Matrix_image[as.matrix(Location_pixels_crop)] = x
  # convert matrix to image
  Matrix_image = as.cimg(Matrix_image-min(Matrix_image, na.rm = TRUE))

  # add colour channels to the image
  Matrix_image = add.color(Matrix_image,simple = TRUE)

  # set red and blue channels to zero to get only green
  R(Matrix_image) <- 0
  B(Matrix_image) <- 0
  #G(Matrix_image) <- 0

  plot((Matrix_image))

}

Plot_channel_cropped(my_cropped_df, 4)

region_dir = paste("/mnt/msi_volume/processed_files/dilution_experiment/", region, sep = "")
print(region_dir)
# if directory does not exist, create it
if( !dir.exists(region_dir)) {
  dir.create(region_dir)
}

hist_dir = paste(region_dir, "/", "histograms", sep = "")
print(hist_dir)
# if directory does not exist, create it
if( !dir.exists(hist_dir)) {
  dir.create(hist_dir)
}

# 11. Histogram ####
for (n in 1:length(peakAnnotation$Name_vector)){

  hist_path = paste(hist_dir,"/", peakAnnotation$Name_vector[n], ".png", sep = "")
  print(hist_path)
  png(file=hist_path, width=1200, height=700)
  hist(log(my_df[,n]), 100, main = paste(region, peakAnnotation$Name_vector[n], sep = "_"), xlab  = "Log Intensity")
  dev.off()

}


Spatial_dir = paste(region_dir, "/", "ion_images",sep = "")
print(Spatial_dir)

# if directory does not exist, create it
if( !dir.exists(Spatial_dir)) {
  dir.create(Spatial_dir)
}

# loop for uncropped images
for (n in 1:ncol(my_df)){
  Spatial_path = paste(Spatial_dir,"/",colnames(my_df)[n], ".png", sep = "")
  print(Spatial_path)
  png(file=Spatial_path, width=1200, height=700)
  Plot_channel(my_df, n)
  dev.off()
}


# for cropped images
for (n in 1:ncol(my_df)){
    Spatial_path = paste(Spatial_dir,"/", region, "/", "_crop2_",colnames(my_df)[n], ".png", sep = "")
    print(Spatial_path)
    png(file=Spatial_path, width=1200, height=700)
    Plot_channel_cropped(my_cropped_df, n)
    dev.off()
}


write.table(mean_df, file = paste("/mnt/msi_volume/processed_files/dilution_experiment/", region, "/", "mean_value.txt", sep = ""), sep="\t", quote = FALSE, row.names = FALSE)

print(paste("/mnt/msi_volume/processed_files/dilution_experiment/", region, "/", "mean_value.txt", sep = ""))

