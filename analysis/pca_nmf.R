#### PCA ####

# Load packages
library(irlba)
library(RcppML)
library(EBImage)
library(imager)

# Vector for plotting images
Plot_vector = function(x, quantile_lim = 0.99) {
  Matrix_image = matrix(NA,ncol = max(Location_pixels$y),nrow=max(Location_pixels$x))

  # get intensity value for 99th percentile most intense pixels
  x_max = quantile(x, probs = quantile_lim)
  # set all pixel values greater than x_max to that of x_max
  x[x>x_max] = x_max
  Matrix_image[as.matrix(Location_pixels)] = x
  # convert matrix to image
  Matrix_image = as.cimg(Matrix_image-min(Matrix_image,na.rm = T))

  # add colour channels to the image
  Matrix_image = add.color(Matrix_image,simple = TRUE)

  # set red and blue channels to zero to get only green
  R(Matrix_image) <- 0
  B(Matrix_image) <- 0


  plot((Matrix_image))#,main=colnames(Final_intensity_matrix[channel_number]))
  text(x = max(Location_pixels$x)*1.02, y = max(Location_pixels$y)*0.1,
       adj = 0,
       cex = 1,
       col = "green")
}

# check dimensions of
dim(Final_intensity_matrix_targeted)

# variance and mean of columns (each marker)
var_temp = apply(Final_intensity_matrix_targeted,MARGIN = 2,FUN = var)
mean_temp = colMeans(Final_intensity_matrix_targeted)

# check linearity
plot(mean_temp,var_temp,log="xy")

# compute PCA (adjust n to set number of principal components)
PCA_simple = prcomp_irlba(Final_intensity_matrix_targeted, n = 6)

# summary of Principal Componants

summary(PCA_simple)

# get contribution of each feature to each principal component
Contribution = PCA_simple$rotation
rownames(Contribution) = colnames(Final_intensity_matrix_targeted)
View(Contribution)

# scree plot
par(las=1,bty="l")
barplot((PCA_simple$sdev)^2/PCA_simple$totalvar*100, ylab="Proportion of variance explained")

# plot principal components
Plot_vector(PCA_simple$x[,1])

# plot the inverse of the principal components
Plot_vector(-PCA_simple$x[,1])


#### NMF ####

# run NMF
nmf_temp = nmf(scale(Final_intensity_matrix_targeted,center = FALSE),k = 10)

# scree plot
barplot(nmf_temp$d)
dim(nmf_temp$w)

# plot factors
Plot_vector(nmf_temp$w[,6])

# contribution of markers to factors
Contribution_NFM = nmf_temp$h
colnames(Contribution_NFM) = colnames(Final_intensity_matrix_targeted)
View(t(Contribution_NFM))




