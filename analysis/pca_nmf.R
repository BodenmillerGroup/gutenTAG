#### PCA ####

# Load packages
library(irlba)
library(RcppML)
library(EBImage)
library(imager)

# Load in functions from helper script
source("/mnt/msi_volume/Rscripts/maldi-processing/code/helper_functions.R")

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
Plot_channel(PCA_simple$x[,1])

# plot the inverse of the principal components
Plot_channel(-PCA_simple$x[,1])


#### NMF ####

# run NMF
nmf_temp = nmf(scale(Final_intensity_matrix_targeted,center = FALSE),k = 10)

# scree plot
barplot(nmf_temp$d)
dim(nmf_temp$w)

# plot factors
Plot_channel(nmf_temp$w, 3)

# contribution of markers to factors
Contribution_NFM = nmf_temp$h
colnames(Contribution_NFM) = colnames(Final_intensity_matrix_targeted)
View(t(Contribution_NFM))




