# 1. Load Packages ####

# Load packages
{
library(irlba)
library(RcppML)
library(EBImage)
library(imager)
library(ggplot2)
  }

# Load in functions from helper script
source("/mnt/msi_volume/Rscripts/maldi-processing/code/helper_functions.R")

# 2. PCA ####

# check dimensions of
dim(Final_intensity_matrix_targeted)

# variance and mean of columns (each marker)
var_temp = apply(Final_intensity_matrix_targeted,MARGIN = 2,FUN = var)
mean_temp = colMeans(Final_intensity_matrix_targeted)

# check linearity
plot(mean_temp,var_temp,log="xy")

# compute PCA (adjust n to set number of principal components)
PCA_simple = prcomp_irlba(Final_intensity_matrix_targeted, n = 4)

# summary of Principal Componants

summary(PCA_simple)

# get contribution of each feature to each principal component
Contribution = PCA_simple$rotation
rownames(Contribution) = colnames(Final_intensity_matrix_targeted)
View(Contribution)


# create df for scree plot
PCA_variance <- data.frame(Variation_Explained = (PCA_simple$sdev)^2/PCA_simple$totalvar*100,
                           Component = colnames(PCA_simple$x))


# ggplot scree plot
scree_plot <- ggplot(data=PCA_variance, aes(x = as.factor(Component), y = Variation_Explained, fill = Component)) +
  geom_bar(stat = "identity") +
  labs(title = "Scree Plot", x = "", y = "Variance Explained") +
  scale_y_continuous(limit = c(0, 100))

#### Save PCA ####

pca_dir = paste(sample_path,"/", "PCA",sep = "")
# if directory does not exist, create it
if( !dir.exists(pca_dir)) {
  dir.create(pca_dir)
}
print(pca_dir)

# Save scree plot
ggsave(filename = "scree_plot.png",
       plot = scree_plot,
       device = "png",
       path = pca_dir)


# plot principal components
Plot_channel(image = PCA_simple$x, channel_number = 3)

# plot the inverse of the principal components
Plot_channel(image = -PCA_simple$x, channel_number = 1)

PCA_df <- data.frame(PCA_simple$x)

# Save PCA images #
for (n in 1:ncol(PCA_df)){

    pca_image_path = paste(pca_dir,"/", colnames(PCA_df)[n], ".png", sep = "")
    print(pca_image_path)
    png(file=pca_image_path, width=1200, height=700)
    Plot_channel(image = PCA_df, channel_number = n)
    dev.off()

}


# Save PCA loadings
write.table(Contribution, file = paste(pca_dir,"/", "PCA_loadings.txt", sep = ""), sep="\t", quote = FALSE, row.names = FALSE)



# 3. NMF ####

# run NMF
nmf_temp = nmf(scale(Final_intensity_matrix_targeted,center = FALSE),k = 10)

# scree plot
barplot(nmf_temp$d, main = "NMF", names = seq(1:10))
dim(nmf_temp$w)

# plot factors
Plot_channel(nmf_temp$w, 10)

# contribution of markers to factors
Contribution_NFM = nmf_temp$h
colnames(Contribution_NFM) = colnames(Final_intensity_matrix_targeted)
Contribution_NFM <- t(Contribution_NFM)


#### Save NMF ####

nmf_dir = paste(sample_path,"/", "NMF",sep = "")
# if directory does not exist, create it
if( !dir.exists(nmf_dir)) {
  dir.create(nmf_dir)
}
print(nmf_dir)

# Save NMF contribution
write.table(Contribution, file = paste(nmf_dir,"/", "NMF_contribution.txt", sep = ""), sep="\t", quote = FALSE, row.names = FALSE)

NMF_df <- data.frame(nmf_temp$w)

# Save NMF images #
for (n in 1:ncol(NMF_df)){

  nmf_image_path = paste(nmf_dir,"/", colnames(NMF_df)[n], ".png", sep = "")
  print(nmf_image_path)
  png(file=nmf_image_path, width=1200, height=700)
  Plot_channel(image = NMF_df, channel_number = n)
  dev.off()

}


