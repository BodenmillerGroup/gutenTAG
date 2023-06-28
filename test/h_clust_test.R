# Hierarchical clustering

# 1. Load Packages ####

library(clValid)
library(dendextend)


# 2. Create distance matrix ####

distance_matrix <- dist(head(Final_intensity_matrix_targeted, n = 100), method = "euclidean")
distance_matrix <- as.matrix(distance_matrix)

# 3. Perform hierarchical clustering and dendrogram ####
h_cluster <- hclust(distance_matrix, method = "ward.D")

dendrogram <- as.dendrogram(h_cluster)

# 4. Plot dendrogram ####
plot(dendrogram)
plot(h_cluster, labels = NULL, hang = -1, cex = 0.5,
     main = "Cluster dendrogram",
     xlab = NULL, ylab = "Height")

# plot coloured dendrogram

color_dendrogram <- color_branches(dendrogram, h = 3)
plot(color_dendrogram)



# 5. Calculate Dunn's index for different cluster number k ####

dunn_vector <- data.frame(length = 50)
dunn_df <- data.frame(matrix(NA, nrow = 49, ncol = 1))
for (i in 2:50){

  cut_init <- as.numeric(cutree(h_cluster, k = i, h = NULL))
  cur_dunn <- dunn(distance_matrix, cut_init)
  dunn_df[i,] <- cur_dunn

}

#### Plot Dunn's index as a function of k number ####
plot(x = 1:50, y = dunn_df$matrix.NA..nrow...49..ncol...1.,
     type = "l", main = "Dunn's Index of H-clustering using k of 2-50",
     xlab = "Number of clusters",
     ylab = "Dunn's Index")


# 6. Silhouette plot ####

# cut tree at a certain height
cut_sil <- as.numeric(cutree(h_cluster, k = 5, h = NULL))

# Silhouette plot
sil <- silhouette(cut_sil, distance_matrix)
plot(sil)


#### Create list of all silhouette plots for different k ####
list_sil <- list()
for (j in 2:10){

  cur_cut <- as.numeric(cutree(h_cluster, k = j, h = NULL))
  cur_sil <- silhouette(cur_cut, distance_matrix)
  list_sil[[j]] <- cur_sil


}

plot(list_sil[[3]])


# Count how many observations are in each cluster

clustered_df <- mutate(head(Final_intensity_matrix_targeted, n = 100), cluster = cut_sil)
count(clustered_df, cluster)




