#### UMAP ####

install.packages("umap")
library(umap)

set.seed(123)


####  ####

head <- head(Final_intensity_matrix_targeted, n=2000)
head_umap <- umap(head, n_components = 1, n_neighbours = 50)

plot(head_umap$layout[,1], head_umap$layout[,2], pch=21, bg = c("red3", "lightgrey"))

