# Leiden clustering 

# Import libraries
library(igraph)

# Use test dataframe for experimenting
test_data <- head(Final_intensity_matrix_targeted, n = 50)
index_col <- seq(test_data$YAP1)
test_data <- cbind(index_col, test_data)
test_data_long <- pivot_longer(test_data, cols = -c(index_col), names_to = "name", values_to = "intensity")
test <- test_data_long[, -1]

# create distance matrix
dissimilarity_matrix <- dist(test_data, method = "euclidean") 

# Convert the dist object to a matrix
dissimilarity_matrix <- as.matrix(dissimilarity_matrix)

# Compute similarity_matrix from dissimilarity_matrix
similarity_matrix <- max(dissimilarity_matrix) - dissimilarity_matrix
# rescale
similarity_matrix <- similarity_matrix / max(similarity_matrix)


# Dataframe to graph
g <- graph_from_adjacency_matrix(similarity_matrix, mode = "undirected")

print(g)
plot(g)


# Toy example 

actors <- data.frame(name=c("Alice", "Bob", "Cecil", "David",
                            "Esmeralda"),
                     age=c(48,33,45,34,21),
                     gender=c("F","M","F","M","F"))
relations <- data.frame(from=c("Bob", "Cecil", "Cecil", "David",
                               "David", "Esmeralda"),
                        to=c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"),
                        same.dept=c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),
                        friendship=c(4,5,5,2,1,1), advice=c(4,5,5,4,2,3))
g <- graph_from_data_frame(relations, directed=FALSE, vertices=actors)
print(g, e=TRUE, v=TRUE)
plot(g)



## The opposite operation
igraph::as_data_frame(g, what="vertices")
igraph::as_data_frame(g, what="edges")


