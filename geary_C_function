# Test script for Geary's C spatial autocorrelation statistic
# Values significantly lower than 1 demonstrate increasing positive spatial autocorrelation, whilst values significantly higher than 1 illustrate increasing negative spatial autocorrelation.

# I). Load librarys
library(Matrix)
library(DelayedMatrixStats)
library(fields)

# II). Create Spatial Weight Matrix

# create a distance matrix based on Location points
distance_matrix = fields::rdist(x1 = Location_points, x2 = NULL, compact = FALSE)

# create spatial weighted matrix (adjacency matrix)
spatial_weight_matrix = distance_matrix
spatial_weight_matrix[spatial_weight_matrix > 1] <- 0
# create sparse spatial weight matrix
spatial_weight_matrix = as(spatial_weight_matrix,"dgCMatrix")

# III). Base Code

geary_vector = c()

for (i in 1:ncol(Final_intensity_matrix)){

  x = matrix(Final_intensity_matrix[,i],ncol = 1)
  x_squared = x^2
  
  # 
  product_temp_1 = sum(spatial_weight_matrix%*%x_squared)
  product_temp_2 = t(x)%*%spatial_weight_matrix%*%x
  
  N = length(x)
  # sum of spatial weight matrix
  W = sum(spatial_weight_matrix)
  Var_X = var(x)*N
  
  # unbiased Geary's C formaula
  Geary_C = as.numeric((N-1)*(product_temp_1-product_temp_2)/(2*Var_X*W))
  
  geary_vector = cbind(geary_vector, Geary_C)
}

# IV). Function for Geary's C. x is input final intensity matrix of metapeaks. weights is the spatial weight matrix

gearyC <- function(x =  Final_intensity_matrix, weights = spatial_weight_matrix){
  
  geary_vector = c()
  for (i in 1:ncol(Final_intensity_matrix)){
    x = matrix(Final_intensity_matrix[,i],ncol = 1)
    x_squared = x^2
    product_temp_1 = 2*sum(spatial_weight_matrix%*%x_squared)
    product_temp_2 = 2*t(x)%*%spatial_weight_matrix%*%x
    N = length(x)
    W = sum(spatial_weight_matrix)
    Var_X = var(x)*N
    Geary_C = as.numeric((N-1)*(product_temp_1-product_temp_2)/(2*Var_X*W))
    geary_vector = cbind(geary_vector, Geary_C)
  }
  
  geary_vector = as.numeric(geary_vector)
  
  return(geary_vector)
}

gearyC(Final_intensity_matrix, spatial_weight_matrix)

geary_vector = as.numeric(geary_vector)
View(data.frame(geary_vector))

# boxplot comparing values of annotated metapeaks vs. untargetted
boxplot(geary_vector~is.na(Correpondance_matrix$Expected_mz_location), 
        main = "Targetted vs. Untargetted Geary's C",
        xlab = Sample_name,
        ylab = "Geary's C score", 
        names = c("Annotated Peaks", "Untargetted Peaks"))
