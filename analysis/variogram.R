library(spatial)
library(gstat)
library(spatstat)
library(sp)

# create single marker intensity df with coordinates for every marker 

Final_intensity_matrix_targeted$x <- Location_pixels$x
Final_intensity_matrix_targeted$y <- Location_pixels$y
GLUT1 <- data.frame(Final_intensity_matrix_targeted$GLUT1)
colnames(GLUT1) <- "glut1"
GLUT1$x <- Location_pixels$x
GLUT1$y <- Location_pixels$y
coordinates(GLUT1) <- ~x+y

# Create variogram
vgram <- variogram(glut1~1, GLUT1)

# Plot variogram
plot(vgram)

# Specify an initial model
initial_model <- vgm(psill = 1, model = "Sph", range = 100, nugget = 0.5)

# Fit the model to your empirical variogram
fitted_model <- fit.variogram(vgram, initial_model)
plot(vgram, model = fitted_model)


