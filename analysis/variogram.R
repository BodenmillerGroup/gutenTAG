library(spatial)
library(gstat)
library(spatstat)
library(sp)

# create single marker intensity df with coordinates for every marker


constructVariogram <- function(df = Final_intensity_matrix_targeted, coords = Location_pixels){

  names(df) <- make.names(names(df))
  names <- colnames(df)
  variograms <- list()

  for (i in 1:ncol(df)){

    marker <- data.frame(df[i])
    colnames(marker) <- names[i]

    marker$x <- coords$x
    marker$y <- coords$y

    coordinates(marker) <- ~x+y

    # construct a formula from the current column name
    formula_string <- paste0(names[i], "~1")
    formula_object <- as.formula(formula_string)

    marker_variogram <- variogram(formula_object, marker)

    #marker_variogram <- variogram(names[i], marker)



    # Store the result in the list
    variograms[[names[i]]] <- marker_variogram



  }

  return(variograms)

  }


vgrams <- constructVariogram(df = Final_intensity_matrix_targeted, coords = Location_pixels)


GLUT1 <- data.frame(Final_intensity_matrix_targeted$GLUT1)
colnames(GLUT1) <- "glut1"
GLUT1$x <- Location_pixels$x
GLUT1$y <- Location_pixels$y
coordinates(GLUT1) <- ~x+y

# Create variogram
vgram <- variogram(glut1~1, GLUT1)



# Plot variogram
plot(vgrams$FN)

# Specify an initial model
initial_model <- vgm(psill = 1, model = "Sph", range = 100, nugget = 0.5)

# Fit the model to your empirical variogram
fitted_model <- fit.variogram(vgrams$GLUT1, initial_model)
plot(vgrams$GLUT1, model = fitted_model)


