library(SynchWave)

class(fitted_models)
List_range = unlist(lapply(fitted_models,FUN = function(x) {x$range[2]}))
View(data.frame(List_range ))

Plot_channel(Final_intensity_matrix_targeted, channel = 1)

Ratio_variance = unlist(lapply(fitted_models,FUN = function(x) {x$psill[1]/x$psill[2]}))
View(data.frame(Ratio_variance ))


Matrix_image = matrix(0,ncol = max(Location_pixels$y),nrow=max(Location_pixels$x))

# intensity vector
x = Final_intensity_matrix_targeted[,21]
# get intensity value for 99th percentile most intense pixels
x_max = quantile(x, probs = 0.99)
# set all pixel values greater than x_max to that of x_max
x[x>x_max] = x_max
Matrix_image[as.matrix(Location_pixels)] = x
# convert matrix to image
Matrix_image = as.cimg(Matrix_image-min(Matrix_image, na.rm = TRUE))
plot(Matrix_image)


mimg <- fftw::FFT(Matrix_image)

m = imager::FFT(Matrix_image)
plot(mimg)


# ChatGPT code ####

# Perform Fourier transform on the 2D image matrix
ft <- fft(Matrix_image)

# Compute the magnitude spectrum
ft_magnitude <- Mod(ft)

# Compute the phase spectrum
ft_phase <- Arg(ft)

# Shift the zero frequency component to the center of the spectrum
ft_magnitude_centered <- fftshift(ft_magnitude)

# Plot the magnitude spectrum
image(log(ft_magnitude_centered), main="Magnitude Spectrum")

image(log(ft_magnitude), main="Magnitude Spectrum")

# Plot the phase spectrum
image(ft_phase, main="Phase Spectrum")



