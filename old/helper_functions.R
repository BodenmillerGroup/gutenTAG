find_peaks <- function(x, ignore_threshold = 0, span = 3, strict = TRUE, na.rm = FALSE){
  # find peaks
  if(is.null(span)) {
    pks <- x == max(x, na.rm = na.rm)
    if (strict && sum(pks) != 1L) {
      pks <- logical(length(x)) # all FALSE
    }
  } else {
    pks <- splus2R::peaks(x = x, span = span, strict = strict)
  }
  # apply threshold to found peaks
  if (abs(ignore_threshold) < 1e-5) {
    pks
  } else {
    range_x <- range(x, na.rm = na.rm, finite = TRUE)
    min_x <- range_x[1]
    max_x <- range_x[2]
    x <- ifelse(!is.finite(x), min_x, x)
    # this can cater for the case when max_x < 0, as with logs
    delta <- max_x - min_x
    top_flag <- ignore_threshold > 0.0
    scaled_threshold <- delta * abs(ignore_threshold)
    if (top_flag) {
      ifelse(x - min_x > scaled_threshold, pks , FALSE)
    } else {
      ifelse(max_x - x > scaled_threshold, pks , FALSE)
    }
  }
}


# Load in functions
{ # estimatenoiseMAD function from Cardinal GitHub
  .estimateNoiseMAD <- function(x, blocks=1, fun=mean, tform=diff) {
    mad <- function(y, na.rm) {
      center <- fun(tform(y), na.rm=na.rm)
      fun(abs(y - center), na.rm=na.rm)
    }
    if ( blocks > 1 ) {
      t <- seq_along(x)
      xint <- split_blocks(x, blocks=blocks)
      tint <- split_blocks(t, blocks=blocks)
      noiseval <- sapply(xint, mad, na.rm=TRUE)
      noiseidx <- sapply(tint, mean, na.rm=TRUE)
      noise <- interp1(noiseidx, noiseval, xi=t, method="linear",
                       extrap=fun(noiseval, na.rm=TRUE))
      noise <- supsmu(x=t, y=noise)$y
    } else {
      noise <- mad(x, na.rm=TRUE)
      noise <- rep(noise, length(x))
    }
    noise
  }
  
  # cardinal MAD noise estimation function. set SNR to 10
  peakPick.mad <- function(x, SNR=10, window=25, blocks=1, fun=mean, tform=diff, ...) {
    noise <- .estimateNoiseMAD(x, blocks=blocks, fun=fun, tform=tform)
    maxs <- locmax(x)
    peaks <- intersect(maxs, which(x / noise >= SNR))
    return(peaks)
  }
  locmax = function(x) {
    Local_max_temp = which(diff(sign(diff(x)))==-2)+1
    return(Local_max_temp)
  }
  
  
  # define function to convert to mz scale
  Convert_to_mz_scale = function(x,range_peaks, N_features) {
    scale_vector = base::seq(range_peaks[1], range_peaks[2], length.out =N_features )
    return(scale_vector[x])
  }
  
  # define function for hierarchical clustering with complete linkage
  HC_single_linkage_function = function(x,Threshold_height=2) {
    sub_clustering = 1
    # if there is more than one element in the region, run the distance function
    if (length(x)>1) {
      dist_matrix = dist(x)
      hc_clustering = hclust(dist_matrix,method = "single")
      sub_clustering = cutree(hc_clustering,h =Threshold_height )
    }
    return(sub_clustering)
  }
  
  # define a modified which max function
  which_max_modified = function(x) {
    if (sum(x)==0) {
      y = NA
    }
    else {
      x[x==0] = NA
      y = which.min(x)
    }
    return(y)
  }
  
  string.to.colors = function (string, colors = NULL){
    if (is.factor(string)) {
      string = as.character(string)
    }
    if (!is.null(colors)) {
      if (length(colors) != length(unique(string))) {
        (break)("The number of colors must be equal to the number of unique elements.")
      }
      else {
        conv = cbind(unique(string), colors)
      }
    }
    else {
      conv = cbind(unique(string), rainbow(length(unique(string))))
    }
    unlist(lapply(string, FUN = function(x) {
      conv[which(conv[, 1] == x), 2]
    }))
  }
  
  Otsu_thresholding = function(x,number_bins=100) {
    List_bin = quantile(x,base::seq(from=0,to=1,length.out=number_bins))
    Intravariance_vector = c()
    for (k in 1:number_bins) {
      threshold_temp = List_bin[k]
      s = length(x[x<threshold_temp])*var(x[x<threshold_temp]) + length(x[x>threshold_temp])*var(x[x>threshold_temp])
      Intravariance_vector= c(Intravariance_vector,s)
    }
    Selected_values = List_bin[which.min(Intravariance_vector)]
    return(Selected_values)
  }
  
  
  # NOTE: if trying to image a channel in Rstudio, can change "Final_intensity_matrix" to "Final_intensity_matrix_targeted" (easier to get channel index)
  Plot_channel = function(channel_number=1, quantile_lim = 0.99) {
    Matrix_image = matrix(0,ncol = max(Location_pixels$y),nrow=max(Location_pixels$x))
    
    x = Final_intensity_matrix[,channel_number]
    # get intensity value for 99th percentile most intense pixels
    x_max = quantile(x, probs = quantile_lim)
    # set all pixel values greater than x_max to that of x_max
    x[x>x_max] = x_max
    Matrix_image[as.matrix(Location_pixels)] = x
    # convert matrix to image
    Matrix_image = as.cimg(Matrix_image-min(Matrix_image))
    
    # add colour channels to the image
    Matrix_image = add.color(Matrix_image,simple = TRUE)
    
    # set red and blue channels to zero to get only green
    R(Matrix_image) <- 0
    B(Matrix_image) <- 0
    
    #ifelse(test = Corresponance_matrix$Annotation_peaks[channel_number] == "NA", yes = print("Untargeted Peak"), no = print(Corresponance_matrix$Annotation_peaks[channel_number]))  
    plot((Matrix_image),main=colnames(Final_intensity_matrix[channel_number]))
  }
  
  
  # overlay two channels
  Plot_two_channel = function(channel_1=1,channel_2=2,quantile_lim = 0.99) {
    Matrix_image_1 = matrix(0,ncol = max(Location_pixels$y),nrow=max(Location_pixels$x))
    # remove hotspots channel 1
    x = Final_intensity_matrix_targeted[,channel_1]
    x_max = quantile(x,probs = quantile_lim)
    x[x>x_max]=x_max
    Matrix_image_1[as.matrix(Location_pixels)] = x
    
    Matrix_image_1 = Matrix_image_1/x_max
    # convert matrix to image
    Matrix_image_1 = as.cimg(Matrix_image_1-min(Matrix_image_1))
    
    # channel 2 as image
    Matrix_image_2 = matrix(0,ncol = max(Location_pixels$y),nrow=max(Location_pixels$x))
    # remove hotspots channel 2
    x = Final_intensity_matrix_targeted[,channel_2]
    x_max = quantile(x,probs = quantile_lim)
    x[x>x_max]=x_max
    
    Matrix_image_2[as.matrix(Location_pixels)] = x
    Matrix_image_2 = Matrix_image_2/x_max
    Matrix_image_2 = as.cimg(Matrix_image_2-min(Matrix_image_2))
    
    Matrix_image_1 = add.color(Matrix_image_1,simple = TRUE)
    Matrix_image_1[,,,3] = Matrix_image_2
    
    plot(Matrix_image_1)
  }
  
  # overlay 3 channels on same plot
  Plot_three_channel = function(channel_1=1,channel_2=2,channel_3=3,quantile_lim = 0.99){
    
    Matrix_image_1 = matrix(0,ncol = max(Location_pixels$y),nrow=max(Location_pixels$x))
    # remove hotspots channel 1
    x = Final_intensity_matrix_targeted[,channel_1]
    x_max = quantile(x,probs = quantile_lim)
    x[x>x_max]=x_max
    Matrix_image_1[as.matrix(Location_pixels)] = x
    Matrix_image_1 = Matrix_image_1/x_max
    # convert matrix to image
    Matrix_image_1 = as.cimg(Matrix_image_1-min(Matrix_image_1))
    
    Matrix_image_2 = matrix(0,ncol = max(Location_pixels$y),nrow=max(Location_pixels$x))
    # remove hotspots channel 2
    x = Final_intensity_matrix_targeted[,channel_2]
    x_max = quantile(x,probs = quantile_lim)
    x[x>x_max]=x_max
    Matrix_image_2[as.matrix(Location_pixels)] = x
    Matrix_image_2 = Matrix_image_2/x_max
    # convert matrix to image
    Matrix_image_2 = as.cimg(Matrix_image_2-min(Matrix_image_2))
    
    Matrix_image_3 = matrix(0,ncol = max(Location_pixels$y),nrow=max(Location_pixels$x))
    # remove hotspots channel 3
    x = Final_intensity_matrix_targeted[,channel_3]
    x_max = quantile(x,probs = quantile_lim)
    x[x>x_max]=x_max
    Matrix_image_3[as.matrix(Location_pixels)] = x
    Matrix_image_3 = Matrix_image_3/x_max
    # convert matrix to image
    Matrix_image_3 = as.cimg(Matrix_image_3-min(Matrix_image_3))
    
    
    Matrix_image_1 = add.color(Matrix_image_1,simple = TRUE)
    Matrix_image_1[,,,2] = Matrix_image_2
    Matrix_image_1[,,,3] = Matrix_image_3
    
    plot(Matrix_image_1)
  }
  
}