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

#Peak detection function

# Inputs are a pre-processed Cardinal MSImagingExperiment object and signal to noise ratio

#peakDetection <- function(preProcessed = preprocessedData, snr = 10){
#  
#  raw_intensity = preProcessed@imageData$data$intensity
#  
#  List_peaks_test = apply(raw_intensity, MARGIN = 2, FUN = function(x) {peakPick.mad(x,SNR = snr,window = 400)})
#  
#  return(List_peaks_test)
#}

#### Development: Full Peak alignment function #### 

# Inputs are a pre-processed Cardinal MSImagingExperiment object for setting parameters, a list of detected peaks in each pixel, a reference list
peakAlignment <- function(preProcessed = preprocessedData, detectedPeaks = detected_peaks, refList = peakAnnotation){
  
  # define parameters
  N_peak_per_pixel = unlist(lapply(detectedPeaks,FUN = length))
  range_peaks = range(preProcessed@featureData@mz)
  N_features = length(preProcessed@featureData@mz)
  N_pixels_total = dim(preProcessed)[2]
  Threshold_Detection = N_pixels_total*0.01
  mz_vector = as.data.frame(preProcessed@featureData)
  Location_pixels = as.data.frame(pData(peakPre))[,c("x","y")]
  
  # how many pixels is each peak in?
  detected_peaks_rescaled = mz_vector[unlist(detectedPeaks),1]
  detected_peaks_unique = unique(unlist(detectedPeaks))
  detected_peaks_unique_rescaled = mz_vector$mz[detected_peaks_unique]
  detected_peaks_unique_rescaled_ordered = detected_peaks_unique_rescaled[order(detected_peaks_unique_rescaled)]
  Freq_peak_rescaled_ordered = table(factor(detected_peaks_rescaled,levels = mz_vector$mz))
  Freq_peak_rescaled_ordered = as.numeric(Freq_peak_rescaled_ordered)
  
  # smooth the intensities with a Gaussian curve with a sigma of 1mz so that isotopic peaks are merged
  mz_index_ratio = (max(mz_vector$mz)-min(mz_vector$mz))/N_features
  sigma_smoothing_isotopic = 1/mz_index_ratio
  list_values = seq(-N_features/2, N_features/2, length.out = N_features)
  Gaussian_vector = exp(-list_values^2/sigma_smoothing_isotopic^2)
  Freq_peak_rescaled_ordered_smooth = stats::convolve(y = Gaussian_vector,x = Freq_peak_rescaled_ordered, conj = TRUE)
  Freq_peak_rescaled_ordered_smooth = Freq_peak_rescaled_ordered_smooth[c((N_features/2):N_features,1:(N_features/2-1))]
  
  # remove values under the pixel detection threshold
  Freq_peak_smoothed_thresholded = Freq_peak_rescaled_ordered_smooth
  Freq_peak_smoothed_thresholded[Freq_peak_smoothed_thresholded<Threshold_Detection] = 0
  
  # We want to find local maxima within a range of 3mz 
  span_local_maxima = 3/mz_index_ratio
  Local_maxima_peak_freq = base::which(find_peaks(Freq_peak_smoothed_thresholded,span = span_local_maxima))
  Local_maxima_peak_freq_reshaped = matrix(rep(0,length(Freq_peak_rescaled_ordered)),ncol = 1)
  Local_maxima_peak_freq_reshaped[Local_maxima_peak_freq]=1:length(Local_maxima_peak_freq)
  
  Freq_peak_smoothed_thresholded = matrix(Freq_peak_smoothed_thresholded,ncol = 1)
  Propagation_selection = propagate(Freq_peak_smoothed_thresholded, seeds = Local_maxima_peak_freq_reshaped, mask = Freq_peak_smoothed_thresholded>Threshold_Detection)
  
  # Find peak centers, peak maxima and peak width
  metapeak_location = aggregate.data.frame(data.frame(Location = mz_vector$mz,Freq = Freq_peak_rescaled_ordered),
                                       by=list(as.numeric(Propagation_selection)),FUN = function(x) {sum(x[[1]]*x[[2]])/sum(x[[2]])})
  metapeak_center = c()
  metapeak_max = c()
  metapeak_delimitation = c()
  for (k in 1:max(Propagation_selection)) {
    
    Location_centered_temp = weighted.mean(x=mz_vector$mz[as.numeric(Propagation_selection)==k],
                                           w =Freq_peak_rescaled_ordered[as.numeric(Propagation_selection)==k] )
    
    metapeak_center = c(metapeak_center,Location_centered_temp)
    
    Location_max_temp  = (mz_vector$mz[as.numeric(Propagation_selection)==k])[which.max(Freq_peak_rescaled_ordered[as.numeric(Propagation_selection)==k])]
    metapeak_max = c(metapeak_max,Location_max_temp)
    
    Cumsum_freq = (cumsum(Freq_peak_rescaled_ordered[as.numeric(Propagation_selection)==k])/sum(Freq_peak_rescaled_ordered[as.numeric(Propagation_selection)==k]))
    
    Beginning_peak = max(which(Cumsum_freq < 0.01))
    End_peak = min(which(Cumsum_freq > 0.99))
    if (is.infinite(Beginning_peak)) {
      Beginning_peak = 1
    }
    
    metapeak_delimitation = rbind(metapeak_delimitation,c(mz_vector$mz[as.numeric(Propagation_selection)==k][Beginning_peak],mz_vector$mz[as.numeric(Propagation_selection)==k][End_peak]))
  }
  
  metapeak_width = metapeak_delimitation[,2]-metapeak_delimitation[,1]
  metapeak_width = metapeak_width*mz_index_ratio
  
  # Create reference list including ion adducts
  FeatureMassHydrogen = refList$FeatureMass
  FeatureMassSodium = (FeatureMassHydrogen - 1) + 23
  FeatureMassAmmonium = (FeatureMassHydrogen - 1) + 18
  refList = cbind(refList, FeatureMassSodium, FeatureMassAmmonium)
  refList = rbind(cbind(refList$Name, refList$FeatureMass),
                  cbind(paste(refList$Name,"Sodium", sep = "_"), refList$FeatureMassSodium),
                  cbind(paste(refList$Name,"Ammonium", sep = "_"), refList$FeatureMassAmmonium))
  refList = as.data.frame(refList)
  colnames(refList) = c("FeatureAnnotation","Location")
  refList$Location = as.numeric(refList$Location)
  
  # cross compare total reference list with metapeaks (find closest peak across lists)
  Mapping_meta = N2R::crossKnn(mA = matrix(metapeak_max, ncol = 1),
                          mB= matrix(refList$Location, ncol = 1), k = 10, indexType = "L2")
  mz_threshold_association = 1
  Mapping_meta[Mapping_meta>mz_threshold_association]=0
  Mapping_meta_cleaned = apply(as.matrix(Mapping_meta),MARGIN = 2,FUN = which_max_modified)
  
  Associated_marker = c()
  for (k in 1:length(metapeak_max)) {
    selected_mz_location = range(mz_vector$mz[as.numeric(Propagation_selection)==k])
    Annotation_assigned = refList$FeatureAnnotation[refList$Location>selected_mz_location[1] & refList$Location<selected_mz_location[2]]
    if (length(Annotation_assigned)==0){
      Annotation_assigned = NA
    }
    if (length(Annotation_assigned)>1){
      Annotation_assigned = paste(Annotation_assigned,collapse="__OR__")
    }
    Associated_marker = c(Associated_marker,Annotation_assigned)
  }
  
  # Create a correspondence matrix
  Correspondence_matrix = data.frame(mz_location = metapeak_max,
                                     Expected_mz_location = refList$Location[Mapping_meta_cleaned],
                                     Annotation_peaks = refList$FeatureAnnotation[Mapping_meta_cleaned])
  rownames(refList) = refList$FeatureAnnotation
  
  # Processed Intensity Matrix
  raw_intensity = preProcessed@imageData$data$intensity
  Final_intensity_matrix = c()
  for (k in 1:nrow(Correspondence_matrix)) {
    selected_mz_location = mz_vector$mz >= metapeak_delimitation[k,1] & mz_vector$mz <= metapeak_delimitation[k,2]
    mz_vector$mz[which(selected_mz_location == T)]
    if (sum(selected_mz_location)>1) {
      intensity_temp = colSums(raw_intensity[selected_mz_location,])
    }
    if (sum(selected_mz_location)==1) {
      intensity_temp = (raw_intensity[selected_mz_location,])
    }
    Final_intensity_matrix = cbind(Final_intensity_matrix,intensity_temp)
    
  }
  
  # Remove untargeted peaks
  Final_intensity_matrix_targeted = Final_intensity_matrix[,!is.na(Correspondence_matrix$Annotation_peaks)]
  colnames(Final_intensity_matrix_targeted) = Correspondence_matrix$Annotation_peaks[!is.na(Correspondence_matrix$Annotation_peaks)]
  Final_intensity_matrix_targeted = as.data.frame(Final_intensity_matrix_targeted)
  
  processed_data <- list("IntensityMatrix" = Final_intensity_matrix_targeted, "CorrespondenceMatrix" = Correspondence_matrix, "SpatialCoordinates" = Location_pixels)
  
  return(processed_data)
  
}


#### peakReshape function ####

peakReshape <- function(pre = peakPre, detected = List_peaks, plot = TRUE){
  
  range_peaks = range(mz(pre))
  N_features = length(mz(pre))
  N_pixels = dim(pre)["Pixels"]
  Threshold_Detection = N_pixels*0.01
  # all m/z values 
  mz_vector = as.data.frame(mz(pre))
  Location_pixels = as.data.frame(pData(pre))[,c("x","y")]
  
  # how many pixels is each peak in?
  detected <- peakData(detected)[["mz"]]
  detected <- unlist(detected)
  detected <- detected[order(detected)]
  # how many pixels is each peak in?
  Freq_peak = table(factor(detected,levels = mz_vector$mz))
  Freq_peak = as.numeric(Freq_peak)
  
  # smooth the intensities with a Gaussian curve with a sigma of 1mz so that isotopic peaks are merged
  mz_index_ratio = (max(mz_vector$mz)-min(mz_vector$mz))/N_features
  sigma_smoothing_isotopic = 1/mz_index_ratio
  # mass range vector centered around zero
  list_values = seq(-N_features/2, N_features/2, length.out = N_features)
  Gaussian_vector = exp(-list_values^2/sigma_smoothing_isotopic^2)
  Freq_peak_smooth = stats::convolve(y = Gaussian_vector,x = Freq_peak, conj = TRUE)
  Freq_peak_smooth = Freq_peak_smooth[c((N_features/2):N_features,1:(N_features/2-1))]
  
  # remove values under the pixel detection threshold
  Freq_peak_smooth_thresholded = Freq_peak_smooth
  Freq_peak_smooth_thresholded[Freq_peak_smooth_thresholded<Threshold_Detection] = 0
  
  # We want to find local maxima within a range of 3mz 
  span_local_maxima = 3/mz_index_ratio
  Local_maxima_peak_freq = base::which(find_peaks(Freq_peak_smooth_thresholded,span = span_local_maxima))
  Local_maxima_peak_freq_reshaped = matrix(rep(0,length(Freq_peak)),ncol = 1)
  Local_maxima_peak_freq_reshaped[Local_maxima_peak_freq]=1:length(Local_maxima_peak_freq)
  
  Freq_peak_smooth_thresholded = matrix(Freq_peak_smooth_thresholded,ncol = 1)
  Propagation_selection = propagate(Freq_peak_smooth_thresholded, seeds = Local_maxima_peak_freq_reshaped, mask = Freq_peak_smooth_thresholded>Threshold_Detection)
  
  # Find peak centers, peak maxima and peak width
  metapeak_location = aggregate.data.frame(data.frame(Location = mz_vector$mz,Freq = Freq_peak),
                                       by=list(as.numeric(Propagation_selection)),FUN = function(x) {sum(x[[1]]*x[[2]])/sum(x[[2]])})
  metapeak_center = c()
  metapeak_max = c()
  metapeak_delimitation = c()
  for (k in 1:max(Propagation_selection)) {
    
    Location_centered_temp = weighted.mean(x=mz_vector$mz[as.numeric(Propagation_selection)==k],
                                           w =Freq_peak[as.numeric(Propagation_selection)==k] )
    
    metapeak_center = c(metapeak_center,Location_centered_temp)
    
    Location_max_temp  = (mz_vector$mz[as.numeric(Propagation_selection)==k])[which.max(Freq_peak[as.numeric(Propagation_selection)==k])]
    metapeak_max = c(metapeak_max,Location_max_temp)
    
    Cumsum_freq = (cumsum(Freq_peak[as.numeric(Propagation_selection)==k])/sum(Freq_peak[as.numeric(Propagation_selection)==k]))
    
    Beginning_peak = max(which(Cumsum_freq < 0.01))
    End_peak = min(which(Cumsum_freq > 0.99))
    if (is.infinite(Beginning_peak)) {
      Beginning_peak = 1
    }
    
    metapeak_delimitation = rbind(metapeak_delimitation,c(mz_vector$mz[as.numeric(Propagation_selection)==k][Beginning_peak],mz_vector$mz[as.numeric(Propagation_selection)==k][End_peak]))
  }
  
  metapeak_width = metapeak_delimitation[,2] - metapeak_delimitation[,1]
  metapeak_width = metapeak_width*mz_index_ratio
  
  if(plot == TRUE){
    plot(mz_vector$mz,Freq_peak_smooth,,xlab="mz",ylab="Number of peaks detected",xlim=c(800,2000),type="l", main = paste("Metapeaks: ", length(metapeak_location$Location)))
    points(mz_vector$mz,Freq_peak,lty=2,col="red",type="l")
    points(mz_vector$mz[Local_maxima_peak_freq],Freq_peak_smooth[Local_maxima_peak_freq],pch=21,bg="red3")
  }
  
  peak_data <- list("peakMax" = metapeak_max, "peakCentre" = metapeak_center, "peakWidth" = metapeak_width, "peakDelim" = metapeak_delimitation)
  
  metapeaks <- list("location" = metapeak_location$Location, "propagation" = Propagation_selection)
  
  reshaped_peaks <- list("peakData" = peak_data, "metapeaks" = metapeaks)
  
  return(reshaped_peaks)
  
}

#### peakAlignment function ####

peakAlignment <- function(pre, detected, metapeaks, ref){
  
  # reference list
  ref = as.data.frame(ref)
  colnames(ref) = c("FeatureAnnotation","Location")
  ref$Location = as.numeric(ref$Location)
  
  mz_vector = as.data.frame(mz(pre))
  Location_pixels = as.data.frame(pData(pre))[,c("x","y")]
  
  # cross compare total reference list with metapeaks (find closest peak across lists)
  Mapping_meta = N2R::crossKnn(mA = matrix(metapeaks$peakData$peakMax, ncol = 1),
                               mB= matrix(ref$Location, ncol = 1), k = 10, indexType = "L2", verbose = FALSE)
  mz_threshold_association = 1
  Mapping_meta[Mapping_meta>mz_threshold_association]=0
  Mapping_meta_cleaned = apply(as.matrix(Mapping_meta),MARGIN = 2,FUN = which_max_modified)
  
  Associated_marker = c()
  for (k in 1:length(metapeaks$peakData$peakMax)) {
    selected_mz_location = range(mz_vector$mz[as.numeric(metapeaks$metapeaks$propagation)==k])
    Annotation_assigned = ref$FeatureAnnotation[ref$Location>selected_mz_location[1] & ref$Location<selected_mz_location[2]]
    if (length(Annotation_assigned)==0){
      Annotation_assigned = NA
    }
    if (length(Annotation_assigned)>1){
      Annotation_assigned = paste(Annotation_assigned,collapse="__OR__")
    }
    Associated_marker = c(Associated_marker,Annotation_assigned)
  }
  
  # Create a correspondence matrix
  Correspondence_matrix = data.frame(mz_location = metapeaks$peakData$peakMax,
                                     Expected_mz_location = ref$Location[Mapping_meta_cleaned],
                                     Annotation_peaks = ref$FeatureAnnotation[Mapping_meta_cleaned])
  rownames(ref) = ref$FeatureAnnotation
  
  # Processed Intensity Matrix
  raw_intensity = pre@imageData$data$intensity
  Final_intensity_matrix = c()
  for (k in 1:nrow(Correspondence_matrix)) {
    selected_mz_location = mz_vector$mz >= metapeaks$peakData$peakDelim[k,1] & mz_vector$mz <= metapeaks$peakData$peakDelim[k,2]
    mz_vector$mz[which(selected_mz_location == T)]
    if (sum(selected_mz_location)>1) {
      intensity_temp = colSums(raw_intensity[selected_mz_location,])
    }
    if (sum(selected_mz_location)==1) {
      intensity_temp = (raw_intensity[selected_mz_location,])
    }
    Final_intensity_matrix = cbind(Final_intensity_matrix,intensity_temp)
    
  }
  
  # Remove untargeted peaks
  Final_intensity_matrix_targeted = Final_intensity_matrix[,!is.na(Correspondence_matrix$Annotation_peaks)]
  colnames(Final_intensity_matrix_targeted) = Correspondence_matrix$Annotation_peaks[!is.na(Correspondence_matrix$Annotation_peaks)]
  Final_intensity_matrix_targeted = as.data.frame(Final_intensity_matrix_targeted)
  
  processed_data <- list("IntensityData" = Final_intensity_matrix_targeted, "CorrespondenceMatrix" = Correspondence_matrix, "SpatialCoordinates" = Location_pixels)
  
  
  
  return(processed_data)
  
}
  
  #### Marker panel cleaning function ####
  
  cleanPanel <- function(panel){
    
    # 1. ensure column names are correct
    
    if (is.numeric(panel$Name) == T){
      colnames(panel) = c("FeatureMass","Name")
    }
    
    # rearrange column order so Name is first column
    relocate(panel, "Name", .before = "FeatureMass")
    
    # 2. sort by mass tag size
    
    panel = dplyr::arrange(panel, panel$FeatureMass)
    
    # 3. clean marker names
    
    for(a in seq_along(panel$Name)){
      # replace spaces between + symbol
      if (grepl("+", panel$Name[a], fixed=T)){
        new_string = gsub(" ", "", panel$Name[a])
        panel$Name[a] <- new_string
      }
      else{
        # split strings with spaces into list with individual strings as elements
        new_string = gsub(" ", "-", panel$Name[a])
        # replace fullstops with underscores
        new_string = gsub(".", "-", new_string, fixed = T)
        # if there is a dash at the end of the name, remove it
        if(endsWith(new_string, "-")){
          new_string = substr(new_string,1, nchar(new_string)-1)
        }
        # replace names in peakAnnotation
        panel$Name[a] <- new_string
      }
      # remove slashes
      if (grepl("/", panel$Name[a], fixed=T)){
        new_string = gsub("/", "", panel$Name[a])
        panel$Name[a] <- new_string
      }
    }
    
    # 4. rearrange column order so Name is first column
    panel <- relocate(panel, "Name", .before = "FeatureMass")
    
    return(panel)
    
  }


  #### Convert_to_mz_scale function ####
  
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
  
  #### which_max_modified function ####
  
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
  
  #### Otsu thresholding ####
  
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
  

  
  #### Plot_mz_channel function ####
  
  Plot_mz_channel = function(dataframe = processedData$IntensityMatrix, coords = processedData$SpatialCoordinates,
                             channel =1, quantile_lim = 0.99, col = "blue") {
    
    coords$x = coords$x - min(coords$x) + 1
    coords$y = coords$y - min(coords$y) + 1
    
    Matrix_image = matrix(0,ncol = max(coords$y), nrow=max(coords$x))
    
    x = dataframe[,channel]
    x_max = quantile(x, probs = quantile_lim)
    x[x>x_max] = x_max
    Matrix_image[as.matrix(coords)] = x
    Matrix_image = as.cimg(Matrix_image-min(Matrix_image))
    
    Matrix_image = add.color(Matrix_image,simple = TRUE)
    
    if (col == "green"){
      R(Matrix_image) <- 0
      B(Matrix_image) <- 0
    }
    
    if (col == "blue"){
      R(Matrix_image) <- 0
      G(Matrix_image) <- 0
    }
    
    if (col == "red"){
      G(Matrix_image) <- 0
      B(Matrix_image) <- 0
    }
    
    if (col == "cyan"){
      R(Matrix_image) <- 0
    }
    
    if (col == "magenta"){
      G(Matrix_image) <- 0
    }
    
    if (col == "yellow"){
      B(Matrix_image) <- 0
    }
    
    if (col == "white"){
    }
    
    plot((Matrix_image))
    
  }
  


  
  #### Plot_two_mz_channel function ####
  
  ## plotting function for comp group meeting
  Plot_two_mz_channel = function(dataframe = processedData$IntensityData, coords = processedData$SpatialCoordinates, 
                                 channel_1=1, channel_2=2, quantile_lim = 0.99) {
    
    coords$x = coords$x - min(coords$x) + 1
    coords$y = coords$y - min(coords$y) + 1
    
    Matrix_image_1 = matrix(0,ncol = max(coords$y),nrow=max(coords$x))
    # remove hotspots channel 1
    x = dataframe[,channel_1]
    x_max = quantile(x,probs = quantile_lim)
    x[x>x_max]=x_max
    Matrix_image_1[as.matrix(coords)] = x
    
    Matrix_image_1 = Matrix_image_1/x_max
    # convert matrix to image
    Matrix_image_1 = as.cimg(Matrix_image_1-min(Matrix_image_1))
    
    # channel 2 as image
    Matrix_image_2 = matrix(0,ncol = max(coords$y),nrow=max(coords$x))
    # remove hotspots channel 2
    x = dataframe[,channel_2]
    x_max = quantile(x,probs = quantile_lim)
    x[x>x_max]=x_max
    
    Matrix_image_2[as.matrix(coords)] = x
    Matrix_image_2 = Matrix_image_2/x_max
    Matrix_image_2 = as.cimg(Matrix_image_2-min(Matrix_image_2))
    
    Matrix_image_1 = add.color(Matrix_image_1,simple = TRUE)
    Matrix_image_1[,,,3] = Matrix_image_2
    
    plot(Matrix_image_1)
    
    text(x = max(coords$x)*1.02, y = max(coords$y)*0.1,
         labels = colnames(dataframe[channel_1]),
         adj = 0,
         cex = 1,
         col = "blue")
    text(x = max(coords$x)*1.02, y = max(coords$y)*0.18,
         adj = 0,
         labels = colnames(dataframe[channel_2]),
         cex = 1,
         col = "gold")
    
  }
  
  

  
}


#### Original Plot_channel functions ####

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
  
  plot((Matrix_image))#,main=colnames(Final_intensity_matrix[channel_number]))
  text(x = max(Location_pixels$x)*1.02, y = max(Location_pixels$y)*0.1,
       labels = colnames(Final_intensity_matrix[channel_number]),
       adj = 0,
       cex = 1,
       col = "green")
}

Plot_channel_rstudio = function(channel_number=1, quantile_lim = 0.99) {
  Matrix_image = matrix(0,ncol = max(Location_pixels$y),nrow=max(Location_pixels$x))
  
  x = Final_intensity_matrix_targeted[,channel_number]
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
  #G(Matrix_image) <- 0
  R(Matrix_image) <- 0
  #B(Matrix_image) <- 0
  
  plot((Matrix_image))
  #text(x = max(Location_pixels$x)*1.02, y = max(Location_pixels$y)*0.1,
  #    labels = colnames(Final_intensity_matrix_targeted[channel_number]),
  #   adj = 0,
  #  cex = 1,
  # col = "green")
}

# overlay two channels
Plot_two_channel = function(channel_1=1,channel_2=2,quantile_lim = 0.99) {
  Matrix_image_1 = matrix(0,ncol = max(Location_pixels$y),nrow=max(Location_pixels$x))
  x = Final_intensity_matrix_targeted[,channel_1]
  x_max = quantile(x,probs = quantile_lim)
  x[x>x_max]=x_max
  Matrix_image_1[as.matrix(Location_pixels)] = x
  
  Matrix_image_1 = Matrix_image_1/x_max
  Matrix_image_1 = as.cimg(Matrix_image_1-min(Matrix_image_1))
  
  Matrix_image_2 = matrix(0,ncol = max(Location_pixels$y),nrow=max(Location_pixels$x))
  x = Final_intensity_matrix_targeted[,channel_2]
  x_max = quantile(x,probs = quantile_lim)
  x[x>x_max]=x_max
  
  Matrix_image_2[as.matrix(Location_pixels)] = x
  Matrix_image_2 = Matrix_image_2/x_max
  Matrix_image_2 = as.cimg(Matrix_image_2-min(Matrix_image_2))
  
  Matrix_image_1 = add.color(Matrix_image_1,simple = TRUE)
  Matrix_image_1[,,,3] = Matrix_image_2
  
  plot(Matrix_image_1)
  
  text(x = max(Location_pixels$x)*1.02, y = max(Location_pixels$y)*0.1,
       labels = colnames(Final_intensity_matrix_targeted[channel_1]),
       adj = 0,
       cex = 1,
       col = "blue")
  text(x = max(Location_pixels$x)*1.02, y = max(Location_pixels$y)*0.18,
       adj = 0,
       labels = colnames(Final_intensity_matrix_targeted[channel_2]),
       cex = 1,
       col = "gold")
  
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