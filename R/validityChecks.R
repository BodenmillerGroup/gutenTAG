#### Integration tests ####

# preProcess ####
.valid.preProcess <- function(x, cores) {

  # break if input is not a MSImagingExperiment
  if (!is(x, "MSImagingExperiment")) {
    stop("'x' should be of class 'MSImagingExperiment'")
  }

  # break if cores argument is anything other than a single number
  if (length(cores) > 1) {
    stop("'cores' should be a single numeric.")
  }

  if (!is.numeric(cores)) {
    stop("'cores' should be a single numeric.")
  }

}


# peakDetection ####
.valid.peakDetection <- function(x, snr, win, cores){

  # break if input is not a MSImagingExperiment
  if (!is(x, "MSImagingExperiment")) {
    stop("'x' should be of class 'MSImagingExperiment'")
  }

  # break if snr is not a single numeric
  if (length(snr) > 1) {
    stop("'snr' should be a single numeric.")
  }

  if (!is.numeric(snr)) {
    stop("'snr' should be a single numeric.")
  }

  # break if win is not a single numeric
  if (length(win) > 1) {
    stop("'win' should be a single numeric.")
  }

  if (!is.numeric(win)) {
    stop("'win' should be a single numeric.")
  }

  # break if cores is not a single numeric
  if (length(cores) > 1) {
    stop("'cores' should be a single numeric.")
  }

  if (!is.numeric(cores)) {
    stop("'cores' should be a single numeric.")
  }

  # break if input is
  if(!isS4(x)){
    stop("'x' should be an MSImagingExperiement object.")
  }

}


# metapeakGeneration ####
.valid.generateMetapeaks <- function(x, threshold){

  # break if input is not a MSProcessedImagingExperiment
  if (!is(x, "MSProcessedImagingExperiment")) {
    stop("'x' should be of class 'MSProcessedImagingExperiment'")
  }

  # break if threshold is not a single numeric
  if (length(threshold) > 1) {
    stop("'threshold' should be a single numeric.")
  }

  if (!is.numeric(threshold)) {
    stop("'threshold' should be a single numeric.")
  }


  # break if threshold is greater than or equal to 1
  if (threshold >= 1) {
    stop("'threshold' should be a decimal.")
  }

}


# estimateMetapeaks ####
.valid.estimateMetapeaks <- function(count_df, seed_mz, detection_threshold) {
  if (!is.data.frame(count_df)) {
    stop("count_df must be a data frame")
  }
  
  if (!all(c("mz", "count") %in% names(count_df))) {
    stop("count_df must contain 'mz' and 'count' columns")
  }
  
  if (!is.numeric(seed_mz)) {
    stop("seed_mz must be a numeric vector")
  }
  
  if (!is.numeric(detection_threshold) || length(detection_threshold) != 1) {
    stop("detection_threshold must be a single numeric value")
  }
  
  return(TRUE)
}


# assignMetapeaks ####
.valid.assignMetapeaks <- function(x, pre, refList, mz_threshold){

  # break if input is not a list
  if (!is(x, "list")) {
    stop("'x' should be a list object. It must explicitly be the output of the 'metapeakGeneration' function.")
  }

  # break if mz_threshold is not a single numeric
  if (!is.numeric(mz_threshold)) {
    stop("'mz_threshold' should be a single numeric.")
  }

  # break if mz_threshold is not a single numeric
  if (length(mz_threshold) > 1) {
    stop("'mz_threshold' should be a single numeric.")
  }

  # break if pre is not a MSImagingExperiment
  if (!is(pre, "MSImagingExperiment")) {
    stop("'x' should be of class 'MSImagingExperiment'")
  }


}


# asCytoImageList ####
.valid.asCytoImageList <- function(x, name){

  # break if input is not a list
  if (!is(x, "list")) {
    stop("'x' should be a list object. It must explicitly be the output of the 'assignMetapeaks' function.")
  }


}


# asSpatialExperiment ####
.valid.asSpatialExperiment <- function(x){

  # break if input is not a list
  if (!is(x, "list")) {
    stop("'x' should be a list object. It must explicitly be the output of the 'assignMetapeaks' function.")
  }

}


# asMSImagingExperiment ####
.valid.asMSImagingExperiment <- function(x){

  # break if input is not a list
  if (!is(x, "list")) {
    stop("'x' should be a list object. It must explicitly be the output of the 'assignMetapeaks' function.")
  }

}


# pca ####
.valid.pca <- function(x, comp, seed, scree, plot){

  # break if input class is not a dataframe
  if (!is.data.frame(x)) {
    stop("'x' should be a dataframe. It should be the intensity dataframe generated from the output of the 'assignMetapeaks' function. Callable via x$IntensityDF, where x is the output of the 'assignMetapeaks' function")
  }

  # break if comp is not a single numeric
  if (!is.numeric(comp)) {
    stop("'comp' should be a single numeric.")
  }

  # break if comp is not a single numeric
  if (length(comp) > 1) {
    stop("'comp' should be a single numeric.")
  }

  # break if seed is not a single numeric
  if (!is.numeric(seed)) {
    stop("'seed' should be a single numeric.")
  }

  # break if seed is not a single numeric
  if (length(seed) > 1) {
    stop("'seed' should be a single numeric.")
  }

  # break if scree is not a bool
  #if (!is(x, "bool") {
  #  stop("'scree' should be a single numeric.")
  #}

}


# computeNMF ####
.valid.computeNMF <- function(x, comp, cntr){

  # break if input class is not a dataframe
  if (!is.data.frame(x)) {
    stop("'x' should be a dataframe. It should be the intensity dataframe generated from the output of the 'assignMetapeaks' function. Callable via x$IntensityDF, where x is the output of the 'assignMetapeaks' function")
  }

  # break if comp is not a single numeric
  if (!is.numeric(comp)) {
    stop("'comp' should be a single numeric.")
  }

  # break if comp is not a single numeric
  if (length(comp) > 1) {
    stop("'comp' should be a single numeric.")
  }


}

# computeGearysC ####
.valid.computeGearysC <- function(x, verbose, update_correspondence){

  # break if intensity dataframe class is not a dataframe
  if (!is.list(x)) {
    stop("'x' should be a dataframe. It should be the output of the 'assignMetapeaks' function.")
  }

  # break if verbose is not bool
  if(!is.logical(verbose)){
    stop("'verbose' should be a boolean (TRUE/FALSE)")
  }

  # break if update_correspondence is not bool
  if(!is.logical(update_correspondence)){
    stop("'update_correspondence' should be a boolean (TRUE/FALSE)")
  }

}

# computeVariogram ####
.valid.computeVariogram <- function(df, coords){

  # break if intensity dataframe class is not a dataframe
  if (!is.data.frame(df)) {
    stop("'x' should be a dataframe. It should be the intensity dataframe generated from the output of the 'assignMetapeaks' function. Callable via x$IntensityDF, where x is the output of the 'assignMetapeaks' function")
  }

  # break if spatial coords class is not a dataframe
  if (!is.data.frame(coords)) {
    stop("'x' should be a dataframe. It should be the intensity dataframe generated from the output of the 'assignMetapeaks' function. Callable via x$IntensityDF, where x is the output of the 'assignMetapeaks' function")
  }


}


# imageChannel ####
.valid.imageChannel <- function(x, coords, channel_number, quantile_lim, interpolate, axes, colna){

  # break if intensity dataframe class is not a dataframe
  #if (!is.list(x)) {
  #  stop("'x' should be a dataframe. It should be the intensity dataframe generated from the output of the 'assignMetapeaks' function. Callable via x$IntensityDF, where x is the output of the 'assignMetapeaks' function")
  #}

  # break if spatial coords class is not a dataframe
  if (!is.data.frame(coords)) {
    stop("'coords' should be a dataframe.")
  }

  # break if channel_number is not a single numeric
  if (!is.numeric(channel_number)) {
    stop("'channel_number' should be a single numeric.")
  }

  # break if channel_number is not a single numeric
  if (length(channel_number) > 1) {
    stop("'channel_number' should be a single numeric.")
  }

  # break if quantile_lim is not a single numeric
  if (!is.numeric(quantile_lim)) {
    stop("'quantile_lim' should be a single numeric.")
  }

  # break if quantile_lim is not a single numeric
  if (length(quantile_lim) > 1) {
    stop("'quantile_lim' should be a single numeric.")
  }

  # break if interpolate is not bool
  if(!is.logical(interpolate)){
    stop("'interpolate' should be a boolean (TRUE/FALSE)")
  }

  # break if axes is not bool
  if(!is.logical(axes)){
    stop("'axes' should be a boolean (TRUE/FALSE)")
  }

  # break if colna is not either 'black' or 'white'
  if(!(colna %in% c("black", "white"))){
    stop("'colna' should be either 'black' or 'white'")
  }


}

# computeGMM ####
.valid.computeGMM <- function(x, hist){

  # break if intensity dataframe class is not a dataframe
  if (!is.list(x)) {
    stop("'x' should be a list. It should be the output of the 'assignMetapeaks' function.")
  }

  # break if hist is not bool
  if(!is.logical(hist)){
    stop("'hist' should be a bool specifying TRUE/FALSE if you want to plot the histograms fitted with the GMM")
  }


}

# computeUMAP ####
.valid.computeUMAP <- function(x, seed){

  # break if input is not a dataframe
  if (!is.data.frame(x)) {
    stop("'x' should be a dataframe. It should be the intensity dataframe generated from the output of the 'assignMetapeaks' function. Callable via x$IntensityDF, where x is the output of the 'assignMetapeaks' function")
  }

  # break if seed is not a single numeric
  if (!is.numeric(seed)) {
    stop("'seed' should be a single numeric.")
  }

  # break if seed is not a single numeric
  if (length(seed) > 1) {
    stop("'seed' should be a single numeric.")
  }

}


# scaleData ####
.valid.scaleData <- function(x, method){

  # break if spatial coords class is not a dataframe
  if (!is.list(x)) {
    stop("'x' should be a list. It should be the output of the 'assignMetapeaks' function.")
  }

  # break if method is not in list of acceptable strings
  if (!method %in% c("corsd", "geary")) {
    stop("method for scaling must be either 'corsd' or 'geary'")
  }

}


# louvainCluster ####
.valid.louvainCluster <- function(x, coords, k, metric, resolution){

  # break if coords class is not a dataframe
  if (!is.data.frame(coords)) {
    stop("'coords' should be a data.frame. It should be the coordinates dataframe generated from the output of the 'assignMetapeaks' function. Callable via x$SpatialCoords, where x is the output of the 'assignMetapeaks' function")
  }

  # break if k is not a single numeric
  if (!is.numeric(k)) {
    stop("'k' should be a single numeric.")
  }

  # break if k is not a single numeric
  if (length(k) > 1) {
    stop("'k' should be a single numeric.")
  }

  # break if resolution is not a single numeric
  if (!is.numeric(resolution)) {
    stop("'resolution' should be a single numeric.")
  }

  # break if resolution is not a single numeric
  if (length(resolution) > 1) {
    stop("'resolution' should be a single numeric.")
  }

}

