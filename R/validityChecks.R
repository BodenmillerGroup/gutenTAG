#### Integration tests ####

# preProcess ####
.valid.preProcess <- function(x, BPPARAM) {

  # break if input is not a MSImagingExperiment
  if (!is(x, "MSImagingExperiment")) {
    stop("'x' should be of class 'MSImagingExperiment'")
  }

  # break if BPPARAM is not a BiocParallelParam object
  if (!is(BPPARAM, "BiocParallelParam")) {
    stop("'BPPARAM' should be a BiocParallelParam object (e.g. BiocParallel::SerialParam()).")
  }

}


# peakDetection ####
.valid.peakDetection <- function(x, snr, win, BPPARAM){

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

  # break if BPPARAM is not a BiocParallelParam object
  if (!is(BPPARAM, "BiocParallelParam")) {
    stop("'BPPARAM' should be a BiocParallelParam object (e.g. BiocParallel::MulticoreParam(4)).")
  }

}


# generateMetapeaks ####
.valid.generateMetapeaks <- function(x, threshold, fixed.limits, sparsity){

  # break if x is not a MSImagingExperiment
  if (!is(x, "MSImagingExperiment")) {
    stop("'x' should be of class 'MSImagingExperiment'")
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
    stop("'threshold' should be a decimal between 0 and 1.")
  }

  if (threshold < 0) {
    stop("'threshold' should be a decimal between 0 and 1.")
  }

  # break if fixed limits are anything other than NULL or a numeric value
  if (!(is.null(fixed.limits) | is.numeric(fixed.limits))){
    stop("fixed.limits must be a single numeric value or null")
  }

  # break if sparsity parameter anything other than NULL or a positive numeric value
  if (!(is.null(sparsity) | is.numeric(sparsity))){
    stop("sparsity must be a single numeric value.")
  }

  #if (sparsity <= 0){
  #  stop("sparsity must be a single numeric value.")
  #}

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

  # break if input is not a list with the required fields
  if (!is(x, "list")) {
    stop("'x' should be a list object. It must be the output of the 'generateMetapeaks' function.")
  }

  if (!all(c("metapeaks", "propagation_selection") %in% names(x))) {
    stop("'x' must contain 'metapeaks' and 'propagation_selection'. It must be the output of the 'generateMetapeaks' function.")
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
#.valid.pca <- function(x, comp, seed, scree, plot){
#
#  # break if input class is not a dataframe
#  if (!is.data.frame(x)) {
#    stop("'x' should be a dataframe. It should be the intensity dataframe generated from the output of the 'assignMetapeaks' function. Callable via x$IntensityDF, where x is the output of the 'assignMetapeaks' function")
#  }
#
#  # break if comp is not a single numeric
#  if (!is.numeric(comp)) {
#    stop("'comp' should be a single numeric.")
#  }
#
#  # break if comp is not a single numeric
#  if (length(comp) > 1) {
#    stop("'comp' should be a single numeric.")
#  }
#
#  # break if seed is not a single numeric
#  if (!is.numeric(seed)) {
#    stop("'seed' should be a single numeric.")
#  }
#
#  # break if seed is not a single numeric
#  if (length(seed) > 1) {
#    stop("'seed' should be a single numeric.")
#  }
#
#  # break if scree is not a bool
#  #if (!is(x, "bool") {
#  #  stop("'scree' should be a single numeric.")
#  #}
#
#}


# computeNMF ####
#.valid.computeNMF <- function(x, comp, cntr){
#
#  # break if input class is not a dataframe
#  if (!is.data.frame(x)) {
#    stop("'x' should be a dataframe. It should be the intensity dataframe generated from the output of the 'assignMetapeaks' function. Callable via x$IntensityDF, where x is the output of the 'assignMetapeaks' function")
#  }
#
#  # break if comp is not a single numeric
#  if (!is.numeric(comp)) {
#    stop("'comp' should be a single numeric.")
#  }
#
#  # break if comp is not a single numeric
#  if (length(comp) > 1) {
#    stop("'comp' should be a single numeric.")
#  }
#
#
#}

# computeGearysC ####
.valid.computeGearysC <- function(x, verbose, update_correspondence){

  # break if x is not a list with the required fields
  required_fields <- c("IntensityDF", "SpatialCoords", "CorrespondenceMatrix")
  if (!is.list(x) || !all(required_fields %in% names(x))) {
    stop("'x' should be a list. It should be the output of the 'assignMetapeaks' function.")
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

# imageChannel ####
.valid.imageChannel <- function(x, coords, channel, quantile_threshold, palette, na_colour){

  # break if coords is not a dataframe when x is a df/matrix
  if ((is.data.frame(x) || is.matrix(x)) && !is.data.frame(coords)) {
    stop("'coords' must be a data.frame when 'x' is a data.frame or matrix.")
  }

  # break if channel is not a single character or positive integer
  if (!(is.character(channel) || is.numeric(channel)) || length(channel) != 1) {
    stop("'channel' should be a single character name or integer index.")
  }

  # break if quantile_threshold is not a single numeric in (0, 1]
  if (!is.numeric(quantile_threshold) || length(quantile_threshold) != 1 ||
      quantile_threshold <= 0 || quantile_threshold > 1) {
    stop("'quantile_threshold' should be a single numeric in (0, 1].")
  }

  # break if palette is not a character string
  if (!is.character(palette) || length(palette) != 1) {
    stop("'palette' should be a single character string (e.g. \"viridis\", \"magma\").")
  }

  # break if na_colour is not a character string
  if (!is.character(na_colour) || length(na_colour) != 1) {
    stop("'na_colour' should be a single character string.")
  }

}

# computeGMM ####
#.valid.computeGMM <- function(x, hist){
#
#  # break if intensity dataframe class is not a dataframe
#  if (!is.list(x)) {
#    stop("'x' should be a list. It should be the output of the 'assignMetapeaks' function.")
#  }
#
#  # break if hist is not bool
#  if(!is.logical(hist)){
#    stop("'hist' should be a bool specifying TRUE/FALSE if you want to plot the histograms fitted with the GMM")
#  }
#
#
#}

# computeUMAP ####
#.valid.computeUMAP <- function(x, seed){
#
#  # break if input is not a dataframe
#  if (!is.data.frame(x)) {
#    stop("'x' should be a dataframe. It should be the intensity dataframe generated from the output of the 'assignMetapeaks' function. Callable via x$IntensityDF, where x is the output of the 'assignMetapeaks' function")
#  }
#
#  # break if seed is not a single numeric
#  if (!is.numeric(seed)) {
#    stop("'seed' should be a single numeric.")
#  }
#
#  # break if seed is not a single numeric
#  if (length(seed) > 1) {
#    stop("'seed' should be a single numeric.")
#  }
#
#}


# scaleData ####
#.valid.scaleData <- function(x, method){
#
#  # break if spatial coords class is not a dataframe
#  if (!is.list(x)) {
#    stop("'x' should be a list. It should be the output of the 'assignMetapeaks' function.")
#  }
#
#  # break if method is not in list of acceptable strings
#  if (!method %in% c("corsd", "geary")) {
#    stop("method for scaling must be either 'corsd' or 'geary'")
#  }
#
#}



