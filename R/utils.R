############################# Reader helpers ###################################

#' @importFrom splus2R  peaks
#' @importFrom dplyr arrange
#' @importFrom dplyr relocate
#' @importFrom stats hclust
#' @importFrom stats cutree
#' @importFrom stats quantile
#'
# find peaks function ####
.find_peaks <- function(x, ignore_threshold = 0, span = 3, strict = TRUE, na.rm = FALSE){
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
      ifelse(x - min_x > scaled_threshold, pks, FALSE)
    } else {
      ifelse(max_x - x > scaled_threshold, pks, FALSE)
    }
  }
}


# Find sample name function ####

# experimental function, could be nice but doesn't currently support getting experiement name

#.sampleNameFinder <- function(path = Path_to_imzml_file){
#
#  Sample_name <- strsplit(path, split = "/", fixed = TRUE)[[1]]
#
#  for (i in length(Sample_name)){
#    if (grepl(".imzML", Sample_name[i], fixed=TRUE)){
#      temp_name <- (Sample_name[i])
#      temp_name <- strsplit(temp_name, split = ".", fixed = TRUE)[[1]]
#      Sample_name <- temp_name[1]
#    }else{
#      print("No .imzML detected")
#    }
#  }
#
#  return(Sample_name)
#}


# Get the name of the experiment and the sample (applicable only for JA's directory structure)
.sampleNameFinder <- function(path = Path_to_imzml_file){

  sample_name <- strsplit(path, split = "/", fixed = TRUE)[[1]]
  h <- 1
  cur_folder <- sample_name[h]
  if(cur_folder == "experiments"){
    print("Correct on first iteration")
    print(paste("Current folder is ", cur_folder, sep = ""))
    h = h+1
  }else{while(cur_folder != "experiments"){
    cur_folder <- sample_name[h]
    print(paste("Current folder is ", cur_folder, sep = ""))
    h = h+1
    if (cur_folder == "experiments"){
      print(paste("The current folder is ", cur_folder, ".", sep = ""))
    }
  }}

  experiment_name <- sample_name[h]
  print(paste("The experiment is ", experiment_name, ".", sep = ""))
  sample_name <- sample_name[length(sample_name)]
  sample_name <- strsplit(sample_name ,split = ".", fixed = TRUE)[[1]]
  sample_name <- sample_name[1]
  print(paste("The sample is ", sample_name, ".", sep = ""))

}



# Marker panel cleaning function ####

.cleanPanel <- function(panel){

  # 1. ensure column names are correct

  if (is.numeric(panel$Name) == T){
    colnames(panel) <- c("FeatureMass","Name")
  }

  # rearrange column order so Name is first column
  dplyr::relocate(panel, "Name", .before = "FeatureMass")

  # 2. sort by mass tag size

  panel <- dplyr::arrange(panel, panel$FeatureMass)

  # 3. clean marker names

  for(a in seq_along(panel$Name)){

    if (grepl("+", panel$Name[a], fixed=TRUE)){
      new_string <- gsub(" ", "", panel$Name[a])
      panel$Name[a] <- new_string

    }
    else{

      # split strings with spaces into list with individual strings as elements
      new_string <- gsub(" ", "", panel$Name[a])
      # replace fullstops with underscores
      new_string <- gsub("-", ".", new_string, fixed = TRUE)
      # if there is a dash at the end of the name, remove it
      if(endsWith(new_string, "-")){

        new_string <- substr(new_string,1, nchar(new_string)-1)

      }
      # replace names in peakAnnotation
      panel$Name[a] <- new_string

    }

    # remove slashes
    if (grepl("/", panel$Name[a], fixed=TRUE)){

      new_string <- gsub("/", "", panel$Name[a])
      panel$Name[a] <- new_string

    }

  }

  # 4. rearrange column order so Name is first column
  panel <- relocate(panel, "Name", .before = "FeatureMass")

  return(panel)

}

# If the x coordinates don't begin at 1, adjust the coordinates
.correctCoordinates <- function(coords){

  if(!min(coords$x) == 1){

    coords$x <- coords$x - (min(coords$x) - 1)
    coords$y <- coords$y - (min(coords$y) - 1)

  }

  return(coords)
}


# Convert_to_mz_scale function ####

.convert_to_mz_scale <- function(x,range_peaks, N_features) {

  scale_vector = base::seq(range_peaks[1], range_peaks[2], length.out = N_features )
  return(scale_vector[x])

}

# define function for hierarchical clustering with complete linkage ####

.hc_single_linkage_function <- function(x, threshold_height = 2) {

  sub_clustering <- 1

  # if there is more than one element in the region, run the distance function
  if (length(x) > 1) {
    dist_matrix <- dist(x)
    hc_clustering <- hclust(dist_matrix, method = "single")
    sub_clustering <- cutree(hc_clustering, h = threshold_height)
  }

  return(sub_clustering)

}


# which_max_modified function ####

.which_max_modified <- function(x) {

  if (sum(x) == 0) {
    y <- NA
  }
  else {

    x[x == 0] <- NA
    y <- which.min(x)

  }

  return(y)

}


# One dimensional otsu thresholding  ####

.otsu_thresholding = function(x, number_bins = 100) {

  list_bin = quantile(x,base::seq(from = 0, to = 1, length.out = number_bins))
  intravariance_vector = c()

  for (k in 1:number_bins) {

    threshold_temp <- list_bin[k]
    s <- length(x[x < threshold_temp]) * var(x[x < threshold_temp]) + length(x[x > threshold_temp]) * var(x[x > threshold_temp])
    intravariance_vector <- c(intravariance_vector, s)

  }

  selected_values <- list_bin[which.min(intravariance_vector)]

  return(selected_values)
}


# Strings to colours ####
.string.to.colors = function (string, colors = NULL)
{
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





