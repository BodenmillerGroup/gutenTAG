#############################  Helper functions ###################################

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
.translate_coordinates <- function(coords){

  if(!min(coords$x) == 1){

    coords$x <- coords$x - (min(coords$x) - 1)
    coords$y <- coords$y - (min(coords$y) - 1)

  }

  return(coords)
}


# which_min_ignore_zero function ####
# returns which.min ignoring 0 values
.which_min_ignore_zero <- function(x) {

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


# utility function for making matrix from flattened dataframe
.curateMatrix <- function(dataframe, channel, coords){

  Matrix_image = matrix(0, ncol = max(coords$y), nrow=max(coords$x)) # initialise matrix of correct shape
  x = dataframe[,channel] # isolate single channel
  Matrix_image[as.matrix(coords)] = x # input channel values into matrix
  return(Matrix_image)

}

# remove duplicated metapeaks for a targeted MSI experiment.
## Specify which direction of the expected mass you expect the mass shift to occur in. removeDuplicates will choose the closest metapeak in that direction.
.removeDuplicates <- function(sample, shift){

  cur_df <- sample$processed$IntensityDF
  cur_correspondence <- sample$processed$CorrespondenceMatrix

  # get the correspondence matrix of just the duplicated values
  dup_tags <- unique(cur_correspondence$expected_mz_location)[which(table(cur_correspondence$expected_mz_location) > 1)]
  dup_correspondence <- cur_correspondence[which(cur_correspondence$expected_mz_location %in% dup_tags), ]

  # loop through all duplicates and choose the correct peak for each
  duplicate_expected_mzs <- unique(dup_correspondence$expected_mz_location)
  correct_peaks <- c()

  for(i in duplicate_expected_mzs){

    idx <- which(duplicate_expected_mzs == i)

    # extract cur pairing
    cur_pairing <- dup_correspondence[dup_correspondence$expected_mz_location == i, ]

    # if you expect mass shift to be to the right of the expected mass
    if(expected_shift == "right"){
      # first select masses that are greater than the expected value
      correct_shifted_mzs <- cur_pairing$mz_location[cur_pairing$mz_location > unique(cur_pairing$expected_mz_location)]
      # choose the closest one (moot if there is only 1, most common case)
      correct_peak <- min(correct_shifted_mzs)
    }

    # if you expect mass shift to be to the right of the expected mass
    if(expected_shift == "left"){
      # first select masses that are less than the expected value
      correct_shifted_mzs <- cur_pairing$mz_location[cur_pairing$mz_location < unique(cur_pairing$expected_mz_location)]
      # choose the closest one (moot if there is only 1, most common case)
      correct_peak <- max(correct_shifted_mzs)
    }

    correct_peaks[idx] <- correct_peak

  }

  # remove all columns in original df that are NOT in the keep list

  ## get non-duplicated tags
  non_dup_tags <- cur_correspondence$mz_location[!cur_correspondence$expected_mz_location %in% dup_tags]
  na_tags_idx <- which(is.na(cur_correspondence$mz_location))

  # these are the mz locations of the tags to keep (non-duplicated peaks plus the correct duplicated ones)
  keep_peaks <- sort(c(correct_peaks, non_dup_tags))
  # get the indexes of peaks to keep -- IMPORTANT: we move to index rather than mz_location so that we can include the NA-valued tags, which would be excluded if we worked in mz_location (since these markers didn't get a metapeak)
  keep_idx <- sort(c(which(cur_correspondence$mz_location %in% keep_peaks), na_tags_idx))

  # correspondence matrix of only the peaks to keep
  updated_correspondence <- cur_correspondence[keep_idx, ]

  # get the marker names only
  ## initialse binary vector (TRUE/FALSE if marker name is duplicate)
  bin_vector <- c()
  for (i in 1:length(colnames(cur_df))){

    ## identify markers with duplicate in the name
    needle <- "duplicate"
    haystack <- colnames(cur_df)[i]
    ## find needle in haystack
    cur_val <- grepl(needle, haystack, fixed = TRUE)
    ## add to binary vector
    bin_vector <- c(bin_vector, cur_val)

  }

  # update marker names in new correspondence matrix
  updated_correspondence$marker <- colnames(cur_df)[!bin_vector]

  # update intensity dataframe
  updated_df <- cur_df[, keep_idx]
  colnames(updated_df) <- colnames(cur_df)[!bin_vector]

  ## Create new sample list
  updated_sample <- list("IntensityDF" = updated_df,
                         "CorrespondenceMatrix" = updated_correspondence,
                         "SpatialCoords" = sample$processed$SpatialCoords)

  return(updated_sample)

}

