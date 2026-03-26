#' @importFrom stats var
NULL

#############################  Helper functions ###################################

# Marker panel cleaning function ####

.cleanPanel <- function(panel){

  # 1. ensure column names are correct

  if (is.numeric(panel$Name)){
    colnames(panel) <- c("FeatureMass","Name")
  }

  # 2. sort by mass tag size

  panel <- dplyr::arrange(panel, panel$FeatureMass)

  # 3. clean marker names
  panel$OriginalName <- panel$Name   # record pre-clean names

  for(a in seq_along(panel$Name)){

    if (grepl("+", panel$Name[a], fixed=TRUE)){
      new_string <- gsub(" ", "", panel$Name[a])
      panel$Name[a] <- new_string

    }
    else{

      # split strings with spaces into list with individual strings as elements
      new_string <- gsub(" ", "", panel$Name[a])
      # replace dashes with dots
      new_string <- gsub("-", ".", new_string, fixed = TRUE)
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
  panel <- relocate(panel, c("Name", "OriginalName"), .before = "FeatureMass")

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


# KNN spatial weight matrix ####
# Build a symmetric sparse n x n weight matrix of squared L2 distances to k nearest
# neighbours. Symmetric: if j is a KNN of i OR i is a KNN of j, both [i,j] and [j,i]
# are set — matching the output of N2R::Knn(indexType = "L2").
.knn_weight_matrix <- function(coords, k) {
  n  <- nrow(coords)
  d2 <- as.matrix(dist(coords))^2
  diag(d2) <- Inf  # exclude self

  # asymmetric KNN adjacency: knn_adj[i,j] = TRUE if j is among k nearest of i
  knn_adj <- matrix(FALSE, n, n)
  for (i in seq_len(n)) {
    knn_adj[i, order(d2[i, ])[seq_len(k)]] <- TRUE
  }

  # symmetrize: include (i,j) if j is KNN of i OR i is KNN of j
  knn_sym <- knn_adj | t(knn_adj)

  # restore diagonal to 0 for value lookup
  diag(d2) <- 0

  idx <- which(knn_sym, arr.ind = TRUE)
  Matrix::sparseMatrix(i = idx[, 1], j = idx[, 2], x = d2[idx], dims = c(n, n))
}


# One dimensional otsu thresholding  ####
## x is logTIC
## number_bins is the number of histogram bins
## safe_var returns 0 for partitions with fewer than 2 elements to avoid NA from var() on length-0 or length-1 vectors
.otsu_thresholding = function(x, number_bins = 100) {

  list_bin = quantile(x,base::seq(from = 0, to = 1, length.out = number_bins))
  intravariance_vector <- numeric(number_bins)

  for (k in seq_len(number_bins)) {

    threshold_temp <- list_bin[k]
    safe_var <- function(v) if (length(v) < 2) 0 else var(v)
    s <- length(x[x < threshold_temp]) * safe_var(x[x < threshold_temp]) +
         length(x[x > threshold_temp]) * safe_var(x[x > threshold_temp])
    intravariance_vector[k] <- s

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

# helper: return the candidate mz with the smallest absolute distance to expected_mz
.choose_closest <- function(cand_mz, expected_mz) {
  cand_mz[which.min(abs(cand_mz - expected_mz))]
}

# remove duplicated metapeaks for a targeted MSI experiment.
# x: output from assignMetapeaks (a list with IntensityDF, CorrespondenceMatrix, SpatialCoords).
# shift: direction of expected mass shift ("left", "right", or "closest").
.removeDuplicates <- function(x, shift){

  if (!shift %in% c("left", "right", "closest")) {
    stop("'shift' must be \"left\", \"right\", or \"closest\".")
  }

  cur_df <- x$IntensityDF
  cur_correspondence <- x$CorrespondenceMatrix

  # get the correspondence matrix of just the duplicated values
  dup_tags <- unique(cur_correspondence$expected_mz_location)[which(table(cur_correspondence$expected_mz_location) > 1)]
  dup_correspondence <- cur_correspondence[which(cur_correspondence$expected_mz_location %in% dup_tags), ]

  # loop through all duplicates and choose the correct peak for each
  duplicate_expected_mzs <- unique(dup_correspondence$expected_mz_location)
  correct_peaks <- numeric(length(duplicate_expected_mzs))


  for(i in duplicate_expected_mzs){

    idx <- which(duplicate_expected_mzs == i)

    # extract cur pairing
    cur_pairing <- dup_correspondence[dup_correspondence$expected_mz_location == i, ]

    if (shift == "closest") {

      correct_peak <- .choose_closest(cur_pairing$mz_location, unique(cur_pairing$expected_mz_location))

    } else {

      compare_fn <- if (shift == "right") `>` else `<`
      select_fn  <- if (shift == "right") min else max

      correct_shifted_mzs <- cur_pairing$mz_location[
        compare_fn(cur_pairing$mz_location, unique(cur_pairing$expected_mz_location))
      ]

      if (length(correct_shifted_mzs) == 0) {
        warning(paste0("No metapeak found to the ", shift, " of expected mz ", i, ". Falling back to closest peak."))
        correct_shifted_mzs <- cur_pairing$mz_location
      }

      correct_peak <- select_fn(correct_shifted_mzs)

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
  bin_vector <- logical(ncol(cur_df))
  for (i in seq_len(ncol(cur_df))){

    ## identify markers with duplicate in the name
    needle <- "duplicate"
    haystack <- colnames(cur_df)[i]
    ## find needle in haystack
    cur_val <- grepl(needle, haystack, fixed = TRUE)
    ## add to binary vector
    bin_vector[i] <- cur_val

  }

  # update marker names in new correspondence matrix
  updated_correspondence$marker <- colnames(cur_df)[!bin_vector]

  # update intensity dataframe
  updated_df <- cur_df[, keep_idx]
  colnames(updated_df) <- colnames(cur_df)[!bin_vector]

  ## Update relevant fields in the original structure
  updated_sample <- x
  updated_sample$IntensityDF <- updated_df
  updated_sample$CorrespondenceMatrix <- updated_correspondence
  updated_sample$FilteredDF <- x$FilteredDF[, colnames(updated_df)]

  return(updated_sample)

}

