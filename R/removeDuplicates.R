#' Remove Duplicated Metapeaks
#'
#' @description When multiple metapeaks map to the same expected marker mass,
#'   this function resolves the ambiguity by keeping a single metapeak per
#'   expected mass. The metapeak to keep is chosen according to the expected
#'   direction of mass shift via the \code{shift} argument. Markers that did not
#'   receive a metapeak (\code{NA} \code{mz_location}) and non-duplicated markers
#'   are retained unchanged.
#'
#' @param x A list. Output of \code{\link{assignMetapeaks}}, containing
#'   \code{IntensityDF}, \code{CorrespondenceMatrix}, \code{SpatialCoords} and
#'   \code{FilteredDF}.
#' @param shift Character; the expected direction of mass shift used to select
#'   the correct metapeak among duplicates. One of \code{"left"} (lower m/z than
#'   expected), \code{"right"} (higher m/z than expected), or \code{"closest"}
#'   (nearest m/z to expected). When \code{"left"} or \code{"right"} is requested
#'   but no candidate exists on that side, the function falls back to the closest
#'   candidate with a warning.
#'
#' @return A list with the same structure as \code{x}, with \code{IntensityDF},
#'   \code{CorrespondenceMatrix} and \code{FilteredDF} updated to retain a single
#'   metapeak per expected marker mass.
#'
#' @export
#'
#' @examples
#' rdata_path <- system.file("extdata/Example_data/Example_processed.Rdata", package = "gutenTAG")
#' load(rdata_path)
#' removeDuplicates(x = results$processed, shift = "closest")
removeDuplicates <- function(x, shift){

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
        warning("No metapeak found to the ", shift, " of expected mz ", i, ". Falling back to closest peak.", call. = FALSE)
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
