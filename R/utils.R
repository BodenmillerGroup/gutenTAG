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

  Matrix_image <- matrix(0, ncol = max(coords$y), nrow = max(coords$x))
  x <- dataframe[, channel]
  Matrix_image[as.matrix(coords)] <- x
  return(Matrix_image)

}

# helper: return the candidate mz with the smallest absolute distance to expected_mz
.choose_closest <- function(cand_mz, expected_mz) {
  cand_mz[which.min(abs(cand_mz - expected_mz))]
}

# NOTE: removeDuplicates() was moved to its own file R/removeDuplicates.R and
# exported. Its private helper .choose_closest remains above.