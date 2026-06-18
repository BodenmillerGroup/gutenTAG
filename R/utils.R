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



# Rook-adjacency spatial weight matrix ####
# MALDI pixels sit on a regular integer lattice, so the k = 4 "rook" neighbourhood
# of a pixel at (x, y) is simply its up/down/left/right grid neighbours:
# (x+1, y), (x-1, y), (x, y+1), (x, y-1). We build the symmetric sparse adjacency
# DIRECTLY from the integer coordinates in O(n) time and memory — never forming a
# pairwise distance matrix. The old implementation materialised a dense n x n
# distance matrix via as.matrix(dist(coords)), which overflows R's 32-bit vector
# index at n >= 65,536 pixels (and OOMs well before that on real samples).
#
# Edges carry UNIT weights (value 1). Geary's C is invariant to global scaling of
# the weight matrix, so on a uniform-spacing grid unit weights are equivalent to
# the squared-distance weights the old code used — but exact and simpler.
#
# Edge/hole behaviour (intentional change): a pixel on a tissue edge or beside a
# hole is connected ONLY to the rook neighbours that actually exist (2 at a corner,
# 3 along an edge). The old force-k-nearest code always picked 4 neighbours, pulling
# in diagonal pixels at boundaries; the rook graph is the geometrically correct
# neighbourhood and slightly shifts edge-pixel Geary's C scores.
#
# `k` is retained for signature compatibility; k = 4 denotes the rook case.
.knn_weight_matrix <- function(coords, k) {
  n <- nrow(coords)
  xs <- coords[, 1]
  ys <- coords[, 2]

  # integer hash key per pixel: x + y * stride, with stride large enough that no
  # two distinct (x, y) collide. match() on these keys maps a coordinate back to
  # its pixel index.
  stride <- max(xs) + 1
  keys <- xs + ys * stride

  # the four rook offsets; build edges by matching each shifted coordinate back to
  # an existing pixel. Non-NA matches are the rook neighbours that actually exist.
  offsets <- list(c(1, 0), c(-1, 0), c(0, 1), c(0, -1))

  from <- integer(0)
  to   <- integer(0)
  for (off in offsets) {
    neigh_keys <- (xs + off[1]) + (ys + off[2]) * stride
    j <- match(neigh_keys, keys)
    hit <- which(!is.na(j))
    from <- c(from, hit)
    to   <- c(to, j[hit])
  }

  # The four-offset scan already yields each undirected edge in both directions
  # (the +1 scan from pixel i produces (i, i+1); the -1 scan from i+1 produces
  # (i+1, i)), so the assembled matrix is symmetric with a zero diagonal.
  Matrix::sparseMatrix(i = from, j = to, x = rep(1, length(from)), dims = c(n, n))
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