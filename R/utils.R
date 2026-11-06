#' @importFrom stats var dist
NULL

#############################  Helper functions ###################################

# Marker panel cleaning function ####

.cleanPanel <- function(panel) {
    # 1. ensure column names are correct

    if (is.numeric(panel$Name)) {
        colnames(panel) <- c("FeatureMass", "Name")
    }

    # 2. sort by mass tag size

    panel <- dplyr::arrange(panel, panel$FeatureMass)

    # 3. clean marker names
    panel$OriginalName <- panel$Name # record pre-clean names

    # remove spaces and slashes; replace dashes with dots unless the name contains "+"
    has_plus <- grepl("+", panel$Name, fixed = TRUE)
    panel$Name <- gsub(" ", "", panel$Name, fixed = TRUE)
    panel$Name[!has_plus] <- gsub("-", ".", panel$Name[!has_plus], fixed = TRUE)
    panel$Name <- gsub("/", "", panel$Name, fixed = TRUE)

    # 4. rearrange column order so Name is first column
    panel <- relocate(panel, c("Name", "OriginalName"), .before = "FeatureMass")

    return(panel)
}

# If the x coordinates don't begin at 1, adjust the coordinates
.translate_coordinates <- function(coords) {
    if (!min(coords$x) == 1) {
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
    n <- nrow(coords)
    d2 <- as.matrix(dist(coords))^2
    diag(d2) <- Inf # exclude self

    # asymmetric KNN adjacency: knn_adj[i,j] = TRUE if j is among k nearest of i
    nn <- vapply(seq_len(n), function(i) order(d2[i, ])[seq_len(k)], integer(k))
    knn_adj <- matrix(FALSE, n, n)
    knn_adj[cbind(rep(seq_len(n), each = k), as.vector(nn))] <- TRUE

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
.otsu_thresholding <- function(x, number_bins = 100) {
    list_bin <- quantile(x, base::seq(from = 0, to = 1, length.out = number_bins))
    safe_var <- function(v) if (length(v) < 2) 0 else var(v)

    intravariance_vector <- vapply(list_bin, function(threshold_temp) {
        length(x[x < threshold_temp]) * safe_var(x[x < threshold_temp]) +
            length(x[x > threshold_temp]) * safe_var(x[x > threshold_temp])
    }, numeric(1))

    return(list_bin[which.min(intravariance_vector)])
}


# utility function for making matrix from flattened dataframe
.curateMatrix <- function(dataframe, channel, coords) {
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
