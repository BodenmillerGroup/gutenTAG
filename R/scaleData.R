#' Scale processed data.
#'
#' @param x Output of getIntensityDF function
#' @param method Choose according to which statistic to scale the data. Currently supported are corrected standard deviation as "corsd" and Geary's C score as "geary".
#'
#' @return Scaled data
#'
#' @importFrom dplyr relocate
#' @importFrom stats sd
#' @importFrom stats var

#' @export
scaleData <- function(x, method = "corsd"){

  # extract dataframe and correspondence matrix
  df <- x$IntensityDF
  correspondence <- x$CorrespondenceMatrix

  # Remove all un-annotated peaks from correspondence matrix
  correspondence_matrix_targeted <- na.omit(correspondence)
  correspondence_matrix_targeted <- dplyr::relocate(correspondence_matrix_targeted, "Annotation_peaks") # change to Annotation_peaks or Marker accordingly

  # specify channel
  scaled_data <- c()
  for (i in 1:ncol(df)){

    channel <- df[, i]
    name <- colnames(df)[i]

    var_x <- stats::var(channel)
    sd_x <- stats::sd(channel)
    x_prime <- channel / sd_x

    # compute alpha based on method
    if(method == "corsd"){
      alpha <- correspondence["Corrected_sd"][i,]
    }

    if(method == "geary"){
      geary <- correspondence["Geary"][i,]
      alpha <- 1 - geary
    }

    # add method for variogram, e.g. distance to 50% of semivariance ???

    # calculate x_prime_prime
    if(alpha > 0){
      x_prime_prime <- x_prime * sqrt(alpha)
    }

    # if alpha is negative
    if(alpha < 0){
      x_prime_prime <- x_prime * -sqrt(abs(alpha))
    }

    # bind columns
    scaled_data <- data.frame(cbind(scaled_data, x_prime_prime))
    colnames(scaled_data)[i] <- name

  }

  return(scaled_data)

}
