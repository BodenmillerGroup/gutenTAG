computeGMM <- function(x, hist = FALSE){

  targeted_intensity <- x$IntensityDF
  untargeted_intensity <- x$Untargeted$UntargetedIntensity

  tic_intensity <- rowSums(untargeted_intensity)

  #pixelsSummarised <- as.data.frame(cbind(tic_intensity, mean_intensity))

  tic_intensity[tic_intensity == 0] <- NA
  tic_intensity <- na.omit(tic_intensity)

  tic_threshold <- .otsu_thresholding(log10(tic_intensity))

  columns <- c("MarkerID", "SeparationScore", "PositivePixels", "MeanNoiseIntensity", "MeanSignalIntensity")
  gmmTable <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(gmmTable) <- columns
  markerID <- colnames(targeted_intensity)

  # need to run this outside the pdf so that variables are saved
  for (k in 1:ncol(targeted_intensity)) {

    # threshold on otsu threshold of TIC
    x <- targeted_intensity[, k][tic_intensity > tic_threshold]

    # no thresholding
    ##x <- targeted_intensity[, k]

    # remove zeros
    x <- x[x > 0]

    # optimal number of Gaussians
    ICL_score <- mclust::mclustICL(data = log10(1 + x), G = 1:2, modelNames = c("V"))
    ICL_score <- c(ICL_score[[1]], ICL_score[[2]])
    #ICL_score <- c(ICL_score[[1]], ICL_score[[2]], ICL_score[[3]], ICL_score[[4]], ICL_score[[5]])

    N_gaussian_curves <- as.numeric(which.max((ICL_score)))
    gmm_model <- Mclust(data = log10(1 + x), G = N_gaussian_curves, modelNames = "V")

    density <- densityMclust(x, plot = FALSE) # same as Mclust but has additional parameter density
    density <- data.frame(density[17]) # convert to dataframe to get max value
    maxDensity <- max(density)

    # Define separability score

    mean1 <- gmm_model$parameters$mean[1]
    mean2 <- gmm_model$parameters$mean[2]

    sd1 <- sqrt(gmm_model$parameters$variance$sigmasq[1])
    sd2 <- sqrt(gmm_model$parameters$variance$sigmasq[2])

    # calculate separation score
    meanDiff <- abs(mean1 - mean2)
    SNR <- mean2 / mean1

    # proportion positive pixels
    positivePixels <- length(which(gmm_model$classification == 2))
    proportionPositive <- positivePixels / length(gmm_model$classification)
    pi2 <- proportionPositive
    pi1 <- 1 - pi2

    denominator <- sqrt((sd1^2 / pi1) + (sd2^2 / pi2))
    separationScore <- meanDiff / denominator

    # Create table for export

    # add markerID and separation score to dataframe
    gmmTable[k, 1] <- as.character(markerID[k])
    gmmTable[k, 2] <- separationScore
    gmmTable[k, 3] <- proportionPositive
    gmmTable[k, 4] <- mean1
    gmmTable[k, 5] <- mean2

    # plot GMM
    if(hist == TRUE){

      hist(log10(x + 1),
           breaks = 50,
           freq = FALSE,
           main = paste(colnames(targeted_intensity[1, ])[k]),
           xlab = paste(colnames(targeted_intensity[1, ])[k],"Intensity"),
           cex.main = 1.5,
           cex.lab = 1.2) # blank xlab for next line
      for (j in 1:N_gaussian_curves) {
        curve(dnorm(x,
                    mean = gmm_model$parameters$mean[j],
                    sd = sqrt(gmm_model$parameter$variance$sigmasq[j])) * gmm_model$parameters$pro[j],
              add = TRUE,
              lwd = 2,
              col = rainbow(N_gaussian_curves)[j])
        text(x = max(gmm_model$data), y = (maxDensity - (0.35 * maxDensity)), # coordinates of the text
             labels = paste("Separation Score: ", round(separationScore, 3)),
             adj = 1,
             cex = 1,
             font = 1)
      }

    }

  }

  gmmTable <- gmmTable %>% relocate(MarkerID)
  gmmTable[is.na(gmmTable)] <- 0

  return(list(gmmTable, gmm_model))

}

