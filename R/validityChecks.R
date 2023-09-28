.valid.preProcess <- function(x, cores) {

  if (!is(x, "MSContinuousImagingExperiment")) {
    stop("'x' should be of class 'MSContinuousImagingExperiment'")
  }

  if (length(cores) > 1) {
    stop("'cores' should be a single numeric.")
  }

  if (!is.numeric(cores)) {
    stop("'cores' should be a single numeric.")
  }

}
