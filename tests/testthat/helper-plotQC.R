# Shared fixtures for plotQC tests.
# Sourced once per testthat session; objects are available in all test files.

.rdata_path <- system.file("extdata/Example_processed.Rdata", package = "gutenTAG")
load(.rdata_path)  # loads 'results' with $processed, $metapeaks, $panel

qc_panel     <- results$panel
qc_metapeaks <- results$metapeaks
qc_processed <- results$processed

qc_processed_geary <- computeGearysC(qc_processed, verbose = FALSE,
                                     update_correspondence = TRUE)
qc_processed_snr   <- computeSNR(qc_processed_geary, update_correspondence = TRUE,
                                 package = "mclust")
