test_that("computeUMAP works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "maldipackage")
  panel_path <- system.file("extdata/ref_list.csv", package = "maldipackage")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)
  peaks <- peakDetection(pre, core = 2)
  metapeaks <- generateMetapeaks(peaks)
  processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)

  # test that input is correct
  expect_type(processed, "list")

  # compute Geary's C
  expect_silent(cur_test <- computeUMAP(processed$IntensityDF, seed = 3))

  # test that output is correct is correct
  expect_equal(cur_test[[1]][1:10], c(6.43675610702485, 5.82524153869599, 6.00166508834809, 0.747123631648719,
                                      -0.472601381130517, -5.86854603607208, -6.08983853179961, -3.97023156005889,
                                      -3.72775032836944, -5.27523997146636))


  # test if character in first argument throws error
  expect_error(computeUMAP("processed", seed = 345))
  expect_error(computeUMAP(processed, seed = "345"))




})
