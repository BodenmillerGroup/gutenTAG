test_that("computeGearysC works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "maldipackage")
  panel_path <- system.file("extdata/ref_list.csv", package = "maldipackage")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)
  peaks <- peakDetection(pre, core = 2)
  metapeaks <- generateMetapeaks(peaks)
  processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)


  # compute Geary's C
  expect_silent(cur_test <- computeGearysC(processed, verbose = FALSE))

  # test that output is correct is correct
  expect_equal(cur_test, c(0.710039425480746, 0.483790051545452, 0.435088747635843, 0.434696662055006,
                           0.360060956047036, 0.301997161578108, 0.556339728384072, 0.759333786571679,
                           0.877932000280288, 0.391559739360737))


  # test if character in first argument throws error
  expect_error(computeGearysC("processed", verbose = FALSE, update_correspondence = FALSE))
  expect_error(computeGearysC(processed, verbose = "FALSE", update_correspondence = FALSE))


})
