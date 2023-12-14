test_that("nmf works",{

  # get input data
  path <- system.file("extdata/Example_data.imzML", package = "maldipackage")
  panel_path <- system.file("extdata/ref_list.csv", package = "maldipackage")
  panel <- readPanel(path = panel_path)
  raw <- readMSIData(path)
  pre <- preProcess(raw, cores = 2)
  peaks <- peakDetection(pre, core = 2)
  metapeaks <- generateMetapeaks(peaks)
  processed <- assignMetapeaks(metapeaks, pre = pre, refList = panel)

  # test that processed data is a list
  expect_type(processed, "list")

  # test that entries into intensity dataframe are correct
  expect_equal(processed$IntensityDF[1:20, ][[1]], c(7.82970978419113, 7.46634943048983, 7.32187611962462, 6.54747422274114,
                                                     7.49722025536448, 3.55251445305979, 5.75086256527666, 3.51621186218398,
                                                     5.27250246930296, 3.86230533109536, 3.44947707629045, 9.10591335737913,
                                                     8.21399509393651, 7.16647250531501, 6.75161308428109, 7.91584977381878,
                                                    6.24461353570202, 9.09662347324635, 6.15357788356246, 6.6490763080287))

  # make CytoImageList object
  cur_test <- computeNMF(processed$IntensityDF[-c(3,4,6)], comp = 5, seed = 123, cntr = FALSE)

  # test that class of s3 object is correct
  expect_type(cur_test, "list")

  # test that feature factor matrix is correct
  expect_equal(cur_test$w[1:5, 1:5], structure(c(0.00754686352746144, 0.00500196233683792, 0.00509466781693867,
                                                 0.00220770799101638, 0.00234366230267723, 0.00794985285653198,
                                                 0.00937889777374931, 0.00818138983906627, 0.00586658147447449,
                                                 0.00373223382282122, 0.000622021022130529, 0.000198790584397358,
                                                 0, 0.00394124723107695, 0.00609020357489004, 0.00723705140133675,
                                                 0.00693503695553305, 0.010405481045203, 0.00659066837779087,
                                                 0.00187423686105082, 0.00033438781408952, 0.000393353280272339,
                                                 0.00138814944463373, 0, 0.000481081266290184), dim = c(5L, 5L)))

  # test that factor sample matrix is correct
  expect_equal(cur_test$h, structure(c(0.0357593781775371, 0.165892672745083, 0.250598574534461,
                                       0, 0.0478804521800722, 0.0897986163349986, 0.135804189600256,
                                       0.199901502792097, 0.00964766909238889, 0.0475035990943403, 0.107432912815836,
                                       0.141270263238217, 0.087400580820735, 0.126993396767021, 0.0219183293839359,
                                       0.150406501214626, 0.0968890909940774, 0, 0.178949059411205,
                                       0.0146650774999818, 0.245993487291453, 0, 0.0356335582948698,
                                       0, 0.0290842156933556, 0.176884432742191, 0.0935126340971572,
                                       0.0208226028205408, 0.0755063550597264, 0.0337796743222781, 0.00858704990853088,
                                       0.0580426026447249, 0.0474168061146664, 0, 0.507915367639716,
                                       0.0233335137350519, 0, 0.135557158969937, 0.390430869339034,
                                       0.191561899644196, 0.07405243643506, 0.0929302038316005, 0.222669215652692,
                                       0.084830381936598, 0.0631053768973551, 0.0877516713447151, 0.215658342848884,
                                       0, 0.133642268394027, 0.042586007644769), dim = c(5L, 10L)))


  # test if arguments are in correct format
  expect_error(computeNMF(x = "test", comp = 5))
  expect_error(computeNMF(x = processed$IntensityDF, comp = "5"))


})
