# removeDuplicates
test_that("removeDuplicates works", {

  # minimal fixture — two candidates for one expected mz (100) plus one
  # non-duplicate marker. The duplicate column uses the ".duplicate" suffix
  # matching real assignMetapeaks output. A non-duplicate column is required
  # because the function uses data.frame subsetting that drops dimensions when
  # only one column is kept.
  make_x <- function(mz_locations) {
    dup_cols  <- c(paste0("m", mz_locations[1]),
                   paste0("m", mz_locations[1], ".duplicate"))
    all_cols  <- c(dup_cols, "m200")
    list(
      IntensityDF = as.data.frame(
        setNames(replicate(3L, rnorm(10), simplify = FALSE), all_cols)
      ),
      CorrespondenceMatrix = data.frame(
        marker                = all_cols,
        expected_mz_location  = c(100, 100, 200),
        mz_location           = c(mz_locations, 200),
        stringsAsFactors      = FALSE
      ),
      FilteredDF = as.data.frame(
        setNames(replicate(3L, rnorm(10), simplify = FALSE), all_cols)
      )
    )
  }

  x_lr <- make_x(c(98, 103))   # one left, one right of expected mz 100

  # "closest" picks the nearer candidate (98, distance 2) over the farther (103, distance 3)
  result_closest <- removeDuplicates(x_lr, shift = "closest")
  expect_equal(ncol(result_closest$IntensityDF), 2L)
  expect_equal(result_closest$CorrespondenceMatrix$mz_location, c(98, 200))

  # "closest" runs without a warning even when all candidates are on one side
  x_right_only <- make_x(c(102, 105))
  expect_silent(removeDuplicates(x_right_only, shift = "closest"))

  x_left_only <- make_x(c(95, 97))
  expect_silent(removeDuplicates(x_left_only, shift = "closest"))

  # "right" still selects the nearest right candidate
  result_right <- removeDuplicates(x_lr, shift = "right")
  expect_equal(result_right$CorrespondenceMatrix$mz_location, c(103, 200))

  # "left" still selects the nearest left candidate
  result_left <- removeDuplicates(x_lr, shift = "left")
  expect_equal(result_left$CorrespondenceMatrix$mz_location, c(98, 200))

  # "right" warns when no candidates exist to the right
  x_left_only2 <- make_x(c(95, 97))
  expect_warning(removeDuplicates(x_left_only2, shift = "right"))

  # invalid shift errors
  expect_error(removeDuplicates(x_lr, shift = "diagonal"))

})
