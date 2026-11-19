test_that("one to five mediator timepoints reproduce professor functions", {
  for (K in 1:5) {
    args <- make_path_case(K)
    for (scale in c("hazard", "sp")) {
      actual <- run_with_raw_draws(args, scale, seed = 7800L + K)
      reference <- run_reference(args, scale, seed = 7800L + K)
      expect_equal(unname(actual$raw), unname(reference$raw), tolerance = 1e-9,
                   info = paste(K, scale))
      expect_identical(actual$rng, reference$rng, info = paste(K, scale))
      expect_equal(unname(actual$raw[, "TE"]),
                   unname(rowSums(actual$raw[, c("DE", "IEM", "IED"), drop = FALSE])),
                   tolerance = 1e-9)
    }
  }
})

test_that("no competing-event models gives zero IED", {
  for (K in c(1L, 3L, 5L)) {
    args <- make_path_case(K, competing = FALSE)
    for (scale in c("hazard", "sp")) {
      actual <- run_with_raw_draws(args, scale, seed = 8100L + K)
      reference <- run_reference(args, scale, seed = 8100L + K)
      expect_equal(unname(actual$raw), unname(reference$raw), tolerance = 1e-9)
      expect_equal(as.numeric(actual$raw[, "IED"]), rep(0, 2))
    }
  }
})
