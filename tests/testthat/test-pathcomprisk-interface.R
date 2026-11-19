test_that("interface documents the current implementation", {
  expect_identical(formals(cmest_pathcomprisk)$ncores, 1L)
  expect_false("seed" %in% names(formals(cmest_pathcomprisk)))
  args <- make_path_case(3)
  out <- run_with_raw_draws(args, "hazard", 111, nboot = 1L)
  expect_identical(names(out$summary), c("DE", "IEM", "IED", "TE"))
  expect_true(all(vapply(out$summary, is.character, logical(1))))
  expect_identical(getvarnames(y ~ x + z)$xvar, c("x", "z"))
  expect_identical(getvarnames(y ~ ., data.frame(y = 1:3, x = 4:6))$xvar, "x")
})

test_that("invalid models and controls are rejected", {
  args <- make_path_case(2)
  args$nboot <- 1L
  bad <- args
  bad$mreg[[1]] <- glm(M_1 ~ E + C, data = args$data)
  expect_error(do.call(cmest_pathcomprisk, bad), "stats::lm")
  bad <- args
  bad$dreg[[1]] <- lm(D1 ~ E + C, data = args$data)
  expect_error(do.call(cmest_pathcomprisk, bad), "binomial glm")
  bad <- args
  bad$yreg <- args$mreg[[1]]
  expect_error(do.call(cmest_pathcomprisk, bad), "timereg::aalen")
  bad <- args
  bad$conf_level <- 1
  expect_error(do.call(cmest_pathcomprisk, bad), "conf_level")
  bad <- args
  bad$data[[bad$time_points[2]]] <- 0
  expect_error(do.call(cmest_pathcomprisk, bad), "strictly increasing")
})
