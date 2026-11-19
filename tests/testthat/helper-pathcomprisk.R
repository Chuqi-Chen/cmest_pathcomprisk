# Small synthetic fixtures exercise code behavior, not statistical true effects.
make_path_case <- function(K = 3L, competing = TRUE) {
  set.seed(6200L + K)
  n <- 600L
  d <- data.frame(id = seq_len(n), E = rnorm(n), C = rbinom(n, 1, .4))
  mvar <- paste0("M_", seq_len(K))
  tp <- paste0("time.since.first.exam_", seq_len(K))
  for (k in seq_len(K)) {
    d[[mvar[k]]] <- .4 + .3 * d$E + .2 * d$C +
      (if (k > 1L) .4 * d[[mvar[k - 1L]]] else 0) + rnorm(n)
    d[[tp[k]]] <- k - 1
  }
  ty <- rexp(n, .09 * exp(.2 * d$E + .1 * d$C + .1 * d[[mvar[1L]]]))
  td <- if (competing) rexp(n, .08 * exp(.15 * d$E + .15 * d$C)) else rep(Inf, n)
  d$time_to_event <- pmin(ty, td, 6)
  d$eventHappened <- as.integer(ty <= td & ty <= 6)
  for (k in seq_len(K)) {
    eligible <- d$time_to_event > k - 1L
    # Cumulative D is only read within the at-risk rows during fitting.
    d[[paste0("D", k)]] <- ifelse(eligible,
      as.integer(td <= ty & td <= (if (k == K) 6 else k)), NA_integer_)
    d[[mvar[k]]][!eligible] <- NA_real_
  }
  refit_case(d, K, competing, n_sim = 0L)
}

refit_case <- function(d, K, competing = TRUE, n_sim = 1000L) {
  mvar <- paste0("M_", seq_len(K))
  tp <- paste0("time.since.first.exam_", seq_len(K))
  M <- lapply(seq_len(K), function(k) {
    f <- stats::reformulate(c("E", if (k > 1L) mvar[k - 1L], "C"), mvar[k])
    stats::lm(f, data = d)
  })
  L <- if (competing) lapply(seq_len(K), function(k) {
    stats::glm(stats::reformulate(c("E", mvar[k], "C"), paste0("D", k)),
               data = d, family = stats::binomial())
  }) else NULL
  if (!"id_unique" %in% names(d)) d$id_unique <- seq_len(nrow(d))
  long <- survival::tmerge(d, d, id = id_unique,
                           endpt = event(time_to_event, eventHappened))
  med <- stats::reshape(d, direction = "long", varying = list(mvar, tp),
                        v.names = c("M", "visit_time"), timevar = "visit",
                        times = seq_len(K), idvar = "id_unique")
  med <- med[order(med$id_unique, med$visit_time), ]
  long <- survival::tmerge(long, med, id = id_unique, M = tdc(visit_time, M))
  const <- timereg::const
  Y <- timereg::aalen(survival::Surv(tstart, tstop, endpt) ~
                       const(E) + const(M) + const(C), data = long,
                     id = NULL, n.sim = n_sim, robust = 1, resample.iid = 1)
  list(dreg = L, mreg = M, mvar = mvar, yreg = Y, avar = "E", a = -.3,
       astar = .7, data = d, time_points = tp, yreg_time = "time_to_event",
       yreg_event = "eventHappened", yreg_mvar = "M")
}

# Observe raw draws through the existing validation call; the estimator body is
# unchanged. The cloned function uses the actual installed namespace as parent.
run_with_raw_draws <- function(args, scale, seed, nboot = 2L, ncores = 1L) {
  f <- cmest_pathcomprisk
  env <- new.env(parent = environment(f))
  original_validator <- get(".mps_validate_bootstrap_results", envir = environment(f))
  env$.mps_validate_bootstrap_results <- function(x) {
    env$raw <- x
    original_validator(x)
  }
  environment(f) <- env
  set.seed(seed)
  out <- do.call(f, c(args, list(estimate = scale, t = 5, nboot = nboot,
                                ncores = ncores)))
  list(summary = out, raw = env$raw, rng = .Random.seed)
}

run_reference <- function(args, scale, seed, nboot = 2L) {
  ref <- new.env(parent = asNamespace("timereg"))
  ref$rmvnorm <- mvtnorm::rmvnorm
  sys.source(testthat::test_path("fixtures", "pathcomprisk-professor-functions.R"), envir = ref)
  if (length(args$mvar) == 1L && scale == "hazard") {
    # The reference skips reshape() for NM=1 and expects the outcome mediator
    # to already have its wide name. Adapt only the column-name lookup.
    original_getvarnames <- ref$getvarnames
    ref$getvarnames <- function(formula, data = NULL) {
      out <- original_getvarnames(formula, data)
      out$xvar[out$xvar == "M"] <- args$mvar[1L]
      out
    }
  }
  set.seed(seed)
  raw <- t(vapply(seq_len(nboot), function(i) {
    d <- args$data[sample.int(nrow(args$data), nrow(args$data), replace = TRUE), ]
    d$id_unique <- seq_len(nrow(d))
    fit <- refit_case(d, length(args$mvar), length(args$dreg) > 0L, n_sim = 1000L)
    f <- if (scale == "hazard") ref$med_longitudinal else ref$med_sp
    x <- list(L = fit$dreg, M = fit$mreg, m = fit$mvar, Y = fit$yreg,
              treat = "E", control.value = args$a, treat.value = args$astar,
              data = fit$data, time_points = fit$time_points)
    if (scale == "sp") x$t <- 5
    as.numeric(do.call(f, x))
  }, numeric(4L)))
  colnames(raw) <- c("DE", "IEM", "IED", "TE")
  list(raw = raw, rng = .Random.seed)
}
