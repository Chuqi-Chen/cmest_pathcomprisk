#' Causal mediation with time-varying mediators and competing risks
#'
#' @description
#' Implements the mediational g-formula for a time-varying mediator, a
#' time-to-event outcome, and a competing event. Effects can be estimated as
#' hazard differences or differences in survival probability. All mediator,
#' competing-event, and outcome models are refitted in every nonparametric
#' bootstrap sample.
#'
#' @param dreg List of fitted binomial `glm` competing-event models, one
#'   per interval. Each model must be fitted in the risk set that is alive and
#'   free of the primary event at the start of that interval; its response is
#'   1 for a new competing event during that interval; it may be numeric or
#'   logical 0/1, or a two-level factor whose second level is the competing
#'   event. The link may be user-chosen, provided counterfactual predictions
#'   remain valid probabilities. The risk-set restriction
#'   must be reproducible from `data` during bootstrap refitting, for example
#'   with `subset = risk_k` where `risk_k` is a column in `data`, or with an
#'   `NA` response outside the risk set; do not fit `dreg` from a detached,
#'   prefiltered data frame. Use `NULL` when competing events are not modeled.
#' @param mreg List of fitted mediator models, one per time point. Each model
#'   must currently be an ordinary univariate `lm`; binary mediator models
#'   are not yet supported. For later
#'   mediators, any risk-set restriction must be reproducible from `data`
#'   during bootstrap refitting, for example by setting the response to `NA`
#'   outside the risk set or by using a risk-set indicator stored in `data`;
#'   do not fit a model from a detached, prefiltered data frame. Following the
#'   authors' reference functions, hazard differences use one coefficient draw
#'   from each bootstrap-refitted model, whereas survival-probability differences
#'   use the fitted point coefficients.
#' @param mvar Character vector giving the wide mediator variable names in
#'   time order.
#' @param yreg Fitted `timereg::aalen` outcome model. In accordance with the
#'   authors' implementation, the cumulative part may contain only the
#'   baseline intercept and all outcome predictors must have constant
#'   coefficients (use `const()` after loading `timereg`; `timereg::const()`
#'   is not recognized as a constant-coefficient term inside the formula).
#'   Data-dependent outcome bases such
#'   as `poly()`, `splines::ns()`, `splines::bs()`, and `scale()` are not
#'   supported; create any such transformed columns before fitting `yreg`.
#'   The outcome formula must list its variables explicitly rather than use
#'   `.`.
#' @param avar Character string naming the exposure variable.
#' @param a Reference exposure value.
#' @param astar Comparison exposure value.
#' @param data Data frame containing the original wide data.
#' @param time_points Character vector giving the wide measurement-time
#'   variable names in time order. Times must be finite, nonnegative, strictly
#'   increasing within subject, and start at 0.
#' @param yreg_time Character string naming the time-to-event variable in
#'   `data`. Follow-up must end at the primary event, competing event, or
#'   censoring time, whichever occurs first.
#' @param yreg_event Character string naming the event indicator in `data`.
#' @param yreg_mvar Character string naming the long mediator in `yreg`.
#' @param estimate Effect scale. The default is `"hazard"`; use `"sp"` for
#'   survival-probability differences, which are reported in percentage points.
#' @param t Positive follow-up horizon used only when `estimate = "sp"`.
#'   Survival probability is evaluated at `t`; this argument is ignored on the
#'   hazard scale.
#' @param peryr Scaling constant used only for hazard differences. The default
#'   reports cases per 100,000 person-years.
#' @param nboot Number of nonparametric bootstrap samples. Defaults to 1000.
#'   A value of 1 is useful only for smoke testing and gives no inferential
#'   confidence interval.
#' @param ncores Number of bootstrap workers when no external parallel backend
#'   is registered. Defaults to 1; use 1 for sequential execution even if an
#'   external backend is registered. An already registered parallel backend
#'   controls its own worker count when `ncores > 1`.
#' @param conf_level Confidence level for the percentile bootstrap interval.
#'
#' @return A list containing formatted direct (`DE`), mediator indirect
#'   (`IEM`), competing-event indirect (`IED`), and total (`TE`) effects.
#' @export
#' @importFrom foreach foreach %dopar%
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats coef formula predict quantile reshape rbinom terms update
#'   vcov
#' @importFrom survival Surv tmerge
#' @importFrom timereg aalen.des2 predict.aalen
cmest_pathcomprisk <- function(dreg = NULL,
                               mreg,
                               mvar,
                               yreg,
                               avar = "E",
                               a = 0,
                               astar = 1,
                               data,
                               time_points,
                               yreg_time,
                               yreg_event,
                               yreg_mvar,
                               estimate = c("hazard", "sp"),
                               t = 20,
                               peryr = 100000,
                               nboot = 1000,
                               conf_level = 0.95,
                               ncores = 1L) {
  estimate <- match.arg(estimate)

  validated <- .mps_validate_pathcomprisk(
    dreg = dreg, mreg = mreg, mvar = mvar, yreg = yreg,
    avar = avar, a = a, astar = astar, data = data,
    time_points = time_points, yreg_time = yreg_time,
    yreg_event = yreg_event, yreg_mvar = yreg_mvar,
    estimate = estimate, t = t, peryr = peryr,
    nboot = nboot, conf_level = conf_level
  )
  ncores <- .mps_validate_ncores(ncores)
  N <- validated$N
  NL <- validated$NL
  NM <- validated$NM
  mreg_formulas <- validated$mreg_formulas
  dreg_formulas <- validated$dreg_formulas
  yreg_formula <- validated$yreg_formula
  yreg_variables <- validated$yreg_variables
  yreg_response_vars <- validated$yreg_response_vars

  draw_mvn <- function(mean, sigma, label) {
    validated_mvn <- .mps_validate_mvn_parameters(mean, sigma, label)
    draw <- drop(mvtnorm::rmvnorm(
      1L, mean = validated_mvn$mean, sigma = validated_mvn$sigma
    ))
    draw <- .mps_validate_mvn_draw(
      draw, validated_mvn$dimension, label
    )
    names(draw) <- names(validated_mvn$mean)
    draw
  }

  set_exposure <- function(newdata, value) {
    .mps_validate_exposure_assignment(newdata, value, avar)
    if (is.factor(newdata[[avar]])) {
      value <- as.character(value)
      newdata[[avar]] <- factor(
        rep(value, nrow(newdata)), levels = levels(newdata[[avar]]),
        ordered = is.ordered(newdata[[avar]])
      )
    } else {
      newdata[[avar]] <- rep(value, nrow(newdata))
    }
    newdata
  }

  set_variable <- function(newdata, variable, value) {
    .mps_validate_variable_assignment(newdata, variable, value)
    if (is.factor(newdata[[variable]])) {
      value <- as.character(value)
      newdata[[variable]] <- factor(
        value, levels = levels(newdata[[variable]])
      )
    } else {
      newdata[[variable]] <- value
    }
    newdata
  }

  predict_mediator <- function(model, newdata, coefficient_draw) {
    .mps_validate_mediator_coefficient_draw(model, coefficient_draw)
    if (estimate == "sp") {
      mu <- stats::predict(model, newdata = newdata)
    } else {
      design <- stats::model.matrix(
        stats::delete.response(stats::terms(model)), data = newdata,
        contrasts.arg = model$contrasts, xlev = model$xlevels
      )
      coefficient_index <- match(names(stats::coef(model)), colnames(design))
      coefficient_index <- .mps_validate_mediator_design(
        coefficient_index, design, coefficient_draw
      )
      mu <- drop(tcrossprod(
        coefficient_draw, design[, coefficient_index, drop = FALSE]
      ))
    }
    mu <- .mps_validate_mediator_prediction(mu, nrow(newdata))
    as.numeric(mu)
  }

  draw_competing_event <- function(model, newdata) {
    probability <- stats::predict(model, newdata = newdata, type = "response")
    probability <- .mps_validate_competing_event_probability(
      probability, nrow(newdata)
    )
    as.integer(stats::rbinom(length(probability), size = 1L, prob = probability))
  }

  `%boot%` <- if (ncores == 1L) foreach::`%do%` else foreach::`%dopar%`

  if (ncores > 1L &&
      (!foreach::getDoParRegistered() ||
       identical(foreach::getDoParName(), "doSEQ"))) {
    boot_cluster <- parallel::makeCluster(ncores)
    on.exit({
      parallel::stopCluster(boot_cluster)
      foreach::registerDoSEQ()
    }, add = TRUE)
    doParallel::registerDoParallel(boot_cluster)
  }
  parallel_boot <- ncores > 1L
  if (parallel_boot) {
    bootstrap_seeds <- sample.int(.Machine$integer.max, nboot)
  }

  boot_res <- foreach::foreach(
    bootstrap_index = seq_len(nboot),
    .combine = "rbind",
    .packages = c("mvtnorm", "stats", "survival", "timereg"),
    .export = c(
      "getvarnames", ".mps_validate_regression_fit",
      ".mps_validate_aalen_fit", ".mps_validate_reference_aalen",
      ".mps_validate_mvn_parameters", ".mps_validate_mvn_draw",
      ".mps_validate_exposure_assignment",
      ".mps_validate_variable_assignment",
      ".mps_validate_mediator_coefficient_draw",
      ".mps_validate_mediator_design",
      ".mps_validate_mediator_prediction",
      ".mps_validate_competing_event_probability",
      ".mps_validate_outcome_columns",
      ".mps_validate_outcome_history",
      ".mps_validate_counterfactual_aalen_design",
      ".mps_validate_outcome_prediction"
    )
  ) %boot% {
    if (parallel_boot) set.seed(bootstrap_seeds[[bootstrap_index]])
    ind <- sample.int(N, size = N, replace = TRUE)
    data_boot <- data[ind, , drop = FALSE]
    data_boot$id_boot <- seq_len(N)

    # Refit every observed-data model in the current nonparametric bootstrap
    # sample, as required by the mediational g-formula algorithm.
    mreg_boot <- mreg
    for (k in seq_len(NM)) {
      mreg_boot[[k]][["call"]][["formula"]] <- mreg_formulas[[k]]
      mreg_boot[[k]] <- stats::update(mreg_boot[[k]], data = data_boot)
      .mps_validate_regression_fit(
        mreg_boot[[k]],
        paste0("mreg[[", k, "]] in bootstrap ", bootstrap_index)
      )
    }

    dreg_boot <- dreg
    if (NL > 0L) {
      for (k in seq_len(NL)) {
        dreg_boot[[k]][["call"]][["formula"]] <- dreg_formulas[[k]]
        dreg_boot[[k]] <- stats::update(dreg_boot[[k]], data = data_boot)
        .mps_validate_regression_fit(
          dreg_boot[[k]],
          paste0("dreg[[", k, "]] in bootstrap ", bootstrap_index)
        )
      }
    }

    # Reconstruct the counting-process data and refit the Aalen outcome model.
    df_process_boot <- survival::tmerge(
      data_boot,
      data_boot,
      id = id_boot,
      endpt = event(data_boot[[yreg_time]], data_boot[[yreg_event]])
    )

    df_tv_boot <- stats::reshape(
      data_boot,
      direction = "long",
      varying = list(mvar, time_points),
      v.names = c(".mps_mediator", ".mps_time"),
      timevar = ".mps_visit",
      times = seq_len(NM),
      idvar = "id_boot"
    )
    df_tv_boot <- df_tv_boot[
      order(df_tv_boot$id_boot, df_tv_boot$.mps_time),
      ,
      drop = FALSE
    ]

    # The wide data can contain a variable with the same name as the long
    # mediator. The time-dependent version must replace it in the outcome fit.
    if (yreg_mvar %in% names(df_process_boot)) {
      df_process_boot[[yreg_mvar]] <- NULL
    }
    df_process_boot <- survival::tmerge(
      df_process_boot,
      df_tv_boot,
      id = id_boot,
      .mps_yreg_mvar = tdc(.mps_time, .mps_mediator)
    )
    names(df_process_boot)[names(df_process_boot) == ".mps_yreg_mvar"] <- yreg_mvar

    # Preserve the response-variable names used by the supplied counting-
    # process outcome model, even when they are not tstart/tstop/endpt.
    df_process_boot[[yreg_response_vars[1L]]] <- df_process_boot$tstart
    df_process_boot[[yreg_response_vars[2L]]] <- df_process_boot$tstop
    df_process_boot[[yreg_response_vars[3L]]] <- df_process_boot$endpt
    yreg_boot <- yreg
    yreg_boot[["call"]][["formula"]] <- yreg_formula
    yreg_boot <- stats::update(
      yreg_boot, data = df_process_boot,
      id = NULL, n.sim = 1000, robust = 1, resample.iid = 1
    )
    .mps_validate_aalen_fit(
      yreg_boot, paste0("yreg in bootstrap ", bootstrap_index)
    )
    .mps_validate_reference_aalen(
      yreg_boot, paste0("yreg in bootstrap ", bootstrap_index)
    )

    MModel <- if (estimate == "hazard") {
      lapply(seq_len(NM), function(k) {
        draw_mvn(
          stats::coef(mreg_boot[[k]]), stats::vcov(mreg_boot[[k]]),
          paste0("mreg[[", k, "]] in bootstrap ", bootstrap_index)
        )
      })
    } else {
      lapply(mreg_boot, stats::coef)
    }
    # The hazard implementation additionally samples the constant outcome
    # coefficients before generating the competing-event paths.
    YModel <- if (estimate == "hazard") {
      draw_mvn(
        yreg_boot$gamma, yreg_boot$robvar.gamma,
        paste0("yreg in bootstrap ", bootstrap_index)
      )
    } else {
      NULL
    }

    PredictD_a <- PredictD_astar <- PredictD_astar_a <-
      matrix(NA_integer_, nrow = N, ncol = NL)
    PredictM_a <- PredictM_astar <- PredictM_a_astar <-
      matrix(NA_real_, nrow = N, ncol = NM)



    pred.data.a.m1 <- set_exposure(data_boot, a)
    pred.data.astar.m1 <- set_exposure(data_boot, astar)
    PredictM_a[, 1L] <- predict_mediator(
      mreg_boot[[1L]], pred.data.a.m1, MModel[[1L]]
    )
    PredictM_astar[, 1L] <-
      predict_mediator(mreg_boot[[1L]], pred.data.astar.m1, MModel[[1L]])
    PredictM_a_astar[, 1L] <- PredictM_a[, 1L]

    # Match the reference functions: med_sp() skips D1 when NL == 1.
    if (NL > 0L && (estimate == "hazard" || NL > 1L)) {
      # Reuse the exposure assignments already made for M1.
      pred.data.a.d1 <- set_variable(
        pred.data.a.m1, mvar[1L], PredictM_a[, 1L]
      )
      pred.data.astar.d1 <- set_variable(
        pred.data.astar.m1, mvar[1L], PredictM_astar[, 1L]
      )
      pred.data.astar.a.d1 <- set_variable(
        pred.data.astar.m1, mvar[1L], PredictM_a_astar[, 1L]
      )
      # Preserve the independent Bernoulli draws and their order from the
      # authors' reference code: mixed path, a path, then astar path.
      PredictD_astar_a[, 1L] <- draw_competing_event(
        dreg_boot[[1L]], pred.data.astar.a.d1
      )
      PredictD_a[, 1L] <- draw_competing_event(
        dreg_boot[[1L]], pred.data.a.d1
      )
      PredictD_astar[, 1L] <- draw_competing_event(
        dreg_boot[[1L]], pred.data.astar.d1
      )
    }

    # Predict later Mi and Di values.
    if (NM > 1L) {
      for (i in 2:NM) {
        if (NL > 0L) {
          alive_a <- which(PredictD_a[, i - 1L] == 0)
          alive_astar <- which(PredictD_astar[, i - 1L] == 0)
          alive_astar_a <- which(PredictD_astar_a[, i - 1L] == 0)
        } else {
          alive_a <- alive_astar <- alive_astar_a <- seq_len(N)
        }

        if (length(alive_a) > 0L) {
          pred.data.a.m <- pred.data.a.m1[alive_a, , drop = FALSE]
          for (j in seq_len(i - 1L)) {
            pred.data.a.m <- set_variable(
              pred.data.a.m, mvar[j], PredictM_a[alive_a, j]
            )
          }
          PredictM_a[alive_a, i] <-
            predict_mediator(mreg_boot[[i]], pred.data.a.m, MModel[[i]])
        }
        if (length(alive_astar) > 0L) {
          pred.data.astar.m <- pred.data.astar.m1[alive_astar, , drop = FALSE]
          for (j in seq_len(i - 1L)) {
            pred.data.astar.m <- set_variable(
              pred.data.astar.m, mvar[j], PredictM_astar[alive_astar, j]
            )
          }
          PredictM_astar[alive_astar, i] <-
            predict_mediator(
              mreg_boot[[i]], pred.data.astar.m, MModel[[i]]
            )
        }
        if (length(alive_astar_a) > 0L) {
          pred.data.a.astar.m <- pred.data.a.m1[alive_astar_a, , drop = FALSE]
          for (j in seq_len(i - 1L)) {
            pred.data.a.astar.m <- set_variable(
              pred.data.a.astar.m, mvar[j],
              PredictM_a_astar[alive_astar_a, j]
            )
          }
          PredictM_a_astar[alive_astar_a, i] <-
            predict_mediator(
              mreg_boot[[i]], pred.data.a.astar.m, MModel[[i]]
            )
        }

        if (NL > 0L && i <= NL) {
          # Subjects who died previously retain NA at later visits. Later mediator and outcome rows are thereby
          # excluded from their corresponding risk sets.
          if (length(alive_a) > 0L) {
            pred.data.a.d <- set_variable(
              pred.data.a.m, mvar[i], PredictM_a[alive_a, i]
            )
            PredictD_a[alive_a, i] <- draw_competing_event(
              dreg_boot[[i]], pred.data.a.d
            )
          }
          if (length(alive_astar) > 0L) {
            pred.data.astar.d <- set_variable(
              pred.data.astar.m, mvar[i], PredictM_astar[alive_astar, i]
            )
            PredictD_astar[alive_astar, i] <- draw_competing_event(
              dreg_boot[[i]], pred.data.astar.d
            )
          }
          if (length(alive_astar_a) > 0L) {
            # The mixed mediator path used a; its death model uses astar.
            pred.data.astar.a.d <- set_variable(
              set_exposure(pred.data.a.astar.m, astar),
              mvar[i], PredictM_a_astar[alive_astar_a, i]
            )
            PredictD_astar_a[alive_astar_a, i] <- draw_competing_event(
              dreg_boot[[i]], pred.data.astar.a.d
            )
          }
        }
      }
    }

    # Construct the four counterfactual wide data sets.
    yvars <- yreg_variables$xvar
    baseline_yvars <- yvars[yvars != yreg_mvar]
    outcome_columns <- unique(c(
      "id_boot", baseline_yvars, mvar, time_points
    ))
    .mps_validate_outcome_columns(outcome_columns, data_boot)

    pred.data.a.y <-
      pred.data.astar.y <-
      pred.data.astar.a.y <-
      pred.data.astar.a.astar.a.y <-
      data_boot[, outcome_columns, drop = FALSE]

    pred.data.a.y <- set_exposure(pred.data.a.y, a)
    pred.data.astar.y <- set_exposure(pred.data.astar.y, astar)
    pred.data.astar.a.y <- set_exposure(pred.data.astar.a.y, astar)
    pred.data.astar.a.astar.a.y <-
      set_exposure(pred.data.astar.a.astar.a.y, astar)
    for (j in seq_len(NM)) {
      pred.data.a.y <- set_variable(
        pred.data.a.y, mvar[j], PredictM_a[, j]
      )
      pred.data.astar.a.y <- set_variable(
        pred.data.astar.a.y, mvar[j], PredictM_a[, j]
      )
      pred.data.astar.y <- set_variable(
        pred.data.astar.y, mvar[j], PredictM_astar[, j]
      )
      pred.data.astar.a.astar.a.y <- set_variable(
        pred.data.astar.a.astar.a.y, mvar[j], PredictM_a_astar[, j]
      )
    }


    # Mediators after a competing event remain NA; aalen.des2() omits those
    # post-event rows in the same way as the reference model.matrix() calls.
    evaluate_outcome <- function(counterfactual_data) {
      long_data <- stats::reshape(
        counterfactual_data,
        direction = "long",
        varying = list(mvar, time_points),
        v.names = c(yreg_mvar, ".mps_time"),
        timevar = ".mps_visit",
        times = seq_len(NM),
        idvar = "id_boot"
      )
      long_data <- long_data[
        order(long_data$id_boot, long_data$.mps_visit),
        ,
        drop = FALSE
      ]

      .mps_validate_outcome_history(long_data, N, NM)

      design_formula <- stats::formula(stats::delete.response(
        stats::terms(yreg_formula)
      ))
      design <- timereg::aalen.des2(
        design_formula, data = long_data, model = "aalen"
      )
      gamma_names <- rownames(yreg_boot$gamma)
      z_design <- .mps_validate_counterfactual_aalen_design(
        design, gamma_names, nrow(yreg_boot$gamma)
      )

      if (estimate == "hazard") {
        prediction <- drop(tcrossprod(YModel, z_design))
      } else {
        prediction <- drop(timereg::predict.aalen(
          object = yreg_boot,
          Z = z_design,
          times = t,
          n.sim = 0,
          se = FALSE,
          uniform = FALSE,
          resample.iid = 0
        )$S0)
        prediction[prediction > 1] <- 1
      }
      .mps_validate_outcome_prediction(prediction, nrow(z_design))
    }

    PredictY_TEDE_2 <- evaluate_outcome(pred.data.a.y)
    PredictY_DEIEM <- evaluate_outcome(pred.data.astar.a.y)
    PredictY_IEMIED <- evaluate_outcome(
      pred.data.astar.a.astar.a.y
    )
    PredictY_IEDTE_1 <- evaluate_outcome(pred.data.astar.y)
    effect_scale <- if (estimate == "hazard") peryr else 100


    DE <- (mean(PredictY_DEIEM) - mean(PredictY_TEDE_2)) * effect_scale
    IEM <- (mean(PredictY_IEDTE_1) - mean(PredictY_IEMIED)) * effect_scale
    IED <- (mean(PredictY_IEMIED) - mean(PredictY_DEIEM)) * effect_scale
    TE <- (mean(PredictY_IEDTE_1) - mean(PredictY_TEDE_2)) * effect_scale

    c(DE = DE, IEM = IEM, IED = IED, TE = TE)
  }

  if (is.null(dim(boot_res))) {
    boot_names <- names(boot_res)
    boot_res <- matrix(
      boot_res, nrow = 1L,
      dimnames = list(NULL, boot_names)
    )
  }
  boot_res <- .mps_validate_bootstrap_results(boot_res)

  alpha <- (1 - conf_level) / 2
  interval_probs <- c(alpha, 0.5, 1 - alpha)

  summarize_effect <- function(x) {
    qs <- stats::quantile(x, probs = interval_probs, names = FALSE)
    paste0(
      round(qs[2], 3), " (", round(qs[1], 3), ", ",
      round(qs[3], 3), ")"
    )
  }

  DE_result <- summarize_effect(boot_res[, "DE"])
  IEM_result <- summarize_effect(boot_res[, "IEM"])
  IED_result <- summarize_effect(boot_res[, "IED"])
  TE_result <- summarize_effect(boot_res[, "TE"])

  list(
    DE = DE_result,
    IEM = IEM_result,
    IED = IED_result,
    TE = TE_result
  )
}


#' Extract variable names from a regression formula
#'
#' This helper supports calls stored by `timereg::aalen` objects.
#'
#' @param formula A formula, fitted model, model call, character vector, or
#'   `NULL`.
#' @param data Optional data used to expand `.` in a formula.
#'
#' @return A list with all variable names, predictor names, and outcome names.
#' @export
getvarnames <- function(formula, data = NULL) {
  if (is.character(formula)) {
    return(list(varnames = formula, xvar = formula, yvar = NULL))
  }
  if (is.null(formula)) {
    return(list(varnames = NULL, xvar = NULL, yvar = NULL))
  }

  formula <- stats::formula(formula)
  yvars <- NULL
  xvars <- varnames <- all.vars(formula[1:2])

  if (length(formula) == 3L) {
    yvars <- xvars
    xvars <- all.vars(formula[-2])
    if ("." %in% xvars) {
      .mps_validate_formula_expansion_data(data)
      expanded_formula <- stats::formula(stats::terms(formula, data = data))
      xvars <- all.vars(expanded_formula[-2])
    }
    varnames <- c(xvars, varnames)
  }

  list(varnames = varnames, xvar = xvars, yvar = yvars)
}


# Nonstandard-evaluation names used by survival::tmerge().
utils::globalVariables(c(
  ".mps_mediator", ".mps_time", "bootstrap_index", "event", "id_boot", "tdc"
))
