# Internal input, fitted-model, and runtime checks for cmest_pathcomprisk().
# Keep these checks separate from the mediation estimand calculations.

.mps_validate_regression_fit <- function(model, label) {
  if (inherits(model, "glm") && !isTRUE(model$converged)) {
    stop(label, " did not converge.")
  }
  coefficients <- stats::coef(model)
  if (length(coefficients) < 1L || any(!is.finite(coefficients)) ||
      (!is.null(model$rank) && model$rank < length(coefficients))) {
    stop(label, " is rank-deficient or has non-finite coefficients.")
  }
  if (inherits(model, "glm") &&
      model$family$family %in% c("binomial", "quasibinomial")) {
    probability <- stats::fitted(model)
    boundary_tolerance <- sqrt(.Machine$double.eps)
    if (length(probability) < 1L || any(!is.finite(probability)) ||
        any(probability <= boundary_tolerance |
            probability >= 1 - boundary_tolerance)) {
      stop(label, " has fitted probabilities at a numerical boundary; ",
           "possible complete or quasi-complete separation.")
    }
    if (!is.null(model$y)) {
      positive_weight <- if (is.null(model$prior.weights)) {
        rep(TRUE, length(model$y))
      } else {
        model$prior.weights > 0
      }
      observed_response <- model$y[positive_weight & !is.na(model$y)]
      if (length(unique(observed_response)) < 2L) {
        stop(label, " must contain both competing-event outcomes.")
      }
    }
  }
  model
}

.mps_validate_aalen_fit <- function(model, label) {
  coefficients <- numeric(0)
  if (is.matrix(model$cum) && ncol(model$cum) > 1L) {
    coefficients <- c(coefficients, model$cum[, -1L, drop = FALSE])
  }
  if (!is.null(model$gamma)) {
    coefficients <- c(coefficients, model$gamma)
  }
  if (length(coefficients) < 1L || any(!is.finite(coefficients))) {
    stop(label, " has non-finite or unavailable Aalen coefficients.")
  }
  model
}

.mps_validate_reference_aalen <- function(model, label) {
  cumulative_names <- if (is.matrix(model$cum)) {
    colnames(model$cum)[-1L]
  } else {
    character(0)
  }
  if (!identical(cumulative_names, "(Intercept)")) {
    stop(label, " must use const() for every outcome predictor; only the ",
         "cumulative baseline intercept may be time-varying.")
  }
  if (!is.matrix(model$gamma) || ncol(model$gamma) != 1L ||
      nrow(model$gamma) < 1L || is.null(rownames(model$gamma)) ||
      anyDuplicated(rownames(model$gamma)) ||
      any(!is.finite(model$gamma))) {
    stop(label, " must contain finite, named constant Aalen coefficients.")
  }
  model
}

.mps_validate_mvn_parameters <- function(mean, sigma, label) {
  mean <- drop(mean)
  sigma <- as.matrix(sigma)
  dimension <- length(mean)
  tolerance <- sqrt(.Machine$double.eps)
  if (dimension < 1L || any(!is.finite(mean)) ||
      !identical(dim(sigma), c(dimension, dimension)) ||
      any(!is.finite(sigma)) ||
      !isTRUE(isSymmetric(
        sigma, tol = tolerance, check.attributes = FALSE
      ))) {
    stop(label, " has an invalid coefficient covariance matrix.")
  }
  eigenvalues <- eigen(sigma, symmetric = TRUE, only.values = TRUE)$values
  eigen_tolerance <- tolerance * max(1, max(abs(eigenvalues)))
  if (min(eigenvalues) < -eigen_tolerance) {
    stop(label, " has a coefficient covariance matrix that is not ",
         "positive semidefinite.")
  }
  list(mean = mean, sigma = sigma, dimension = dimension)
}

.mps_validate_mvn_draw <- function(draw, dimension, label) {
  if (length(draw) != dimension || any(!is.finite(draw))) {
    stop(label, " produced an invalid coefficient draw.")
  }
  draw
}

.mps_validate_exposure_assignment <- function(newdata, value, avar) {
  if (is.factor(newdata[[avar]])) {
    if (length(value) != 1L ||
        !as.character(value) %in% levels(newdata[[avar]])) {
      stop("a and astar must be valid levels of the factor exposure.")
    }
  } else if (length(value) != 1L || is.na(value)) {
    stop("a and astar must each be one nonmissing exposure value.")
  }
  invisible(TRUE)
}

.mps_validate_variable_assignment <- function(newdata, variable, value) {
  if (is.factor(newdata[[variable]]) &&
      any(!as.character(value) %in% levels(newdata[[variable]]))) {
    stop("A simulated mediator value is not a valid factor level.")
  }
  invisible(TRUE)
}

.mps_validate_mediator_coefficient_draw <- function(model, coefficient_draw) {
  if (length(coefficient_draw) != length(stats::coef(model)) ||
      any(!is.finite(coefficient_draw))) {
    stop("Internal mediator coefficient draws are invalid.")
  }
  invisible(TRUE)
}

.mps_validate_mediator_design <- function(
    coefficient_index, design, coefficient_draw) {
  if (anyNA(coefficient_index) || anyDuplicated(coefficient_index) ||
      ncol(design) != length(coefficient_draw)) {
    stop("A mediator design matrix does not match its fitted model.")
  }
  coefficient_index
}

.mps_validate_mediator_prediction <- function(prediction, expected_length) {
  if (length(prediction) != expected_length || any(!is.finite(prediction))) {
    stop("A mediator model produced non-finite or incomplete predictions.")
  }
  prediction
}

.mps_validate_competing_event_probability <- function(
    probability, expected_length) {
  if (length(probability) != expected_length ||
      any(!is.finite(probability)) ||
      any(probability < 0 | probability > 1)) {
    stop("A competing-event model must predict finite probabilities in [0, 1].")
  }
  probability
}

.mps_validate_outcome_columns <- function(outcome_columns, data) {
  if (anyNA(match(outcome_columns, names(data)))) {
    stop("All non-mediator outcome predictors must be available in data.")
  }
  invisible(TRUE)
}

.mps_validate_outcome_history <- function(long_data, N, NM) {
  expected_id <- rep(seq_len(N), each = NM)
  expected_visit <- rep(seq_len(NM), times = N)
  if (!identical(as.integer(long_data$id_boot), expected_id) ||
      !identical(as.integer(long_data$.mps_visit), expected_visit)) {
    stop("Internal outcome histories are not ordered by subject and visit.")
  }
  invisible(TRUE)
}

.mps_validate_counterfactual_aalen_design <- function(
    design, gamma_names, gamma_rows) {
  z_index <- match(gamma_names, design$covnamesZ)
  if (anyNA(z_index) || anyDuplicated(z_index) ||
      length(design$covnamesZ) != length(gamma_names) ||
      !identical(design$covnamesX, "(Intercept)")) {
    stop("The counterfactual Aalen design does not match yreg.")
  }
  z_design <- design$Z[, z_index, drop = FALSE]
  if (nrow(z_design) < 1L || ncol(z_design) != gamma_rows ||
      any(!is.finite(z_design))) {
    stop("The counterfactual Aalen design is empty or non-finite.")
  }
  z_design
}

.mps_validate_outcome_prediction <- function(prediction, expected_length) {
  if (length(prediction) != expected_length || any(!is.finite(prediction))) {
    stop("The outcome model produced invalid counterfactual predictions.")
  }
  prediction
}

.mps_validate_bootstrap_results <- function(boot_res) {
  if (!identical(colnames(boot_res), c("DE", "IEM", "IED", "TE")) ||
      any(!is.finite(boot_res))) {
    stop("One or more bootstrap iterations produced invalid effect estimates.")
  }
  boot_res
}

.mps_validate_formula_expansion_data <- function(data) {
  if (is.null(data) || length(data) == 0L) {
    stop("'.' appears in formula but data was not supplied.")
  }
  invisible(TRUE)
}

.mps_validate_ncores <- function(ncores) {
  if (length(ncores) != 1L || !is.numeric(ncores) ||
      !is.finite(ncores) || ncores < 1L ||
      ncores > .Machine$integer.max || ncores != floor(ncores)) {
    stop("ncores must be a positive integer.")
  }
  as.integer(ncores)
}

.mps_validate_pathcomprisk <- function(
    dreg, mreg, mvar, yreg, avar, a, astar, data, time_points,
    yreg_time, yreg_event, yreg_mvar, estimate, t, peryr,
    nboot, conf_level) {
  if (!is.data.frame(data) || nrow(data) < 1L) {
    stop("data must be a nonempty data frame.")
  }
  if (!is.list(mreg) || (!is.null(dreg) && !is.list(dreg))) {
    stop("mreg must be a list and dreg must be NULL or a list.")
  }
  N <- nrow(data)
  NL <- length(dreg)
  NM <- length(mreg)

  if (NM < 1L) {
    stop("mreg must contain at least one fitted mediator model.")
  }
  if (!is.character(mvar) || !is.character(time_points) ||
      length(mvar) != NM || length(time_points) != NM ||
      anyDuplicated(mvar) || anyDuplicated(time_points) ||
      length(intersect(mvar, time_points)) > 0L) {
    stop("mreg, mvar, and time_points must have the same length.")
  }
  if (!NL %in% c(0L, NM)) {
    stop("dreg must be NULL or contain one model for each mediator time point.")
  }
  if (!is.character(avar) || length(avar) != 1L || !avar %in% names(data)) {
    stop("avar must name one exposure variable in data.")
  }
  if (!is.numeric(data[[avar]]) && !is.factor(data[[avar]])) {
    stop("avar must name a numeric or factor exposure variable.")
  }
  internal_names <- c(
    "id_boot", ".mps_mediator", ".mps_time", ".mps_visit",
    ".mps_yreg_mvar"
  )
  if (any(internal_names %in% names(data))) {
    stop("data contains a variable name reserved for internal processing.")
  }
  if (!is.character(yreg_time) || length(yreg_time) != 1L ||
      !is.character(yreg_event) || length(yreg_event) != 1L) {
    stop("yreg_time and yreg_event must each name one variable in data.")
  }
  required_data_vars <- c(mvar, time_points, yreg_time, yreg_event)
  if (anyNA(match(required_data_vars, names(data)))) {
    stop("mvar, time_points, yreg_time, and yreg_event must name variables in data.")
  }
  for (k in seq_len(NM)) {
    if (!identical(class(mreg[[k]]), "lm")) {
      stop("mreg[[", k, "]] must currently be a fitted stats::lm model; ",
           "binary mediator models are not yet supported.")
    }
  }
  for (k in seq_len(NL)) {
    if (!inherits(dreg[[k]], "glm") ||
        !identical(dreg[[k]]$family$family, "binomial")) {
      stop("dreg[[", k, "]] must be a fitted binomial glm model for a ",
           "binary competing event.")
    }
  }
  if (!inherits(yreg, "aalen")) {
    stop("yreg must be a fitted timereg::aalen model.")
  }
  recover_model_formula <- function(model, label) {
    model_formula <- tryCatch(
      stats::formula(model),
      error = function(e) NULL
    )
    if (!inherits(model_formula, "formula")) {
      stop(label, " does not contain a recoverable model formula.")
    }
    model_formula
  }
  mreg_formulas <- lapply(seq_len(NM), function(k) {
    recover_model_formula(mreg[[k]], paste0("mreg[[", k, "]]"))
  })
  if (any(vapply(mreg, function(model) {
    length(attr(stats::terms(model), "offset")) > 0L || !is.null(model$offset)
  }, logical(1)))) {
    stop("Mediator models with offsets are not supported by the authors' method.")
  }
  mreg_responses <- vapply(seq_len(NM), function(k) {
    response <- mreg_formulas[[k]][[2L]]
    if (is.name(response)) as.character(response) else NA_character_
  }, character(1))
  if (anyNA(mreg_responses) || !identical(unname(mreg_responses), mvar)) {
    stop("The response of mreg[[k]] must be the untransformed variable mvar[k].")
  }
  dreg_formulas <- lapply(seq_len(NL), function(k) {
    recover_model_formula(dreg[[k]], paste0("dreg[[", k, "]]"))
  })
  validate_model_data_scope <- function(model, model_formula, label) {
    predictor_names <- all.vars(model_formula[[3L]])
    if ("." %in% predictor_names ||
        anyNA(match(predictor_names, names(data)))) {
      stop(label, " must list predictors explicitly and use only variables in data.")
    }
    subset_expression <- model$call$subset
    if (!is.null(subset_expression) &&
        !identical(subset_expression, quote(NULL)) &&
        !identical(subset_expression, TRUE)) {
      subset_names <- all.vars(subset_expression)
      if (length(subset_names) < 1L ||
          anyNA(match(subset_names, names(data)))) {
        stop(label, " subset must be defined only by variables stored in data.")
      }
    }
  }
  for (k in seq_len(NM)) {
    validate_model_data_scope(
      mreg[[k]], mreg_formulas[[k]], paste0("mreg[[", k, "]] ")
    )
  }
  for (k in seq_len(NL)) {
    validate_model_data_scope(
      dreg[[k]], dreg_formulas[[k]], paste0("dreg[[", k, "]] ")
    )
  }
  dreg_responses <- vapply(seq_len(NL), function(k) {
    response <- dreg_formulas[[k]][[2L]]
    if (is.name(response)) as.character(response) else NA_character_
  }, character(1))
  for (k in seq_len(NL)) {
    response <- dreg_responses[k]
    if (is.na(response) || !response %in% names(data)) {
      stop("dreg[[", k, "]] response must be an untransformed binary ",
           "variable in data.")
    }
    value <- data[[response]]
    numeric_binary <- (is.numeric(value) || is.logical(value)) &&
      all(is.na(value) | value %in% c(0, 1))
    factor_binary <- is.factor(value) && nlevels(value) == 2L
    if (!numeric_binary && !factor_binary) {
      stop("dreg[[", k, "]] response must be numeric/logical 0/1 or a ",
           "two-level factor (the second level is the competing event); ",
           "NA is allowed outside its risk set.")
    }
  }
  for (k in seq_len(NM)) {
    .mps_validate_regression_fit(mreg[[k]], paste0("mreg[[", k, "]]"))
  }
  for (k in seq_len(NL)) {
    .mps_validate_regression_fit(dreg[[k]], paste0("dreg[[", k, "]]"))
  }
  if (!all(vapply(data[, mvar, drop = FALSE], is.numeric, logical(1)))) {
    stop("The mediator variables named by mvar must be numeric.")
  }
  yreg_formula <- attr(yreg, "Formula", exact = TRUE)
  if (!inherits(yreg_formula, "formula")) {
    yreg_formula <- tryCatch(
      stats::formula(yreg),
      error = function(e) NULL
    )
  }
  if (!inherits(yreg_formula, "formula")) {
    stop("The formula stored in yreg could not be recovered.")
  }
  yreg_rhs <- yreg_formula[[3L]]
  if ("." %in% all.vars(yreg_rhs)) {
    stop("yreg must list every outcome predictor explicitly; '~ .' is not supported.")
  }
  call_name <- function(call) {
    if (is.name(call)) {
      return(as.character(call))
    }
    if (is.call(call) && length(call) >= 3L &&
        as.character(call[[1L]]) %in% c("::", ":::")) {
      return(as.character(call[[3L]]))
    }
    ""
  }
  contains_data_dependent_basis <- function(expression) {
    if (!is.call(expression)) {
      return(FALSE)
    }
    if (call_name(expression[[1L]]) %in% c("poly", "ns", "bs", "scale")) {
      return(TRUE)
    }
    arguments <- as.list(expression)[-1L]
    length(arguments) > 0L && any(vapply(
      arguments, contains_data_dependent_basis, logical(1)
    ))
  }
  if (contains_data_dependent_basis(yreg_rhs)) {
    stop("yreg cannot contain poly(), ns(), bs(), or scale(); create the ",
         "transformed outcome predictors in data before fitting yreg.")
  }
  .mps_validate_aalen_fit(yreg, "yreg")
  yreg_variables <- getvarnames(yreg_formula)
  if (!is.character(yreg_mvar) || length(yreg_mvar) != 1L ||
      yreg_mvar %in% internal_names ||
      !yreg_mvar %in% yreg_variables$xvar) {
    stop("yreg_mvar must name a non-reserved long mediator in the outcome model.")
  }
  yreg_response_vars <- yreg_variables$yvar
  yreg_response <- yreg_formula[[2L]]
  if (!is.call(yreg_response) ||
      call_name(yreg_response[[1L]]) != "Surv" ||
      length(yreg_response) != 4L ||
      !all(vapply(as.list(yreg_response)[-1L], is.name, logical(1))) ||
      length(yreg_response_vars) != 3L ||
      anyDuplicated(yreg_response_vars)) {
    stop("yreg must use a counting-process Surv(start, stop, event) response.")
  }
  .mps_validate_reference_aalen(yreg, "yreg")
  if (length(nboot) != 1L || !is.numeric(nboot) || !is.finite(nboot) ||
      nboot < 1 || nboot > .Machine$integer.max || nboot != floor(nboot)) {
    stop("nboot must be a positive integer.")
  }
  if (length(conf_level) != 1L || !is.numeric(conf_level) ||
      !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be strictly between 0 and 1.")
  }
  if (estimate == "hazard" &&
      (length(peryr) != 1L || !is.numeric(peryr) || !is.finite(peryr) ||
       peryr <= 0)) {
    stop("peryr must be a positive number for hazard differences.")
  }
  if (estimate == "sp" &&
      (length(t) != 1L || !is.numeric(t) || !is.finite(t) || t <= 0)) {
    stop("t must be one positive numeric follow-up time for survival probability.")
  }
  validate_exposure_value <- function(value, label) {
    if (length(value) != 1L || is.na(value)) {
      stop(label, " must be one nonmissing exposure value.")
    }
    if (is.factor(data[[avar]])) {
      if (!as.character(value) %in% levels(data[[avar]])) {
        stop(label, " must be a valid level of the factor exposure.")
      }
    } else if (is.numeric(data[[avar]]) &&
               (!is.numeric(value) || !is.finite(value))) {
      stop(label, " must be one finite numeric exposure value.")
    }
  }
  validate_exposure_value(a, "a")
  validate_exposure_value(astar, "astar")
  time_matrix <- as.matrix(data[, time_points, drop = FALSE])
  if (!is.numeric(time_matrix) || any(!is.finite(time_matrix)) ||
      any(time_matrix < 0)) {
    stop("time_points must contain finite nonnegative numeric values.")
  }
  if (NM > 1L && any(apply(time_matrix, 1L, function(x) any(diff(x) <= 0)))) {
    stop("time_points must be strictly increasing within each subject.")
  }
  if (any(time_matrix[, 1L] != 0)) {
    stop("The first mediator measurement time must be 0 for every subject.")
  }
  if (!is.numeric(data[[yreg_time]]) || any(!is.finite(data[[yreg_time]])) ||
      any(data[[yreg_time]] <= 0)) {
    stop("yreg_time must contain finite strictly positive numeric values.")
  }
  if ((!is.numeric(data[[yreg_event]]) &&
       !is.logical(data[[yreg_event]])) ||
      anyNA(data[[yreg_event]]) ||
      any(!data[[yreg_event]] %in% c(0, 1))) {
    stop("yreg_event must be a complete binary event indicator.")
  }
  if (estimate == "sp" && t > max(data[[yreg_time]])) {
    stop("t cannot exceed the largest observed follow-up time.")
  }

  list(
    N = N, NL = NL, NM = NM,
    mreg_formulas = mreg_formulas, dreg_formulas = dreg_formulas,
    yreg_formula = yreg_formula, yreg_variables = yreg_variables,
    yreg_response_vars = yreg_response_vars
  )
}
