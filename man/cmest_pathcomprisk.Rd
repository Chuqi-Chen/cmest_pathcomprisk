\name{cmest_pathcomprisk}
\alias{cmest_pathcomprisk}
\title{Estimate Path-Specific Effects with Time-Varying Mediators and Competing Risks}
\description{
Implements the authors' mediational g-formula for a longitudinal continuous
mediator, a time-to-event primary outcome, and an optional time-varying
competing-event process. The function returns path-specific direct, mediator,
and competing-event effects, together with their total effect, on either the
additive-hazard or survival-probability scale. Uncertainty is evaluated by a
subject-level nonparametric bootstrap in which every supplied regression
model is refitted.
}
\usage{
cmest_pathcomprisk(
  dreg = NULL,
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
  ncores = 1L
)
}
\arguments{
  \item{dreg}{
    Either \code{NULL}, when no competing event is modeled, or a list of
    fitted \code{stats::glm} objects in measurement-time order. The list must
    have exactly one element for every element of \code{mreg} and
    \code{mvar}.

    Every element must use \code{family = binomial(...)}. Any binomial link
    supported by \code{predict.glm(type = "response")} may be used. Model
    \code{k} must describe the occurrence of the new competing event in the
    interval associated with mediator visit \code{k}, among individuals who
    are still eligible for that interval. Its response must be an
    untransformed variable in \code{data}, coded as numeric or logical 0/1,
    or as a two-level factor whose second level represents the competing
    event. \code{NA} responses may be used outside the applicable risk set.

    The formula must list predictors explicitly; \code{.} is not supported.
    Every predictor and every variable used in \code{subset} must be stored in
    \code{data}, because the model is refitted after every bootstrap resample.
    A risk-set restriction should therefore be represented by a variable in
    \code{data}, or by an \code{NA} response outside the risk set. A model
    fitted from a detached, prefiltered data frame cannot in general be
    reconstructed correctly and should not be supplied.
  }
  \item{mreg}{
    Nonempty list of fitted mediator models in measurement-time order. The
    current implementation accepts only objects whose class is exactly
    \code{"lm"}; generalized, binary, ordinal, multivariate, mixed, and other
    mediator models are not currently supported. Model \code{k} must have
    the untransformed response \code{mvar[k]}. The response variables named
    by \code{mvar} must be numeric. Offsets are not supported.

    Formulas must list all predictors explicitly and all predictors must be
    columns of \code{data}. Any \code{subset} condition must depend only on
    columns in \code{data}. As for \code{dreg}, restrictions for later visits
    must be reproducible during bootstrap refitting and must not rely on a
    detached data frame. The initial fits and all bootstrap fits must be full
    rank and have finite coefficients.
  }
  \item{mvar}{
    Character vector of unique wide-format mediator column names, ordered
    from the first to the last measurement. It must have length
    \code{length(mreg)}; \code{mvar[k]} must be the response of
    \code{mreg[[k]]}. All named columns must exist in \code{data} and be
    numeric.
  }
  \item{yreg}{
    A fitted \code{timereg::aalen} additive-hazards model. Its response must
    be a three-variable counting-process response of the form
    \code{survival::Surv(start, stop, event)}. The cumulative part of the
    Aalen model may contain only the baseline intercept. Every outcome
    predictor must have a constant coefficient and must therefore be wrapped
    in \code{const()} in the supplied model formula, after loading
    \pkg{timereg}. Using \code{timereg::const()} inside the formula is not
    recognized as a constant-coefficient term by \code{timereg::aalen()}.

    The formula must list predictors explicitly and must include the long
    mediator named by \code{yreg_mvar}. All other outcome predictors must be
    available as columns of \code{data}. \code{poly()}, \code{scale()},
    \code{splines::ns()}, and \code{splines::bs()} are not supported inside
    the outcome formula. If such a basis is required, calculate its columns
    in \code{data} before fitting the model and refer to those columns
    explicitly. Constant coefficients must be finite, and their names must be
    present, unique, and reproducible in the counterfactual Aalen design
    matrix.

    The supplied object defines and validates the outcome-model
    specification. Within each bootstrap sample the function rebuilds the
    counting-process data from \code{data}, \code{mvar}, \code{time_points},
    \code{yreg_time}, and \code{yreg_event}, then refits this formula with
    \code{id = NULL}, \code{n.sim = 1000}, \code{robust = 1}, and
    \code{resample.iid = 1}, matching the authors' implementation.
  }
  \item{avar}{
    Length-one character string naming the baseline exposure column in
    \code{data}. Numeric and factor exposures are supported. The exposure
    should be measured before the mediator and outcome processes represented
    by the fitted models.
  }
  \item{a}{
    Length-one, nonmissing reference exposure value. If \code{avar} is
    numeric, \code{a} must be finite and numeric. If it is a factor,
    \code{a} must identify an existing factor level.
  }
  \item{astar}{
    Length-one, nonmissing comparison exposure value. The same type and
    factor-level requirements as for \code{a} apply. Effect signs represent
    movement from the reference value \code{a} to the comparison value
    \code{astar} along the path specified by each estimand.
  }
  \item{data}{
    A nonempty wide-format \code{data.frame} with one row per individual.
    Rows, rather than preconstructed long-format records, are sampled with
    replacement in the nonparametric bootstrap. The data frame must contain
    the exposure, all wide mediator and measurement-time variables, primary
    outcome follow-up and event variables, the responses and predictors of
    the mediator and competing-event models, all non-mediator outcome
    predictors, and any variables used in model \code{subset} expressions.
    The counting-process response variables in the supplied \code{yreg}
    object do not need to be columns of this wide data frame because they are
    reconstructed internally.

    The following names are reserved for internal processing and must not
    occur in \code{data}: \code{id_boot}, \code{.mps_mediator},
    \code{.mps_time}, \code{.mps_visit}, and \code{.mps_yreg_mvar}.
  }
  \item{time_points}{
    Character vector naming the wide-format mediator measurement-time
    columns, in the same order and with the same length as \code{mvar} and
    \code{mreg}. Names must be unique and must not overlap \code{mvar}.
    Values must be numeric, finite, nonnegative, complete, and strictly
    increasing within every individual. The first measurement time must be
    exactly 0. Measurement times and \code{yreg_time} must use a common time
    unit.
  }
  \item{yreg_time}{
    Length-one character string naming the primary-outcome follow-up time in
    \code{data}. Values must be numeric, finite, complete, and strictly
    positive. For the intended competing-risk interpretation, follow-up must
    end at the earliest of the primary event, competing event, or censoring.
  }
  \item{yreg_event}{
    Length-one character string naming the primary-event indicator in
    \code{data}. It must be a complete numeric or logical 0/1 variable, where
    1 denotes occurrence of the primary event at \code{yreg_time}.
  }
  \item{yreg_mvar}{
    Length-one character string naming the long-format mediator term in the
    \code{yreg} formula. It must occur among the outcome predictors and must
    not be one of the internally reserved names listed under \code{data}.
    The function creates this long variable from \code{mvar} and
    \code{time_points} separately within every bootstrap sample and
    counterfactual path.
  }
  \item{estimate}{
    Character string selecting the effect scale. \code{"hazard"}, the
    default, returns differences in the constant-coefficient linear predictor
    of the additive-hazards model. \code{"sp"} returns differences in
    predicted survival probability at time \code{t}. Standard
    \code{match.arg()} matching is used.
  }
  \item{t}{
    Length-one finite positive follow-up horizon used only when
    \code{estimate = "sp"}. It must use the same time unit as
    \code{yreg_time} and cannot exceed the largest observed value of
    \code{data[[yreg_time]]}. It is ignored for hazard differences.
  }
  \item{peryr}{
    Length-one finite positive multiplier used only when
    \code{estimate = "hazard"}. Every hazard difference is multiplied by
    this value. The default, 100000, reports a difference per 100,000 units
    of follow-up time; it is interpretable as per 100,000 person-years only
    when the input time variables are measured in years. It is ignored for
    survival-probability differences.
  }
  \item{nboot}{
    Positive integer giving the number of subject-level nonparametric
    bootstrap replicates. The default is 1000. Larger values reduce Monte
    Carlo error but increase run time. \code{nboot = 1} is suitable only for
    a smoke test: its reported lower limit, median, and upper limit are all
    the single bootstrap value and are not an inferential confidence
    interval.
  }
  \item{conf_level}{
    Length-one finite number strictly between 0 and 1 giving the level of the
    equal-tailed percentile bootstrap interval. The default is 0.95, which
    uses the 0.025 and 0.975 empirical quantiles.
  }
  \item{ncores}{
    Positive integer giving the number of bootstrap workers when no external
    parallel backend is registered. The default is 1. Use \code{ncores = 1}
    for sequential execution even when an external backend is registered.
    If a parallel \pkg{foreach} backend is already registered and
    \code{ncores > 1}, that backend determines the actual worker count.
    \pkg{doParallel} is needed when the function creates its own workers.
  }
}
\details{
The function generalizes the authors' concrete analysis code while retaining
its model-refitting, counterfactual-history generation, coefficient-sampling,
and effect-calculation steps. Supplying fitted models gives users flexibility
over covariates and supported formula terms; it does not cause the original
fitted coefficients to be reused without bootstrap refitting.

The mediator is represented in wide form in \code{data}. During outcome-model
refitting and counterfactual evaluation, the function reshapes the mediator
history to counting-process long form. A mediator value is treated as current
from its recorded measurement time until the next measurement. When a
competing event is simulated, subsequent mediators on that path remain
missing, and their post-event long-format records are omitted from the Aalen
design in the same manner as in the authors' reference code.
}
\section{Counterfactual paths and effect definitions}{
Within each bootstrap replicate, the function evaluates four mean
counterfactual outcome predictions. For concise notation, let:
\describe{
  \item{\eqn{\mu_0}}{
    the mean prediction under exposure \code{a}, with the mediator and, when
    present, competing-event histories generated under the reference path;
  }
  \item{\eqn{\mu_1}}{
    the mean prediction after changing the exposure in the outcome model to
    \code{astar}, while retaining the reference mediator and competing-event
    histories;
  }
  \item{\eqn{\mu_2}}{
    the mean prediction after additionally changing the exposure used for the
    competing-event process to \code{astar}, while the mediator process is
    still generated under \code{a}; and
  }
  \item{\eqn{\mu_3}}{
    the mean prediction under the full comparison path, with exposure,
    mediator, and competing-event processes generated under \code{astar}.
  }
}

The unscaled effects in every bootstrap replicate are
\deqn{DE = \mu_1 - \mu_0,}
\deqn{IED = \mu_2 - \mu_1,}
\deqn{IEM = \mu_3 - \mu_2,}
\deqn{TE = \mu_3 - \mu_0.}
Consequently, before display rounding, \eqn{TE = DE + IED + IEM} within each
replicate. \code{DE} is the path-specific direct effect, \code{IED} is the
indirect effect through the competing-event process, and \code{IEM} is the
indirect effect through the longitudinal mediator process. When
\code{dreg = NULL}, no competing-event path is simulated and \code{IED} is
zero by construction.

For \code{estimate = "hazard"}, a positive value denotes a higher additive
hazard linear predictor along the changed path. Effects are multiplied by
\code{peryr}. The cumulative baseline intercept is omitted from the authors' hazard
calculation; this function reports the constant-coefficient linear-predictor
contrast.

For \code{estimate = "sp"}, a positive value denotes greater predicted
survival probability at \code{t} along the changed path. Differences are
multiplied by 100 and are therefore reported in percentage points, not as
proportions. As in the implemented reference method, predicted survival
values above 1 are truncated to 1 before contrasts are formed.
}
\section{Bootstrap and model simulation}{
For each of the \code{nboot} replicates, the function performs the following
operations:
\enumerate{
  \item Samples \code{nrow(data)} individuals with replacement and assigns a
  new internal subject identifier.
  \item Refits every mediator model and, when supplied, every competing-event
  model in the resampled wide data.
  \item Reconstructs counting-process outcome data and refits the Aalen model.
  \item Generates the mediator histories from conditional mean predictions.
  No mediator residual is drawn. For hazard effects, one multivariate-normal
  coefficient vector is drawn from each refitted mediator model using its
  fitted coefficient vector and covariance matrix. For survival-probability
  effects, the refitted mediator point coefficients are used without this
  additional draw.
  \item When \code{dreg} is supplied, generates competing-event histories by
  independent Bernoulli draws from the refitted models' response-scale
  probabilities, in the path order of the authors' reference code.
  \item For hazard effects, draws one multivariate-normal vector for the
  constant Aalen coefficients from \code{yreg_boot$gamma} and
  \code{yreg_boot$robvar.gamma}. For survival-probability effects, obtains
  survival from \code{timereg::predict.aalen()} at \code{t} without this
  additional coefficient draw.
  \item Forms \code{DE}, \code{IED}, \code{IEM}, and \code{TE} from the four
  counterfactual means described above.
}

The displayed point summary is the median of the bootstrap effect values; it
is not a separate plug-in estimate computed from the original, unresampled
data. The interval is the equal-tailed percentile interval with probabilities
\code{(1 - conf_level) / 2} and \code{1 - (1 - conf_level) / 2}, computed by
\code{stats::quantile()} using its default quantile type. This is a percentile
interval, not a normal, basic, studentized, or BCa interval.
}
\section{Parallel execution and reproducibility}{
The function has no seed argument. Call \code{set.seed()} before the function
if reproducible results are needed. With \code{ncores = 1}, bootstrap
replicates run sequentially using the caller's random-number stream. With
\code{ncores > 1} and no registered parallel backend (or only a registered
\code{doSEQ} backend), the function creates a \pkg{parallel} cluster with
\code{ncores} workers through \pkg{doParallel} and stops it on return or
error. An already registered parallel \pkg{foreach} backend is used instead.

For parallel execution, the function draws one seed per bootstrap replicate
from the caller's random-number stream and sets that replicate's seed inside
its worker. Thus repeated parallel calls with the same input, seed, worker
setup, and package versions are reproducible across worker counts. Parallel
and sequential calls do not generally produce identical draws under the same
initial seed. For comparison with the authors' sequential code, use
\code{ncores = 1} and match the input rows and row order, model
specifications, arguments, seed, R version, and relevant package versions.
}
\section{Validation, missing values, and bootstrap failures}{
Input structure and initial fitted objects are checked before bootstrap
execution. The mediator and competing-event fits must have finite,
full-rank coefficients. Binomial models must have converged, contain both
outcome values among positive-weight observations, and avoid fitted
probabilities at numerical boundaries. Covariance matrices used for
multivariate-normal draws must be finite, symmetric, and positive
semidefinite. The Aalen model must provide finite, named constant
coefficients and a compatible counterfactual design matrix.

The primary-event and time-point variables have the completeness requirements
listed under their arguments. Other missing values are governed by the
supplied model calls and their \code{na.action}; however, every active
counterfactual history must ultimately yield finite and complete mediator,
competing-event, and outcome predictions.

All models are checked again after each bootstrap refit. A resample can, for
example, become singular, lose one competing-event response category, exhibit
separation, or produce an invalid prediction even when the original fit is
valid. Such a replicate stops the call with an error identifying the affected
model or bootstrap iteration. Invalid replicates are not silently discarded,
replaced, or retried, because doing so would change the authors' method.
}
\value{
A named list, in the order \code{DE}, \code{IEM}, \code{IED}, and \code{TE}.
Each component is a length-one character string of the form
\code{"median (lower, upper)"}. Values are rounded to at most three decimal
places, so trailing zeroes may not be printed.

\describe{
  \item{DE}{Path-specific direct effect.}
  \item{IEM}{Indirect effect through the longitudinal mediator process.}
  \item{IED}{Indirect effect through the competing-event process.}
  \item{TE}{Total effect.}
}

Hazard-scale values include the \code{peryr} multiplier. Survival-scale
values are percentage-point differences. The function returns formatted
character summaries rather than the underlying replicate-level numeric
matrix. It does not calculate or return a mediation ratio such as
\code{IEM / TE}.
}
\note{
The function implements the supplied model specifications and the authors'
algorithm; it does not by itself establish causal identification. Scientific
interpretation additionally requires the identification, consistency,
positivity, temporal-ordering, and model-correctness assumptions appropriate
to the path-specific estimands and study design.
}
\examples{
\dontrun{
## Required wide layout: one row per individual. This example assumes dat
## already contains E, C, M_1--M_3, D1--D3, time_to_event, eventHappened,
## and the three mediator measurement-time columns.
mvar <- c("M_1", "M_2", "M_3")
time_points <- c("time_1", "time_2", "time_3")

## Continuous mediator models in visit order.
mreg <- list(
  stats::lm(M_1 ~ E + C, data = dat),
  stats::lm(M_2 ~ E + M_1 + C, data = dat),
  stats::lm(M_3 ~ E + M_2 + C, data = dat)
)

## Binary competing-event models in the same order. If risk-set restrictions
## are needed, encode them in dat and use subset= or NA responses.
dreg <- list(
  stats::glm(D1 ~ E + M_1 + C, data = dat,
             family = stats::binomial("logit")),
  stats::glm(D2 ~ E + M_2 + C, data = dat,
             family = stats::binomial("logit")),
  stats::glm(D3 ~ E + M_3 + C, data = dat,
             family = stats::binomial("logit"))
)

## Construct the observed counting-process data used to fit yreg.
dat$id_for_yreg <- seq_len(nrow(dat))
outcome_long <- survival::tmerge(
  dat, dat, id = id_for_yreg,
  endpt = event(time_to_event, eventHappened)
)
mediator_long <- stats::reshape(
  dat,
  direction = "long",
  varying = list(mvar, time_points),
  v.names = c("M", "mediator_time"),
  timevar = "visit",
  times = seq_along(mvar),
  idvar = "id_for_yreg"
)
mediator_long <- mediator_long[
  order(mediator_long$id_for_yreg, mediator_long$mediator_time),
]
outcome_long <- survival::tmerge(
  outcome_long, mediator_long, id = id_for_yreg,
  M = tdc(mediator_time, M)
)

## All outcome predictors must use constant Aalen coefficients.
library(timereg)
yreg <- timereg::aalen(
  survival::Surv(tstart, tstop, endpt) ~
    const(E) + const(M) + const(C),
  data = outcome_long,
  id = NULL,
  n.sim = 0,
  robust = 1,
  resample.iid = 1
)

exposure_values <- unname(stats::quantile(dat$E, c(0.25, 0.75)))

## Use one core for comparison with the reference code.
set.seed(111)
fit <- cmest_pathcomprisk(
  dreg = dreg,
  mreg = mreg,
  mvar = mvar,
  yreg = yreg,
  avar = "E",
  a = exposure_values[1],
  astar = exposure_values[2],
  data = dat,
  time_points = time_points,
  yreg_time = "time_to_event",
  yreg_event = "eventHappened",
  yreg_mvar = "M",
  estimate = "hazard",
  peryr = 100000,
  nboot = 1000,
  conf_level = 0.95,
  ncores = 1
)
}
}
\seealso{
\code{\link[stats]{lm}}, \code{\link[stats]{glm}},
\code{\link[stats]{binomial}}, \code{\link[survival]{Surv}},
\code{\link[survival]{tmerge}}, \code{\link[timereg]{aalen}},
\code{\link[timereg]{predict.aalen}},
\code{\link[foreach]{registerDoSEQ}},
\code{\link[parallel]{makeCluster}}, and
\code{\link[base]{set.seed}}.
}
