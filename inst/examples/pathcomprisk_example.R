# Runnable after installing this development fork of CMAverse.
# ---- setup ----
library(survival)
library(timereg)
library(CMAverse)

# ---- data ----
# Synthetic data; one row per subject; follow-up is in years.
dat <- read.csv(system.file("extdata", "pathcomprisk_example.csv",
                            package = "CMAverse"), na.strings = "NA")
mvar <- c("M_1", "M_2", "M_3")
time_points <- c("time_1", "time_2", "time_3")

# ---- models ----
mreg <- list(
  lm(M_1 ~ E + C, data = dat),
  lm(M_2 ~ E + M_1 + C, data = dat, subset = risk_2),
  lm(M_3 ~ E + M_2 + C, data = dat, subset = risk_3)
)
dreg <- list(
  glm(D1 ~ E + M_1 + C, data = dat, family = binomial("logit")),
  glm(D2 ~ E + M_2 + C, data = dat,
      subset = risk_2, family = binomial("logit")),
  glm(D3 ~ E + M_3 + C, data = dat,
      subset = risk_3, family = binomial("logit"))
)

# ---- outcome ----
outcome_long <- tmerge(
  dat, dat, id = id,
  endpt = event(time_to_event, eventHappened)
)
mediator_long <- reshape(
  dat, direction = "long",
  varying = list(mvar, time_points),
  v.names = c("M", "mediator_time"),
  timevar = "visit", times = seq_along(mvar), idvar = "id"
)
mediator_long <- mediator_long[
  order(mediator_long$id, mediator_long$mediator_time),
]
outcome_long <- tmerge(
  outcome_long, mediator_long, id = id,
  M = tdc(mediator_time, M)
)
set.seed(2026)
yreg <- aalen(
  Surv(tstart, tstop, endpt) ~ const(E) + const(M) + const(C),
  data = outcome_long, n.sim = 0, robust = 1, resample.iid = 1
)

# ---- estimate ----
# 20 replicates make a short walkthrough, not a stable final inference.
# For a full analysis, begin with nboot = 1000 and assess stability.
# set.seed() is optional and belongs in the analysis, not in the command.
set.seed(2026)
fit <- cmest_pathcomprisk(
  data = dat, avar = "E", a = 1.434303, astar = 2.569418,
  mreg = mreg, mvar = mvar, dreg = dreg, yreg = yreg,
  time_points = time_points,
  yreg_time = "time_to_event", yreg_event = "eventHappened",
  yreg_mvar = "M", estimate = "hazard", peryr = 100000,
  nboot = 20, conf_level = 0.95, ncores = 4
)

# ---- output ----
print(fit)

# Optional alternative: repeat the same call with estimate = "sp", t = 20.
# That scale reports predicted survival differences in percentage points.
