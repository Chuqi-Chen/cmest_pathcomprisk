# cmest_pathcomprisk

`cmest_pathcomprisk` estimates direct, mediator, competing-event, and total
path-specific effects for time-varying mediators and a time-to-event outcome.
It is available in [this CMAverse development fork](https://github.com/Chuqi-Chen/cmest_pathcomprisk),
under Linda Valeri's supervision. The method is described by
[Domingo-Relloso et al. (2026)](https://doi.org/10.1002/sim.70425).

Install the development fork and load CMAverse:

```r
remotes::install_github("Chuqi-Chen/cmest_pathcomprisk", ref = "master",
                        build_vignettes = FALSE, upgrade = "never")
library(CMAverse)
```

After fitting the mediator (`M1`-`M3`), competing-event (`L1`-`L3`), and
Aalen outcome (`Y`) models, run:

```r
nboot <- 1000L
ncores <- 4L
fit <- cmest_pathcomprisk(
  dreg = list(L1, L2, L3), mreg = list(M1, M2, M3), mvar = mvar,
  yreg = Y, avar = "E", a = 0, astar = 1,
  data = dat, time_points = time_points,
  yreg_time = "time_to_event", yreg_event = "eventHappened", yreg_mvar = "M",
  estimate = "hazard", peryr = 100000, nboot = nboot, ncores = ncores
)
```

Here `a = 0` and `astar = 1` are the function's default exposure values;
choose values appropriate to your data. The function's default `ncores` is
one; this example requests four. For model and data requirements, see
`?cmest_pathcomprisk`. Call `set.seed()` before fitting if you want reproducible
results. The command returns `DE`, `IEM`, `IED`, and `TE`. For a runnable
example with included synthetic data:

```r
source(system.file("examples", "pathcomprisk_example.R", package = "CMAverse"))
```

[![Command tests](https://github.com/Chuqi-Chen/cmest_pathcomprisk/actions/workflows/pathcomprisk-check.yaml/badge.svg)](https://github.com/Chuqi-Chen/cmest_pathcomprisk/actions/workflows/pathcomprisk-check.yaml)

