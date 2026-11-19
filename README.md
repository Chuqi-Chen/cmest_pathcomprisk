# cmest_pathcomprisk() in this CMAverse fork

## Development fork: cmest_pathcomprisk()

[![Command tests](https://github.com/Chuqi-Chen/cmest_pathcomprisk/actions/workflows/pathcomprisk-check.yaml/badge.svg)](https://github.com/Chuqi-Chen/cmest_pathcomprisk/actions/workflows/pathcomprisk-check.yaml)

This is Chuqi Chen's development fork of CMAverse, adding the proposed
`cmest_pathcomprisk()` command under Linda Valeri's supervision. The package
name remains **CMAverse**. This command is not yet part of the upstream release.
The statistical method is by Domingo-Relloso, Zhang, Wang, and colleagues
(2026), [doi:10.1002/sim.70425](https://doi.org/10.1002/sim.70425).

Install this fork, then use a new R session:

```r
install.packages("remotes")  # once, if needed
remotes::install_github("Chuqi-Chen/cmest_pathcomprisk", ref = "master",
                        build_vignettes = FALSE, upgrade = "never")
library(CMAverse)
?cmest_pathcomprisk
```

After fitting your models, call the command directly:

```r
nboot <- 1000L
ncores <- 4L
fit <- cmest_pathcomprisk(
  dreg = list(L1, L2, L3), mreg = list(M1, M2, M3), mvar = mvar,
  yreg = Y, avar = "E", a = 1.434303, astar = 2.569418,
  data = dat, time_points = time_points,
  yreg_time = "time_to_event", yreg_event = "eventHappened", yreg_mvar = "M",
  estimate = "hazard", peryr = 100000, nboot = nboot, ncores = ncores
)
```

`M1`--`M3` are fitted continuous-mediator `lm` models; `L1`--`L3` are fitted
binomial competing-event models; `Y` is a fitted `timereg::aalen` model.
See `?cmest_pathcomprisk` for data, risk-set and formula requirements.
Call `set.seed()` before the analysis if you want reproducibility.
The output contains `DE`, `IEM`, `IED`, and `TE`; there is no `Q`.

A complete runnable example uses bundled synthetic data and 20 bootstrap
replicates for demonstration:

```r
source(system.file("examples", "pathcomprisk_example.R", package = "CMAverse"))
```

Use a larger bootstrap count for inference and assess stability. Matching the
professor's functions tests implementation fidelity; it does not prove
unbiasedness under all simulation models. See [validation notes](validation/pathcomprisk/README.md).


