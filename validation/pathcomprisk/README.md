# Provenance and scope of validation

## Research inputs

The two files in `R/` are byte-for-byte copies of the current standalone
`cmest_pathcomprisk.R` and `cmest_pathcomprisk_validation.R` in the research
folder. Their hashes are in `source-manifest.csv`. Packaging does not change
the estimator body. The detailed Rd was copied from the research `man/` folder;
its stale default of four workers was corrected to the actual default of one.
The duplicate old helper and obsolete `pathcomprisk()` alias were removed
from the old package skeleton.

The professor's method functions are parsed from
`professor_reference/mediation_path_specific.R`, without executing its
application, and stored in the test fixture. the contribution repository extraction script
rebuilds that fixture and its source hash. The original code takes precedence
over paper/code discrepancies, as requested. For the one-timepoint hazard
reference, a test-only column-name lookup maps long `M` to wide `M_1`, because
that function skips its reshape when there is only one mediator. No scientific
formula or stochastic calculation is changed by that naming adapter.

The synthetic example comes from the first three-timepoint simulated dataset
used in the research project. It was prepared as a 3,000-person example with
mediators missing outside the event-free risk sets and interval-specific
incident competing-event indicators. Its creator script is recorded in
the contribution repository example preparation script. It contains no real study participants. It is a
usage example.

The paper is Domingo-Relloso A, Zhang Y, Wang Z, et al. (2026),
"A Path-Specific Effect Approach to Mediation Analysis With Time-Varying
Mediators and Time-to-Event Outcomes Accounting for Competing Risks,"
Statistics in Medicine 45:e70425, doi:10.1002/sim.70425. The CodeReview comments
and two-, three-, and five-timepoint analyses informed the interface and
limitations. Patient data, application materials, PDFs, and 250-run simulation
outputs are not required to execute or integrate this command.

## What the tests check

- Installation, namespace loading, documentation usage, and expected return API.
- Rejection of unsupported mediator, competing-event, and outcome models.
- One to five timepoints on hazard and survival scales: raw bootstrap effects
  and random-number states against the professor's functions.
- Cases without a competing-event model.
- Reproducibility of parallel runs and preservation of a caller-owned cluster.

Bootstrap samples, complete refits (including Aalen `n.sim = 1000`), and
random draws are aligned in reference comparisons. Matching only an initial
seed across arbitrary sequential and parallel implementations does not align
all random draws. Short tests use a few bootstrap samples to check algorithms;
this is not a choice of bootstrap count for inference.

## Reference convention

The one-timepoint survival behavior (the reference skips its first death draw)
is deliberately retained as a source convention.

The earlier temporary-package checks are recorded in `local-check.txt`.
Checks against the complete fork are recorded separately in `fork-check.txt`.
