# Namespace declaration for the top-level NSE registration in the estimator.
#' @importFrom utils globalVariables
NULL

# NSE symbols used by the inherited foreach/dplyr multistate implementation.
utils::globalVariables(c("index", "RD", "SD", "TE", "time", "Haz"))
