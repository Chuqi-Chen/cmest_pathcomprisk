\name{getvarnames}
\alias{getvarnames}
\title{Extract Variable Names from a Regression Formula}
\usage{getvarnames(formula, data = NULL)}
\arguments{
\item{formula}{A formula, fitted model, model call, character vector, or \code{NULL}.}
\item{data}{Optional data used to expand \code{.} in a formula.}
}
\value{A list with \code{varnames}, \code{xvar} (predictors), and \code{yvar} (responses).}
\description{Extracts the variable names used by the mediation implementation, including formulas stored in Aalen model calls.}
\examples{
getvarnames(y ~ x + z)
getvarnames(y ~ ., data = data.frame(y = 1:3, x = 4:6))
}
