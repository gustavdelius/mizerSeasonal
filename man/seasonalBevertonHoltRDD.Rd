% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seasonal.R
\name{seasonalBevertonHoltRDD}
\alias{seasonalBevertonHoltRDD}
\title{Beverton Holt function to calculate density-dependent reproduction rate}
\usage{
seasonalBevertonHoltRDD(rdi, params, t, ...)
}
\arguments{
\item{rdi}{Vector of density-independent reproduction rates
\eqn{R_{di}}{R_di} for all species.}

\item{params}{A params object that must contain \code{other_params$R_max}}

\item{t}{The time at which to calculate RDD}

\item{...}{Unused}
}
\value{
Vector of density-dependent reproduction rates.
}
\description{
Takes the density-independent rates \eqn{R_{di}(t)}{R_di(t)} of egg
production and returns reduced, density-dependent reproduction rates
\eqn{R_{dd}(t)}{R_dd(t)} given as
\deqn{R_{dd}(t) = R_{di}(t)
\frac{R_{max}(t)}{R_{di}(t) + R_{max}(t)}}{R_dd(t) = R_di(t) R_max(t)/(R_di(t) + R_max(t))}
where \eqn{R_{max}(t)}{R_max(t)} are the maximum possible reproduction rates
that must be specified in an array (time x species) saved in the
\code{other_params$R_max} slot of the \code{MizerParams} object. This is simply a
time-dependent version of \code{mizer::BevertonHoltRDD()}. (All quantities in the
above equation are species-specific but we dropped the species index for
simplicity.)
}
\seealso{
Other functions calculating density-dependent reproduction rates: 
\code{\link{seasonalVonMisesRDD}()}
}
\concept{functions calculating density-dependent reproduction rates}
