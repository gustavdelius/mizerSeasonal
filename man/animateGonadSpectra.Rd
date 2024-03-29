% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/animate.R
\name{animateGonadSpectra}
\alias{animateGonadSpectra}
\title{Animation of the gonadic mass spectra}
\usage{
animateGonadSpectra(
  sim,
  species = NULL,
  time_range,
  wlim = c(NA, NA),
  ylim = c(NA, NA),
  power = 0
)
}
\arguments{
\item{sim}{A MizerSim object}

\item{species}{Name or vector of names of the species to be plotted. By
default all species are plotted.}

\item{time_range}{The time range to animate over. Either a vector of values
or a vector of min and max time. Default is the entire time range of the
simulation.}

\item{wlim}{A numeric vector of length two providing lower and upper limits
for the w axis. Use NA to refer to the existing minimum or maximum.}

\item{ylim}{A numeric vector of length two providing lower and upper limits
for the y axis. Use NA to refer to the existing minimum or maximum. Any
values below 1e-20 are always cut off.}

\item{power}{The gonadic mass density is plotted as the number density times
the per-capita gonadic mass times the weight raised to \code{power}. The
default \code{power = o} gives the gonadic mass density, whereas \code{power =
  1} gives the gonadic mass density with respect to logarithmic size bins.}
}
\value{
A plotly object
}
\description{
Animation of the gonadic mass spectra
}
\seealso{
Other plotting functions: 
\code{\link{plotRDD}()},
\code{\link{plotRDI}()}
}
\concept{plotting functions}
