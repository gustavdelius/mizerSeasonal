% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plotGonadsVsTime}
\alias{plotGonadsVsTime}
\title{Plot per-capita gonadic mass through time at a given size}
\usage{
plotGonadsVsTime(sim, time_range, sizes = sim@params@species_params$w_max)
}
\arguments{
\item{sim}{A MizerSim object}

\item{time_range}{The time range over which to plot}

\item{sizes}{A vector of sizes for each species at which to plot the
gonadic mass of an individual of that size.}
}
\description{
Plot per-capita gonadic mass through time at a given size
}
