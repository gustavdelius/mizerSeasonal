# Get time series of summary functions from simulation

Given a function that can calculate a vector of quantities of interest
at a single time for all species, this function returns an array (time x
species) with the values at all the time steps saved in a simulation.

## Usage

``` r
getTimeseries(sim, func = getRDI, ...)
```

## Arguments

- sim:

  A MizerSim object

- func:

  The function calculating the quantities at a single time step

- ...:

  Unused

## Value

A matrix (time x species)
