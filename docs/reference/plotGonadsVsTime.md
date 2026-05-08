# Plot per-capita gonadic mass through time at a given size

This function plots the gonadic mass of an individual of a given size
through time for all species.

## Usage

``` r
plotGonadsVsTime(sim, time_range, sizes = sim@params@species_params$w_max)
```

## Arguments

- sim:

  A MizerSim object

- time_range:

  The time range over which to plot

- sizes:

  A vector of sizes for each species at which to plot the gonadic mass
  of an individual of that size.

## See also

Other plotting functions:
[`animateGonadSpectra()`](https://gustavdelius.github.io/mizerSeasonal/reference/animateGonadSpectra.md),
[`plotRDD()`](https://gustavdelius.github.io/mizerSeasonal/reference/plotRDD.md),
[`plotRDI()`](https://gustavdelius.github.io/mizerSeasonal/reference/plotRDI.md)
