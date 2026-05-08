# Animation of the gonadic mass spectra

This function creates an animation of the gonadic mass spectra of the
species in the simulation. The gonadic mass spectrum is the size
distribution of the total gonadic mass of each species.

## Usage

``` r
animateGonadSpectra(
  sim,
  species = NULL,
  time_range,
  wlim = c(NA, NA),
  ylim = c(NA, NA),
  power = 0
)
```

## Arguments

- sim:

  A MizerSim object

- species:

  Name or vector of names of the species to be plotted. By default all
  species are plotted.

- time_range:

  The time range to animate over. Either a vector of values or a vector
  of min and max time. Default is the entire time range of the
  simulation.

- wlim:

  A numeric vector of length two providing lower and upper limits for
  the w axis. Use NA to refer to the existing minimum or maximum.

- ylim:

  A numeric vector of length two providing lower and upper limits for
  the y axis. Use NA to refer to the existing minimum or maximum. Any
  values below 1e-20 are always cut off.

- power:

  The gonadic mass density is plotted as the number density times the
  per-capita gonadic mass times the weight raised to `power`. The
  default `power = o` gives the gonadic mass density, whereas
  `power = 1` gives the gonadic mass density with respect to logarithmic
  size bins.

## Value

A plotly object

## See also

Other plotting functions:
[`plotGonadsVsTime()`](https://gustavdelius.github.io/mizerSeasonal/reference/plotGonadsVsTime.md),
[`plotRDD()`](https://gustavdelius.github.io/mizerSeasonal/reference/plotRDD.md),
[`plotRDI()`](https://gustavdelius.github.io/mizerSeasonal/reference/plotRDI.md)
