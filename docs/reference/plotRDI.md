# Plot density-independent reproduction rate over time

Plot density-independent reproduction rate over time

## Usage

``` r
plotRDI(
  sim,
  sim2,
  species = NULL,
  total = FALSE,
  log = FALSE,
  highlight = NULL,
  return_data = FALSE,
  ...
)
```

## Arguments

- sim:

  An object of class MizerSim

- sim2:

  An optional second object of class MizerSim. If this is provided its
  RDIs will be shown on the same plot in bolder lines.

- species:

  The species to be selected. Optional. By default all target species
  are selected. A vector of species names, or a numeric vector with the
  species indices, or a logical vector indicating for each species
  whether it is to be selected (TRUE) or not.

- total:

  A boolean value that determines whether the total over all species in
  the system is plotted as well. Note that even if the plot only shows a
  selection of species, the total is including all species. Default is
  FALSE.

- log:

  Boolean whether RDI should be plotted on a logarithmic axis. Defaults
  to true.

- highlight:

  Name or vector of names of the species to be highlighted.

- return_data:

  A boolean value that determines whether the formatted data used for
  the plot is returned instead of the plot itself. Default value is
  FALSE

- ...:

  Other arguments (currently unused)

## Value

A ggplot2 object, unless return_data = TRUE, in which case a data frame
with the three variables 'Year', 'RDI', 'Species' is returned.

## See also

Other plotting functions:
[`animateGonadSpectra()`](https://gustavdelius.github.io/mizerSeasonal/reference/animateGonadSpectra.md),
[`plotGonadsVsTime()`](https://gustavdelius.github.io/mizerSeasonal/reference/plotGonadsVsTime.md),
[`plotRDD()`](https://gustavdelius.github.io/mizerSeasonal/reference/plotRDD.md)
