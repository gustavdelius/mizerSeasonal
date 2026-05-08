# Set seasonal reproduction

This returns a new model in which the reproduction rate varies
throughout the year. See details below.

## Usage

``` r
setSeasonalReproduction(
  params,
  release_func = "seasonalVonMisesRelease",
  RDD = "seasonalBevertonHoltRDD",
  include_gonads = TRUE
)
```

## Arguments

- params:

  A MizerParams object

- release_func:

  Name of the function giving the time-dependent mass-specific release
  rate. This function should take a time and a MizerParams object as
  arguments and return a vector with the reproduction rates for all
  species.

- RDD:

  Name of the function for calculating the density-dependent
  reproduction rate RDD.

- include_gonads:

  Boolean. If TRUE (default) then the gonadic mass is included in the
  prey encounter rate.

## Value

A MizerParams object with seasonal reproduction

## Details

This function does not change the rate at which fish invest energy into
reproduction. It however changes what happens with this investment.
Rather than being released immediately to produce offspring, it is used
to accumulate gonadic mass. You specify a function that gives the
time-dependent mass-specific rate \\r(t)\\ at which this gonadic mass is
then released for reproduction.

The package provides several candidate functions for calculating release
rates:
[`seasonalVonMisesRelease()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalVonMisesRelease.md),
[`seasonalBetaHazardRelease()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalBetaHazardRelease.md),
[`seasonalBetaRelease()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalBetaRelease.md)
and
[`seasonalGaussianRelease()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalGaussianRelease.md),
and you can use these as templates for writing your own release
functions.

You also specify a non-linear function that calculates the
density-dependent reproduction rate \\R\_{dd}\\ from the
density-independent rate \\R\_{di}\\ of egg production.
