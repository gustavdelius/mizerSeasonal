# Beverton Holt function to calculate density-dependent reproduction rate

Takes the density-independent rates \\R\_{di}(t)\\ of egg production and
returns reduced, density-dependent reproduction rates \\R\_{dd}(t)\\
given as \$\$R\_{dd}(t) = R\_{di}(t) \frac{R\_{max}(t)}{R\_{di}(t) +
R\_{max}(t)}\$\$ where \\R\_{max}(t)\\ are the maximum possible
reproduction rates that must be specified in an array (time x species)
saved in the `other_params$R_max` slot of the `MizerParams` object. This
is simply a time-dependent version of
[`mizer::BevertonHoltRDD()`](https://sizespectrum.org/mizer/reference/BevertonHoltRDD.html).
(All quantities in the above equation are species-specific but we
dropped the species index for simplicity.)

## Usage

``` r
seasonalBevertonHoltRDD(rdi, params, t, ...)
```

## Arguments

- rdi:

  Vector of density-independent reproduction rates \\R\_{di}\\ for all
  species.

- params:

  A params object that must contain `other_params$R_max`

- t:

  The time at which to calculate RDD

- ...:

  Unused

## Value

Vector of density-dependent reproduction rates.

## See also

Other functions calculating density-dependent reproduction rates:
[`seasonalVonMisesRDD()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalVonMisesRDD.md)
