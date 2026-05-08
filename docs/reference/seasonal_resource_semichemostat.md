# Seasonal semichemostat resource dynamics

This implements the standard semichemostat dynamics for the resource
(see
[`mizer::resource_semichemostat()`](https://sizespectrum.org/mizer/reference/resource_semichemostat.html))
but with a time-dependent carrying capacity. The carrying capacity is
given by \$\$c_R(w,t) = (1 + \text{maxR} \cdot \text{vonMises}(t, \mu,
\kappa)) c_R(w)\$\$ where \\c_R(w)\\ is the standard carrying capacity,
\\\mu\\ is the mean of the von-Mises distribution, and \\\kappa\\ is the
concentration parameter of the von-Mises distribution. These parameters
must be supplied in the `rp$maxR`, `rp$mu` and `rp$kappa` slots of the
`resource_params` slot of the `params` object.

## Usage

``` r
seasonal_resource_semichemostat(
  params,
  n,
  n_pp,
  n_other,
  rates,
  t,
  dt,
  resource_rate,
  resource_capacity,
  ...
)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
  object

- n:

  A matrix of species abundances (species x size)

- n_pp:

  A vector of the resource abundance by size

- n_other:

  A list with the abundances of other components

- rates:

  A list of rates as returned by
  [`mizerRates()`](https://sizespectrum.org/mizer/reference/mizerRates.html)

- t:

  The current time

- dt:

  Time step

- resource_rate:

  Resource replenishment rate

- resource_capacity:

  Resource carrying capacity

- ...:

  Unused

## Value

Vector containing resource spectrum at next timestep
