# von-Mises distributed reproduction rate independent of abundance

This function calculates the mass-specific reproduction rate of a
species at a given time. The reproduction rate is given by \$\$r(t) =
r_0 \frac{\exp(\kappa \cos(2\pi(t - \mu)))}{2\pi I_0(\kappa)}\$\$ where
\\r_0\\ is the amplitude of the reproduction rate, \\\kappa\\ is the
concentration parameter of the von-Mises distribution, and \\\mu\\ is
the mean of the von-Mises distribution. These parameters must be
supplied in the slots `rdd_vonMises_r0`, `rdd_vonMises_kappa` and
`rdd_vonMises_mu` of the `species_params` data frame of the `params`
object.

## Usage

``` r
seasonalVonMisesRDD(params, t, ...)
```

## Arguments

- params:

  A MizerParams object

- t:

  The time at which to calculate RDD

- ...:

  Unused

## Value

A vector of species-specific reproduction rates

## Details

Note that this reproduction rate is independent of the rate at which the
species release gonadic mass for reproduction and is thus only to be
used as an aid during the process of setting up the model. It should be
replaced by another reproduction rate function before simulating the
dynamics of the model.

## See also

Other functions calculating density-dependent reproduction rates:
[`seasonalBevertonHoltRDD()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalBevertonHoltRDD.md)
