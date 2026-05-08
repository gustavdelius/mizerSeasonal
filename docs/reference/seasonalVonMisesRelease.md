# von-Mises distributed gonad release rate

Uses the formula \$\$r(w, t) = r(t) = r_0 \frac{\exp(\kappa
\cos(2\pi(t - \mu)))}{2\pi I_0(\kappa)}\$\$ where the parameters are
taken from the `vonMises_r0`, `vonMises_kappa` and `vonMises_mu` columns
in the `species_params` data frame in `params`.

## Usage

``` r
seasonalVonMisesRelease(t, params, ...)
```

## Arguments

- t:

  The time at which to calculate the release rate

- params:

  A MizerParams object

- ...:

  Unused

## Value

A vector of species-specific release rates at time t

## See also

Other release functions:
[`seasonalBetaHazardRelease()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalBetaHazardRelease.md),
[`seasonalBetaRelease()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalBetaRelease.md),
[`seasonalGaussianRelease()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalGaussianRelease.md)
