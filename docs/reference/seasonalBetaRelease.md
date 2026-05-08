# Beta distributed gonad release rate

Uses the formula \$\$r(w, t) = r(t) = f(t-\lfloor t \rfloor)\$\$ where
\\f(t)\\ is the probability density function of the beta distribution
(see see [`?dbeta`](https://rdrr.io/r/stats/Beta.html)) with parameters
`shape1 = beta_a` and `shape2 = beta_b`. Because mizer measures time in
years, \\t-\lfloor t \rfloor\\ gives the time within the year and so
\\r(t)\\ is a periodic function with the period of one year.

## Usage

``` r
seasonalBetaRelease(t, params)
```

## Arguments

- t:

  The time at which to calculate the reproduction rate

- params:

  A MizerParams object

## Value

A vector of species-specific release rates at time t

## See also

Other release functions:
[`seasonalBetaHazardRelease()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalBetaHazardRelease.md),
[`seasonalGaussianRelease()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalGaussianRelease.md),
[`seasonalVonMisesRelease()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalVonMisesRelease.md)
