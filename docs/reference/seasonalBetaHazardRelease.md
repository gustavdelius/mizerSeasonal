# Beta hazard mass-specific gonad release rate

Uses the formula \$\$r(w, t) = r(t) = \frac{f(t-\lfloor t \rfloor)}{1 -
F(t-\lfloor t \rfloor)}\$\$ where \\f(t)\\ and \\F(t)\\ are the
probability density function and the cumulative distribution function of
the beta distribution (see see
[`?dbeta`](https://rdrr.io/r/stats/Beta.html)) with parameters
`shape1 = beta_a` and `shape2 = beta_b`. Because mizer measures time in
years, \\t-\lfloor t \rfloor\\ gives the time within the year and so
\\r(t)\\ is a periodic function with the period of one year.

## Usage

``` r
seasonalBetaHazardRelease(t, params)
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
[`seasonalBetaRelease()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalBetaRelease.md),
[`seasonalGaussianRelease()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalGaussianRelease.md),
[`seasonalVonMisesRelease()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalVonMisesRelease.md)
