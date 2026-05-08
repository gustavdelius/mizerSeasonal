# Gaussian mass-specific gonad release rate

Uses the formula \$\$r(w, t) = r(t) = r_0
\exp{\left(-\dfrac{(t-\lfloor{t}\rfloor-t_0)^2}{2\sigma^2}\right)}\$\$
where parameters \\r_0, \sigma\\ and \\t_0\\ are given by new columns
`sr_r0`, `sr_sigma` and `sr_t0` in the species parameter data frame.
Because mizer measures time in years, \\t-\lfloor t \rfloor\\ gives the
time within the year and so \\r(t)\\ is a periodic function with the
period of one year.

## Usage

``` r
seasonalGaussianRelease(t, params)
```

## Arguments

- t:

  The time at which to calculate the release rate

- params:

  A MizerParams object

## Value

A vector of species-specific release rates at time t

## See also

Other release functions:
[`seasonalBetaHazardRelease()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalBetaHazardRelease.md),
[`seasonalBetaRelease()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalBetaRelease.md),
[`seasonalVonMisesRelease()`](https://gustavdelius.github.io/mizerSeasonal/reference/seasonalVonMisesRelease.md)
