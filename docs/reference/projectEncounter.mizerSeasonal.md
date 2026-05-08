# Seasonal encounter rate including gonadic mass of prey via S3 dispatch

Method for `projectEncounter()` that adds the extra encounter arising
from the gonadic mass carried by prey fish. It first calls
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html) to obtain the
standard encounter (using somatic mass only), then appends the
contribution from the gonadic mass of prey fish. This is only done when
`include_gonads = TRUE` was passed to
[`setSeasonalReproduction()`](https://gustavdelius.github.io/mizerSeasonal/reference/setSeasonalReproduction.md).

## Usage

``` r
# S3 method for class 'mizerSeasonal'
projectEncounter(params, n, n_pp, n_other, t = 0, ...)
```

## Arguments

- params:

  A MizerParams object of class `mizerSeasonal`.

- n:

  A matrix of species abundances (species x size).

- n_pp:

  A vector of the resource abundance by size.

- n_other:

  A list of other model components, including `gonads`.

- t:

  The current time. Unused.

- ...:

  Passed to [`NextMethod()`](https://rdrr.io/r/base/UseMethod.html).

## Value

A named two-dimensional array (predator species x predator size) with
the encounter rates.
