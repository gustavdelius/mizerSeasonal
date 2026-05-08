# Get encounter rate that includes gonadic mass of prey

This is doing the same as the `mizerEncounter()` function in core mizer
except that the prey mass is the sum of its somatic mass \\w_p\\ and its
gonadic mass \\q(w_p)\\. This function is automatically registered as
the model's encounter rate function when
[`setSeasonalReproduction()`](https://gustavdelius.github.io/mizerSeasonal/reference/setSeasonalReproduction.md)
is called with `include_gonads = TRUE`.

## Usage

``` r
seasonalEncounter(params, n, n_pp, n_other, t, ...)
```

## Arguments

- params:

  A MizerParams object

- n:

  A matrix of species abundances (species x size).

- n_pp:

  A vector of the resource abundance by size

- n_other:

  A list of abundances for other dynamical components of the ecosystem

- t:

  The time for which to do the calculation

- ...:

  Unused

## Value

A named two dimensional array (predator species x predator size) with
the encounter rates.
