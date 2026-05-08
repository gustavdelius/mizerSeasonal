# Seasonal density-independent reproduction rate via S3 dispatch

Method for `projectRDI()` that replaces the standard energy-based RDI
with the gonadic-release calculation used by this extension. Eggs are
produced by the timed release of accumulated gonadic mass rather than by
direct energy investment, so the standard RDI formula does not apply.

## Usage

``` r
# S3 method for class 'mizerSeasonal'
projectRDI(params, n, n_pp, n_other, t = 0, e_growth, mort, e_repro, ...)
```

## Arguments

- params:

  A MizerParams object of class `mizerSeasonal`.

- n:

  A matrix of species abundances (species x size).

- n_pp:

  A vector of the resource abundance by size. Unused.

- n_other:

  A list of other model components, including `gonads`.

- t:

  The current time.

- e_growth, mort, e_repro:

  Unused; present for generic compatibility.

- ...:

  Unused

## Value

A numeric vector with the rate of egg production for each species.
