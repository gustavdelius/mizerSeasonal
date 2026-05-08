# Get density-independent rate of seasonal reproduction

This function calculates the density-independent rate of offspring
production for each species. The rate of offspring production is given
by \$\$R\_{di}(t) =\frac{\epsilon}{2w_0}\int N(w,t)q(w,t)r(w,t)dw\$\$
where \\q(w,t)\\ is the gonadic mass of an individual, \\r(w,t)\\ is the
mass-specific gonad release rate, \\N(w,t)\\ is the abundance density,
\\\epsilon\\ is the reproductive efficiency and \\w_0\\ is the weight of
an offspring.

## Usage

``` r
seasonalRDI(params, n, n_other, t, dt = 0.1, ...)
```

## Arguments

- params:

  MizerParams object

- n:

  Species abundances at current time step

- n_other:

  Other model components at current time step. This will include in
  particular the gonadic biomass in `n_other$gonads`.

- t:

  The current time

- dt:

  The time step size

- ...:

  Unused

## Value

A numeric vector with the rate of egg production for each species.
