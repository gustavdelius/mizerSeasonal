# Update the gonadic mass using a seasonal reproduction rate

Update the gonadic mass using a seasonal reproduction rate

## Usage

``` r
gonadDynamics(params, n_other, rates, t, dt, ...)
```

## Arguments

- params:

  MizerParams object

- n_other:

  Other model components

- rates:

  Previously calculated rates, including in particular e_repro

- t:

  The current time

- dt:

  The time step size

- ...:

  Unused

## Value

Array (species x size) with the current gonadic mass of an individual.
