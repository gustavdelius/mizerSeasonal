#' Update the gonadic mass to the next time step
#'
#' If it is time for the reproductive pulse, the gonadic mass is set to zero,
#' otherwise it is increased by the energy available for reproduction.
#'
#' @param params MizerParams object
#' @param n_other Other model components
#' @param rates Previously calculated rates, including in particular e_repro
#' @param t The current time
#' @param dt The time step size
#' @param ... Unused
#'
#' @return Array (species x size) with the current gonadic mass of an individual.
#' @export
gonadPulsedDynamics <- function(params, n_other, rates, t, dt, ...) {
    if (params@other_params$gonads >= t - trunc(t) &&
        params@other_params$gonads < t - trunc(t) + dt) {
        dims <- dim(n_other$gonads)
        zero_array <- array(0, dims)
        dimnames(zero_array) <- dimnames(n_other$gonads)
        return(zero_array)
    }
    n_other$gonads + rates$e_repro * dt
}

#' Get density-dependent rate of pulsed reproduction
#'
#' The rate of reproduction is zero except during one timestep, during which
#' all gonadic mass is released.
#'
#' @param params MizerParams object
#' @param n Species abundances at current time step
#' @param n_other Other model components at current time step. This will include
#'   in particular the gonadic biomass in `n_other$gonads`.
#' @param t The current time
#' @param dt The time step size
#' @param ... Unused
#'
#' @return A numeric vector with the rate of egg production for each species.
#' @export
pulsedRDI <- function(params, n, n_other, t, dt = 0.1, ...) {
    # Is it time for the reproductive pulse?
    if (params@other_params$gonads >= t - trunc(t) &&
        params@other_params$gonads < t - trunc(t) + dt) {
        # Calculate total mass from per capita gonadic mass
        total <- drop((n_other$gonads * n) %*% params@dw)
        # Assume sex_ratio = 0.5.
        # Divide by dt to convert to a rate, given that all the energy
        # will be released in a single time step
        rdi <- 0.5 * (total * params@species_params$erepro) /
            params@w[params@w_min_idx] / dt
        return(rdi)
    }
    # else return 0 vector
    rep(0, nrow(params@species_params))
}

#' Set pulsed reproduction
#'
#' This returns a new model in which the reproduction is set to take place
#' only within a single time step each year. During the rest of the year the energy
#' available for reproduction is accumulated as gonadic mass, all of which is
#' realeased during the time step into which the `pulse_time` falls.
#'
#' @param params A MizerParams object
#' @param pulse_time A number between 0 and 1 that specifies the time within
#'   each year at which reproduction should take place.
#'
#' @return A MizerParams object with the updated model
#' @export
setPulsedReproduction <- function(params, pulse_time) {
    # start with zero gonadic mass
    initial <- initialN(params) # to get the right dimensions
    initial[] <- 0
    # add gonads component and register new RDI function
    setComponent(params, component = "gonads",
                 initial_value = initial,
                 dynamics_fun = "gonadPulsedDynamics",
                 component_params = pulse_time) |>
        setRateFunction("RDI", "pulsedRDI")
}
