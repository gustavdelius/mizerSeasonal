#' Set reproduction with vonMises distribution
#'
#' This returns a new model in which the reproduction rate varies throughout
#' the year so that
#'
#' @param params A MizerParams object
#' @param pulse_time A number between 0 and 1 that specifies the time within
#'   each year at which reproduction should take place.
#'
#' @return A MizerParams object with von Mises reproduction rate
#' @export
setVonMisesReproduction <- function(params, pulse_time) {
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
