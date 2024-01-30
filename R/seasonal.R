#' Set seasonal reproduction
#'
#' This returns a new model in which the reproduction rate varies throughout
#' the year
#'
#' @param params A MizerParams object
#' @param repro_func Name of the function giving the time-dependent
#'     mass-specific reproduction rate.
#'
#' @return A MizerParams object with seasonal reproduction
#' @export
setSeasonalReproduction <- function(params, repro_func = "repro_gaussian") {
    # start with zero gonadic mass
    initial <- initialN(params) # to get the right dimensions
    initial[] <- 0

    # TODO: check that `repro_func` is valid and that all necessary parameters
    #  are contained in the species parameters
    other_params(params)$repro_func <- repro_func

    # add gonads component and register new RDI function
    setComponent(params, component = "gonads",
                 initial_value = initial,
                 dynamics_fun = "gonadDynamics") |>
        setRateFunction("RDI", "seasonalRDI") |>
        # Turn of density dependence for now
        setRateFunction("seasonalRDD", "seasonalBevertonHoltRDD")
}


#' Update the gonadic mass using a seasonal reproduction rate
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
gonadDynamics <- function(params, n_other, rates, t, dt, ...) {
    # Handy things
    no_sp <- nrow(params@species_params) # number of species
    no_w <- length(params@w) # number of fish size bins
    repro_func <- get0(other_params(params)$repro_func)
    r <- repro_func(t, params)
    # Matrices for solver
    # a_{ij} = - g_i(w_j) / dw_j dt
    a <- sweep(-rates$e_growth * dt, 2, params@dw, "/")
    # b_{ij} = 1 + g_i(w_j) / dw_j dt + r_i(w_j) dt
    b <- 1 + sweep(rates$e_growth * dt, 2, params@dw, "/") + r * dt
    # S_{ij} <- q_i(w_j) + rates$e_repro dt
    s <- n_other$gonads + rates$e_repro * dt
    # Update q
    # for (i in 1:no_sp) # number of species assumed small, so no need to
    #                      vectorize this loop over species
    #     for (j in (params@w_min_idx[i]+1):no_w)
    #         q[i,j] <- (s[i,j] - a[i,j]*q[i,j-1]) / b[i,j]
    # This is implemented via Rcpp
    mizer:::inner_project_loop(no_sp = no_sp, no_w = no_w,
                       n = n_other$gonads,
                       A = a, B = b, S = s,
                       w_min_idx = params@w_min_idx)
}

#' Get density-independent rate of seasonal reproduction
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
seasonalRDI <- function(params, n, n_other, t, dt = 0.1, ...) {
    repro_func <- get0(other_params(params)$repro_func)
    r <- repro_func(t, params)
    total <- drop((sweep(n_other$gonads, 1, r, "*") * n) %*% params@dw)
    # Assume sex_ratio = 0.5.
    0.5 * (total * params@species_params$erepro) /
        params@w[params@w_min_idx]
}

#' Beverton Holt function to calculate density-dependent reproduction rate
#'
#' Takes the density-independent rates \eqn{R_{di}}{R_di} of egg production (as
#' calculated by [getRDI()]) and returns
#' reduced, density-dependent reproduction rates \eqn{R_{dd}}{R_dd} given as
#' \deqn{R_{dd} = R_{di}
#' \frac{R_{max}}{R_{di} + R_{max}}}{R_dd = R_di R_max/(R_di + R_max)} where
#' \eqn{R_{max}}{R_max} are the maximum possible reproduction rates that must be
#' specified in a column in the species parameter dataframe.
#' (All quantities in the above equation are species-specific but we dropped
#' the species index for simplicity.)
#'
#' This is only one example of a density-dependence. You can write your own
#' function based on this example, returning different density-dependent
#' reproduction rates. Three other examples provided are [RickerRDD()],
#' [SheperdRDD()], [noRDD()] and [constantRDD()]. For more explanation see
#' [setReproduction()].
#'
#' @param rdi Vector of density-independent reproduction rates
#'   \eqn{R_{di}}{R_di} for all species.
#' @param species_params A species parameter dataframe. Must contain a column
#'   `R_max` holding the maximum reproduction rate \eqn{R_{max}}{R_max} for each
#'   species.
#' @param ... Unused
#'
#' @return Vector of density-dependent reproduction rates.
#' @export
#' @family functions calculating density-dependent reproduction rate
seasonalBevertonHoltRDD <- function(rdi, params, t, ...) {
    if (!("r_max" %in% names(other_params))) {
        stop("The r_max paramter is missing.")
    }
    return(rdi / (1 + rdi / other_params$r_max))
}

#' Seasonal version of mizerRates()
#'
#' This is needed because we need to pass the time argument to the RDD function
#'
#' @param params A \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size).
#' @param n_pp A vector of the resource abundance by size
#' @param n_other A list of abundances for other dynamical components of the
#'   ecosystem
#' @param t The time for which to do the calculation (Not used by standard
#'   mizer rate functions but useful for extensions with time-dependent
#'   parameters.)
#' @param effort The effort for each fishing gear
#' @param rates_fns Named list of the functions to call to calculate the rates.
#'   Note that this list holds the functions themselves, not their names.
#' @param ... Unused
#' @return List of rates.
#' @export
seasonalRates <- function(params, n, n_pp, n_other,
                          t = 0, effort, rates_fns, ...) {
    r <- mizerRates(params = params, n = n, n_pp = n_pp, n_other = n_other,
                    t = t, effort = effort, rates_fns = rates_fns)
    # Overwrite rdd with call to seasonal RDD
    r$rdd <- rates_fns$seasonalRDD(rdi = r$rdi, params = params, t = t, ...)
    return(r)
}

#' Gaussian mass-specific reproduction rate
#'
#' Uses the formula
#' \deqn{r(w, t) = r(t) = r_0 \exp{\left(-\dfrac{(t-\lfloor{t}\rfloor-t_0)^2}{2\sigma^2}\right)}}{r(w,t) = r(t)= r_0 \exp(-(t-|_t_|-t_0)^2/(2\sigma^2))}
#' where parameters \eqn{r_0, \sigma} and \eqn{t_0} are given by new columns
#' `sr_r0`, `sr_sigma` and `sr_t0` in the species parameter data frame.
#' Because mizer measures time in years, \eqn{t-\lfloor t \rfloor}{t-|_t_|} gives the time
#' within the year and so \eqn{r(t)} is a periodic function with
#' the period of one year.
#'
#' @param t The time at which to calculate the reproduction rate
#' @param params A MizerParams object
#'
#' @return A vector of species-specific reproduction rates
#' @export
repro_gaussian <- function(t, params) {
    sp <- params@species_params
    sp$sr_r0 * exp(-(t - trunc(t) - sp$sr_t0)^2/(2*sp$sr_sigma^2))
}
