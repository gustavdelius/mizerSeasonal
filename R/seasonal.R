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
        setRateFunction("RDD", "noRDD")
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
