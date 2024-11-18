#' Set seasonal reproduction
#'
#' This returns a new model in which the reproduction rate varies throughout
#' the year
#'
#' @param params A MizerParams object
#' @param release_func Name of the function giving the time-dependent
#'     mass-specific reproduction rate.
#' @param RDD Name of the function for calculating the density-dependent 
#'     reproduction rate RDD.
#' @param include_gonads Boolean. If TRUE (default) then the gonadic mass is
#'     included in the prey encounter rate.
#'
#' @return A MizerParams object with seasonal reproduction
#' @export
setSeasonalReproduction <- function(params, release_func = "repro_vonMises",
                                    RDD = "seasonalBevertonHoltRDD",
                                    include_gonads = TRUE) {
    # start with zero gonadic mass
    initial <- initialN(params) # to get the right dimensions
    initial[] <- 0

    # TODO: check that `release_func` is valid and that all necessary parameters
    #  are contained in the species parameters
    other_params(params)$release_func <- release_func

    # add gonads component and register new RDI function
    p <- setComponent(params, component = "gonads",
                      initial_value = initial,
                      dynamics_fun = "gonadDynamics") |>
        setRateFunction("RDI", "seasonalRDI") |>
        setRateFunction("RDD", RDD)
    # If requested, use Encounter function that includes the gonadic mass
    if (include_gonads) {
        p <- setRateFunction(p, "Encounter", "seasonalEncounter")
    }
    
    return(p)
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
    release_func <- get0(other_params(params)$release_func)
    r <- release_func(t, params)
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

#' Get encounter rate that includes gonadic mass of prey
#' 
#' This is doing the same as the `mizerEncounter()` function in core mizer
#' except that the prey mass is the sum of its somatic mass \eqn{w_p} and
#' its gonadic mass \eqn{q(w_p)}.
#' 
#' @param params A \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size).
#' @param n_pp A vector of the resource abundance by size
#' @param n_other A list of abundances for other dynamical components of the
#'   ecosystem
#' @param t The time for which to do the calculation
#' @param ... Unused
#'   
#' @return A named two dimensional array (predator species x predator size) with
#'   the encounter rates.
#' @export
#' @family mizer rate functions
seasonalEncounter <- function(params, n, n_pp, n_other, t, ...) {
    
    # idx_sp are the index values of params@w_full such that
    # params@w_full[idx_sp] = params@w
    idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
    
    # If the the user has set a custom pred_kernel we can not use fft.
    # In this case we use the code from mizer version 0.3
    if (!is.null(comment(params@pred_kernel))) {
        # n_eff_prey is the total prey abundance by size exposed to each
        # predator (prey not broken into species - here we are just working out
        # how much a predator eats - not which species are being eaten - that is
        # in the mortality calculation
        # \sum_j \theta_{ij} N_j(w_p) wt_p dw_p where
        # wt_p = w_p + q_j(w_p) is the total mass, including gonads, of the prey
        wt <- sweep(n_other$gonads, 2, params@w, "+")
        n_eff_prey <- sweep(params@interaction %*% (n * wt), 2, 
                            params@dw, "*", check.margin = FALSE) 
        # pred_kernel is predator species x predator size x prey size
        # So multiply 3rd dimension of pred_kernel by the prey biomass density
        # Then sum over 3rd dimension to get consumption rate of each predator by 
        # predator size
        # This line is a bottle neck
        phi_prey_species <- rowSums(sweep(
            params@pred_kernel[, , idx_sp, drop = FALSE],
            c(1, 3), n_eff_prey, "*", check.margin = FALSE), dims = 2)
        # Eating the background
        # This line is a bottle neck
        phi_prey_background <- params@species_params$interaction_resource *
            rowSums(sweep(
                params@pred_kernel, 3, params@dw_full * params@w_full * n_pp,
                "*", check.margin = FALSE), dims = 2)
        encounter <- params@search_vol * (phi_prey_species + phi_prey_background)
    } else {
        # resource biomass
        prey <- outer(params@species_params$interaction_resource, n_pp)
        prey <- sweep(prey, 2, params@w_full, "*")
        # add fish biomass including gonads
        wt <- sweep(n_other$gonads, 2, params@w, "+")
        prey[, idx_sp] <- prey[, idx_sp] + params@interaction %*% (n * wt)
        # multiply everything by dw_full
        prey <- sweep(prey, 2,  params@dw_full, "*")
        # The vector prey equals everything inside integral (3.4) except the feeding
        # kernel phi_i(w_p/w).
        # Eq (3.4) is then a convolution integral in terms of prey[w_p] and phi[w_p/w].
        # We approximate the integral by the trapezoidal method. Using the
        # convolution theorem we can evaluate the resulting sum via fast fourier
        # transform.
        # mvfft() does a Fourier transform of each column of its argument, but
        # we need the Fourier transforms of each row, so we need to apply mvfft()
        # to the transposed matrices and then transpose again at the end.
        avail_energy <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_e) * 
                                             mvfft(base::t(prey)),
                                         inverse = TRUE))) / length(params@w_full)
        # Only keep the bit for fish sizes
        avail_energy <- avail_energy[, idx_sp, drop = FALSE]
        # Due to numerical errors we might get negative or very small entries that
        # should be 0
        avail_energy[avail_energy < 1e-18] <- 0
        
        encounter <- params@search_vol * avail_energy
    }
    
    # Add contributions from other components
    for (i in seq_along(params@other_encounter)) {
        encounter <- encounter + 
            do.call(params@other_encounter[[i]], 
                    list(params = params,
                         n = n, n_pp = n_pp, n_other = n_other,
                         component = names(params@other_encounter)[[i]], ...))
    }
    
    # Add external encounter
    return(encounter + params@ext_encounter)
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
    release_func <- get0(other_params(params)$release_func)
    r <- release_func(t, params)
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
#' @param params A params object that must contain `other_params$R_max`
#' @param t The time at which to calculate RDD
#' @param ... Unused
#'
#' @return Vector of density-dependent reproduction rates.
#' @export
#' @family functions calculating density-dependent reproduction rate
seasonalBevertonHoltRDD <- function(rdi, params, t, ...) {
    t_name = as.character(round(t - floor(t), 5))
    if (!(t_name %in% dimnames(other_params(params)$r_max)$time)) {
        stop("r_max not defined at time ", t_name)
    }
    r_max <- other_params(params)$r_max[t_name, ]
    if (is.null(r_max)) {
        stop("r_max is NULL at time ", t_name)
    }
    return(rdi / (1 + rdi / r_max))
}

#' von-Mises distributed reproduction rate independent of abundance
#'
#' @param params A MizerParams object
#' @param t The time at which to calculate RDD
#' @param ... Unused
#'
#' @return A vector of species-specific reproduction rates
#' @export
seasonalVonMisesRDD <- function(params, t, ...) {
    sp <- params@species_params
    kappa <- sp$rdd_vonMises_kappa
    mu <- sp$rdd_vonMises_mu
    H <- sp$rdd_vonMises_r0 * exp(kappa * cos(2*pi*(t - mu))) /
        (2 * pi * besselI(kappa, nu = 0))
    return(H)
}

#' Seasonal semichemostat resource dynamics
#' 
#' @param params A [MizerParams] object
#' @param n A matrix of species abundances (species x size)
#' @param n_pp A vector of the resource abundance by size
#' @param n_other A list with the abundances of other components
#' @param rates A list of rates as returned by [mizerRates()]
#' @param t The current time
#' @param dt Time step
#' @param resource_rate Resource replenishment rate
#' @param resource_capacity Resource carrying capacity
#' @param ... Unused
#'
#' @return Vector containing resource spectrum at next timestep
#' @export
#' @family resource dynamics
seasonal_resource_semichemostat <- function(params, n, n_pp, n_other,
                                            rates, t, dt,
                                            resource_rate, resource_capacity, 
                                            ...) {
    # The resource capacity is now time-dependent
    resource_capacity <- (1 + resource_vonMises(params = params, t = t)) *
        resource_capacity
    # We use the exact solution under the assumption of constant mortality
    # during timestep
    mur <- resource_rate + rates$resource_mort
    n_steady <- resource_rate * resource_capacity / mur
    n_pp_new <- n_steady + (n_pp - n_steady) * exp(-mur * dt)

    # Here is an alternative expression that looks as if it might be more
    # precise when the sum of the rates is small due to the use of expm1.
    # However the above has the advantage of preserving the steady state
    # n_steady exactly.
    # n_pp_new <- n_pp * exp(-mur * dt) + n_steady * expm1(-mur * dt)

    # if growth rate and death rate are zero then the above would give NaN
    # whereas the value should simply not change
    sel <- !is.finite(n_pp_new)
    n_pp_new[sel] <- n_pp[sel]

    n_pp_new
}

pulsed_rate <- function(params, t) {
    # Constant velocity v in log size
    v <- params@resource_params$pulse_speed
    kappa <- params@resource_params$pulse_kappa
    mu <- params@resource_params$pulse_mu
    r0 <- params@resource_params$pulse_r0
    w_full <- params@w_full
    n <- params@resource_params$n
    r <- r0 * w_full ^ (n - 1) *
        vonMises(log(w_full) - v * t, mu = mu, kappa = kappa)
}

vonMises <- function(x, mu, kappa) {
    exp(kappa * cos(2*pi*(x - mu))) / (2 * pi * besselI(kappa, nu = 0))
}

resource_vonMises <- function(t, params, ...){
    new_t <- t - floor(t)
    kappa <- params@resource_params$rp$kappa
    mu <- params@resource_params$rp$mu
    params@resource_params$rp$maxR * exp(kappa * cos(2 * pi * (new_t - mu))) / 
        (besselI(kappa,nu=0))
}
