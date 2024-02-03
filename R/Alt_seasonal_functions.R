#' von-Mises distributed release rate
#'
#' Uses the formula
#' \deqn{r(t) = r_0 \frac{\exp(\kappa \cos(2\pi(t - \mu)))}{2\pi I_0(\kappa)}}{r(t) = r_0 exp(\kappa cos(2\pi(t - \mu))) / (2\pi I_0(\kappa))}
#' where the parameters are taken from the `vonMises_r0`, `vonMises_kappa` and 
#' `vonMises_mu` columns in the `species_params` data frame in `params`.
#'
#' @param t The time at which to calculate the release rate
#' @param params A MizerParams object
#'
#' @return A vector of species-specific release rates at time t
#' @export
seasonalVonMisesRelease <- function(t, params, ...) {
    sp <- params@species_params
    kappa <- sp$vonMises_kappa
    mu <- sp$vonMises_mu
    H <- sp$vonMises_r0 * exp(kappa * cos(2 * pi * (t - mu))) /
        (2 * pi * besselI(kappa, nu = 0))
    return(H)
}

#' Beta hazard mass-specific release rate
#'
#' @param t The time at which to calculate the reproduction rate
#' @param params A MizerParams object
#'
#' @return A vector of species-specific release rates at time t
#' @export
seasonalBetaHazardRelease <- function(t,params){
  new_t <- t - floor(t)
  H <- stats::dbeta(new_t, params@species_params$beta_a,
                    params@species_params$beta_b) / 
      (1 - stats::pbeta(new_t, params@species_params$beta_a,
                      params@species_params$beta_b))
  return(H)
}

#' Beta distributed release rate
#'
#' @param t The time at which to calculate the reproduction rate
#' @param params A MizerParams object
#'
#' @return A vector of species-specific release rates at time t
#' @export
seasonalBetaRelease <- function(t,params){
  new_t <- t - floor(t)
  H <- params@species_params$beta_r * stats::dbeta(new_t,params@species_params$beta_a,params@species_params$beta_b)
  return(H)
}
