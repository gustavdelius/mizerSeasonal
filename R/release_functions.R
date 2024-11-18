#' von-Mises distributed gonad release rate
#'
#' Uses the formula
#' \deqn{r(t) = r_0 \frac{\exp(\kappa \cos(2\pi(t - \mu)))}{2\pi I_0(\kappa)}}{r(t) = r_0 exp(\kappa cos(2\pi(t - \mu))) / (2\pi I_0(\kappa))}
#' where the parameters are taken from the `vonMises_r0`, `vonMises_kappa` and 
#' `vonMises_mu` columns in the `species_params` data frame in `params`.
#'
#' @param t The time at which to calculate the release rate
#' @param params A MizerParams object
#' @param ... Unused
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

#' Beta hazard mass-specific gonad release rate
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

#' Beta distributed gonad release rate
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

#' Gaussian mass-specific gonad release rate
#'
#' Uses the formula
#' \deqn{r(w, t) = r(t) = r_0 \exp{\left(-\dfrac{(t-\lfloor{t}\rfloor-t_0)^2}{2\sigma^2}\right)}}{r(w,t) = r(t)= r_0 \exp(-(t-|_t_|-t_0)^2/(2\sigma^2))}
#' where parameters \eqn{r_0, \sigma} and \eqn{t_0} are given by new columns
#' `sr_r0`, `sr_sigma` and `sr_t0` in the species parameter data frame.
#' Because mizer measures time in years, \eqn{t-\lfloor t \rfloor}{t-|_t_|} gives the time
#' within the year and so \eqn{r(t)} is a periodic function with
#' the period of one year.
#'
#' @param t The time at which to calculate the release rate
#' @param params A MizerParams object
#'
#' @return A vector of species-specific release rates at time t
#' @export
seasonalGaussianRelease <- function(t, params) {
    sp <- params@species_params
    sp$sr_r0 * exp(-(t - trunc(t) - sp$sr_t0)^2/(2*sp$sr_sigma^2))
}