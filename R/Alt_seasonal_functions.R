#' Beta hazard mass-specific reproduction rate
#' 
#' Uses the formula
#' 
#' @param t The time at which to calculate the reproduction rate
#' @param params A MizerParams object
#' 
#' @return A vector of species-specific reproduction rates
#' @export
repro_beta_hazard <- function(t,params){
  new_t <- t - floor(t)
  H <- stats::dbeta(new_t,params@species_params$beta_a,params@species_params$beta_b) / (1-stats::pbeta(new_t,params@species_params$beta_a,params@species_params$beta_b))
  return(H)
}

#' Beta distributed reproduction rate
#' 
#' Uses the formula
#' 
#' @param t The time at which to calculate the reproduction rate
#' @param params A MizerParams object
#' 
#' @return A vector of species-specific reproduction rates
#' @export
repro_beta <- function(t,params){
  new_t <- t - floor(t)
  H <- params@species_params$beta_r * stats::dbeta(new_t,params@species_params$beta_a,params@species_params$beta_b)
  return(H)
}

#' von-Mises distributed reproduction rate
#' 
#' Uses the formula
#' 
#' @param t The time at which to calculate the reproduction rate
#' @param params A MizerParams object
#' 
#' @return A vector of species-specific reproduction rates
#' @export
repro_vonMises <- function(t,params,...){
  sp <- params@species_params
  new_t <- t - floor(t)
  kappa <- sp$vonMises_k
  mu <- sp$vonMises_mu
  H <- sp$vonMises_r * exp(kappa * cos(new_t - mu)) / (2 * pi * besselI(kappa,nu=0))
  return(H)
}



#' not a real function. Used to follow gonadic mass through time
plotGonads <- function(sim,
                       time_range,
                       prop_size=1,
                       ...){
  def.par <- graphics::par(no.readonly=TRUE) #old pars
  
  if (missing(time_range)) {
    time_range  <- as.numeric(dimnames(sim@n)$time)
  }
  max_size <- sapply(sim@params@species_params$w_max * prop_size,function(x,w){which.max(x <= w)},w=sim@params@w)
  time_elements <- get_time_elements(sim, time_range)
  
  ###
  
  q_max <- matrix(0,length(time_elements),dim(sim@n)[2])
  qf <- sim@n_other[time_elements, "gonads"]
  
  qs_mat<-matrix(unlist(lapply(qf,function(x){return(diag(x[,max_size]))})),length(time_elements),dim(sim@n)[2],byrow=T)
  
  rownames(qs_mat) <- (as.numeric(dimnames(sim@n)$time))
  
  graphics::par(mfrow=c(4,3))
  graphics::par(oma=c(1,1,1,3))
  graphics::par(mar=c(2,4,2,0))
  for(i in 1:dim(sim@n)[2]){
    plot((as.numeric(dimnames(sim@n)$time)),qs_mat[,i],main=i,...,type="l")
  }
  graphics::par(def.par)
}
