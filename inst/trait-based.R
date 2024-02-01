# Set up Samik's model
install_github("sizespectrum/mizer", ref = "t-dependent_RDD")

library(mizerSeasonal)
p <- newTraitParams()

sp <- p@species_params
sp$rdd_vonMises_k <- 5
sp$rdd_vonMises_mu <- 0.3
sp$rdd_vonMises_r <- getRDD(p)

# Make some guesses for r(t)
# Parameters for Von Mises r(t)
sp$vonMises_k <- sp$rdd_vonMises_k
sp$vonMises_mu <- sp$rdd_vonMises_mu
sp$vonMises_r <- 100

species_params(p) <- sp

# Set seasonal investment into reproduction
# # Using gaussian
# p <- setSeasonalReproduction(p)
# Using Von Mises
p <- setSeasonalReproduction(p, repro_func = "repro_vonMises")

# Set RDD to observations
p <- setRateFunction(p, "RDD", "seasonalVonMisesRDD")

# Force resource to stay at current level
old_resource_dynamics <- p@resource_dynamics
p@resource_dynamics <- "resource_constant"

sim <- projectToSteady(p, distance_func = distanceMaxRelRDI,
                      t_per = 1, t_max = 50, dt = 0.01, tol = 0.01,
                      return_sim = TRUE)
plotSpectra(sim)
plotBiomass(sim)

ps <- setInitialValues(p, sim)
sim <- project(ps, t_max = 20, dt = 0.01, t_save = 0.1)
plotBiomass(sim)


# Plot it for just 2 years at higher time resolution
ps <- setInitialValues(ps, sim)
sim2 <- project(ps, t_max = 2, dt = 0.01, t_save = 0.01)
species = 12
plotRDI(sim2, species = species)
plotRDD(sim2, species = species)
plotGonads(sim2,sizes=sim2@params@species_params$w_max, prop_size = 0.8)

# For animation store it with lower time resolution
sim3 <- project(ps, t_max = 2, dt = 0.01, t_save = 0.1)
animateGonadSpectra(sim3)
animateSpectra(sim3, power = 2)

Times <- 20
feed_bit <-matrix(0,Times * 10 +1 ,12)
for (i in 1:12){
    w_max_idx <- sum(ps@w < ps@species_params$w_mat[i])
  feed_bit[,i] <- sapply(1:(Times * 10 + 1),function(t){getEncounter(ps,n=sim@n[t,,],n_pp=sim@n_pp[t,],n_other=sim@n_other[t,])[i,w_max_idx]})
}

graphics::par(mfrow=c(4,3))
graphics::par(oma=c(1,1,1,3))
graphics::par(mar=c(2,4,2,0))
for(i in 1:12){
    plot(seq(0,Times +1,length.out=Times * 10 + 1),feed_bit[,i],type="l",main=ps@species_params$species[i])
}


# Turn on resource dynamics

# We choose C = 2*N and r = average \mu
getResourceMort_annual <- function(sim,
                                   time_range,
                                   ...){
    params <- sim@params
    no_sp <- nrow(params@species_params)
    times <- getTimes(sim)
    time_indices <- seq_along(times)
    value <- array(0, dim = c(length(times), length(sim@params@w_full)),
                   dimnames = dimnames(sim@n_pp)[1:2])
    for (i in time_indices) {
        n <- sim@n[i, , ]
        dim(n) <- dim(params@initial_n)
        n_other <- sim@n_other[i, ]
        names(n_other) <- dimnames(sim@n_other)[[2]]
        value[i, ] <- getResourceMort(params,
                                      n = n,
                                      n_pp = sim@n_pp[i, ],
                                      n_other = n_other,
                                      t = times[i])
    }
    value
}
mu <- colMeans(getResourceMort_annual(sim2))
pr <- ps
pr@rr_pp <- mu / 200000000000000000
pr@cc_pp <- 200000000000000001 * pr@initial_n_pp
pr@resource_dynamics <- "resource_semichemostat"

simr <- project(pr, t_max = 20, dt = 0.01, t_save = 0.1)
plotBiomass(simr)

pr <- setInitialValues(pr, simr)
sim2r <- project(pr, t_max = 2, dt = 0.01, t_save = 0.01)
plotBiomass(sim2r)
plotGonads(sim2r)
plotRDD(sim2r)
plot.ts(sim2r@n_pp[, 10*(1:10)])

sim2rr <- project(pr, t_max = 3, dt = 0.01, t_save = 0.1)
animateGonadSpectra(sim2rr)
animateSpectra(sim2rr, power = 2)

# Pulsed resource ----
calc_rp <- function(params,spe,init,kappa,maxR){
    ret <- list()
    Nws <- length(params@w_full) - length(params@w)
    LW <- log(params@w_full)
    speed <- spe / diff(range(LW))
    const <- init - speed * LW[1]
    tmp <-  speed * LW + const
    ret$mu <- tmp - floor(tmp)
    ret$kappa <- kappa
    ret$ maxR<- maxR
    params@resource_params$rp <- ret
    return(params)
}

################
resource_vonMises <- function(t,params,...){
    new_t <- t - floor(t)
    kappa <- params@resource_params$rp$kappa
    mu <- params@resource_params$rp$mu
    params@resource_params$rp$maxR * exp(kappa * cos(2 * pi * (new_t - mu))) / (besselI(kappa,nu=0))
}

pp <- calc_rp(ps, spe=2,init=0.1,kappa=5,maxR=0)
pp@resource_dynamics <- "seasonal_resource_semichemostat"
pp@cc_pp <- 2 * pp@initial_n_pp

simpp <- project(pp, t_max = 20, dt = 0.01, t_save = 0.1)
plotBiomass(simpp)
plotSpectra(simpp)

# For animation store it with lower time resolution
sim3 <- project(pp, t_max = 2, dt = 0.01, t_save = 0.1)
animateGonadSpectra(sim3)
animateSpectra(sim3, power = 2)

# Turn on reproduction
prr <- pr
repro_level <- 0.5
rdd <- getRDD(prr)
rdi <- getRDI(prr)
prr@other_params$r_max <- rdd * rdi / (rdi - rdd)
