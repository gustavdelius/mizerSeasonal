# Set up Samik's model
install_github("sizespectrum/mizer", ref = "t-dependent_RDD")

library(mizerSeasonal)
p <- NS_params

sp <- p@species_params
# Von Mises parameters for RDD from Samik's paper
order <- c(1, 2, 3, 5, 4, 8, 7, 6, 9, 10, 11, 12)
sp$rdd_vonMises_k <- c(3.6047, 2.9994, 1.944, 1.141, 0.40493, 0.795, 4.973,1.4257, 4.0495, 3.6567, 5.4732, 1.951)[order]
sp$rdd_vonMises_mu <- c(0.5477, 0.0643, 0.3574, 0.3825, 0.8576, 0.5021,0.4856, 0.3716, 0.2245, 0.4181, 0.3123, 0.3333)[order]
sp$rdd_vonMises_r <- getRDD(p)

# Make some guesses for r(t)
# Parameters for Gaussian r(t)
sp$sr_t0 <- sp$rdd_vonMises_mu
sp$sr_sigma <- 1 / sqrt(sp$rdd_vonMises_k) / 10
sp$sr_r0 <- 50
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

sim <- projectToSteady(ps, distance_func = distanceMaxRelRDI,
                      t_per = 1, t_max = 50, dt = 0.01, tol = 0.01,
                      return_sim = TRUE)
plotSpectra(sim)
plotBiomass(sim)

ps <- setInitialValues(ps, sim)
sim <- project(ps, t_max = 20, dt = 0.01, t_save = 0.1)
plotBiomass(sim)
plotRDI(sim)

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
animateSpectra(sim3, power = 0)

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
plot.ts(log(sim2r@n_pp[, 10*(1:10)]), log="y")


sim2rr <- project(pr, t_max = 2, dt = 0.01, t_save = 0.1)
animateGonadSpectra(sim2rr)
animateSpectra(sim2rr, power = 2)
# Turn on reproduction
