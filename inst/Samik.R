# Set up Samik's model
library(mizerSeasonal)
p <- NS_params

sp <- p@species_params
# Von Mises parameters for RDD from Samik's paper
sp$vonMises_k <- c(3.6047, 2.9994, 1.944, 1.141, 0.40493, 0.795, 4.973,1.4257, 4.0495, 3.6567, 5.4732, 1.951)
sp$vonMises_mu <- c(0.5477, 0.0643, 0.3574, 0.3825, 0.8576, 0.5021,0.4856, 0.3716, 0.2245, 0.4181, 0.3123, 0.3333)
sp$vonMises_r <- getRDD(p) 

# Make some guesses for r(t)
sp$sr_t0 <- sp$vonMises_mu
sp$sr_sigma <- 1 / sqrt(sp$vonMises_k) / 10
sp$sr_r0 <- 50

species_params(p) <- sp

# Set seasonal investment into reproduction
p <- setSeasonalReproduction(p)
# Set RDD to observations
p <- setRateFunction(p, "RDD", "seasonalVonMisesRDD")

# Force resource to stay at current level
old_resource_dynamics <- p@resource_dynamics
p@resource_dynamics <- "resource_constant"

ps <- projectToSteady(p,
                      distance_func = distanceMaxRelRDI,
                      t_per = 1,
                      t_max = 50,
                      dt = 0.01)
plotSpectra(ps)

sim <- project(ps, t_max = 20, dt = 0.01, t_save = 0.1)
plotBiomass(sim)
plotRDI(sim)

ps <- setInitialValues(ps, sim)
sim2 <- project(ps, t_max = 2, dt = 0.01, t_save = 0.01)
plotRDI(sim2, species = 5)

sim3 <- project(ps, t_max = 2, dt = 0.01, t_save = 0.1)
animateGonadSpectra(sim3)
