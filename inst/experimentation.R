p <- newTraitParams()
p <- setBevertonHolt(p, reproduction_level = 0)
p <- steady(p)
p@species_params$sr_r0 <- 5
p@species_params$sr_sigma <- sqrt(1/2000)
p@species_params$sr_t0 <- 0.2
p <- setSeasonalReproduction(p)

sim <- project(p, t_max = 10, t_save = 0.05, dt = 0.01)
p <- setInitialValues(p, sim)

plotRDI(sim, log = TRUE)
plotRDI(sim, log = FALSE)
plotSpectra(sim)
animateSpectra(sim, power = 2, time_range = c(9,10))
animateGonadSpectra(sim, time_range = c(9,10))
