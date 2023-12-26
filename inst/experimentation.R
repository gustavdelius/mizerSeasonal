p <- newTraitParams()
p <- setBevertonHolt(p, reproduction_level = 0)
p <- steady(p)
p@species_params$sr_r0 <- 10
p@species_params$sr_exp <- 1e3
p@species_params$sr_t0 <- 0.2
p <- setSeasonalReproduction(p)

sim <- project(p, t_max = 10, t_save = 0.05, dt = 0.01)
p <- setInitialValues(p, sim)

plotRDI(sim, log = TRUE)
plotRDI(sim, log = FALSE)
plotRDD(sim, log = TRUE)
plotSpectra(sim)
animateSpectra(sim, power = 2, time_range = c(9,10))
animateGonadSpectra(sim)

t <- seq(0, 1, by = 0.01)
r <- sapply(t, seasonalRepro, params = p)
plot(t, log(r[1, ]), type = "l")
plot(t, r[1, ], type = "l")


plot(getTimes(sim), unlist(sim@n_other[, "gonads"]), type = "l")
