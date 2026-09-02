suppressPackageStartupMessages(library(mizer))

# Single-species params with all seasonal columns needed for tests
base_params <- suppressMessages(newSingleSpeciesParams())

sp <- species_params(base_params)

sp$vonMises_r0 <- 10
sp$vonMises_kappa <- 2
sp$vonMises_mu <- 0.25

sp$beta_a <- 2
sp$beta_b <- 5
sp$beta_r <- 1

sp$sr_r0 <- 10
sp$sr_sigma <- 0.1
sp$sr_t0 <- 0.25

sp$rdd_vonMises_r0 <- 100
sp$rdd_vonMises_kappa <- 2
sp$rdd_vonMises_mu <- 0.25

species_params(base_params) <- sp

# Seasonal params using seasonalVonMisesRDD (avoids r_max dependency)
seasonal_params <- suppressMessages(
    setSeasonalReproduction(
        base_params,
        release_func = "seasonalVonMisesRelease",
        RDD = "seasonalVonMisesRDD"
    )
)

# Short simulation for plot tests (2 years, saves every 0.5 years)
seasonal_sim <- suppressMessages(
    project(seasonal_params, t_max = 2, dt = 0.1, t_save = 0.5)
)
