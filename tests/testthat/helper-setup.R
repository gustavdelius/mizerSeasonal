suppressPackageStartupMessages(library(mizer))

# Single-species params with all seasonal columns needed for tests
base_params <- suppressMessages(newSingleSpeciesParams())

base_params@species_params$vonMises_r0 <- 10
base_params@species_params$vonMises_kappa <- 2
base_params@species_params$vonMises_mu <- 0.25

base_params@species_params$beta_a <- 2
base_params@species_params$beta_b <- 5
base_params@species_params$beta_r <- 1

base_params@species_params$sr_r0 <- 10
base_params@species_params$sr_sigma <- 0.1
base_params@species_params$sr_t0 <- 0.25

base_params@species_params$rdd_vonMises_r0 <- 100
base_params@species_params$rdd_vonMises_kappa <- 2
base_params@species_params$rdd_vonMises_mu <- 0.25

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
