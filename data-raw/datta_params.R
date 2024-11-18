# We extract the parameters from the `baseModel` object created with the
# script `Seasonality paper results.R` in the code from the paper by Datta &
# Blanchard (2016), that can be downloaded from 
# <https://figshare.com/s/e75f29d4cc9b94ae393b>.
# So before running this code, you must have that `baseModel` object in your
# environment.

library(dplyr)
library(mizerSeasonal)
library(mizerExperimental)

param <- baseModel$param

# Extract species parameters ----
sp <- param$species
sp <- dplyr::rename(sp, w_inf = Winf, w_mat = Wmat, w_min = Wmin, 
             erepro = eRepro, R_max = R0)
sp$w_max <- sp$w_inf

# Extract gear parameters ----
# Some checks on the selectivity parameters
sel_params <- param$sel_params
all(sel_params$species == rep(sp$species, 4))
all(sel_params$gear == sel_params$species)
all(sel_params$param_value[sel_params$param_name == "a"] == sp$a)
all(sel_params$param_value[sel_params$param_name == "b"] == sp$b)
all(rownames(param$Q) == sp$species)
# Model uses one sigmoid_length gear for each species
gp <- data.frame(
    species = sp$species,
    gear = sp$species,
    sel_func = "sigmoid_length",
    catchability = diag(param$Q),
    l50 = param$sel_params$param_value[param$sel_params$param_name == "L50"],
    l25 = param$sel_params$param_value[param$sel_params$param_name == "L25"]
)

# MizerParams ----
p <- newMultispeciesParams(
    species_params = sp,
    gear_params = gp,
    interaction = param$theta,
    no_w = param$ngrid,
    min_w = param$w0,
    max_w = param$wMax,
    min_w_pp = baseModel$wFull[1],
    n = param$n,
    p = param$p,
    lambda = param$lambda,
    w_pp_cutoff = param$wPPcut,
    resource_rate = param$rPP,
    resource_capacity = param$kap,
    kappa = param$kap,
    z0pre = param$Z0pre,
    z0exp = param$Z0exp)

# The base model uses a constant effort of 1 for all species
all(baseModel$effort == 1)
# Therefore we set that as the initial effort in the MizerParams object
initial_effort(p) <- baseModel$effort[1, ]

# Initial abundances ----
initialN(p) <- baseModel$N[1, , ]
# The code by Datta & Blanchard (2016) uses
# wider size-bins for the resource spectrum than for the fish spectrum. The
# modern mizer code no longer supports this. So this MizerParams object uses
# 180 size bins for the full spectra instead of 130.
length(w_full(p))
length(baseModel$wFull)
# This means that we can not simply copy over the initial resource abundances.
# We have to interpolate them at our new bin boundaries.
initialNResource(p) <- approx(baseModel$wFull, baseModel$nPP[1, ], 
                              w_full(p), rule = 2)$y

# Compare ----
# First we check that we understand the different size grids
waldo::compare(baseModel$w, w(p))
idx_fish_old <- (length(baseModel$wFull) - length(baseModel$w) + 1) : length(baseModel$wFull)
idx_fish_new <- (length(w_full(p)) - length(w(p)) + 1) : length(w_full(p))
waldo::compare(baseModel$wFull[idx_fish_old], p@w_full[idx_fish_new])
# We now check that the rates in the MizerParams object agree with those in the
# baseModel object
waldo::compare(baseModel$psi, p@psi, tolerance = 1e-6, ignore_attr = TRUE)
waldo::compare(baseModel$IntakeMax, intake_max(p), ignore_attr = TRUE)
waldo::compare(baseModel$Z0, species_params(p)$z0)
waldo::compare(baseModel$SearchVol, search_vol(p), tolerance = 1e-14, ignore_attr = TRUE)
waldo::compare(baseModel$StdMetab, metab(p), tolerance = 1e-14, ignore_attr = TRUE)           
waldo::compare(baseModel$selectivity, aperm(p@selectivity, c(2,3,1)),
               tolerance = 1e-14, ignore_attr = TRUE)
all(baseModel$Aktivity == 0)
all.equal(getFMort(p), baseModel$F[1, , ], 
          check.attributes = FALSE)
# Because of the different size grids for the resource, for the pred kernel 
# we can make the comparison only where the prey are fish
all.equal(baseModel$predkernel[1, , idx_fish_old], 
          pred_kernel(p)[1, , idx_fish_new], 
          check.attributes = FALSE)

# The feeding level is off by a little bit, presumably because of the slightly
# different resource abundance due to the different size grids
all.equal(getFeedingLevel(p), baseModel$f[1, , ], 
          check.attributes = FALSE, tolerance = 0.006)
# We can look at this graphically for individual species
i <- 1
plot(w(p), baseModel$f[1, i, ], type = "l", log = "xy")
lines(w(p), getFeedingLevel(p)[i, ], col = "red")
# We do not need to be too bothered by such differences by less than 1%

all.equal(getM2(p), baseModel$M2[1, , ], 
          check.attributes = FALSE, tolerance = 1e-4)

all.equal(getResourceMort(p)[idx_fish_new], 
          baseModel$M2background[1, idx_fish_old], 
          check.attributes = FALSE, tolerance = 1e-4)
all.equal(getERepro(p), baseModel$eSpawning[1, , ], 
          check.attributes = FALSE, tolerance = 0.007)

all.equal(getRDI(p), baseModel$RDI[1, ], 
          check.attributes = FALSE, tolerance = 0.014)
all.equal(getRDD(p), baseModel$RDD[1, ], 
          check.attributes = FALSE, tolerance = 0.014)

# check resource graphically
plot(baseModel$wFull, baseModel$rrPP, type = "l", log = "xy")
lines(w_full(p), resource_rate(p), col = "red")

plot(baseModel$wFull, baseModel$NinfPP, type = "l", log = "xy")
lines(w_full(p), resource_capacity(p), col = "red")

# Simulation ----
sim <- project(p, t_max = 500)
plotBiomass(sim)
# We see that the initial state of the Datta and Blanchard model is far from
# steady state. In the first few years the biomass of the fish species changes
# by up to 10^6%! It is therefore not surprising that if we compare the dynamics
# between the two models we find differences.
sim2 <- sim
sim2@n[2:501, , ] <- baseModel$N[(1:500) * 52, , ]
plotlyBiomassRelative(sim, sim2)

# Both models however reach a steady state quite quickly and the
# steady states are quite similar.
# Create params object with mizer steady state
ps <- setInitialValues(p, sim)
# Crate params object with Datta & Blanchard steady state
final_time_idx <- dim(baseModel$N)[1]
ps_datta <- p
initialN(ps_datta) <- baseModel$N[final_time_idx, , ]
initialNResource(ps_datta) <- approx(baseModel$wFull, baseModel$nPP[final_time_idx, ], 
                                     w_full(p), rule = 2)$y
# Compare the steady states
plotSpectra2(ps, ps_datta)
plotSpectraRelative(ps, ps_datta)

# We'll now make this MizerParams object available in the package.
datta_params <- 
    setMetadata(ps,
                title = "Base model from Datta & Blanchard (2016)",
                description = "This is a re-implementation of the base model from Datta & Blanchard (2016) using the mizer package. The initial state is set to the steady state.")
usethis::use_data(datta_params, overwrite = TRUE)
