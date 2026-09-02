# We extract the parameters from the trimmed `baseModel` object in
# `data-raw/baseModel.RData`. The full object was created with the script
# `Seasonality paper results.R` in the code from the paper by Datta & Blanchard
# (2016), available from <https://figshare.com/s/e75f29d4cc9b94ae393b>.

library(dplyr)
library(mizerSeasonal)
library(mizerExperimental)

load("data-raw/baseModel.RData")
load("data-raw/legacy_targets.RData")
param <- baseModel$param

# Extract species parameters ----
sp <- param$species
sp <- dplyr::rename(sp, w_inf = Winf, w_mat = Wmat, w_min = Wmin,
             erepro = eRepro, R_max = R0, l_inf = Linf, l_mat = Lmat)
sp$w_max <- sp$w_inf
# The length paramters do not match their weight parameters
all.equal(l2w(sp$l_inf, sp), sp$w_inf)
all.equal(l2w(sp$l_mat, sp), sp$w_mat)
# So we will remove them
sp$l_inf <- NULL
sp$l_mat <- NULL

# Fix dimnames on interaction matrix ----
interaction <- param$theta
rownames(interaction) <- sp$species
colnames(interaction) <- sp$species

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
create_params <- function(second_order_w = FALSE) {
    newMultispeciesParams(
        species_params = sp,
        gear_params = gp,
        interaction = interaction,
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
        z0exp = param$Z0exp,
        second_order_w = second_order_w)
}

p <- create_params()
p_second_order <- create_params(second_order_w = TRUE)

# The base model uses a constant effort of 1 for all species
all(baseModel$effort == 1)
# Therefore we set that as the initial effort in the MizerParams object
initial_effort(p) <- baseModel$effort[1, ]
initial_effort(p_second_order) <- baseModel$effort[1, ]

# Initial abundances ----
# The trimmed baseModel retains only the initial abundance snapshots needed to
# reproduce the original MizerParams objects, not the full simulated arrays.
initialN(p) <- baseModel$initialN
initialN(p_second_order) <- baseModel$initialN

# The code by Datta & Blanchard (2016) uses wider size-bins for the resource
# spectrum than for the fish spectrum. Modern mizer no longer supports this.
# The new MizerParams objects therefore use 180 size bins for the full spectra
# instead of 130. newMultispeciesParams() has already initialised the resource
# power law on the new grid, using point values for the first-order object and
# bin averages for the second-order object. Do not interpolate the abundance
# density from the old grid.
length(w_full(p))
length(baseModel$wFull)

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
all(baseModel$Activity == 0)
# Because of the different size grids for the resource, for the pred kernel
# we can make the comparison only where the prey are fish
all.equal(baseModel$predkernel[1, , idx_fish_old],
          pred_kernel(p)[1, , idx_fish_new],
          check.attributes = FALSE)

# Compare rates and summaries without reconstructing the legacy resource on the
# modern grid. The legacy targets contain the values calculated on the original
# grid at the initial and final states.
compare_values <- function(current, target, ...) {
    all.equal(unclass(current), unclass(target),
              check.attributes = FALSE, ...)
}

compare_legacy_snapshot <- function(params, target) {
    list(
        fish_density = compare_values(initialN(params), target$N),
        resource_density = compare_values(
            initialNResource(params)[idx_fish_new],
            target$nPP[idx_fish_old]),
        fishing_mortality = compare_values(
            getFMort(params), target$fishing_mortality),
        feeding_level = compare_values(
            getFeedingLevel(params), target$f, tolerance = 0.006),
        encounter = compare_values(
            getEncounter(params), target$encounter, tolerance = 0.006),
        somatic_growth = compare_values(
            getEGrowth(params), target$somatic_growth, tolerance = 0.007),
        predation_mortality = compare_values(
            getPredMort(params), target$predation_mortality,
            tolerance = 1e-4),
        resource_mortality = compare_values(
            getResourceMort(params)[idx_fish_new],
            target$resource_mortality[idx_fish_old],
            tolerance = 1e-4),
        reproductive_investment = compare_values(
            getERepro(params), target$reproductive_investment,
            tolerance = 0.007),
        recruitment_density_independent = compare_values(
            getRDI(params), target$recruitment_density_independent,
            tolerance = 0.014),
        recruitment_density_dependent = compare_values(
            getRDD(params), target$recruitment_density_dependent,
            tolerance = 0.014),
        biomass = compare_values(getBiomass(params), target$biomass),
        yield = compare_values(getYield(params), target$yield)
    )
}

legacy_initial_checks <- compare_legacy_snapshot(p, legacy_targets$initial)
legacy_initial_checks
# We see that there are discrepancies, resulting from the different
# resource bin sizes.
i <- 1
plot(w(p), legacy_targets$initial$f[i, ], type = "l", log = "xy")
lines(w(p), getFeedingLevel(p)[i, ], col = "red")

# check resource graphically
plot(baseModel$wFull, baseModel$rrPP, type = "l", log = "xy")
lines(w_full(p), resource_rate(p), col = "red")

plot(baseModel$wFull, baseModel$NinfPP, type = "l", log = "xy")
lines(w_full(p), resource_capacity(p), col = "red")

# Simulation ----
sim <- projectUntilSettled(p, t_max = 500)
plotHover(getSteadyResidual(finalParams(sim)))
sim_second_order <- projectUntilSettled(p_second_order, t_max = 200, dt = 0.001,
                                        method = "tr_bdf2")
plotHover(getSteadyResidual(finalParams(sim_second_order)))

plotBiomass(sim)
# We see that the initial state of the Datta and Blanchard model is far from
# steady state. In the first few years the biomass of the fish species changes
# by up to 10^6%! It is therefore not surprising that if we compare the dynamics
# between the two models we find differences.
# The compact legacy targets deliberately retain only endpoint snapshots, so
# the historical year-by-year biomass comparison is not reproduced here.

# Both models however reach a steady state quite quickly and the
# steady states are quite similar.
# Create params object with mizer steady state
ps <- finalParams(sim)
ps_second_order <- finalParams(sim_second_order)
ps_datta <- p
initialN(ps_datta) <- legacy_targets$final$N

# Compare the fish spectra on their shared grid. The resource spectra are
# plotted separately below because their grids differ.
plotSpectra2(ps, ps_datta, "Modern", "Legacy", resource = FALSE)
plotSpectraRelative(ps, ps_datta, resource = FALSE)

plot(legacy_targets$grid$w_full, legacy_targets$final$nPP,
     type = "l", log = "xy", xlab = "Resource size", ylab = "Density")
lines(w_full(ps), initialNResource(ps), col = "red")

legacy_final_checks <- compare_legacy_snapshot(ps, legacy_targets$final)
legacy_final_checks

# We'll now make this MizerParams object available in the package.
datta_params <-
    setMetadata(ps,
                title = "Base model from Datta & Blanchard (2016)",
                description = "This is a re-implementation of the base model from Datta & Blanchard (2016) using the mizer package. The initial state is set to the steady state.")
usethis::use_data(datta_params, overwrite = TRUE)

datta_params_second_order <-
    setMetadata(
        ps_second_order,
        title = "Base model from Datta & Blanchard (2016), second-order scheme",
        description = paste(
            "This is a re-implementation of the base model from Datta &",
            "Blanchard (2016) using the second-order size scheme in the mizer",
            "package. The initial state is set to the steady state."
        ))
usethis::use_data(datta_params_second_order, overwrite = TRUE)
