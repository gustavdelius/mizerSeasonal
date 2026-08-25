# Re-run the non-seasonal Datta & Blanchard (2016) model and retain compact
# initial- and final-state calibration targets. The original source bundle is
# available from https://figshare.com/s/e75f29d4cc9b94ae393b.
# Usage:
# Rscript data-raw/run_legacy_model.R [legacy-code-directory] [output-file]

args <- commandArgs(trailingOnly = TRUE)
legacy_dir <- if (length(args) >= 1) {
    args[[1]]
} else {
    "data-raw/legacy"
}
output_file <- if (length(args) >= 2) {
    args[[2]]
} else {
    "data-raw/legacy_targets.RData"
}

legacy_dir <- normalizePath(legacy_dir, mustWork = TRUE)
source_files <- file.path(
    legacy_dir,
    c(
        "paramNorthSeaModel.R",
        "SizeBasedModel.r",
        "SelectivityFuncs.r"
    )
)
input_files <- file.path(
    legacy_dir,
    c(
        "nsea_params.csv",
        "interactionmatrix_Schoener_twostage2D.RData"
    )
)
stopifnot(all(file.exists(source_files)), all(file.exists(input_files)))

source(source_files[[1]])
source(source_files[[2]])
source(source_files[[3]])
suppressPackageStartupMessages(library(plyr))
load(input_files[[2]])

make_params <- function(tmax, isave) {
    params <- paramNSModel(
        fileSpecies = input_files[[1]],
        theta = theta,
        kap = 1e11
    )
    params$tmax <- tmax
    params$dt <- 1 / 52
    params$isave <- isave
    params
}

# Globals read by the legacy Project() function. These values reproduce the
# paper's non-seasonal base case.
ss_frac <- 0.7
p_peak <- 0
p_peakt <- 0.4
r_peak <- numeric(12)
r_peakt <- c(
    0.54765, 0.064324, 0.3574, 0.85759, 0.38245, 0.37164,
    0.48564, 0.50213, 0.22454, 0.41806, 0.31236, 0.33333
)

extract_state <- function(model) {
    list(
        N = model$N[dim(model$N)[[1]], , ],
        nPP = model$nPP[dim(model$nPP)[[1]], ]
    )
}

state_model <- function(template, state) {
    model <- template
    model$N[dim(model$N)[[1]], , ] <- state$N
    model$nPP[dim(model$nPP)[[1]], ] <- state$nPP
    model
}

evaluate_state <- function(template, state, time) {
    # Project() calculates rates before advancing the state. The first rates
    # from a probe therefore apply to exactly the supplied N and nPP.
    probe_params <- template$param
    # The legacy code uses `for (i in 2:itimemax)` when constructing a noise
    # vector, so it cannot run with itimemax = 1 on current R. Run two steps
    # and retain only the first step's rates.
    probe_params$tmax <- 2 * probe_params$dt
    probe_params$isave <- 1
    probe <- Setup(
        probe_params,
        ContinueCalculation = TRUE,
        initialcommunity = state_model(template, state)
    )
    probe <- Project(probe)

    f <- probe$f[1, , ]
    fishing_mortality <- probe$F[1, , ]
    weight_by_bin <- matrix(
        probe$w * probe$dw,
        nrow = probe$param$nspp,
        ncol = probe$param$ngrid,
        byrow = TRUE
    )
    encounter <- f * probe$IntakeMax / (1 - f)

    list(
        time = time,
        N = state$N,
        nPP = state$nPP,
        f = f,
        encounter = encounter,
        somatic_growth = probe$gg[1, , ],
        fishing_mortality = fishing_mortality,
        predation_mortality = probe$M2[1, , ],
        resource_mortality = probe$M2background[1, ],
        reproductive_investment = probe$eSpawning[1, , ],
        recruitment_density_independent = probe$RDI[1, ],
        recruitment_density_dependent = probe$RDD[1, ],
        biomass = rowSums(state$N * weight_by_bin),
        yield = rowSums(fishing_mortality * state$N * weight_by_bin)
    )
}

set.seed(1)

# The compact baseModel object retains the original first weekly state. Verify
# that the source bundle still reproduces it before running the long projection.
compact_env <- new.env(parent = emptyenv())
load("data-raw/baseModel.RData", envir = compact_env)
compact_base_model <- compact_env$baseModel

weekly_params <- make_params(tmax = 2 / 52, isave = 1)
weekly_template <- Setup(weekly_params)
first_week <- Project(weekly_template)
initial_state <- list(
    N = first_week$N[1, , ],
    nPP = first_week$nPP[1, ]
)
stopifnot(
    isTRUE(all.equal(initial_state$N, compact_base_model$initialN,
                     tolerance = 0)),
    isTRUE(all.equal(initial_state$nPP,
                     compact_base_model$initialNResource,
                     tolerance = 0))
)

# Save only the final snapshot. isave only controls storage; the integration
# still uses the original weekly step for all 26,000 updates.
legacy_params <- make_params(tmax = 500, isave = 26000)
legacy_template <- Setup(legacy_params)
legacy_run <- Project(legacy_template)
final_state <- extract_state(legacy_run)

legacy_targets <- list(
    initial = evaluate_state(legacy_template, initial_state, time = 1 / 52),
    final = evaluate_state(legacy_template, final_state, time = 500),
    grid = list(
        species = legacy_params$species$species,
        w = legacy_template$w,
        dw = legacy_template$dw,
        w_full = legacy_template$wFull,
        dw_full = legacy_template$dwFull
    ),
    run = list(
        model = "Non-seasonal base case from Datta & Blanchard (2016)",
        t_max = 500,
        dt = 1 / 52,
        integration_steps = 26000L,
        storage_interval = 500,
        rate_evaluation = paste(
            "Rates were evaluated with one non-seasonal legacy step started",
            "from the N and nPP stored in the same snapshot; the probe's",
            "post-step densities were discarded."
        ),
        source_files = setNames(
            unname(tools::md5sum(c(source_files, input_files))),
            basename(c(source_files, input_files))
        ),
        r_version = R.version.string,
        plyr_version = as.character(packageVersion("plyr")),
        created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
)

required_fields <- c(
    "N", "nPP", "f", "encounter", "somatic_growth", "biomass",
    "recruitment_density_independent",
    "recruitment_density_dependent", "yield"
)
stopifnot(
    all(required_fields %in% names(legacy_targets$initial)),
    all(required_fields %in% names(legacy_targets$final)),
    all(is.finite(unlist(legacy_targets[c("initial", "final")],
                         recursive = TRUE, use.names = FALSE)))
)

save(legacy_targets, file = output_file, compress = "xz")

summary <- data.frame(
    snapshot = c("initial", "final"),
    time = c(legacy_targets$initial$time, legacy_targets$final$time),
    fish_density_dim = c(
        paste(dim(legacy_targets$initial$N), collapse = " x "),
        paste(dim(legacy_targets$final$N), collapse = " x ")
    ),
    resource_bins = c(
        length(legacy_targets$initial$nPP),
        length(legacy_targets$final$nPP)
    ),
    total_biomass = c(
        sum(legacy_targets$initial$biomass),
        sum(legacy_targets$final$biomass)
    ),
    total_yield = c(
        sum(legacy_targets$initial$yield),
        sum(legacy_targets$final$yield)
    )
)
print(summary)
