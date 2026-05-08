## setSeasonalReproduction ----

test_that("setSeasonalReproduction returns a MizerParams object", {
    result <- suppressMessages(
        setSeasonalReproduction(base_params,
                                release_func = "seasonalVonMisesRelease",
                                RDD = "seasonalVonMisesRDD")
    )
    expect_s4_class(result, "MizerParams")
})

test_that("setSeasonalReproduction adds a gonads component", {
    expect_true("gonads" %in% names(seasonal_params@other_dynamics))
})

test_that("setSeasonalReproduction registers gonadDynamics as component dynamics", {
    expect_equal(seasonal_params@other_dynamics$gonads, "gonadDynamics")
})

test_that("setSeasonalReproduction leaves RDI at default so S3 dispatch takes over", {
    expect_equal(seasonal_params@rates_funcs$RDI, "mizerRDI")
})

test_that("setSeasonalReproduction registers the specified RDD function", {
    result <- suppressMessages(
        setSeasonalReproduction(base_params,
                                release_func = "seasonalVonMisesRelease",
                                RDD = "seasonalVonMisesRDD")
    )
    expect_equal(result@rates_funcs$RDD, "seasonalVonMisesRDD")
})

test_that("setSeasonalReproduction stores release_func name in other_params", {
    result <- suppressMessages(
        setSeasonalReproduction(base_params,
                                release_func = "seasonalGaussianRelease",
                                RDD = "seasonalVonMisesRDD")
    )
    expect_equal(other_params(result)$release_func, "seasonalGaussianRelease")
})

test_that("setSeasonalReproduction leaves Encounter at default so S3 dispatch takes over", {
    result <- suppressMessages(
        setSeasonalReproduction(base_params,
                                release_func = "seasonalVonMisesRelease",
                                RDD = "seasonalVonMisesRDD",
                                include_gonads = TRUE)
    )
    expect_equal(result@rates_funcs$Encounter, "mizerEncounter")
})

test_that("setSeasonalReproduction stores include_gonads = FALSE in other_params", {
    result <- suppressMessages(
        setSeasonalReproduction(base_params,
                                release_func = "seasonalVonMisesRelease",
                                RDD = "seasonalVonMisesRDD",
                                include_gonads = FALSE)
    )
    expect_false(other_params(result)$include_gonads)
})

test_that("setSeasonalReproduction initialises gonads to zero", {
    initial_gonads <- seasonal_params@initial_n_other$gonads
    expect_true(all(initial_gonads == 0))
})


## gonadDynamics ----

test_that("gonadDynamics returns an array with species x size dimensions", {
    no_sp <- nrow(seasonal_params@species_params)
    no_w <- length(seasonal_params@w)

    gonads <- matrix(0, nrow = no_sp, ncol = no_w,
                     dimnames = dimnames(seasonal_params@initial_n))
    rates <- list(
        e_growth = matrix(0.1, nrow = no_sp, ncol = no_w,
                          dimnames = dimnames(seasonal_params@initial_n)),
        e_repro  = matrix(0.01, nrow = no_sp, ncol = no_w,
                          dimnames = dimnames(seasonal_params@initial_n))
    )

    result <- gonadDynamics(params = seasonal_params,
                            n_other = list(gonads = gonads),
                            rates = rates,
                            t = 0, dt = 0.1)
    expect_equal(dim(result), c(no_sp, no_w))
})

test_that("gonadDynamics returns non-negative values", {
    no_sp <- nrow(seasonal_params@species_params)
    no_w <- length(seasonal_params@w)

    gonads <- matrix(0, nrow = no_sp, ncol = no_w,
                     dimnames = dimnames(seasonal_params@initial_n))
    rates <- list(
        e_growth = matrix(0.1, nrow = no_sp, ncol = no_w,
                          dimnames = dimnames(seasonal_params@initial_n)),
        e_repro  = matrix(0.01, nrow = no_sp, ncol = no_w,
                          dimnames = dimnames(seasonal_params@initial_n))
    )

    result <- gonadDynamics(params = seasonal_params,
                            n_other = list(gonads = gonads),
                            rates = rates,
                            t = 0, dt = 0.1)
    expect_true(all(result >= 0))
})

test_that("gonadDynamics returns zero when e_repro is zero and initial gonads are zero", {
    no_sp <- nrow(seasonal_params@species_params)
    no_w <- length(seasonal_params@w)

    gonads <- matrix(0, nrow = no_sp, ncol = no_w,
                     dimnames = dimnames(seasonal_params@initial_n))
    rates <- list(
        e_growth = matrix(0.1, nrow = no_sp, ncol = no_w,
                          dimnames = dimnames(seasonal_params@initial_n)),
        e_repro  = matrix(0, nrow = no_sp, ncol = no_w,
                          dimnames = dimnames(seasonal_params@initial_n))
    )

    result <- gonadDynamics(params = seasonal_params,
                            n_other = list(gonads = gonads),
                            rates = rates,
                            t = 0, dt = 0.1)
    expect_true(all(result == 0))
})


## projectRDI.mizerSeasonal ----

test_that("projectRDI dispatches to seasonal method and returns a numeric vector", {
    no_sp <- nrow(seasonal_params@species_params)
    no_w <- length(seasonal_params@w)
    gonads <- matrix(1, nrow = no_sp, ncol = no_w,
                     dimnames = dimnames(seasonal_params@initial_n))
    e_repro  <- getERepro(seasonal_params)
    e_growth <- getEGrowth(seasonal_params)
    mort     <- getMort(seasonal_params)

    result <- projectRDI(params = seasonal_params,
                         n = seasonal_params@initial_n,
                         n_pp = seasonal_params@initial_n_pp,
                         n_other = list(gonads = gonads),
                         t = 0,
                         e_growth = e_growth, mort = mort, e_repro = e_repro)
    expect_type(result, "double")
    expect_length(result, no_sp)
    expect_true(all(result >= 0))
})

test_that("projectRDI returns zero when gonadic mass is zero", {
    no_sp <- nrow(seasonal_params@species_params)
    no_w <- length(seasonal_params@w)
    gonads <- matrix(0, nrow = no_sp, ncol = no_w,
                     dimnames = dimnames(seasonal_params@initial_n))
    e_repro  <- getERepro(seasonal_params)
    e_growth <- getEGrowth(seasonal_params)
    mort     <- getMort(seasonal_params)

    result <- projectRDI(params = seasonal_params,
                         n = seasonal_params@initial_n,
                         n_pp = seasonal_params@initial_n_pp,
                         n_other = list(gonads = gonads),
                         t = 0,
                         e_growth = e_growth, mort = mort, e_repro = e_repro)
    expect_equal(result, setNames(rep(0, no_sp),
                                  rownames(seasonal_params@species_params)))
})

test_that("projectRDI scales linearly with gonadic mass", {
    no_sp <- nrow(seasonal_params@species_params)
    no_w <- length(seasonal_params@w)
    gonads1 <- matrix(1, nrow = no_sp, ncol = no_w,
                      dimnames = dimnames(seasonal_params@initial_n))
    gonads2 <- matrix(2, nrow = no_sp, ncol = no_w,
                      dimnames = dimnames(seasonal_params@initial_n))
    e_repro  <- getERepro(seasonal_params)
    e_growth <- getEGrowth(seasonal_params)
    mort     <- getMort(seasonal_params)

    rdi1 <- projectRDI(params = seasonal_params,
                       n = seasonal_params@initial_n,
                       n_pp = seasonal_params@initial_n_pp,
                       n_other = list(gonads = gonads1),
                       t = 0,
                       e_growth = e_growth, mort = mort, e_repro = e_repro)
    rdi2 <- projectRDI(params = seasonal_params,
                       n = seasonal_params@initial_n,
                       n_pp = seasonal_params@initial_n_pp,
                       n_other = list(gonads = gonads2),
                       t = 0,
                       e_growth = e_growth, mort = mort, e_repro = e_repro)
    expect_equal(rdi2, 2 * rdi1)
})


## seasonalBevertonHoltRDD ----

test_that("seasonalBevertonHoltRDD applies the Beverton-Holt formula", {
    p <- seasonal_params
    sp_name <- rownames(seasonal_params@species_params)
    r_max_mat <- matrix(100, nrow = 1, ncol = 1,
                        dimnames = list(time = "0", species = sp_name))
    other_params(p)$r_max <- r_max_mat

    rdi <- setNames(50, sp_name)
    result <- seasonalBevertonHoltRDD(rdi = rdi, params = p, t = 0)

    expected <- rdi * 100 / (rdi + 100)  # = 33.333...
    expect_equal(result, expected)
})

test_that("seasonalBevertonHoltRDD returns rdi when r_max is Inf", {
    p <- seasonal_params
    sp_name <- rownames(seasonal_params@species_params)
    r_max_mat <- matrix(Inf, nrow = 1, ncol = 1,
                        dimnames = list(time = "0", species = sp_name))
    other_params(p)$r_max <- r_max_mat

    rdi <- setNames(50, sp_name)
    result <- seasonalBevertonHoltRDD(rdi = rdi, params = p, t = 0)
    expect_equal(result, rdi)
})

test_that("seasonalBevertonHoltRDD result is never greater than r_max", {
    p <- seasonal_params
    sp_name <- rownames(seasonal_params@species_params)
    r_max_val <- 100
    r_max_mat <- matrix(r_max_val, nrow = 1, ncol = 1,
                        dimnames = list(time = "0", species = sp_name))
    other_params(p)$r_max <- r_max_mat

    for (rdi_val in c(1, 50, 100, 500, 1e6)) {
        rdi <- setNames(rdi_val, sp_name)
        result <- seasonalBevertonHoltRDD(rdi = rdi, params = p, t = 0)
        expect_lte(result, r_max_val)
    }
})

test_that("seasonalBevertonHoltRDD errors when r_max not defined at time t", {
    p <- seasonal_params
    sp_name <- rownames(seasonal_params@species_params)
    r_max_mat <- matrix(100, nrow = 1, ncol = 1,
                        dimnames = list(time = "0.5", species = sp_name))
    other_params(p)$r_max <- r_max_mat

    rdi <- setNames(50, sp_name)
    expect_error(seasonalBevertonHoltRDD(rdi = rdi, params = p, t = 0),
                 "r_max not defined")
})


## seasonalVonMisesRDD ----

test_that("seasonalVonMisesRDD returns a numeric vector of length = no. species", {
    result <- seasonalVonMisesRDD(params = base_params, t = 0.25)
    expect_type(result, "double")
    expect_length(result, nrow(base_params@species_params))
    expect_true(all(result > 0))
})

test_that("seasonalVonMisesRDD is maximal at t = mu", {
    sp <- base_params@species_params
    result <- seasonalVonMisesRDD(params = base_params, t = sp$rdd_vonMises_mu)
    expected <- sp$rdd_vonMises_r0 * exp(sp$rdd_vonMises_kappa) /
        (2 * pi * besselI(sp$rdd_vonMises_kappa, nu = 0))
    expect_equal(result, expected)
})

test_that("seasonalVonMisesRDD is periodic with period 1", {
    t <- 0.3
    expect_equal(seasonalVonMisesRDD(params = base_params, t = t),
                 seasonalVonMisesRDD(params = base_params, t = t + 1))
})


## projectEncounter.mizerSeasonal ----

test_that("projectEncounter dispatches to seasonal method and returns correct dimensions", {
    no_sp <- nrow(seasonal_params@species_params)
    no_w <- length(seasonal_params@w)
    gonads <- matrix(0, nrow = no_sp, ncol = no_w,
                     dimnames = dimnames(seasonal_params@initial_n))

    result <- projectEncounter(params = seasonal_params,
                               n = seasonal_params@initial_n,
                               n_pp = seasonal_params@initial_n_pp,
                               n_other = list(gonads = gonads),
                               t = 0)
    expect_equal(dim(result), c(no_sp, no_w))
    expect_true(all(result >= 0))
})

test_that("projectEncounter with zero gonads equals mizerEncounter", {
    no_sp <- nrow(seasonal_params@species_params)
    no_w <- length(seasonal_params@w)
    gonads <- matrix(0, nrow = no_sp, ncol = no_w,
                     dimnames = dimnames(seasonal_params@initial_n))

    enc_seasonal <- projectEncounter(params = seasonal_params,
                                     n = seasonal_params@initial_n,
                                     n_pp = seasonal_params@initial_n_pp,
                                     n_other = list(gonads = gonads),
                                     t = 0)
    enc_standard <- mizerEncounter(params = seasonal_params,
                                   n = seasonal_params@initial_n,
                                   n_pp = seasonal_params@initial_n_pp,
                                   n_other = list(gonads = gonads),
                                   t = 0)
    expect_equal(enc_seasonal, enc_standard)
})

test_that("projectEncounter is higher than mizerEncounter when gonads are positive", {
    # newSingleSpeciesParams() has interaction[1,1] = 0 (no self-predation),
    # so we enable it to allow gonadic mass of fish prey to enter encounter.
    p <- seasonal_params
    p@interaction[1, 1] <- 1
    n <- p@initial_n
    gonads <- n * 0.1
    n_other <- list(gonads = gonads)

    enc_seasonal <- projectEncounter(params = p,
                                     n = n,
                                     n_pp = p@initial_n_pp,
                                     n_other = n_other,
                                     t = 0)
    enc_standard <- mizerEncounter(params = p,
                                   n = n,
                                   n_pp = p@initial_n_pp,
                                   n_other = n_other,
                                   t = 0)
    expect_true(all(enc_seasonal >= enc_standard))
    expect_true(any(enc_seasonal > enc_standard))
})


## seasonal_resource_semichemostat ----

test_that("seasonal_resource_semichemostat returns a vector of length = no. resource bins", {
    p <- base_params
    p@resource_params$rp <- list(kappa = 2, mu = 0.25, maxR = 0.5)
    no_w_full <- length(p@w_full)

    result <- seasonal_resource_semichemostat(
        params = p,
        n = p@initial_n,
        n_pp = p@initial_n_pp,
        n_other = list(),
        rates = list(resource_mort = rep(0.1, no_w_full)),
        t = 0, dt = 0.1,
        resource_rate = p@rr_pp,
        resource_capacity = p@cc_pp
    )
    expect_type(result, "double")
    expect_length(result, no_w_full)
    expect_true(all(result >= 0))
})

test_that("seasonal_resource_semichemostat gives higher capacity at t = mu than t = mu + 0.5", {
    p <- base_params
    mu <- 0.25
    p@resource_params$rp <- list(kappa = 2, mu = mu, maxR = 1)
    no_w_full <- length(p@w_full)
    rates <- list(resource_mort = rep(0.1, no_w_full))

    result_peak <- seasonal_resource_semichemostat(
        params = p, n = p@initial_n, n_pp = p@initial_n_pp,
        n_other = list(), rates = rates, t = mu, dt = 0.1,
        resource_rate = p@rr_pp, resource_capacity = p@cc_pp
    )
    result_trough <- seasonal_resource_semichemostat(
        params = p, n = p@initial_n, n_pp = p@initial_n_pp,
        n_other = list(), rates = rates, t = mu + 0.5, dt = 0.1,
        resource_rate = p@rr_pp, resource_capacity = p@cc_pp
    )
    # At peak (t = mu), capacity is amplified more, so steady state n_pp is higher
    expect_true(sum(result_peak) > sum(result_trough))
})

test_that("seasonal_resource_semichemostat is periodic with period 1", {
    p <- base_params
    p@resource_params$rp <- list(kappa = 2, mu = 0.25, maxR = 0.5)
    no_w_full <- length(p@w_full)
    rates <- list(resource_mort = rep(0.1, no_w_full))

    result_t <- seasonal_resource_semichemostat(
        params = p, n = p@initial_n, n_pp = p@initial_n_pp,
        n_other = list(), rates = rates, t = 0.3, dt = 0.1,
        resource_rate = p@rr_pp, resource_capacity = p@cc_pp
    )
    result_t1 <- seasonal_resource_semichemostat(
        params = p, n = p@initial_n, n_pp = p@initial_n_pp,
        n_other = list(), rates = rates, t = 1.3, dt = 0.1,
        resource_rate = p@rr_pp, resource_capacity = p@cc_pp
    )
    expect_equal(result_t, result_t1)
})
