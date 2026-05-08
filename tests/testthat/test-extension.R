local_params <- function(...) {
    sp <- NS_species_params
    # columns for seasonalVonMisesRelease (release function)
    sp$vonMises_r0    <- 1e10
    sp$vonMises_kappa <- 2
    sp$vonMises_mu    <- 0.25
    # columns for seasonalVonMisesRDD (density-dependent RDD function)
    sp$rdd_vonMises_r0    <- 1e10
    sp$rdd_vonMises_kappa <- 2
    sp$rdd_vonMises_mu    <- 0.25
    p <- newMultispeciesParams(sp, no_w = 20)
    setSeasonalReproduction(p, RDD = "seasonalVonMisesRDD", ...)
}

test_that("setSeasonalReproduction returns a mizerSeasonal object", {
    p <- local_params()
    expect_s4_class(p, "mizerSeasonal")
    expect_s4_class(p, "MizerParams")
})

test_that("extensions slot is populated after setSeasonalReproduction", {
    p <- local_params()
    expect_true("mizerSeasonal" %in% names(p@extensions))
})

test_that("RDI rate function stays at default (for S3 dispatch)", {
    p <- local_params()
    expect_equal(p@rates_funcs$RDI, "mizerRDI")
})

test_that("Encounter rate function stays at default (for S3 dispatch)", {
    p <- local_params()
    expect_equal(p@rates_funcs$Encounter, "mizerEncounter")
})

test_that("RDD rate function is set to the supplied function", {
    p <- local_params()
    expect_equal(p@rates_funcs$RDD, "seasonalVonMisesRDD")
})

test_that("gonads component is registered", {
    p <- local_params()
    expect_true("gonads" %in% names(p@other_dynamics))
})

test_that("include_gonads flag is stored in other_params", {
    p_on  <- local_params(include_gonads = TRUE)
    p_off <- local_params(include_gonads = FALSE)
    expect_true(other_params(p_on)$include_gonads)
    expect_false(other_params(p_off)$include_gonads)
})

test_that("projectRDI.mizerSeasonal is dispatched and uses gonadic mass", {
    p <- local_params()
    n       <- initialN(p)
    n_pp    <- initialNResource(p)
    e_repro <- getERepro(p)
    e_growth <- getEGrowth(p)
    mort    <- getMort(p)
    # With zero gonadic mass, seasonal RDI should be zero
    rdi_zero <- projectRDI(p, n = n, n_pp = n_pp,
                           n_other = list(gonads = n * 0),
                           t = 0, e_growth = e_growth, mort = mort,
                           e_repro = e_repro)
    expect_true(all(rdi_zero == 0))
    # With positive gonadic mass, seasonal RDI should be positive
    rdi_pos <- projectRDI(p, n = n, n_pp = n_pp,
                          n_other = list(gonads = n * 0.1),
                          t = 0.25, e_growth = e_growth, mort = mort,
                          e_repro = e_repro)
    expect_true(all(rdi_pos >= 0))
    expect_true(any(rdi_pos > 0))
})

test_that("projectEncounter.mizerSeasonal adds gonad contribution when include_gonads = TRUE", {
    p_on  <- local_params(include_gonads = TRUE)
    p_off <- local_params(include_gonads = FALSE)
    n       <- initialN(p_on)
    n_pp    <- initialNResource(p_on)
    gonads  <- n * 0.1
    n_other <- list(gonads = gonads)
    enc_on  <- projectEncounter(p_on,  n = n, n_pp = n_pp, n_other = n_other, t = 0)
    enc_off <- projectEncounter(p_off, n = n, n_pp = n_pp, n_other = n_other, t = 0)
    expect_true(all(enc_on >= enc_off))
    expect_true(any(enc_on > enc_off))
})

test_that("project() runs without error on a mizerSeasonal model", {
    p <- local_params()
    expect_no_error(project(p, t_max = 0.1, dt = 0.01, t_save = 0.1,
                            progress_bar = FALSE))
})

test_that("project() output is a mizerSeasonalSim object", {
    p <- local_params()
    sim <- project(p, t_max = 0.1, dt = 0.01, t_save = 0.1,
                   progress_bar = FALSE)
    expect_s4_class(sim, "mizerSeasonalSim")
})
