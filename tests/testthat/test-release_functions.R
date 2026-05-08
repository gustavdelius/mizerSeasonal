test_that("seasonalVonMisesRelease returns a numeric vector of length = no. species", {
    result <- seasonalVonMisesRelease(0.5, base_params)
    expect_type(result, "double")
    expect_length(result, nrow(base_params@species_params))
    expect_true(all(result > 0))
})

test_that("seasonalVonMisesRelease is maximal at t = mu", {
    sp <- base_params@species_params
    result <- seasonalVonMisesRelease(sp$vonMises_mu, base_params)
    expected <- sp$vonMises_r0 * exp(sp$vonMises_kappa) /
        (2 * pi * besselI(sp$vonMises_kappa, nu = 0))
    expect_equal(result, expected)
})

test_that("seasonalVonMisesRelease is minimal at t = mu + 0.5", {
    sp <- base_params@species_params
    result <- seasonalVonMisesRelease(sp$vonMises_mu + 0.5, base_params)
    expected <- sp$vonMises_r0 * exp(-sp$vonMises_kappa) /
        (2 * pi * besselI(sp$vonMises_kappa, nu = 0))
    expect_equal(result, expected)
})

test_that("seasonalVonMisesRelease is periodic with period 1", {
    t <- 0.3
    expect_equal(seasonalVonMisesRelease(t, base_params),
                 seasonalVonMisesRelease(t + 1, base_params))
    expect_equal(seasonalVonMisesRelease(t, base_params),
                 seasonalVonMisesRelease(t + 2, base_params))
})

test_that("seasonalVonMisesRelease rate at mu is greater than at all other times", {
    sp <- base_params@species_params
    rate_at_mu <- seasonalVonMisesRelease(sp$vonMises_mu, base_params)
    times <- seq(0, 0.99, by = 0.05)
    times <- times[times != sp$vonMises_mu]
    other_rates <- sapply(times, function(t) seasonalVonMisesRelease(t, base_params))
    expect_true(all(rate_at_mu > other_rates))
})

test_that("seasonalBetaHazardRelease returns a numeric vector of length = no. species", {
    result <- seasonalBetaHazardRelease(0.5, base_params)
    expect_type(result, "double")
    expect_length(result, nrow(base_params@species_params))
    expect_true(all(is.finite(result)))
    expect_true(all(result > 0))
})

test_that("seasonalBetaHazardRelease matches the hazard rate of the beta distribution", {
    t <- 0.3
    new_t <- t - floor(t)
    sp <- base_params@species_params
    expected <- stats::dbeta(new_t, sp$beta_a, sp$beta_b) /
        (1 - stats::pbeta(new_t, sp$beta_a, sp$beta_b))
    expect_equal(seasonalBetaHazardRelease(t, base_params), expected)
})

test_that("seasonalBetaHazardRelease is periodic with period 1", {
    t <- 0.4
    expect_equal(seasonalBetaHazardRelease(t, base_params),
                 seasonalBetaHazardRelease(t + 1, base_params))
    expect_equal(seasonalBetaHazardRelease(t, base_params),
                 seasonalBetaHazardRelease(t + 2, base_params))
})

test_that("seasonalBetaRelease returns a numeric vector of length = no. species", {
    result <- seasonalBetaRelease(0.3, base_params)
    expect_type(result, "double")
    expect_length(result, nrow(base_params@species_params))
    expect_true(all(result >= 0))
})

test_that("seasonalBetaRelease equals beta_r * dbeta(t - floor(t), beta_a, beta_b)", {
    t <- 0.35
    new_t <- t - floor(t)
    sp <- base_params@species_params
    expected <- sp$beta_r * stats::dbeta(new_t, sp$beta_a, sp$beta_b)
    expect_equal(seasonalBetaRelease(t, base_params), expected)
})

test_that("seasonalBetaRelease is zero at integer times (for beta_a > 1)", {
    # dbeta(0, 2, 5) = 0 since beta_a = 2 > 1
    expect_equal(seasonalBetaRelease(0, base_params), 0)
    expect_equal(seasonalBetaRelease(1, base_params), 0)
    expect_equal(seasonalBetaRelease(2, base_params), 0)
})

test_that("seasonalBetaRelease is periodic with period 1", {
    t <- 0.35
    expect_equal(seasonalBetaRelease(t, base_params),
                 seasonalBetaRelease(t + 1, base_params))
    expect_equal(seasonalBetaRelease(t, base_params),
                 seasonalBetaRelease(t + 2, base_params))
})

test_that("seasonalGaussianRelease returns a numeric vector of length = no. species", {
    result <- seasonalGaussianRelease(0.3, base_params)
    expect_type(result, "double")
    expect_length(result, nrow(base_params@species_params))
    expect_true(all(result > 0))
})

test_that("seasonalGaussianRelease is maximal at t = sr_t0", {
    sp <- base_params@species_params
    result <- seasonalGaussianRelease(sp$sr_t0, base_params)
    expect_equal(result, sp$sr_r0)
})

test_that("seasonalGaussianRelease follows Gaussian formula", {
    t <- 0.3
    sp <- base_params@species_params
    new_t <- t - trunc(t)
    expected <- sp$sr_r0 * exp(-(new_t - sp$sr_t0)^2 / (2 * sp$sr_sigma^2))
    expect_equal(seasonalGaussianRelease(t, base_params), expected)
})

test_that("seasonalGaussianRelease drops to exp(-0.5) * r0 at t = sr_t0 + sr_sigma", {
    sp <- base_params@species_params
    result <- seasonalGaussianRelease(sp$sr_t0 + sp$sr_sigma, base_params)
    expect_equal(result, sp$sr_r0 * exp(-0.5))
})

test_that("seasonalGaussianRelease is periodic with period 1", {
    t <- 0.3
    expect_equal(seasonalGaussianRelease(t, base_params),
                 seasonalGaussianRelease(t + 1, base_params))
    expect_equal(seasonalGaussianRelease(t, base_params),
                 seasonalGaussianRelease(t + 2, base_params))
})
