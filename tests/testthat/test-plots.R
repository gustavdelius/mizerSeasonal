## getTimeseries ----

test_that("getTimeseries returns a matrix with time x species dimensions", {
    result <- getTimeseries(seasonal_sim)
    times <- getTimes(seasonal_sim)
    no_sp <- nrow(seasonal_sim@params@species_params)

    expect_true(is.matrix(result))
    expect_equal(nrow(result), length(times))
    expect_equal(ncol(result), no_sp)
})

test_that("getTimeseries returns non-negative values for RDI", {
    result <- getTimeseries(seasonal_sim)
    expect_true(all(result >= 0))
})

test_that("getTimeseries returns non-negative values for RDD", {
    result <- getTimeseries(seasonal_sim, func = getRDD)
    expect_true(all(result >= 0))
})

test_that("getTimeseries preserves dimnames from the simulation", {
    result <- getTimeseries(seasonal_sim)
    expect_equal(dimnames(result)[[1]], dimnames(seasonal_sim@n)[[1]])
    expect_equal(dimnames(result)[[2]], dimnames(seasonal_sim@n)[[2]])
})


## plotRDI ----

test_that("plotRDI with return_data = TRUE returns a data frame", {
    result <- plotRDI(seasonal_sim, return_data = TRUE)
    expect_s3_class(result, "data.frame")
})

test_that("plotRDI data frame has columns Year, RDI, Species", {
    result <- plotRDI(seasonal_sim, return_data = TRUE)
    expect_true(all(c("Year", "RDI", "Species") %in% names(result)))
})

test_that("plotRDI data frame has non-negative RDI values", {
    result <- plotRDI(seasonal_sim, return_data = TRUE)
    expect_true(all(result$RDI >= 0))
})

test_that("plotRDI without return_data returns a ggplot object", {
    result <- plotRDI(seasonal_sim, return_data = FALSE)
    expect_s3_class(result, "gg")
})

test_that("plotRDI with total = TRUE includes a Total row", {
    result <- plotRDI(seasonal_sim, total = TRUE, return_data = TRUE)
    expect_true("Total" %in% as.character(result$Species))
})

test_that("plotRDI errors on non-MizerSim input", {
    expect_error(plotRDI("not a sim"))
})


## plotRDD ----

test_that("plotRDD with return_data = TRUE returns a data frame", {
    result <- plotRDD(seasonal_sim, return_data = TRUE)
    expect_s3_class(result, "data.frame")
})

test_that("plotRDD data frame has columns Year, RDD, Species", {
    result <- plotRDD(seasonal_sim, return_data = TRUE)
    expect_true(all(c("Year", "RDD", "Species") %in% names(result)))
})

test_that("plotRDD data frame has non-negative RDD values", {
    result <- plotRDD(seasonal_sim, return_data = TRUE)
    expect_true(all(result$RDD >= 0))
})

test_that("plotRDD without return_data returns a ggplot object", {
    result <- plotRDD(seasonal_sim, return_data = FALSE)
    expect_s3_class(result, "gg")
})

test_that("plotRDD with total = TRUE includes a Total row", {
    result <- plotRDD(seasonal_sim, total = TRUE, return_data = TRUE)
    expect_true("Total" %in% as.character(result$Species))
})

test_that("plotRDD errors on non-MizerSim input", {
    expect_error(plotRDD("not a sim"))
})

test_that("plotRDI and plotRDD both return non-empty data frames for same sim", {
    rdi_data <- plotRDI(seasonal_sim, return_data = TRUE)
    rdd_data <- plotRDD(seasonal_sim, return_data = TRUE)
    expect_gt(nrow(rdi_data), 0)
    expect_gt(nrow(rdd_data), 0)
})

test_that("plotRDI with two simulations returns combined data frame", {
    result <- plotRDI(seasonal_sim, seasonal_sim, return_data = TRUE)
    expect_s3_class(result, "data.frame")
    expect_true("Simulation" %in% names(result))
})

test_that("plotRDD with two simulations returns combined data frame", {
    result <- plotRDD(seasonal_sim, seasonal_sim, return_data = TRUE)
    expect_s3_class(result, "data.frame")
    expect_true("Simulation" %in% names(result))
})


## plotGonadsVsTime ----

test_that("plotGonadsVsTime runs without error", {
    expect_no_error(plotGonadsVsTime(seasonal_sim))
})

test_that("plotGonadsVsTime runs without error for a specified time range", {
    expect_no_error(plotGonadsVsTime(seasonal_sim, time_range = c(0.5, 2)))
})


## animateGonadSpectra ----

test_that("animateGonadSpectra returns a plotly object", {
    result <- animateGonadSpectra(seasonal_sim)
    expect_s3_class(result, "plotly")
})

test_that("animateGonadSpectra works with a specified time_range", {
    result <- animateGonadSpectra(seasonal_sim, time_range = c(0.5, 2))
    expect_s3_class(result, "plotly")
})

test_that("animateGonadSpectra works with power = 1", {
    result <- animateGonadSpectra(seasonal_sim, power = 1)
    expect_s3_class(result, "plotly")
})
