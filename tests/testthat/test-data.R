test_that("Datta params objects are valid for the installed mizer version", {
    params_objects <- list(datta_params, datta_params_second_order)

    for (params in params_objects) {
        validated <- mizer::validParams(params, info_level = 0)

        expect_identical(validated, params)
    }
})

test_that("datta_params_second_order uses the second-order size scheme", {
    expect_identical(
        mizer::second_order_w(datta_params_second_order),
        list(flux = "van_leer", bin_average = TRUE)
    )
})
