test_that("datta_params is valid for the installed mizer version", {
    validated <- mizer::validParams(datta_params, info_level = 0)

    expect_identical(validated, datta_params)
})
