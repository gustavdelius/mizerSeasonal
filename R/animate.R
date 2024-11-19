#' Animation of the gonadic mass spectra
#'
#' This function creates an animation of the gonadic mass spectra of the species
#' in the simulation. The gonadic mass spectrum is the size distribution of the
#' total gonadic mass of each species.
#'
#' @param sim A MizerSim object
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param time_range The time range to animate over. Either a vector of values
#'   or a vector of min and max time. Default is the entire time range of the
#'   simulation.
#' @param wlim A numeric vector of length two providing lower and upper limits
#'   for the w axis. Use NA to refer to the existing minimum or maximum.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis. Use NA to refer to the existing minimum or maximum. Any
#'   values below 1e-20 are always cut off.
#' @param power The gonadic mass density is plotted as the number density times
#'   the per-capita gonadic mass times the weight raised to \code{power}. The
#'   default \code{power = o} gives the gonadic mass density, whereas \code{power =
#'   1} gives the gonadic mass density with respect to logarithmic size bins.
#'
#' @return A plotly object
#' @export
#' @family plotting functions
animateGonadSpectra <- function(sim,
                                species = NULL,
                                time_range,
                                wlim = c(NA, NA),
                                ylim = c(NA, NA),
                                power = 0) {
    assert_that(is.number(power),
                length(wlim) == 2,
                length(ylim) == 2)

    species <- valid_species_arg(sim, species)
    if (missing(time_range)) {
        time_range  <- as.numeric(dimnames(sim@n)$time)
    }
    time_elements <- get_time_elements(sim, time_range)

    qf <- list_of_matrices_to_df(sim@n_other[time_elements, "gonads"])
    nf <- melt(sim@n[time_elements,
                     as.character(dimnames(sim@n)$sp) %in% species,
                     , drop = FALSE])
    nf <- merge(qf, nf, by = c("time", "sp", "w"))
    nf$w <- as.numeric(nf$w)
    nf$Q <- nf$value * nf$q

    # Deal with power argument ----
    if (power %in% c(0, 1)) {
        y_label <- c("Gonad density [1/g]", "Gonad density")[power + 1]
    } else {
        y_label <- paste0("Gonad density * w^", power)
    }
    nf$Q <- nf$Q * nf$w^power

    # Impose limits ----
    if (is.na(wlim[1])) wlim[1] <- min(sim@params@w) / 100
    if (is.na(wlim[2])) wlim[2] <- max(sim@params@w_full)
    if (is.na(ylim[1])) ylim[1] <- 10^-100
    if (is.na(ylim[2])) ylim[2] <- Inf
    nf <- nf %>%
        filter(Q >= ylim[1],
               Q <= ylim[2],
               w >= wlim[1],
               w <= wlim[2])

    nf %>%
        plotly::plot_ly() %>%
        plotly::add_lines(x = ~w, y = ~Q,
                          color = ~sp, colors = sim@params@linecolour,
                          frame = ~time,
                          line = list(simplify = FALSE)) %>%
        plotly::layout(xaxis = list(type = "log", exponentformat = "power",
                                    title = "Size [g]"),
                       yaxis = list(type = "log", exponentformat = "power",
                                    title = y_label))
}

# Function to transform matrix to dataframe
matrix_to_df <- function(mat, time) {
    sp <- rownames(mat)
    w <- colnames(mat)

    data.frame(
        time = rep(time, each = length(sp) * length(w)),
        sp = rep(sp, times = length(w)),
        w = rep(w, each = length(sp)),
        q = as.vector(mat)
    )
}

# Function to transform list of matrices to dataframe
# by applying the function to each matrix and combining the results
list_of_matrices_to_df <- function(mat_list) {
    do.call(rbind, lapply(names(mat_list), function(time) matrix_to_df(mat_list[[time]], time)))
}
