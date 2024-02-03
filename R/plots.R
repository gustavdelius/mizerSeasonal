#' Plot density-independent reproduction rate over time
#'
#' @param sim A MizerSim object
#' @param sim An object of class MizerSim
#' @param sim2 An optional second object of class MizerSim. If this is provided
#'   its RDIs will be shown on the same plot in bolder lines.
#' @param species The species to be selected. Optional. By default all target
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#' @param total A boolean value that determines whether the total over all
#'   species in the system is plotted as well. Note that even if the plot only
#'   shows a selection of species, the total is including all species. Default
#'   is FALSE.
#' @param log Boolean whether RDI should be plotted on a logarithmic axis.
#'   Defaults to true.
#' @param highlight	Name or vector of names of the species to be highlighted.
#' @param return_data A boolean value that determines whether the formatted data
#'   used for the plot is returned instead of the plot itself. Default value is
#'   FALSE
#' @param ... Other arguments (currently unused)
#'
#' @return A ggplot2 object, unless return_data = TRUE, in which case a data
#'   frame with the three variables 'Year', 'RDI', 'Species' is returned.
#' @export
#' @family plotting functions
plotRDI <- function(sim, sim2,
                    species = NULL,
                    total = FALSE, log = FALSE,
                    highlight = NULL, return_data = FALSE,
                    ...) {
    assert_that(is(sim, "MizerSim"),
                is.flag(total),
                is.flag(log),
                is.flag(return_data))
    params <- sim@params
    species <- valid_species_arg(sim, species, error_on_empty = TRUE)
    if (missing(sim2)) {
        y <- getTimeseries(sim, ...)
        y_total <- rowSums(y)
        y <- y[, (as.character(dimnames(y)[[2]]) %in% species),
               drop = FALSE]
        if (total) {
            # Include total
            y <- cbind(y, "Total" = y_total)
        }
        plot_dat <- reshape2::melt(y, varnames = c("Year", "Species"),
                                   value.name = "RDI")
        plot_dat <- subset(plot_dat, plot_dat$RDI > 0)
        # plotDataFrame() needs the columns in a particular order
        plot_dat <- plot_dat[, c(1, 3, 2)]
        
        if (nrow(plot_dat) == 0) {
            warning("There is no RDI to include.")
        }
        if (return_data) return(plot_dat)
        
        plotDataFrame(plot_dat, params,
                      ylab = "RDI [1/year]",
                      ytrans = ifelse(log, "log10", "identity"),
                      highlight = highlight)
    } else {
        # We need to combine two plots
        if (!all(dimnames(sim@n)$time == dimnames(sim2@n)$time)) {
            stop("The two simulations do not have the same times")
        }
        ym <- plotRDI(sim, species = species,
                      total = total, log = log,
                      highlight = highlight, return_data = TRUE, ...)
        ym2 <- plotRDI(sim2, species = species,
                       total = total, log = log,
                       highlight = highlight, return_data = TRUE, ...)
        ym$Simulation <- rep(1, nrow(ym)) # We don't use recycling because that
        # fails when there are zero rows.
        ym2$Simulation <- rep(2, nrow(ym2))
        ym <- rbind(ym, ym2)
        
        if (return_data) return(ym)
        
        plotDataFrame(ym, params,
                      ylab = "RDI [g/year]",
                      ytrans = ifelse(log, "log10", "identity"),
                      highlight = highlight, wrap_var = "Simulation")
    }
}

#' Plot density-dependent reproduction rate over time
#'
#' @param sim A MizerSim object
#' @param sim An object of class MizerSim
#' @param sim2 An optional second object of class MizerSim. If this is provided
#'   its RDDs will be shown on the same plot in bolder lines.
#' @param species The species to be selected. Optional. By default all target
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#' @param total A boolean value that determines whether the total over all
#'   species in the system is plotted as well. Note that even if the plot only
#'   shows a selection of species, the total is including all species. Default
#'   is FALSE.
#' @param log Boolean whether RDD should be plotted on a logarithmic axis.
#'   Defaults to true.
#' @param highlight	Name or vector of names of the species to be highlighted.
#' @param return_data A boolean value that determines whether the formatted data
#'   used for the plot is returned instead of the plot itself. Default value is
#'   FALSE
#' @param ... Other arguments (currently unused)
#'
#' @return A ggplot2 object, unless return_data = TRUE, in which case a data
#'   frame with the three variables 'Year', 'RDD', 'Species' is returned.
#' @export
#' @family plotting functions
plotRDD <- function(sim, sim2,
                    species = NULL,
                    total = FALSE, log = FALSE,
                    highlight = NULL, return_data = FALSE,
                    ...) {
    assert_that(is(sim, "MizerSim"),
                is.flag(total),
                is.flag(log),
                is.flag(return_data))
    params <- sim@params
    species <- valid_species_arg(sim, species, error_on_empty = TRUE)
    if (missing(sim2)) {
        y <- getTimeseries(sim, getRDD, ...)
        y_total <- rowSums(y)
        y <- y[, (as.character(dimnames(y)[[2]]) %in% species),
               drop = FALSE]
        if (total) {
            # Include total
            y <- cbind(y, "Total" = y_total)
        }
        plot_dat <- reshape2::melt(y, varnames = c("Year", "Species"),
                                   value.name = "RDD")
        plot_dat <- subset(plot_dat, plot_dat$RDD > 0)
        # plotDataFrame() needs the columns in a particular order
        plot_dat <- plot_dat[, c(1, 3, 2)]
        
        if (nrow(plot_dat) == 0) {
            warning("There is no RDD to include.")
        }
        if (return_data) return(plot_dat)
        
        plotDataFrame(plot_dat, params,
                      ylab = "RDD [1/year]",
                      ytrans = ifelse(log, "log10", "identity"),
                      highlight = highlight)
    } else {
        # We need to combine two plots
        if (!all(dimnames(sim@n)$time == dimnames(sim2@n)$time)) {
            stop("The two simulations do not have the same times")
        }
        ym <- plotRDD(sim, species = species,
                      total = total, log = log,
                      highlight = highlight, return_data = TRUE, ...)
        ym2 <- plotRDD(sim2, species = species,
                       total = total, log = log,
                       highlight = highlight, return_data = TRUE, ...)
        ym$Simulation <- rep(1, nrow(ym)) # We don't use recycling because that
        # fails when there are zero rows.
        ym2$Simulation <- rep(2, nrow(ym2))
        ym <- rbind(ym, ym2)
        
        if (return_data) return(ym)
        
        plotDataFrame(ym, params,
                      ylab = "RDD [g/year]",
                      ytrans = ifelse(log, "log10", "identity"),
                      highlight = highlight, wrap_var = "Simulation")
    }
}

#' Get time series of summary functions from simulation
#' 
#' Given a function that can calculate a vector of quantities of interest at a 
#' single time for all species, this function returns an array (time x species)
#' with the values at all the time steps saved in a simulation.
#' 
#' @param sim A MizerSim object
#' @param func The function calculating the quantities at a single time step
#' 
#' @return A matrix (time x species)
#' @export
getTimeseries <- function(sim, func = getRDI, ...) {
    params <- sim@params
    no_sp <- nrow(params@species_params)
    times <- getTimes(sim)
    time_indices <- seq_along(times)
    value <- array(0, dim = c(length(times), no_sp),
                   dimnames = dimnames(sim@n)[1:2])
    for (i in time_indices) {
        n <- sim@n[i, , ]
        dim(n) <- dim(params@initial_n)
        n_other <- sim@n_other[i, ]
        names(n_other) <- dimnames(sim@n_other)[[2]]
        value[i, ] <- func(params, 
                           n = n,
                           n_pp = sim@n_pp[i, ],
                           n_other = n_other,
                           t = times[i])
    }
    value
}

#' Plot per-capita gonadic mass through time at a given size
#' 
#' @param sim A MizerSim object
#' @param time_range The time range over which to plot
#' @param sizes A vector of sizes for each species at which to plot the 
#'     gonadic mass of an individual of that size.
plotGonadsVsTime <- function(sim,
                       time_range,
                       sizes=sim@params@species_params$w_max){
    def.par <- graphics::par(no.readonly = TRUE) #old pars
    if (missing(time_range)) {
        time_range  <- as.numeric(dimnames(sim@n)$time)
    }
    max_size <- sapply(sizes * (1 - 1e-6),
                       function(x,w){which.max(x <= w)}, w = sim@params@w)
    time_elements <- get_time_elements(sim, time_range)
    
    ###
    
    q_max <- matrix(0,length(time_elements), dim(sim@n)[2])
    qf <- sim@n_other[time_elements, "gonads"]
    
    qs_mat <- matrix(unlist(lapply(qf, function(x){return(diag(x[,max_size]))})),length(time_elements),dim(sim@n)[2],byrow=T)
    
    rownames(qs_mat) <- (as.numeric(dimnames(sim@n)$time))
    
    graphics::par(mfrow = c(4,3))
    graphics::par(oma = c(1,1,1,3))
    graphics::par(mar = c(2,4,2,0))
    for (i in 1:dim(sim@n)[2]) {
        plot((as.numeric(dimnames(sim@n)$time)), qs_mat[,i], 
             main = i, type = "l")
    }
    graphics::par(def.par)
}
