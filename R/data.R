#' Parameters used by Datta & Blanchard (2016)
#' 
#' This MizerParams object is an approximation of the model used by Datta &
#' Blanchard (2016), which allows to reproduce their results and to extend it
#' to the full gonadic mass dynamics implemented in mizerSeasonal.
#' 
#' The paper "The effects of seasonal processes on size spectrum dynamics" by
#' Samik Datta and Julia Blanchard, published in Can. J. Fish. Aquat. Sci. 73, 
#' 598–610 (2016), was the first paper to investigate seasonal processes in
#' the context of size spectrum dynamics. The authors used an early version of 
#' the mizer code that they then modified. Their code is available at
#' <https://figshare.com/s/e75f29d4cc9b94ae393b>.
#' 
#' This MizerParams object is based on the parameters used in the paper. 
#' However the choice of size bins for the numerical implementation of the
#' mizer equations is different from the one used in the paper. The paper used
#' wider size-bins for the resource spectrum than for the fish spectrum. The
#' modern mizer code no longer supports this. So this MizerParams object uses
#' the same size bins for the resource and fish spectra. Like the model by
#' Datta and Blanchard, we use 100 size bins for the fish spectra but we use
#' 180 size bins for the full spectra instead of 130.
#'
#' `datta_params` uses mizer's default first-order size scheme.
#' `datta_params_second_order` was created with `second_order_w = TRUE`, which
#' uses the second-order van Leer flux scheme and bin-averaged quantities.
#'
#' The initial state of each MizerParams object is set to the steady state of
#' its model.
#' 
#' The script that created this MizerParams object is available at
#' <https://github.com/gustavdelius/mizerSeasonal/blob/main/data-raw/datta_params.R>.
#' 
#' @format A [mizer::MizerParams] object. It is stored as an S3 object, so R's
#'   standard lazy-loading preserves its class vector and metadata; no load hook
#'   or active binding is needed. It is a plain mizer model: call
#'   [setSeasonalReproduction()] on it to obtain a `mizerSeasonal` object.
#' @source Datta, S. & Blanchard, J. L. "The effects of seasonal processes on
#'  size spectrum dynamics". Canadian Journal of Fisheries and Aquatic Sciences
#'  (2016). <https://cdnsciencepub.com/doi/full/10.1139/cjfas-2015-0468>
"datta_params"

#' Second-order parameters used by Datta & Blanchard (2016)
#'
#' A variant of [datta_params] created with `second_order_w = TRUE`. It uses
#' the second-order van Leer flux scheme and bin-averaged quantities, and its
#' initial state is set to the steady state of the second-order model.
#'
#' @format A [mizer::MizerParams] object, stored as an S3 object like
#'   [datta_params].
#' @source Datta, S. & Blanchard, J. L. "The effects of seasonal processes on
#'  size spectrum dynamics". Canadian Journal of Fisheries and Aquatic Sciences
#'  (2016). <https://cdnsciencepub.com/doi/full/10.1139/cjfas-2015-0468>
"datta_params_second_order"
