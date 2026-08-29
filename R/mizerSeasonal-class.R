#' mizerSeasonal extension classes
#'
#' S3 extension classes for [mizer::MizerParams] and [mizer::MizerSim] that
#' enable S3 dispatch for this package's methods.
#'
#' The class names are ordinary entries in the object's S3 class vector, so no
#' class declaration is needed. All extension-specific data lives in
#' `other_params(params)` or in the `gonads` component (see
#' [mizer::setComponent()]).
#'
#' Objects of class `mizerSeasonal` are created by [setSeasonalReproduction()],
#' which records the extension on the object with [mizer::recordExtension()] and
#' then calls [mizer::coerceToExtensionClass()]. Objects of class
#' `mizerSeasonalSim` are returned automatically by [mizer::project()] when it
#' is called on a `mizerSeasonal` params object.
#'
#' @name mizerSeasonal-class
#' @keywords internal
NULL
