# mizerSeasonal 0.0.0.9004

* Raised the minimum supported mizer version to 3.3.0 and updated the package
  to use mizer-managed dynamic marker classes and versioned extension metadata.
* Upgraded the packaged `datta_params` model for mizer 3.3.0 and added
  regression tests for its validity.
* Added `datta_params_second_order`, a separately constructed version of the
  Datta & Blanchard (2016) model using `second_order_w = TRUE` and initialised
  at the steady state of the second-order model.
* Added `data_raw/baseModel.Rda`, a trimmed-down version of the parameter
  object created with the code shared on figshare used by the 
  Datta & Blanchard (2016) paper.

# mizerSeasonal 0.0.0.9003

* Converted seasonal reproduction to mizer's composable extension mechanism.
  Seasonal reproduction and encounter rates now dispatch through the
  `mizerSeasonal` marker class, and the encounter method builds on other mizer
  extensions with `NextMethod()`.
* Added comprehensive tests for extension registration, seasonal dynamics,
  gonad-release functions, and plotting helpers.
* Fixed `plotGonadsVsTime()` so that selected sizes and time ranges are handled
  correctly.
* Added a vignette describing the package architecture and expanded the
  theoretical gonadic-mass vignette to cover growth diffusion, finite-volume
  bin averages, and second-order flux schemes.
* Added repository guidance and workflows for coding agents working with mizer
  extensions.

# mizerSeasonal 0.0.0.9002

* Updated the package for mizer 2.5.1.
* Added and documented `datta_params`, a `MizerParams` reconstruction of the
  base model used by Datta & Blanchard (2016), with its initial state set to the
  mizer steady state.
* Corrected indices in the gonadic-mass derivation and substantially expanded
  the function, getting-started, and theoretical documentation.
* Converted the theoretical vignette to Quarto and refreshed the package
  website.
* Adopted four-space indentation for mizer code in the lint configuration.

# mizerSeasonal 0.0.0.9001

* Implemented seasonal reproduction by accumulating gonadic mass and releasing
  it according to a configurable, time-dependent mass-specific release rate.
* Added von Mises, beta, beta-hazard, and Gaussian gonad-release functions, and
  corrected the expression used for the von Mises distribution.
* Added time-dependent density-dependent reproduction functions, including
  seasonal Beverton-Holt and von Mises forms, together with support for finding
  a steady state while reproduction and resource inputs are held independent
  of abundance.
* Added seasonal semichemostat resource dynamics for modelling resource blooms.
* Included the gonadic mass of prey in encounter rates by default, with an
  option to disable it in `setSeasonalReproduction()`.
* Added `animateGonadSpectra()`, `plotRDI()`, `plotRDD()`,
  `plotGonadsVsTime()`, and `getTimeseries()` for inspecting simulations.
* Added getting-started examples based on the North Sea model and a vignette
  deriving the gonadic-mass dynamics.
* Published the first pkgdown website and added spelling and build checks.

# mizerSeasonal 0.0.0.9000

* Created the package and began the first implementation of pulsed and seasonal
  reproduction.
* Added the gonad component and its dynamics for accumulating reproductive
  investment between spawning events.
* Added the initial package documentation and build configuration.
