# Project notes

- The package now targets mizer 3.3.0; its previous declared minimum was
  2.5.4.9125.
- The package uses the native R pipe and therefore declares R >= 4.1.
- Since mizer 3.2, columns extracted from `species_params` carry species names.
  Release-function tests should preserve and expect those names rather than
  compare against unnamed scalars.
- The package is a dispatching extension. Marker classes are created
  dynamically by mizer, and `setSeasonalReproduction()` records only
  `mizerSeasonal` with `recordExtension()` before coercing the object.
- The bundled `datta_params` object was created with mizer 2.5.3. It was
  upgraded with `validParams()` to mizer 3.3.0; the 12-species model, its
  100 fish-size bins, 180 full-spectrum bins, and initial fish and resource
  abundances were preserved.
