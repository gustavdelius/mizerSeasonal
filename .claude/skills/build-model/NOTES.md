# Project notes

- `data-raw/baseModel.RData` is a trimmed legacy Datta & Blanchard model. It
  deliberately omits all projected time-series arrays but retains
  `initialN` (12 x 100) and `initialNResource` (length 130), the first weekly
  legacy state needed to reproduce the modern model initialization.
- With installed mizer 3.3.0 on 2026-08-25, rebuilding from that compact state
  reproduced the packaged first-order final fish and resource abundances
  exactly. The rebuilt second-order endpoint was not bitwise identical to the
  packaged object: species biomasses differed by at most 0.0141% and resource
  biomass by 0.0000100%. The discrepancy also occurs with the original legacy
  initialization and is therefore not caused by trimming `baseModel`.
- Rebuilding under current mizer also adds bookkeeping fields to the
  first-order species and resource parameter tables and changes the ordering of
  `rates_funcs`, while grids and all numerical rate arrays remain identical.
