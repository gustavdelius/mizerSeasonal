# Legacy Datta–Blanchard model sources

These are the unchanged parts of the original R code bundle needed by
`data-raw/run_legacy_model.R` to reproduce the non-seasonal base model from
Datta & Blanchard (2016). The bundle was shared at
<https://figshare.com/s/e75f29d4cc9b94ae393b>.

- `paramNorthSeaModel.R` constructs the legacy North Sea parameter list.
- `SizeBasedModel.r` provides the legacy `Setup()` and `Project()` functions.
- `SelectivityFuncs.r` provides the fishing selectivity functions used during
  setup.
- `nsea_params.csv` contains the species parameters.
- `interactionmatrix_Schoener_twostage2D.RData` contains the interaction
  matrix.

The compatibility wrapper, output selection, and validation are kept in
`data-raw/run_legacy_model.R`; the original files in this directory should stay
byte-for-byte unchanged. From the repository root, regenerate the compact
targets with:

```sh
Rscript --vanilla data-raw/run_legacy_model.R
```

MD5 checksums of the copied files:

```text
f79522fe937ababded72057197de8a24  paramNorthSeaModel.R
b58f9562ad67f810aa234b6a2b779664  SizeBasedModel.r
3bffc39fa9e948177902fc71ceee7720  SelectivityFuncs.r
3ecf78bed2b0b415a60e17ebe2780dad  nsea_params.csv
57a6aacaf75dbd06be4f4b56392d4c4d  interactionmatrix_Schoener_twostage2D.RData
```
