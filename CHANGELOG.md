# Changelog

## v0.1.4 (2026-11)
Print-ready figures; statistical results unchanged.
- Figures are sized so that text is at least ~7 pt when placed at journal page width;
  medium-weight Inter font (falls back to Arial / DejaVu Sans), bold titles, thicker lines.
- Violin panels: extra headroom so the adjusted p-value never overlaps the data.
- Volcano plot: clearer gap between the left labels and the y-axis title.

## v0.1.3 (2026-11)
Figure readability release; statistical results are identical to v0.1.1/v0.1.2.
- Readable feature names in all figures (e.g. `Cells_AreaShape_FormFactor` ->
  "Circularity (form factor)", `Nuclei_Intensity_StdIntensity_CorrTub` ->
  "Nuc_Intensity_Std [Tub]"); full names remain in all output tables.
- Larger fonts and bold titles; legend inside the PCA plot.
- Volcano plot: labels for the strongest features outside the plot with leader lines
  (higher in treatment on the right, lower on the left), effect-size name in the axis label.
- Violin panel: the adjusted p-value of each comparison is printed in the panel
  (two-group runs); y-axis label follows the normalisation used.

## v0.1.2 (2026-11)
Usability release prompted by an independent test of the Bio-protocol instructions.
Statistical results are identical to v0.1.1.
- `morphostat run --help` now describes every option and shows its default.
- Clear error messages (no Python traceback) when the control label, group, batch or
  metadata column is not found; the message lists the labels/columns that exist.
- REPORT.md now lists the metadata columns, every column excluded automatically, and
  all analysed features; the violin-panel heading says whether features were chosen
  with `--panel`; PNG files are listed.
- `cleaning_report.txt` lists the automatically excluded (non-feature) columns.
- `run_manifest.json` records every setting plus the MorphoStat and Python versions.
- `normalized_features.csv` keeps the `__source_file` column (origin of each row).
- `run_validation.py` prints which article figure/table each output corresponds to.

## v0.1.1 (2026-11)
- Added `run_validation.py`: one-command reproduction of the BBBC021 validation
  reported in the Bio-protocol article (figures, tables, direction check, environment record).
- README: corrected the clone folder name (`cd morphostat`), the path to
  `prepare_bbbc021.py`, added virtual-environment and Windows instructions,
  a validation section, and a support section.
- Fix: CellProfiler `*_Location_*` columns (pixel coordinates of the brightest point
  in an object) are now excluded from the feature set, like other positional columns.
- Fix: Fiji/ImageJ Results-table columns that are not morphology features (the
  blank-named row-number column, X/Y, XM/YM, BX/BY, FeretX/FeretY, Slice, Ch, Frame)
  are now excluded automatically; covered by a new test in `test_smoke.py`.
- These two exclusions are the only changes to the analysis code; the BBBC021
  validation results are unaffected by the Fiji change.

## v0.1.0 (2026-05)
- Initial release.
