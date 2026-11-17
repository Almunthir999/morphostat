# Changelog

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
