# Hexalink v1.0.1 — Screening correction and Imo verification

This is a maintenance release of the **existing harmonic-flux Hexalink framework**. Its dynamic routing formulation is unchanged. The separate algebraic screening workflow was corrected and verified against `reprojected_elevation_2.tif`.

## Changes
- Normalize 1/2/4 terrain weights by valid-DEM-pixel area, retaining the prescribed runoff mean.
- Remove the artificial nonzero screening depth floor; S1 and S2 remain zero.
- Correct flat-top hexagon coverage of valid DEM pixels.
- Include synthetic screening regression tests and verified screening GeoTIFFs and Figure 6.
- Update repository documentation and reproduction instructions.

## Reproduced Imo screening results
- Valid DEM support: 4,980,592 pixels (30 m), area 4,482.5328 km²; 100% covered.
- S1/S2 mean depths: 0 m / 0 m.
- S3 mean: 0.075 m, range 0.04415954–0.17663816 m, runoff volume 336,189,960 m³.
- S4 mean: 0.175 m, range 0.10303893–0.41215571 m, runoff volume 784,443,240 m³.
- Independent raster-readback maximum absolute volume residual: 5.96 × 10⁻⁸ m³.

## Scope
These are equivalent algebraic screening depths, **not hydraulic flood simulations**. The existing harmonic-flux routing engine and its historical Imo results are **not modified or rerun** by this release.
