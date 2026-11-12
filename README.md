# Hexalink

A mass-conservative hexagonal finite-volume framework for hydrologic routing and rapid flood susceptibility screening, demonstrated on the Imo River Basin, Nigeria.

The repository provides the routing code, a separate scenario-based screening workflow, numerical tests, and example outputs.

## Contents

| File | Purpose | Description |
|---|---|---|
| `routing_engine.py` | Dynamic routing engine, six-neighbor edge fluxes, donor limiting and four-store accounting; includes bundled checks A–D | Routing and storage accounting |
| `run_headline.py` | Main 77,434-cell Imo harmonic-flux routing demonstration | Dynamic-routing example |
| `run_resolution.py` | Existing 75/100/300 m resolution-sensitivity tests | Resolution study |
| `run_sensitivity.py` | Existing one-at-a-time parameter-sensitivity experiments | Parameter sensitivity |
| `run_imo_scenarios.py` | **Area-conservative algebraic screening** for four rainfall scenarios | Algebraic screening |
| `extract_metrics.py` | Summary metrics for output rasters | Output summaries |
| `tests/test_screening_mass_balance.py` | Synthetic regression tests for the screening implementation | Screening regression checks |
| `environment.yml` | Conda environment | — |
| `reprojected_elevation_2.tif` | SRTM-derived 30 m Imo DEM in EPSG:32632 | DEM input |

## Requirements and setup

Use Python ≥3.9 and the dependencies in `environment.yml` (NumPy, Rasterio, SciPy, Shapely and GeoPandas). The Imo DEM is derived from public NASA/USGS SRTM data and is stored in the main repository.

```bash
conda env create -f environment.yml
conda activate hexalink
```

## Run the computations

```bash
python routing_engine.py
python run_headline.py
python run_resolution.py 75 result_R75.json
python run_sensitivity.py REF result_REF.json
python run_imo_scenarios.py --dem reprojected_elevation_2.tif --output-dir screening_outputs
python -m unittest discover -s tests -v
```

`run_imo_scenarios.py` generates the output files `imo_flood_{moderate|heavy|extreme|catastrophic}.tif` for compatibility. **The TIFF values are algebraic equivalent screening depths, not simulated flood depths or inundation predictions.** Its four scenarios impose one-hour rainfalls of 50/100/200/300 mm h⁻¹, prescribed infiltration capacity 120 mm h⁻¹ and illustrative ET subtraction 5 mm h⁻¹.

For each valid hexagon, DEM-elevation-class weights of 1, 2 and 4 are normalized by their **valid-DEM-pixel-area-weighted mean**. A proper flat-top hexagonal grid covers all valid raster pixels, and no nonzero depth floor or output clipping is applied. Thus S1 and S2 return exactly zero, while S3 and S4 have valid-area-weighted means of 0.075 m and 0.175 m, respectively. `screening_metrics.csv` and `screening_metrics.json` record scenario volumes and residuals. The screening workflow was executed on the Imo DEM (4,980,592 valid 30 m pixels; 4,482.5328 km²); the maximum screening depths are 0.17663816 m (S3) and 0.41215571 m (S4), and mapped runoff volumes are 336,189,960 m³ (S3) and 784,443,240 m³ (S4). The exported GeoTIFF readback closes to within 5.96 × 10⁻⁸ m³ in absolute volume.

## Model limitations

The dynamic-routing engine retains the manuscript's harmonic face-averaging formulation;The screening workflow is independent of the dynamic-routing engine and does not alter its governing equations. The four prognostic water stores use conceptual vadose, groundwater and channel process representations. The routine `run_imo_scenarios.py` does not simulate lateral flow, travel time, or inundation dynamics.

The screening regression tests assess the screening implementation and conservation of prescribed runoff volume. They do **not** constitute physical validation or demonstrate predictive flood accuracy.

## Citation and releases

Use the versioned GitHub release and its archival DOI (when available) to identify the exact software version used in an analysis. Numerical results should be associated with the code revision and input datasets that generated them.

## License

MIT License; retain `LICENSE` from the main repository.
