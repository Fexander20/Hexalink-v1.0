# Hexalink v1.0

A mass-conservative hexagonal finite-volume framework for hydrologic routing and rapid flood
susceptibility screening, demonstrated on the Imo River Basin, Nigeria.

This repository accompanies the manuscript *"Hexalink: A Mass-Conservative Hexagonal
Finite-Volume Framework for Hydrologic Routing and Rapid Flood Screening"* (submitted to
*Advances in Water Resources*).

## Contents

| File | Purpose | Manuscript reference |
|---|---|---|
| `routing_engine.py` | The core dynamic routing engine: true hexagonal six-neighbor adjacency, antisymmetric edge fluxes, aggregate donor-based positivity limiter, and explicit multi-store (surface/vadose/groundwater/channel) mass accounting. Includes four bundled self-tests (A–D). | Sections 3–6 |
| `run_imo_scenarios.py` | The rapid screening (topographic accumulation) workflow: four standardized storm scenarios (Moderate/Heavy/Extreme/Catastrophic) over the Imo Basin. | Section 7.4, Table 7 |
| `extract_metrics.py` | Extracts summary statistics (area, depth percentiles, coverage) from output flood-depth rasters. | Section 7 |
| `run_resolution.py` | Reproduces the spatial-resolution sensitivity experiment (75 m / 100 m / 300 m grids). | Section 6.7, Table 4 |
| `run_sensitivity.py` | Reproduces the one-at-a-time parameter sensitivity screening (Manning's n, K_s, rainfall intensity; 8 scenarios). | Section 7.7, Table 9 |
| `run_headline.py` | Reproduces the main Imo Basin demonstration (R = 150 m, 77,434 cells). | Section 7.3, Table 6 |
| `environment.yml` | Conda environment specification (Python dependencies). | — |

## Requirements

- Python ≥ 3.9
- See `environment.yml` for the full dependency list (numpy, rasterio, scipy, shapely, geopandas)
- The Imo River Basin digital elevation model, `reprojected_elevation_2.tif` (30 m SRTM,
  reprojected to EPSG:32632), is included in this repository. It is derived from the public
  SRTM dataset distributed by NASA/USGS (available via USGS EarthExplorer,
  https://earthexplorer.usgs.gov); the copy included here is the exact reprojected and
  clipped raster used to produce every result in the manuscript.

## Setup

```bash
conda env create -f environment.yml
conda activate hexalink
```

## Usage

Run the bundled self-tests (no DEM required):

```bash
python3 routing_engine.py
```

Reproduce a specific manuscript result (DEM required, see Requirements above):

```bash
python3 run_headline.py                          # Table 6
python3 run_resolution.py 75  result_R75.json     # Table 4 / Section 6.7 (repeat for 100, 300)
python3 run_sensitivity.py REF result_REF.json    # Table 9 (repeat for N003, N008, KS25, KS75, P050, P100, P150)
python3 run_imo_scenarios.py                      # Table 7
```

## Known limitations

See Section 8.5 of the manuscript for a full discussion. Briefly: the vadose zone is a
single-bucket conceptual store (not a Richards-equation solver), groundwater is a linear
reservoir, and channel storage does not route flow between neighboring cells (only the
designated outlet cell discharges from the domain). The bundled test suite includes one
integration smoke test (Test D) that exercises the real routing engine class directly, but
Tests A–C are standalone reimplementations of the reference math and do not themselves
exercise the production class — see the docstring in `routing_engine.py` for details.

## License

[To be added by the authors — see the step-by-step notes provided alongside this file.]

## Citation

If you use this code, please cite the accompanying manuscript. Citation details will be
added here once the paper is published (or as a preprint, if posted first).
