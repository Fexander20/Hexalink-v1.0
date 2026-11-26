# Hexalink v2.0.0

**Hexagonal finite-volume surface routing and terrain-based flood screening**

Hexalink is a Python research implementation with two distinct components:

1. **Surface routing:** explicit mass-conservative transfer among regular hexagonal finite volumes under a nonlinear, hydraulic-head-driven flux formulation.
2. **Rapid screening:** an independent algebraic terrain-priority calculation that distributes a specified effective-runoff volume. Screening-equivalent depth is **not** simulated inundation depth.

The release contains the authoritative routing engine, numerical verification programs, St. John basin demonstration scripts, flood-screening scripts, Global Flood Database (GFD) assessment scripts, and selected numerical results. It does **not** claim calibrated historical flood prediction or operational forecasting skill.

## Repository structure

| Path | Purpose |
| --- | --- |
| `routing_engine.py` | Canonical surface-routing implementation |
| `02_NUMERICAL_TESTS/` | Analytical nonlinear diffusion, independent Cartesian reference, timestep and equilibrium checks |
| `03_STJOHN_ROUTING/` | Controlled synthetic routing experiment over a 300 m hex mesh in Liberia |
| `04_SCREENING/` | Separate volume-accounting terrain-priority screening |
| `05_GFD_PROVENANCE/` | Google Earth Engine export and event-based spatial-ranking evaluation |
| `06_RESULTS/` | Numerical results and aggregate statistics |
| `DATA_AVAILABILITY.md` | Third-party source information, exclusions, and reproducibility limits |
| `CITATION.cff` | Software citation metadata |

## Setup

Recommended: Python 3.11 or later and an isolated virtual environment. The reference environment used Python 3.13.5. Install the dependencies:

```bash
python -m pip install -r requirements.txt
```

### Numerical verification

From the repository root:

```bash
python 02_NUMERICAL_TESTS/verify_hexalink_numerics.py
python 02_NUMERICAL_TESTS/rebuild_independent_cartesian_reference.py
```

The first command evaluates the analytical nonlinear-diffusion benchmark, sloping-bed comparison, timestep sensitivity, and lake-at-rest condition using the **same root routing engine**. It writes numerical results and figures in `02_NUMERICAL_TESTS/`. The second reconstructs the independent 2 m Cartesian comparison field and compares it with the distributed reference. Neither command needs a basin DEM.

### St. John synthetic surface routing

```bash
python 03_STJOHN_ROUTING/run_stjohn_surface_routing.py
```

Requires a locally supplied `03_STJOHN_ROUTING/hexalink_stjohn_config_300m.npz`. The DEM-derived mesh is not distributed; see `DATA_AVAILABILITY.md`. The script checks its expected geometry before running. As the mesh-construction procedure is not included, the public repository **does not enable independent end-to-end reproduction** of this case from the original DEM alone.

### Terrain-based screening

```bash
python 04_SCREENING/run_stjohn_screening.py --dem PATH_TO_DEM --mesh PATH_TO_MESH --out PATH_TO_OUTPUT
```

The supplied DEM must be a compatible 100 m GeoTIFF in EPSG:32629. The local mesh must contain the attributes expected by the script. `python 04_SCREENING/run_stjohn_screening.py --help` displays the options.

### GFD assessment

The scripts in `05_GFD_PROVENANCE/` document event DFO_3982 (2–7 September 2012). Export legally accessed observation data through Google Earth Engine, then use `python 05_GFD_PROVENANCE/evaluate_gfd_observation.py --help` to supply the appropriate local inputs. The evaluation is a **spatial ranking comparison** for the independent screening pathway, not validation of dynamic routing against observations.

## Reported results

| Check | Result |
| --- | --- |
| Analytical nonlinear diffusion normalized L1 error (R = 40, 20, 10 m) | 2.5618%, 0.9780%, 0.2229% |
| Sloping-bed depth difference versus independent 2 m Cartesian reference | 6.1001% at specified configuration |
| Lake-at-rest maximum depth change | 0.0 m |
| St. John surface-routing example | 68,292 full hexagons, 203,076 internal edges |
| St. John synthetic rain pulse | 300 mm/h for 30 min over 1,945 cells; 60 min modeled |
| St. John recorded mass-closure fraction | 4.587 × 10⁻¹⁵ |
| Screening versus GFD DFO_3982 | ROC-AUC 0.7038; 70 flood-positive of 216,906 valid pixels |

Numerical checks are **verification**, not evidence of predictive performance. The event-based screening assessment has very few mapped flood positives and cannot establish generalizable accuracy. The synthetic pulse is not a historical flood reconstruction.

The St. John DEM catchment mask area is approximately **15,744.41 km²**; its 68,292 uncut full-hexagon computational domain covers approximately **15,968.50 km²**. These areas are not interchangeable. A donor-positivity safeguard is included in the routing algorithm; it does not imply unconditional stability.

## Data, rights and citation

Only original software and selected shareable tabular/aggregate results are distributed. The repository deliberately excludes third-party DEMs, the DEM-derived mesh, GFD raster imagery and derived map GeoTIFFs. Consult `DATA_AVAILABILITY.md` to acquire inputs independently and respect source-specific terms. MIT applies to the original **software**, not to third-party data.

Software author and copyright holder: **Ifeanyi Chidozie Oli**. See `CITATION.cff`. Cite the exact archived release DOI shown on Zenodo after publication. Manuscript authorship is a separate matter.

**Research use notice:** this code has been numerically checked in the reported cases, but is not a validated operational flood forecasting system. No warranty is provided; see `LICENSE`.
