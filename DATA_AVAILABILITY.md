# Data availability and third-party rights

The public release contains source code, independently generated numerical test inputs/results, and selected aggregate summaries. **No third-party terrain, satellite observation raster, or terrain-derived hex mesh is distributed.** The MIT licence covers original software only and does not replace the conditions on external datasets.

| Input | Acquisition / requirements | Redistribution in this release |
| --- | --- | --- |
| MERIT DEM elevation source for St. John | [MERIT DEM official page](https://global-hydrodynamics.github.io/MERIT_DEM/); follow its registration and licensing requirements | Not included |
| St. John 100 m DEM, EPSG:32629 | Locally derived from appropriately licensed elevation data; expected screening filename `StJohn_MtFinlay_Catchment_Elevation_100m_UTM29N.tif` | Not included |
| St. John hex mesh, circumradius 300 m | Locally supplied `hexalink_stjohn_config_300m.npz`; arrays include `x_m`, `y_m`, `elevation_merit_m`, `edge_i`, `edge_j`, `active_area_m2` and `cell_id` as required by scripts | Not included |
| GFD DFO_3982 flood observations (2–7 Sep 2012) | [Google Earth Engine Global Flood Database](https://developers.google.com/earth-engine/datasets/catalog/GLOBAL_FLOOD_DB_MODIS_EVENTS_V1); see scripts in `05_GFD_PROVENANCE/` | Not included (CC BY-NC 4.0 source) |
| Screening GeoTIFFs | Generate locally with licensed DEM and compatible mesh | Not included |

**Reproducibility boundary.** The numerical verification programs run without third-party geographic inputs. The St. John routing and screening scripts require the specific locally supplied DEM-derived mesh. An independently audited mesh-construction program and downloadable licensed mesh are not part of this release; therefore the geospatial case cannot be fully reproduced from a DEM download using this repository alone. Associated summary metrics and the observation evaluation code are provided for transparency.

**Rights.** MERIT DEM uses source-specific CC BY-NC 4.0 / ODbL terms and redistribution provisions, including restrictions on distributing the full source dataset. GFD is licensed CC BY-NC 4.0. Commercial use of the MIT-licensed Hexalink code does not automatically authorize commercial use of either dataset. Obtain inputs at their source and comply with applicable attribution and licence requirements.

**Sources and attribution:**
- Yamazaki et al., MERIT DEM: [official product and licensing information](https://global-hydrodynamics.github.io/MERIT_DEM/).
- Tellman et al. (2021), Global Flood Database, Google Earth Engine collection `GLOBAL_FLOOD_DB/MODIS_EVENTS/V1`: [catalog, description and terms](https://developers.google.com/earth-engine/datasets/catalog/GLOBAL_FLOOD_DB_MODIS_EVENTS_V1).

No claim is made here that third-party data rights have been cleared for public redistribution.
