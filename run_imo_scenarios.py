"""Hexalink v1.0 — area-conservative rapid screening for the Imo River Basin.

This script implements the separate *algebraic* screening pathway described in
Section 5.6 and Section 7.4 of the Environmental Modelling & Software manuscript
ENVSOFT-D-26-02733. It does not calculate routed hydraulic inundation depths.

Four one-hour scenarios use prescribed intensities of 50, 100, 200 and 300 mm/h;
infiltration capacity 120 mm/h and illustrative ET subtraction 5 mm/h.

The valid-DEM-pixel-area-weighted mean of the 1/2/4 terrain weights is used to
normalise mapped depths. There is no positive display-depth floor or clipping.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import rasterio
from rasterio.features import rasterize
from shapely.geometry import Polygon, box


SCENARIOS = {"moderate": 50.0, "heavy": 100.0, "extreme": 200.0, "catastrophic": 300.0}


def create_hexagonal_grid(bounds, hex_radius=150):
    """Flat-top regular hexagons using the same geometric lattice as routing.

    Columns are 1.5R apart and alternate columns have a vertical offset of
    sqrt(3)R/2. Six vertices are at 0°, 60°, ... , 300° about the centre.
    Include exterior padding so boundary pixels have a containing hexagon.
    """
    if hex_radius <= 0:
        raise ValueError("hex_radius must be positive")
    minx, miny, maxx, maxy = map(float, bounds)
    if maxx <= minx or maxy <= miny:
        raise ValueError("Invalid DEM bounds")
    R = float(hex_radius)
    horiz = 1.5 * R
    vert = np.sqrt(3.0) * R
    a = np.arange(6, dtype=np.float64) * np.pi / 3.0
    dx = R * np.cos(a)
    dy = R * np.sin(a)
    domain = box(minx, miny, maxx, maxy)
    polys = []
    ncols = int(np.ceil((maxx - minx) / horiz)) + 4
    nrows = int(np.ceil((maxy - miny) / vert)) + 4
    for col in range(-2, ncols):
        x = minx + col * horiz
        for row in range(-2, nrows):
            y = miny + (row + (col % 2) / 2.0) * vert
            polygon = Polygon(zip(x + dx, y + dy))
            if polygon.intersects(domain):
                polys.append(polygon)
    return polys, R


def screening_depth_by_hex(hex_elev, counts, runoff_m,
                           low_pct=50, channel_pct=10,
                           low_mult=2.0, channel_mult=4.0):
    """Return normalized hex depths (m), with invalid hexes set to NaN.

    Pixel counts are proportional to valid DEM area for a constant-area raster.
    The arithmetic of volume normalization is performed in float64.
    """
    if runoff_m < 0 or not np.isfinite(runoff_m):
        raise ValueError("runoff_m must be finite and nonnegative")
    if low_mult <= 0 or channel_mult <= 0:
        raise ValueError("Terrain multipliers must be positive")
    counts = np.asarray(counts, dtype=np.float64)
    hex_elev = np.asarray(hex_elev, dtype=np.float64)
    valid = (counts > 0) & np.isfinite(hex_elev)
    if not np.any(valid):
        raise ValueError("No valid DEM-supported hexagons")
    depths = np.full(hex_elev.shape, np.nan, dtype=np.float64)
    if runoff_m == 0.0:
        depths[valid] = 0.0
        return depths, 1.0

    low_cutoff = float(np.percentile(hex_elev[valid], low_pct))
    channel_cutoff = float(np.percentile(hex_elev[valid], channel_pct))
    weights = np.ones(hex_elev.shape, dtype=np.float64)
    weights[valid & (hex_elev < low_cutoff)] = low_mult
    weights[valid & (hex_elev < channel_cutoff)] = channel_mult
    weighted_mean = float(np.dot(counts[valid], weights[valid]) / counts[valid].sum())
    if not np.isfinite(weighted_mean) or weighted_mean <= 0:
        raise ValueError("Invalid area-weighted terrain normalization")
    depths[valid] = runoff_m * weights[valid] / weighted_mean
    return depths, weighted_mean


def run_imo_scenarios(dem_path=None, hex_radius=150, duration_min=60,
                      infiltration_mmhr=120, et_mmhr=5,
                      low_pct=50, channel_pct=10, low_mult=2.0,
                      channel_mult=4.0, output_dir="."):
    """Generate four volume-conserving screening rasters and a diagnostic CSV.

    The output GeoTIFF names retain their legacy names for downstream
    compatibility; the rasters are equivalent *screening* depths, not floods.
    """
    if duration_min <= 0 or infiltration_mmhr < 0 or et_mmhr < 0:
        raise ValueError("Duration must be positive; loss rates nonnegative")
    if dem_path is None:
        candidates = sorted(Path(".").glob("*reprojected*elevation*.tif"))
        if not candidates:
            candidates = sorted(Path(".").glob("*elevation*reprojected*.tif"))
        if not candidates:
            raise FileNotFoundError("Please supply the Imo DEM using --dem")
        dem_path = candidates[0]
    dem_path = Path(dem_path)
    if not dem_path.is_file():
        raise FileNotFoundError(f"DEM not found: {dem_path}")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with rasterio.open(dem_path) as src:
        dem = np.asarray(src.read(1, masked=True).astype(np.float64).filled(np.nan), dtype=np.float64)
        transform = src.transform
        crs = src.crs
        bounds = src.bounds
        original_profile = src.profile.copy()
    valid_mask = np.isfinite(dem)
    valid_pixel_count = int(valid_mask.sum())
    if valid_pixel_count == 0:
        raise ValueError("DEM contains no finite valid pixels")
    if crs is None or not crs.is_projected or "met" not in crs.linear_units.lower():
        raise ValueError("Use a meter-based projected DEM; manuscript uses EPSG:32632")
    cell_area = abs(transform.a * transform.e - transform.b * transform.d)
    if cell_area <= 0:
        raise ValueError("DEM pixel area is not positive")

    hexagons, _ = create_hexagonal_grid(bounds, hex_radius=hex_radius)
    if len(hexagons) >= np.iinfo(np.int32).max:
        raise ValueError("Too many hexagons for raster labels")
    hex_ids = rasterize(
        ((p, i + 1) for i, p in enumerate(hexagons)),
        out_shape=dem.shape, transform=transform, fill=0, dtype=np.int32,
    )
    # The published screening script left some valid raster pixels unassigned.
    # Reject that case rather than silently undercount the DEM-supported area.
    uncovered = int(np.count_nonzero(valid_mask & (hex_ids == 0)))
    if uncovered:
        raise ValueError(f"Screening grid leaves {uncovered} valid DEM pixels uncovered")

    valid_ids = hex_ids[valid_mask]
    counts = np.bincount(valid_ids, minlength=len(hexagons) + 1).astype(np.float64)
    sums = np.bincount(valid_ids, weights=dem[valid_mask], minlength=len(hexagons) + 1)
    elev = np.full(len(hexagons) + 1, np.nan, dtype=np.float64)
    keep = counts > 0
    elev[keep] = sums[keep] / counts[keep]
    results = []

    for name, rain_mmhr in SCENARIOS.items():
        runoff_mm = max(0.0, (rain_mmhr - infiltration_mmhr - et_mmhr) * duration_min / 60.0)
        base_depth_m = runoff_mm / 1000.0
        depths, normalization = screening_depth_by_hex(
            elev, counts, base_depth_m, low_pct=low_pct, channel_pct=channel_pct,
            low_mult=low_mult, channel_mult=channel_mult,
        )
        # Use one shared pixel/hex assignment for DEM aggregation, weighting,
        # output and volume integration. No new rasterization of polygons.
        image = np.full(dem.shape, np.nan, dtype=np.float64)
        image[valid_mask] = depths[hex_ids[valid_mask]]
        expected_volume = float(base_depth_m * valid_pixel_count * cell_area)
        mapped_volume = float(np.sum(image[valid_mask], dtype=np.float64) * cell_area)
        residual = mapped_volume - expected_volume
        rel_percent = 0.0 if expected_volume == 0 else abs(residual) / expected_volume * 100.0
        if not np.isclose(mapped_volume, expected_volume, rtol=1e-12, atol=1e-7):
            raise AssertionError(f"Volume mismatch in {name}: {residual} m3")

        out_file = output_dir / f"imo_flood_{name}.tif"  # legacy path, screening quantity
        profile = original_profile.copy()
        profile.update(driver="GTiff", count=1, dtype="float64", nodata=np.nan, compress="deflate")
        with rasterio.open(out_file, "w", **profile) as dest:
            dest.write(image, 1)
            dest.set_band_description(1, "Equivalent screening depth (m), not hydraulic flood depth")
            dest.update_tags(
                quantity="algebraic equivalent screening depth",
                hydrodynamic_simulation="false",
                normalization="valid DEM pixel area-weighted 1/2/4 weights",
                journal="Environmental Modelling & Software",
            )
        metric = {
            "scenario": name,
            "rainfall_mmhr": rain_mmhr,
            "runoff_mm": runoff_mm,
            "valid_pixels": valid_pixel_count,
            "valid_area_m2": valid_pixel_count * cell_area,
            "area_weighted_mean_depth_m": float(np.mean(image[valid_mask])),
            "min_depth_m": float(np.min(image[valid_mask])),
            "max_depth_m": float(np.max(image[valid_mask])),
            "expected_volume_m3": expected_volume,
            "mapped_volume_m3": mapped_volume,
            "volume_residual_m3": residual,
            "relative_volume_error_percent": rel_percent,
            "weight_normalization": normalization,
            "output": str(out_file),
        }
        results.append(metric)
        print(f"{name}: runoff={runoff_mm:.1f} mm | mean={metric['area_weighted_mean_depth_m']:.6f} m "
              f"| residual={residual:.3e} m³ | {out_file}")

    with (output_dir / "screening_metrics.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    with (output_dir / "screening_metrics.json").open("w", encoding="utf-8") as stream:
        json.dump(results, stream, indent=2)
    return results


def main():
    parser = argparse.ArgumentParser(description="Hexalink volume-conserving algebraic screening")
    parser.add_argument("--dem", default="reprojected_elevation_2.tif", help="Metric projected DEM")
    parser.add_argument("--output-dir", default=".")
    parser.add_argument("--hex-radius", type=float, default=150.0)
    args = parser.parse_args()
    run_imo_scenarios(dem_path=args.dem, output_dir=args.output_dir, hex_radius=args.hex_radius)


if __name__ == "__main__":
    main()
