"""
Hexalink v1.0 - Imo River Basin Multi-Scenario Simulation
Manuscript Reference: Section 6.2, Table 4
Parameters: infiltration=120 mm/hr, et=5 mm/hr, low_mult=2.0, channel_mult=4.0
Output: imo_flood_{moderate|heavy|extreme|catastrophic}.tif
"""
import numpy as np
import rasterio
from rasterio.features import rasterize
from scipy.ndimage import mean as ndimage_mean
import os, glob, warnings
from shapely.geometry import Polygon, box

warnings.filterwarnings('ignore', category=RuntimeWarning, module='scipy')

def create_hexagonal_grid(bounds, hex_radius=150):
    minx, miny, maxx, maxy = bounds
    horiz_spacing = 1.5 * hex_radius
    vert_spacing = np.sqrt(3) * hex_radius
    hexagons = []
    cols = int(np.ceil((maxx - minx) / horiz_spacing)) + 1
    rows = int(np.ceil((maxy - miny) / vert_spacing)) + 1
    for row in range(rows):
        for col in range(cols):
            x_offset = (row % 2) * horiz_spacing / 2
            x = minx + col * horiz_spacing + x_offset - hex_radius
            y = miny + row * vert_spacing - hex_radius
            angles = np.linspace(0, 2*np.pi, 7)[:-1] + np.pi/6
            vertices = [(x + hex_radius * np.cos(a), y + hex_radius * np.sin(a)) for a in angles]
            hexagon = Polygon(vertices)
            if hexagon.intersects(box(minx, miny, maxx, maxy)):
                hexagons.append(hexagon)
    return hexagons, hex_radius

def run_imo_scenarios(dem_path=None, hex_radius=150, duration_min=60,
                      infiltration_mmhr=120, et_mmhr=5,
                      low_pct=50, channel_pct=10, low_mult=2.0, channel_mult=4.0):
    
    print("Hexalink: Imo River Basin Multi-Scenario Simulation (Tropical Humid)")
    print("=" * 70)

    # 1. Locate DEM
    if dem_path is None:
        tifs = glob.glob("*.tif")
        dem_path = next((f for f in tifs if "elevation" in f.lower() and "bani" not in f.lower() and "reprojected" in f.lower()), None)
    if not dem_path or not os.path.exists(dem_path):
        raise FileNotFoundError("Imo DEM not found.")
    print(f"Using DEM: {dem_path}")

    # 2. Load DEM
    with rasterio.open(dem_path) as src:
        dem = src.read(1).astype(np.float32)
        transform = src.transform
        crs = src.crs
        rows, cols = src.shape
        bounds = src.bounds

    dem[dem == src.nodata] = np.nan if src.nodata is not None else dem
    valid_mask = ~np.isnan(dem)
    print(f"DEM: {rows}x{cols} | Valid pixels: {np.sum(valid_mask):,}")

    # 3. Build Hex Grid & Aggregate (Run ONCE)
    print(f"\nBuilding hexagonal grid (radius={hex_radius}m)...")
    hexagons, R = create_hexagonal_grid(bounds, hex_radius=hex_radius)
    n_hex = len(hexagons)
    print(f"Generated {n_hex:,} hexagons")

    print("Aggregating DEM to hexagons...")
    hex_shapes = [(p, i+1) for i, p in enumerate(hexagons)]
    hex_ids = rasterize(hex_shapes, out_shape=(rows, cols), transform=transform, fill=0, dtype=np.int32)

    ids_flat = hex_ids.ravel()
    dem_flat = dem.ravel()
    valid_flat = valid_mask.ravel()
    mask = (ids_flat > 0) & valid_flat
    valid_ids = ids_flat[mask]
    valid_dems = dem_flat[mask]
    
    hex_sum = np.zeros(n_hex + 1, dtype=np.float64)
    hex_count = np.zeros(n_hex + 1, dtype=np.float64)
    np.add.at(hex_sum, valid_ids, valid_dems)
    np.add.at(hex_count, valid_ids, 1)

    hex_elev = np.full(n_hex, np.nan, dtype=np.float32)
    nonzero = hex_count[1:] > 0
    hex_elev[nonzero] = (hex_sum[1:][nonzero] / hex_count[1:][nonzero]).astype(np.float32)
    valid_hex = ~np.isnan(hex_elev)
    print(f"Valid hexagons: {np.sum(valid_hex):,}")

    # 4. Run Scenarios
    scenarios = {"moderate": 50, "heavy": 100, "extreme": 200, "catastrophic": 300}
    for name, rain_rate in scenarios.items():
        print(f"\n{'─'*40} Processing {name.upper()} ({rain_rate} mm/hr) {'─'*40}")
        
        total_rain = rain_rate * (duration_min / 60.0)
        total_loss = (infiltration_mmhr + et_mmhr) * (duration_min / 60.0)
        runoff_mm = max(0.0, total_rain - total_loss)
        base_depth = runoff_mm / 1000.0

        elev_valid = hex_elev[valid_hex]
        p_low = np.percentile(elev_valid, low_pct)
        p_chan = np.percentile(elev_valid, channel_pct)

        depths = np.full(n_hex, np.nan, dtype=np.float32)
        depths[valid_hex] = base_depth
        
        depths[valid_hex & (hex_elev < p_low)] = base_depth * low_mult
        depths[valid_hex & (hex_elev < p_chan)] = base_depth * channel_mult
        depths = np.clip(depths, 0.01, 1.0)

        # Rasterize
        shapes = [(p, d) for p, d in zip(hexagons, depths) if not np.isnan(d)]
        out_raster = rasterize(shapes, out_shape=(rows, cols), transform=transform, fill=np.nan, dtype=np.float32)

        out_file = f"imo_flood_{name}.tif"
        with rasterio.open(out_file, 'w', driver='GTiff', height=rows, width=cols,
                           count=1, dtype=np.float32, crs=crs, transform=transform, nodata=np.nan) as dst:
            dst.write(out_raster, 1)

        print(f"Max: {np.nanmax(depths):.3f} m | Mean: {np.nanmean(depths):.3f} m | Saved: {out_file}")

if __name__ == "__main__":
    run_imo_scenarios()