"""Area-conservative terrain ranking for rapid flood-priority screening.
Output is algebraic screening-equivalent depth, not hydraulic inundation.
"""
import numpy as np
SCENARIOS = {"moderate":50.0,"heavy":100.0,"extreme":200.0,"catastrophic":300.0}
def screening_depth_by_hex(hex_elev, counts, runoff_m, low_pct=50, channel_pct=10, low_mult=2.0, channel_mult=4.0):
    """Return normalized hex depths (m), with invalid hexes set to NaN.

    Pixel counts are proportional to valid DEM area for a constant-area raster.
    The arithmetic of volume normalization is performed in float64.
    """
    if runoff_m < 0 or not np.isfinite(runoff_m):
        raise ValueError('runoff_m must be finite and nonnegative')
    if low_mult <= 0 or channel_mult <= 0:
        raise ValueError('Terrain multipliers must be positive')
    counts = np.asarray(counts, dtype=np.float64)
    hex_elev = np.asarray(hex_elev, dtype=np.float64)
    valid = (counts > 0) & np.isfinite(hex_elev)
    if not np.any(valid):
        raise ValueError('No valid DEM-supported hexagons')
    depths = np.full(hex_elev.shape, np.nan, dtype=np.float64)
    if runoff_m == 0.0:
        depths[valid] = 0.0
        return (depths, 1.0)
    low_cutoff = float(np.percentile(hex_elev[valid], low_pct))
    channel_cutoff = float(np.percentile(hex_elev[valid], channel_pct))
    weights = np.ones(hex_elev.shape, dtype=np.float64)
    weights[valid & (hex_elev < low_cutoff)] = low_mult
    weights[valid & (hex_elev < channel_cutoff)] = channel_mult
    weighted_mean = float(np.dot(counts[valid], weights[valid]) / counts[valid].sum())
    if not np.isfinite(weighted_mean) or weighted_mean <= 0:
        raise ValueError('Invalid area-weighted terrain normalization')
    depths[valid] = runoff_m * weights[valid] / weighted_mean
    return (depths, weighted_mean)
