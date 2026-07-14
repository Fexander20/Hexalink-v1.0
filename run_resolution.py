"""
run_resolution.py - Spatial-resolution sensitivity experiment (Hexalink manuscript Section 6.7, Table 4)

Runs the full dynamic routing engine on the Imo DEM at three hexagon radii (75 m, 100 m, 300 m)
under identical forcing (300 mm/hr, 60 min, dt = 2.5 s) and reports the depth/discharge metrics
used to compute the percentage differences reported in Table 4's "Spatial sensitivity" row and
discussed in Section 6.7. The 75 m run is treated as the fine reference; 100 m and 300 m are
compared against it.

Usage:
    python3 run_resolution.py <R_meters> <output_json_path>
    e.g. python3 run_resolution.py 75  result_R75.json
         python3 run_resolution.py 100 result_R100.json
         python3 run_resolution.py 300 result_R300.json

Requires reprojected_elevation_2.tif (the Imo DEM) in the working directory.
"""
import sys, time, json
import numpy as np
from routing_engine import HexalinkRoutingEngine

def run_at_resolution(R, out_json):
    params = {
        'n_manning': 0.05,
        'Ks': 50.0,           # mm/hr
        'ET_pot': 5.0/24.0,   # mm/hr  (5 mm/day)
        'theta_s': 0.40,
        'theta_r': 0.08,
        'vadose_depth': 0.5,
        'alpha_gw': 0.005,
    }

    t0 = time.perf_counter()
    eng = HexalinkRoutingEngine('reprojected_elevation_2.tif', params, hex_radius=R)

    depths, stats = eng.run(precip_rate_mmhr=300.0, duration_min=60.0, dt_sec=2.5,
                             output_path=f'depth_R{R}.tif', print_interval=15)
    total_time = time.perf_counter() - t0

    valid = depths[~np.isnan(depths)]
    result = {
        'R': R,
        'N_cells': int(stats['N_cells']),
        'N_edges': int(stats['N_edges']),
        'max_depth_m': float(np.max(valid)),
        'mean_depth_m': float(np.mean(valid)),
        'p90_depth_m': float(np.percentile(valid, 90)),
        'p95_depth_m': float(np.percentile(valid, 95)),
        'p99_depth_m': float(np.percentile(valid, 99)),
        'M_in': stats['M_in'],
        'M_out': stats['M_out'],
        'M_ET': stats['M_ET'],
        'rel_error_pct': stats['rel_error_pct'],
        'runtime_s': stats['runtime_s'],
        'total_wallclock_s': total_time,
    }
    with open(out_json, 'w') as f:
        json.dump(result, f, indent=2)
    print("\n\nRESULT:", json.dumps(result, indent=2))

if __name__ == '__main__':
    R = int(sys.argv[1])
    out_json = sys.argv[2]
    run_at_resolution(R, out_json)
