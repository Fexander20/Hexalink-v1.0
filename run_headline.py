"""
run_headline.py - Full Imo Basin dynamic routing demonstration (Hexalink manuscript
Section 7.3, Table 6 - the paper's central result)

Runs the full dynamic routing engine on the Imo DEM at R = 150 m under the standardized
storm (300 mm/hr, 60 min, dt = 10 s, Manning's n = 0.05, Ks = 50 mm/hr, ET_pot = 5 mm/hr)
and reports the depth statistics, mass balance, and runtime reported in Table 6. Included
for direct reproducibility of the paper's headline demonstration.

Usage:
    python3 run_headline.py

Requires reprojected_elevation_2.tif (the Imo DEM) in the working directory.
"""
import sys, time, json
import numpy as np
from routing_engine import HexalinkRoutingEngine

def run_headline():
    params = {
        'n_manning': 0.05,
        'Ks': 50.0,           # mm/hr
        'ET_pot': 5.0,        # mm/hr (matches manuscript §7.2: ET_pot = 5 mm/hr, not 5mm/day here)
        'theta_s': 0.40,
        'theta_r': 0.08,
        'vadose_depth': 0.5,
        'alpha_gw': 0.005,
    }

    t0 = time.perf_counter()
    eng = HexalinkRoutingEngine('reprojected_elevation_2.tif', params, hex_radius=150)

    depths, stats = eng.run(precip_rate_mmhr=300.0, duration_min=60.0, dt_sec=10.0,
                             output_path='depth_R150_headline.tif', print_interval=15)
    total_time = time.perf_counter() - t0

    valid = depths[~np.isnan(depths)]
    positive = valid[valid > 0]
    result = {
        'R': 150,
        'N_cells': int(stats['N_cells']),
        'N_edges': int(stats['N_edges']),
        'outlet_elevation_m': float(eng.z[eng.outlet]),
        'max_depth_m': float(np.max(valid)),
        'mean_positive_depth_m': float(np.mean(positive)),
        'p95_depth_m': float(np.percentile(positive, 95)),
        'p99_depth_m': float(np.percentile(positive, 99)),
        'M_in': stats['M_in'],
        'M_out': stats['M_out'],
        'M_ET': stats['M_ET'],
        'M_deficit': stats['M_deficit'],
        'rel_error_pct': stats['rel_error_pct'],
        'runtime_s': stats['runtime_s'],
        'total_wallclock_s': total_time,
    }
    with open('result_R150_headline.json', 'w') as f:
        json.dump(result, f, indent=2)
    print("\n\nRESULT:", json.dumps(result, indent=2))

if __name__ == '__main__':
    run_headline()
