"""
run_sensitivity.py - One-at-a-time parameter sensitivity screening (Hexalink manuscript
Section 7.7, Table 9)

Runs the full dynamic routing engine on the Imo DEM at R = 100 m under a reference parameter
set (Manning's n = 0.05, Ks = 50 mm/hr, rainfall = 300 mm/hr, ET_pot = 5 mm/day) and seven
one-at-a-time perturbations (Manning's n = 0.03/0.08, Ks = 25/75 mm/hr, rainfall =
50/100/150 mm/hr), reporting the depth and outlet-volume metrics used to compute the
percentage changes in Table 9. All eight scenarios must be run to reproduce the full table;
each is independent and can be run in any order.

Usage:
    python3 run_sensitivity.py <scenario_id> <output_json_path>
    Valid scenario_id values: REF, N003, N008, KS25, KS75, P050, P100, P150
    e.g. python3 run_sensitivity.py REF  result_REF.json

Requires reprojected_elevation_2.tif (the Imo DEM) in the working directory.
"""
import sys, time, json
import numpy as np
from routing_engine import HexalinkRoutingEngine

SCENARIOS = {
    'REF':  dict(n_manning=0.05, Ks=50.0, rainfall=300.0),
    'N003': dict(n_manning=0.03, Ks=50.0, rainfall=300.0),
    'N008': dict(n_manning=0.08, Ks=50.0, rainfall=300.0),
    'KS25': dict(n_manning=0.05, Ks=25.0, rainfall=300.0),
    'KS75': dict(n_manning=0.05, Ks=75.0, rainfall=300.0),
    'P050': dict(n_manning=0.05, Ks=50.0, rainfall=50.0),
    'P100': dict(n_manning=0.05, Ks=50.0, rainfall=100.0),
    'P150': dict(n_manning=0.05, Ks=50.0, rainfall=150.0),
}

def run_scenario(scenario_id, out_json):
    s = SCENARIOS[scenario_id]
    params = {
        'n_manning': s['n_manning'],
        'Ks': s['Ks'],
        'ET_pot': 5.0/24.0,   # 5 mm/day -> mm/hr
        'theta_s': 0.40,
        'theta_r': 0.08,
        'vadose_depth': 0.5,
        'alpha_gw': 0.005,
    }

    t0 = time.perf_counter()
    eng = HexalinkRoutingEngine('reprojected_elevation_2.tif', params, hex_radius=100)
    depths, stats = eng.run(precip_rate_mmhr=s['rainfall'], duration_min=60.0, dt_sec=2.5,
                             output_path=f'depth_{scenario_id}.tif', print_interval=30)
    total_time = time.perf_counter() - t0

    valid = depths[~np.isnan(depths)]
    positive = valid[valid > 0]
    if len(positive) == 0:
        max_d = mean_d = p95_d = p99_d = 0.0
    else:
        max_d = float(np.max(positive))
        mean_d = float(np.mean(positive))
        p95_d = float(np.percentile(positive, 95))
        p99_d = float(np.percentile(positive, 99))

    result = {
        'scenario_id': scenario_id,
        'n_manning': s['n_manning'], 'Ks': s['Ks'], 'rainfall_mmhr': s['rainfall'],
        'N_cells': int(stats['N_cells']), 'N_edges': int(stats['N_edges']),
        'max_depth_m': max_d, 'mean_depth_m': mean_d, 'p95_depth_m': p95_d, 'p99_depth_m': p99_d,
        'M_in': stats['M_in'], 'M_out': stats['M_out'], 'M_ET': stats['M_ET'],
        'rel_error_pct': stats['rel_error_pct'], 'runtime_s': stats['runtime_s'],
        'total_wallclock_s': total_time,
    }
    with open(out_json, 'w') as f:
        json.dump(result, f, indent=2)
    print("\nRESULT:", json.dumps(result, indent=2))

if __name__ == '__main__':
    run_scenario(sys.argv[1], sys.argv[2])
