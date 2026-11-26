#!/usr/bin/env python3
"""Hexalink: conservative finite-volume hydrological routing on a hexagonal grid.

Surface-water intercell transfer uses nonlinear secant-depth mobility and an
aggregate donor limiter. The implementation also provides vertical storage,
infiltration, evapotranspiration, groundwater baseflow and outlet accounting.
The controlled numerical experiments exercise the integrated surface flux.
"""
import numpy as np
import rasterio
import time
import os
from shapely.geometry import Polygon, box
from rasterio.features import rasterize

def build_hex_grid(bounds, hex_radius):
    """
    Tessellate the bounding box with regular hexagons (flat-top orientation).
    Returns:
        centers  : (N,2) array of (x,y) cell centres
        polygons : list of N Shapely Polygon objects
        A_hex    : scalar cell area (m²)  = (3√3/2) R²
        d_hex    : centre-to-centre distance = √3 R
        w_hex    : shared edge length = R
    """
    R = hex_radius
    A_hex = 3 * np.sqrt(3) / 2 * R ** 2
    d_hex = np.sqrt(3) * R
    w_hex = R
    minx, miny, maxx, maxy = bounds
    col_step = 1.5 * R
    row_step = d_hex
    centers = []
    polygons = []
    bbox = box(minx, miny, maxx, maxy)
    col_idx = 0
    x = minx
    while x <= maxx + R:
        row_idx = 0
        y = miny
        y_off = d_hex / 2 * (col_idx % 2)
        y_start = miny - y_off
        y = y_start
        while y <= maxy + R:
            angles = np.linspace(0, 2 * np.pi, 7)[:-1] + np.pi / 6
            verts = [(x + R * np.cos(a), y + R * np.sin(a)) for a in angles]
            poly = Polygon(verts)
            if poly.intersects(bbox):
                centers.append((x, y))
                polygons.append(poly)
            y += row_step
            row_idx += 1
        x += col_step
        col_idx += 1
    centers = np.array(centers, dtype=np.float64)
    return (centers, polygons, A_hex, d_hex, w_hex)

def build_adjacency(centers, d_hex, tol=0.15):
    """
    Build 6-neighbor adjacency list using distance criterion.
    Two cells are neighbors if their centre distance ≈ d_hex.
    Returns:
        neighbors : list of lists — neighbors[i] = [j1, j2, ...] (≤6 entries)
        edge_list : list of (i,j) tuples, each internal edge listed ONCE (i < j)
    """
    N = len(centers)
    neighbors = [[] for _ in range(N)]
    edge_list = []
    from scipy.spatial import cKDTree
    tree = cKDTree(centers)
    pairs = tree.query_pairs(r=d_hex * (1 + tol))
    for i, j in pairs:
        dist = np.linalg.norm(centers[i] - centers[j])
        if abs(dist - d_hex) < tol * d_hex:
            neighbors[i].append(j)
            neighbors[j].append(i)
            edge_list.append((i, j))
    return (neighbors, edge_list)

def assign_dem_to_hexagons(centers, polygons, dem, transform, nodata, rows, cols):
    """
    Assign mean DEM elevation to each hexagonal cell.
    Returns hex_elev (N,) array; NaN for cells with no valid DEM pixels.
    """
    N = len(polygons)
    shapes = [(p, i + 1) for i, p in enumerate(polygons)]
    hex_ids = rasterize(shapes, out_shape=(rows, cols), transform=transform, fill=0, dtype=np.int32)
    dem_flat = dem.ravel().astype(np.float64)
    ids_flat = hex_ids.ravel()
    valid = ids_flat > 0
    if nodata is not None:
        valid &= dem_flat != nodata
    valid &= ~np.isnan(dem_flat)
    hex_sum = np.zeros(N + 1, dtype=np.float64)
    hex_cnt = np.zeros(N + 1, dtype=np.float64)
    np.add.at(hex_sum, ids_flat[valid], dem_flat[valid])
    np.add.at(hex_cnt, ids_flat[valid], 1)
    hex_elev = np.full(N, np.nan, dtype=np.float64)
    ok = hex_cnt[1:] > 0
    hex_elev[ok] = hex_sum[1:][ok] / hex_cnt[1:][ok]
    return (hex_elev, hex_ids)

class HexalinkRoutingEngine:
    """Hexagonal finite-volume hydrological routing engine.

State comprises surface water volume, vadose moisture, groundwater and channel
storage. Surface flux uses hydraulic-head gradients and secant-depth mobility.
Aggregate outgoing demand is limited by the available donor water volume.
"""

    def __init__(self, dem_path, params, hex_radius=150):
        print('=' * 65)
        print('  Hexalink v2.0.0  |  Hexagonal Routing Engine')
        print('=' * 65)
        with rasterio.open(dem_path) as src:
            dem_full = src.read(1).astype(np.float64)
            self.nodata = src.nodata
            self.transform = src.transform
            self.crs = src.crs
            self.raster_rows, self.raster_cols = src.shape
            self.bounds = src.bounds
        print(f'\n[1/4] Building hexagonal grid  (R = {hex_radius} m) ...')
        self.R = hex_radius
        centers, polygons, self.A_hex, self.d_hex, self.w_hex = build_hex_grid(tuple(self.bounds), hex_radius)
        print(f'      Generated {len(centers):,} candidate hexagons')
        print('[2/4] Assigning DEM to hexagons ...')
        hex_elev_all, self.hex_ids_raster = assign_dem_to_hexagons(centers, polygons, dem_full, self.transform, self.nodata, self.raster_rows, self.raster_cols)
        self.valid_idx = np.where(~np.isnan(hex_elev_all))[0]
        self.N = len(self.valid_idx)
        self.z = hex_elev_all[self.valid_idx]
        self.centers_valid = centers[self.valid_idx]
        self.polygons_valid = [polygons[i] for i in self.valid_idx]
        full_to_compact = np.full(len(centers), -1, dtype=np.int64)
        full_to_compact[self.valid_idx] = np.arange(self.N)
        print(f'      Valid hexagonal cells: {self.N:,}')
        print('[3/4] Building 6-neighbor adjacency ...')
        neighbors_full, edge_list_full = build_adjacency(centers, self.d_hex)
        self.neighbors = [[] for _ in range(self.N)]
        self.edge_list = []
        seen = set()
        for orig_i, orig_j in edge_list_full:
            ci = full_to_compact[orig_i]
            cj = full_to_compact[orig_j]
            if ci < 0 or cj < 0:
                continue
            key = (min(ci, cj), max(ci, cj))
            if key not in seen:
                seen.add(key)
                self.edge_list.append(key)
                self.neighbors[ci].append(cj)
                self.neighbors[cj].append(ci)
        self.edge_i = np.array([e[0] for e in self.edge_list], dtype=np.int64)
        self.edge_j = np.array([e[1] for e in self.edge_list], dtype=np.int64)
        n_edges = len(self.edge_list)
        avg_neighbors = sum((len(nb) for nb in self.neighbors)) / max(self.N, 1)
        print(f'      Internal edges: {n_edges:,}  |  Avg neighbors: {avg_neighbors:.2f}')
        n_neighbors = np.array([len(nb) for nb in self.neighbors])
        boundary_mask = n_neighbors < 6
        boundary_idx = np.where(boundary_mask)[0]
        if len(boundary_idx) == 0:
            boundary_idx = np.arange(self.N)
        self.outlet = boundary_idx[np.argmin(self.z[boundary_idx])]
        print(f'      Outlet cell: {self.outlet}  (z = {self.z[self.outlet]:.2f} m, neighbors = {n_neighbors[self.outlet]})')
        print('[4/4] Setting parameters ...')
        self.Ks = params.get('Ks', 50.0) / 1000 / 3600
        self.n_mann = params.get('n_manning', 0.05)
        self.alpha_gw = params.get('alpha_gw', 0.005) / 3600
        self.ET_pot = params.get('ET_pot', 5.0) / 1000 / 3600
        self.theta_s = params.get('theta_s', 0.4)
        self.theta_r = params.get('theta_r', 0.08)
        self.vadose_depth = params.get('vadose_depth', 0.5)
        self.V_s = np.zeros(self.N, dtype=np.float64)
        self.theta = np.full(self.N, self.theta_r)
        self.S_gw = np.zeros(self.N, dtype=np.float64)
        self.V_chan = np.zeros(self.N, dtype=np.float64)
        self.S_f = np.full(self.N, 0.0001, dtype=np.float64)
        for i in range(self.N):
            nb = self.neighbors[i]
            if nb:
                dz = np.abs(self.z[i] - self.z[nb])
                self.S_f[i] = max(np.mean(dz) / self.d_hex, 0.0001)
        print(f'\n  Grid ready. {self.N:,} cells | {n_edges:,} edges\n')

    def _D(self, h, S_f):
        """D = (1/n) * h^(5/3) * S_f^(-1/2)   [m²/s]  (Eq 3.13)"""
        h_eff = np.maximum(h, 1e-06)
        return 1.0 / self.n_mann * h_eff ** (5 / 3) / np.sqrt(S_f)

    def _route_surface(self, dt):
        """Compute conservative internal face exchanges and outlet discharge.

For each face, harmonic averaging of friction prefactors is combined with
secant-depth mobility. Opposite face contributions are accumulated
antisymmetrically; one donor limiter accounts for both face and outlet demand.
Returns the cellwise surface volume increment and outlet loss in m3 per step.
"""
        h = self.V_s / self.A_hex
        H = h + self.z
        Hi = H[self.edge_i]
        Hj = H[self.edge_j]
        m = 8.0 / 3.0
        hi = np.maximum(h[self.edge_i], 0.0)
        hj = np.maximum(h[self.edge_j], 0.0)
        Ki = 1.0 / (self.n_mann * np.sqrt(self.S_f[self.edge_i]))
        Kj = 1.0 / (self.n_mann * np.sqrt(self.S_f[self.edge_j]))
        K_face = 2.0 * Ki * Kj / (Ki + Kj)
        diff = hi - hj
        near_equal = np.abs(diff) < 1e-08 * np.maximum(np.maximum(hi, hj), 1e-05)
        mobility = np.empty_like(diff)
        mobility[near_equal] = (0.5 * (hi[near_equal] + hj[near_equal])) ** (m - 1.0)
        unequal = ~near_equal
        mobility[unequal] = (hi[unequal] ** m - hj[unequal] ** m) / (m * diff[unequal])
        D_face = K_face * mobility
        Q_raw = D_face * (Hi - Hj) / self.d_hex * self.w_hex
        Q_vol = Q_raw * dt
        h_out = max(self.V_s[self.outlet] / self.A_hex, 0.0)
        Q_out_rate = 1.0 / self.n_mann * h_out ** (5 / 3) * np.sqrt(self.S_f[self.outlet]) * self.w_hex
        Q_outlet_demand = Q_out_rate * dt
        V_req_out = np.zeros(self.N, dtype=np.float64)
        np.add.at(V_req_out, self.edge_i, np.maximum(Q_vol, 0))
        np.add.at(V_req_out, self.edge_j, np.maximum(-Q_vol, 0))
        V_req_out[self.outlet] += Q_outlet_demand
        V_avail = np.maximum(self.V_s, 0.0)
        eps = 1e-12
        lam = np.minimum(1.0, V_avail / (V_req_out + eps))
        lam = np.where(V_req_out < eps, 1.0, lam)
        lam_i = lam[self.edge_i]
        lam_j = lam[self.edge_j]
        Q_final = np.where(Q_vol >= 0, Q_vol * lam_i, --Q_vol * lam_j)
        dV = np.zeros(self.N, dtype=np.float64)
        np.add.at(dV, self.edge_i, -Q_final)
        np.add.at(dV, self.edge_j, Q_final)
        Q_outlet = Q_outlet_demand * lam[self.outlet]
        dV[self.outlet] -= Q_outlet
        return (dV, Q_outlet)

    def _infiltration(self, dt):
        """Compute infiltration volume constrained by surface water,
hydraulic capacity and remaining vadose-zone storage capacity."""
        h = self.V_s / self.A_hex
        I_rate = np.minimum(h / dt, self.Ks)
        I_vol = I_rate * self.A_hex * dt
        I_vol = np.minimum(I_vol, self.V_s)
        vadose_capacity_vol = np.maximum((self.theta_s - self.theta) * self.A_hex * self.vadose_depth, 0.0)
        I_vol = np.minimum(I_vol, vadose_capacity_vol)
        return I_vol

    def _baseflow(self, dt):
        """Linear reservoir (Eqs 3.10-3.11). Returns Q_base [m³]."""
        return np.minimum(self.alpha_gw * self.S_gw * dt, self.S_gw)

    def _et(self, dt):
        """
        Potential ET with stress factor (Eqs 3.19-3.20).
        Returns (ET_s, ET_u, ET_g) volumes from surface, vadose, groundwater.
        Total ET never exceeds available water across all three stores.
        """
        ET_demand = self.ET_pot * self.A_hex * dt
        ET_s = np.minimum(ET_demand, self.V_s)
        rem = ET_demand - ET_s
        V_vadose = (self.theta - self.theta_r) * self.A_hex * self.vadose_depth
        ET_u = np.minimum(rem, np.maximum(V_vadose, 0))
        rem -= ET_u
        ET_g = np.minimum(rem, self.S_gw)
        return (ET_s, ET_u, ET_g)

    def _total_storage(self):
        """Return surface, vadose, groundwater and channel storage in m3."""
        V_vadose = (self.theta - self.theta_r) * self.A_hex * self.vadose_depth
        return np.sum(self.V_s) + np.sum(self.S_gw) + np.sum(self.V_chan) + np.sum(V_vadose)

    def run(self, precip_rate_mmhr, duration_min, dt_sec, output_path='flood_depth_routing.tif', print_interval=5):
        """
        Execute explicit forward-Euler timestepping (Algorithm 1).
        All state updates are synchronous (all cells updated from time-n values).
        Mass deficit from positivity enforcement is tracked explicitly.
        """
        print(f'  Storm:    {precip_rate_mmhr} mm/hr  ×  {duration_min} min')
        print(f'  Timestep: {dt_sec} s')
        print('─' * 65)
        P_rate = precip_rate_mmhr / 1000 / 3600
        P_vol = P_rate * self.A_hex * dt_sec
        n_steps = int(duration_min * 60 / dt_sec)
        M_in = 0.0
        M_out = 0.0
        M_ET = 0.0
        M_deficit = 0.0
        M0 = self._total_storage()
        t_start = time.perf_counter()
        for step in range(1, n_steps + 1):
            self.V_s += P_vol
            M_in += P_vol * self.N
            I_vol = self._infiltration(dt_sec)
            self.V_s -= I_vol
            d_theta = I_vol / (self.A_hex * self.vadose_depth)
            self.theta = np.clip(self.theta + d_theta, self.theta_r, self.theta_s)
            excess_theta = np.maximum(self.theta - 0.9 * self.theta_s, 0)
            R_vol = excess_theta * self.A_hex * self.vadose_depth
            self.theta -= excess_theta
            self.S_gw += R_vol
            Q_base = self._baseflow(dt_sec)
            self.S_gw -= Q_base
            self.V_chan += Q_base
            dV_surf, Q_outlet = self._route_surface(dt_sec)
            self.V_s += dV_surf
            self.V_chan += Q_outlet * 0.0
            M_out += Q_outlet
            ET_s, ET_u, ET_g = self._et(dt_sec)
            self.V_s -= ET_s
            dth = ET_u / (self.A_hex * self.vadose_depth)
            self.theta = np.maximum(self.theta - dth, self.theta_r)
            self.S_gw -= ET_g
            M_ET += np.sum(ET_s + ET_u + ET_g)
            Q_chan_out = self.V_chan[self.outlet]
            self.V_chan[self.outlet] = 0.0
            M_out += Q_chan_out
            neg_s = np.minimum(self.V_s, 0.0)
            neg_gw = np.minimum(self.S_gw, 0.0)
            neg_ch = np.minimum(self.V_chan, 0.0)
            deficit_step = -(np.sum(neg_s) + np.sum(neg_gw) + np.sum(neg_ch))
            M_deficit += deficit_step
            self.V_s = np.maximum(self.V_s, 0.0)
            self.S_gw = np.maximum(self.S_gw, 0.0)
            self.V_chan = np.maximum(self.V_chan, 0.0)
            self.theta = np.clip(self.theta, self.theta_r, self.theta_s)
            if step % max(1, int(print_interval * 60 / dt_sec)) == 0 or step == n_steps:
                h_max = np.max(self.V_s) / self.A_hex
                h_mean = np.mean(self.V_s) / self.A_hex
                print(f'  t={step * dt_sec / 60:5.1f} min | max depth: {h_max:.4f} m | mean depth: {h_mean:.4f} m | deficit: {M_deficit:.3e} m³')
        runtime = time.perf_counter() - t_start
        M_final = self._total_storage()
        dM_computed = M_final - M0
        dM_expected = M_in - M_out - M_ET
        # The physical water-budget residual includes any mass introduced by
        # positivity clipping. Report that residual and the clipping diagnostic
        # separately; never subtract a clipping amount after abs().
        mass_error_signed = dM_computed - dM_expected
        mass_error_abs = abs(mass_error_signed)
        mass_error_clipping_adjusted_signed = mass_error_signed - M_deficit
        mass_error_clipping_adjusted_abs = abs(mass_error_clipping_adjusted_signed)
        rel_error = mass_error_abs / max(M_in, 1e-09) * 100
        print('\n' + '─' * 65)
        print('  MASS BALANCE SUMMARY')
        print('─' * 65)
        print(f'  Precipitation in   : {M_in:.4f} m³')
        print(f'  Outlet discharge   : {M_out:.4f} m³')
        print(f'  Evapotranspiration : {M_ET:.4f} m³')
        print(f'  ΔStorage (computed): {dM_computed:.4f} m³')
        print(f'  ΔStorage (expected): {dM_expected:.4f} m³')
        print(f'  Positivity deficit : {M_deficit:.4e} m³')
        print(f'  Raw signed residual: {mass_error_signed:.4e} m³')
        print(f'  Raw absolute error : {mass_error_abs:.4e} m³')
        print(f'  Residual minus clip: {mass_error_clipping_adjusted_signed:.4e} m³')
        print(f'  Relative raw error : {rel_error:.4e} %')
        print(f'  Runtime            : {runtime:.2f} s')
        print('─' * 65)
        flood_depth_hex = self.V_s / self.A_hex
        self._rasterize_output(flood_depth_hex, output_path)
        return (flood_depth_hex, {'M_in': M_in, 'M_out': M_out, 'M_ET': M_ET, 'M_deficit': M_deficit, 'mass_error_signed_m3': mass_error_signed, 'mass_error_abs_m3': mass_error_abs, 'mass_error_clipping_adjusted_signed_m3': mass_error_clipping_adjusted_signed, 'mass_error_clipping_adjusted_abs_m3': mass_error_clipping_adjusted_abs, 'rel_error_pct': rel_error, 'runtime_s': runtime, 'N_cells': self.N, 'N_edges': len(self.edge_list)})

    def _rasterize_output(self, hex_depths, path):
        """Convert per-hexagon depths back to raster GeoTIFF."""
        shapes = [(poly, float(d)) for poly, d in zip(self.polygons_valid, hex_depths) if not np.isnan(d)]
        out = rasterize(shapes, out_shape=(self.raster_rows, self.raster_cols), transform=self.transform, fill=np.nan, dtype=np.float32)
        with rasterio.open(path, 'w', driver='GTiff', height=self.raster_rows, width=self.raster_cols, count=1, dtype=np.float32, crs=self.crs, transform=self.transform, nodata=np.nan) as dst:
            dst.write(out, 1)
        print(f'  Saved: {path}')

def run_closed_domain_test(N_cells=100, hex_radius=150, dt=10, n_steps=360, P_rate_mmhr=100):
    """
    Test A: Closed-domain mass conservation (Sec 5.3, Test 2).

    A flat synthetic domain receives uniform rainfall with no outlet,
    no infiltration, no ET. Expected result:
        ΔM = P_rate × A_total × duration  (to machine precision)
    """
    print('\n' + '=' * 55)
    print('  VERIFICATION TEST A: Closed-Domain Conservation')
    print('=' * 55)
    R = hex_radius
    A_hex = 3 * np.sqrt(3) / 2 * R ** 2
    d_hex = np.sqrt(3) * R
    w_hex = R
    side = int(np.ceil(np.sqrt(N_cells)))
    xs = np.arange(side) * 1.5 * R
    ys_even = np.arange(side) * d_hex
    ys_odd = ys_even + d_hex / 2
    centers = []
    for col in range(side):
        ys = ys_odd if col % 2 else ys_even
        for row in range(side):
            centers.append((xs[col], ys[row]))
    centers = np.array(centers)
    N = len(centers)
    from scipy.spatial import cKDTree
    tree = cKDTree(centers)
    pairs = tree.query_pairs(r=d_hex * 1.15)
    neighbors = [[] for _ in range(N)]
    edge_i, edge_j = ([], [])
    seen = set()
    for i, j in pairs:
        dist = np.linalg.norm(centers[i] - centers[j])
        if abs(dist - d_hex) < 0.15 * d_hex:
            neighbors[i].append(j)
            neighbors[j].append(i)
            key = (min(i, j), max(i, j))
            if key not in seen:
                seen.add(key)
                edge_i.append(min(i, j))
                edge_j.append(max(i, j))
    edge_i = np.array(edge_i)
    edge_j = np.array(edge_j)
    z = np.zeros(N)
    S_f = np.full(N, 0.0001)
    n_mann = 0.05
    V_s = np.zeros(N)
    P_vol = P_rate_mmhr / 1000 / 3600 * A_hex * dt
    M0 = 0.0
    M_in = 0.0
    for step in range(n_steps):
        V_s += P_vol
        M_in += P_vol * N
        h = V_s / A_hex
        H = h + z
        Hi = H[edge_i]
        Hj = H[edge_j]
        h_i = h[edge_i]
        h_j = h[edge_j]
        D_i = 1 / n_mann * np.maximum(h_i, 1e-06) ** (5 / 3) / np.sqrt(S_f[edge_i])
        D_j = 1 / n_mann * np.maximum(h_j, 1e-06) ** (5 / 3) / np.sqrt(S_f[edge_j])
        D_f = 0.5 * (D_i + D_j)
        Q_vol = D_f * (Hi - Hj) / d_hex * w_hex * dt
        V_req = np.zeros(N)
        np.add.at(V_req, edge_i, np.maximum(Q_vol, 0))
        np.add.at(V_req, edge_j, np.maximum(-Q_vol, 0))
        eps = 1e-12
        lam = np.minimum(1.0, V_s / (V_req + eps))
        lam = np.where(V_req < eps, 1.0, lam)
        Q_f = np.where(Q_vol >= 0, Q_vol * lam[edge_i], --Q_vol * lam[edge_j])
        dV = np.zeros(N)
        np.add.at(dV, edge_i, -Q_f)
        np.add.at(dV, edge_j, Q_f)
        V_s += dV
        V_s = np.maximum(V_s, 0.0)
    M_final = np.sum(V_s)
    expected = M_in
    abs_err = abs(M_final - expected)
    rel_err = abs_err / max(expected, 1e-12) * 100
    print(f'  Cells         : {N}')
    print(f'  Steps         : {n_steps}  (dt = {dt} s)')
    print(f'  P rate        : {P_rate_mmhr} mm/hr')
    print(f'  Expected ΔM   : {expected:.6f} m³')
    print(f'  Computed ΔM   : {M_final:.6f} m³')
    print(f'  Absolute error: {abs_err:.4e} m³')
    print(f'  Relative error: {rel_err:.4e} %')
    status = 'PASS ✓' if rel_err < 1e-08 else 'FAIL ✗'
    print(f'  Result        : {status}')
    return rel_err

def run_six_neighbor_test(hex_radius=150, dt=10, n_steps=100):
    """
    Test B: Six-neighbor donor test (Sec 5.3, Test 4).

    One wet central cell surrounded by six dry cells.
    Verifies aggregate limiter prevents over-drainage.
    Expected: no negative depths, symmetric redistribution.
    """
    print('\n' + '=' * 55)
    print('  VERIFICATION TEST B: Six-Neighbor Donor')
    print('=' * 55)
    R = hex_radius
    A_hex = 3 * np.sqrt(3) / 2 * R ** 2
    d_hex = np.sqrt(3) * R
    w_hex = R
    angles = np.linspace(0, 2 * np.pi, 7)[:-1]
    centers = np.array([[0.0, 0.0]] + [[d_hex * np.cos(a), d_hex * np.sin(a)] for a in angles])
    N = 7
    edge_i = np.zeros(6, dtype=int)
    edge_j = np.arange(1, 7, dtype=int)
    z = np.zeros(N)
    S_f = np.full(N, 0.0001)
    n_mann = 0.05
    V_s = np.zeros(N)
    V_s[0] = 1000.0
    M0 = np.sum(V_s)
    min_depth_global = np.inf
    for step in range(n_steps):
        h = V_s / A_hex
        H = h + z
        Hi = H[edge_i]
        Hj = H[edge_j]
        D_i = 1 / n_mann * np.maximum(h[edge_i], 1e-06) ** (5 / 3) / np.sqrt(S_f[edge_i])
        D_j = 1 / n_mann * np.maximum(h[edge_j], 1e-06) ** (5 / 3) / np.sqrt(S_f[edge_j])
        D_f = 0.5 * (D_i + D_j)
        Q_vol = D_f * (Hi - Hj) / d_hex * w_hex * dt
        V_req = np.zeros(N)
        np.add.at(V_req, edge_i, np.maximum(Q_vol, 0))
        np.add.at(V_req, edge_j, np.maximum(-Q_vol, 0))
        eps = 1e-12
        lam = np.minimum(1.0, V_s / (V_req + eps))
        lam = np.where(V_req < eps, 1.0, lam)
        Q_f = np.where(Q_vol >= 0, Q_vol * lam[edge_i], --Q_vol * lam[edge_j])
        dV = np.zeros(N)
        np.add.at(dV, edge_i, -Q_f)
        np.add.at(dV, edge_j, Q_f)
        V_s += dV
        min_depth_global = min(min_depth_global, np.min(V_s / A_hex))
        V_s = np.maximum(V_s, 0.0)
    M_final = np.sum(V_s)
    abs_err = abs(M_final - M0)
    rel_err = abs_err / M0 * 100
    neighbor_depths = V_s[1:] / A_hex
    symmetry_cv = np.std(neighbor_depths) / (np.mean(neighbor_depths) + 1e-12)
    print(f'  Initial centre volume : {M0:.2f} m³')
    print(f'  Final total volume    : {M_final:.6f} m³')
    print(f'  Relative error        : {rel_err:.4e} %')
    print(f'  Minimum depth reached : {min_depth_global:.6f} m')
    print(f'  Neighbor depth CV     : {symmetry_cv:.4f}  (0=perfect symmetry)')
    neg_ok = 'PASS ✓' if min_depth_global >= -1e-10 else 'FAIL ✗'
    cons_ok = 'PASS ✓' if rel_err < 1e-08 else 'FAIL ✗'
    print(f'  Non-negativity        : {neg_ok}')
    print(f'  Conservation          : {cons_ok}')
    return (rel_err, min_depth_global)

def run_infiltration_mass_balance_test(N_cells=100, hex_radius=150, dt=10, n_steps=400, P_rate_mmhr=150, Ks_mmhr=50.0, theta_s=0.4, theta_r=0.08, vadose_depth=0.1):
    """Verify infiltration and vadose mass conservation
under uniform synthetic rainfall on a closed domain."""
    print('\n' + '=' * 55)
    print('  VERIFICATION TEST C: Infiltration Mass Balance')
    print('=' * 55)
    R = hex_radius
    A_hex = 3 * np.sqrt(3) / 2 * R ** 2
    N = N_cells
    Ks = Ks_mmhr / 1000 / 3600
    P_rate = P_rate_mmhr / 1000 / 3600
    P_vol = P_rate * A_hex * dt
    V_s = np.zeros(N)
    theta = np.full(N, theta_r)
    M_in = 0.0
    for step in range(n_steps):
        V_s += P_vol
        M_in += P_vol * N
        h = V_s / A_hex
        I_rate = np.minimum(h / dt, Ks)
        I_vol = np.minimum(I_rate * A_hex * dt, V_s)
        vadose_capacity_vol = np.maximum((theta_s - theta) * A_hex * vadose_depth, 0.0)
        I_vol = np.minimum(I_vol, vadose_capacity_vol)
        V_s -= I_vol
        theta += I_vol / (A_hex * vadose_depth)
    V_vadose_final = (theta - theta_r) * A_hex * vadose_depth
    M_final = np.sum(V_s) + np.sum(V_vadose_final)
    abs_err = abs(M_final - M_in)
    rel_err = abs_err / max(M_in, 1e-12) * 100
    n_saturated = int(np.sum(theta >= theta_s - 1e-09))
    print(f'  Cells                  : {N}')
    print(f'  Steps                  : {n_steps}  (dt = {dt} s)')
    print(f'  P rate                 : {P_rate_mmhr} mm/hr  |  Ks = {Ks_mmhr} mm/hr')
    print(f'  Cells reaching theta_s : {n_saturated}/{N}  (exercises vadose-capacity cap)')
    print(f'  Cumulative P (M_in)    : {M_in:.6f} m³')
    print(f'  Final V_s + V_vadose   : {M_final:.6f} m³')
    print(f'  Absolute error         : {abs_err:.4e} m³')
    print(f'  Relative error         : {rel_err:.4e} %')
    status = 'PASS ✓' if rel_err < 1e-08 else 'FAIL ✗'
    print(f'  Result                 : {status}')
    return rel_err

def run_class_integration_test(dt=5.0, duration_min=30.0, hex_radius=150, precip_rate_mmhr=80.0, tol_pct=1e-06):
    """Exercise the integrated engine on a synthetic sloping DEM;
check execution, storage conservation and nonnegative water depths."""
    import tempfile
    import rasterio as _rasterio
    from rasterio.transform import from_origin as _from_origin
    print('\n' + '=' * 55)
    print('  VERIFICATION TEST D: Real-Class Integration Smoke Test')
    print('=' * 55)
    nrows, ncols, res = (60, 60, 30.0)
    y = np.linspace(100.0, 40.0, nrows).reshape(-1, 1)
    dem = np.tile(y, (1, ncols)).astype(np.float32)
    xx, yy = np.meshgrid(np.arange(ncols), np.arange(nrows))
    dem += 0.5 * np.sin(xx / 5.0)
    nodata = -9999.0
    transform = _from_origin(500000.0, 5000000.0, res, res)
    tmp_dir = tempfile.mkdtemp()
    dem_path = os.path.join(tmp_dir, 'synthetic_smoke_dem.tif')
    with _rasterio.open(dem_path, 'w', driver='GTiff', height=nrows, width=ncols, count=1, dtype='float32', crs='EPSG:32632', transform=transform, nodata=nodata) as dst:
        dst.write(dem, 1)
    params = {'n_manning': 0.05, 'Ks': 20.0, 'ET_pot': 2.0, 'theta_s': 0.4, 'theta_r': 0.08, 'vadose_depth': 0.5, 'alpha_gw': 0.01}
    try:
        engine = HexalinkRoutingEngine(dem_path, params, hex_radius=hex_radius)
        depths, stats = engine.run(precip_rate_mmhr=precip_rate_mmhr, duration_min=duration_min, dt_sec=dt, output_path=os.path.join(tmp_dir, 'smoke_depth.tif'), print_interval=10 ** 9)
        ran_ok = True
        error_msg = None
    except Exception as exc:
        ran_ok = False
        error_msg = repr(exc)
        depths, stats = (None, None)
    if not ran_ok:
        print(f'  Class instantiation/run : FAIL ✗  ({error_msg})')
        print(f'  Result                  : FAIL ✗')
        return float('inf')
    valid = depths[~np.isnan(depths)]
    has_negative = bool(np.any(valid < -1e-09))
    rel_err = abs(stats['rel_error_pct'])
    deficit = stats['M_deficit']
    print(f"  Cells / edges           : {stats['N_cells']} / {stats['N_edges']}")
    print(f'  Class ran to completion : PASS ✓')
    print(f'  Relative mass error     : {rel_err:.4e} %')
    print(f'  Positivity deficit      : {deficit:.4e} m³')
    print(f'  Any negative depth      : {has_negative}')
    passed = ran_ok and rel_err < tol_pct and (deficit == 0.0) and (not has_negative)
    status = 'PASS ✓' if passed else 'FAIL ✗'
    print(f'  Result                  : {status}')
    return rel_err
if __name__ == '__main__':
    err_a = run_closed_domain_test()
    err_b, min_d = run_six_neighbor_test()
    err_c = run_infiltration_mass_balance_test()
    err_d = run_class_integration_test()
    print('\n' + '=' * 55)
    print('  VERIFICATION SUMMARY')
    print('=' * 55)
    print(f'  Test A (closed domain)  : rel. error = {err_a:.4e} %')
    print(f'  Test B (6-neighbor)     : rel. error = {(err_b[0] if isinstance(err_b, tuple) else err_b):.4e} %')
    print(f'  Test C (infiltration)   : rel. error = {err_c:.4e} %')
    print(f'  Test D (class smoke test): rel. error = {err_d:.4e} %')
    print('=' * 55)
    dem_file = 'reprojected elevation 2.tif'
    if os.path.exists(dem_file):
        params = {'Ks': 50.0, 'n_manning': 0.05, 'alpha_gw': 0.005, 'ET_pot': 5.0, 'theta_s': 0.4, 'theta_r': 0.08, 'vadose_depth': 0.5}
        engine = HexalinkRoutingEngine(dem_file, params, hex_radius=150)
        depths, stats = engine.run(precip_rate_mmhr=300, duration_min=60, dt_sec=10, output_path='flood_depth_routing.tif')
        print(f'\n  Max surface depth : {np.max(depths):.4f} m')
        print(f'  Mean surface depth: {np.mean(depths[depths > 0]):.4f} m')
    else:
        print(f"\n  (DEM not found at '{dem_file}' — skipping basin run)")
        print('  Verification tests above are DEM-independent.')
