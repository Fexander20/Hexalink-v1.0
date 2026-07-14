#!/usr/bin/env python3
"""
Hexalink v1.0 - Full Six-Node Diffusive-Kinematic Routing Engine
True hexagonal 6-neighbor adjacency via precomputed adjacency list.
Aggregate donor-based positivity limiter (Eqs 5.4-5.8).
Exact mass deficit tracking with no silent clipping.
Manuscript References: Sections 3.2-3.8, 4.1-4.2, 5.1-5.4, Algorithm 1

v1.0 changelog (five code fixes applied to this file):
  Fix #1 - _D() (Eq. 3.13 diffusion/conductance coefficient) previously
           divided by self.S_f, a per-CELL array, while being called
           with edge-indexed depth values (h[edge_i], h[edge_j]) — a
           shape mismatch that crashed as soon as N_edges != N_cells,
           i.e. on any real grid. S_f is now passed in explicitly and
           indexed at the same edge endpoints as h. This bug was never
           caught by the bundled Test A/B/C because none of them
           instantiate HexalinkRoutingEngine or call _route_surface —
           see note in the class docstring below.
  Fix #2 - Vadose storage (theta) now included in M0/M_final via
           _total_storage(); it was previously invisible to the mass
           balance check.
  Fix #3 - Infiltration is now capped by remaining vadose capacity
           (theta_s - theta), not just by available surface water, so
           excess rainfall that cannot physically infiltrate stays on
           the surface instead of being silently discarded.
  Fix #4 - Outlet boundary discharge demand is folded into the donor
           aggregate limiter (V_req_out) before lambda is computed, so
           the outlet cell cannot be overdrawn by internal edges and
           the outlet face simultaneously.
  Fix #5 - D_face now uses the harmonic mean (2*Di*Dj/(Di+Dj)), matching
           the docstring/Sec 3.6; it previously used an arithmetic mean.
  Added  - Test C: infiltration + vadose mass-balance closure test.
  Added  - Test D: real-class integration smoke test. Tests A-C never
           instantiate HexalinkRoutingEngine itself (see Fix #1 note and
           the class docstring); Test D closes that specific gap with a
           fast synthetic-DEM run that actually exercises the class.
"""

import numpy as np
import rasterio
import time
import os
from shapely.geometry import Polygon, box
from rasterio.features import rasterize


# ─────────────────────────────────────────────────────────────────
# SECTION 1: Hexagonal Grid Construction
# ─────────────────────────────────────────────────────────────────

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
    A_hex = (3 * np.sqrt(3) / 2) * R**2          # Eq 3.1
    d_hex = np.sqrt(3) * R                         # Eq 3.2
    w_hex = R                                      # Eq 3.3

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
        # Offset every other column
        y_off = (d_hex / 2) * (col_idx % 2)
        y_start = miny - y_off
        y = y_start
        while y <= maxy + R:
            # Hexagon vertices (pointy-top)
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
    return centers, polygons, A_hex, d_hex, w_hex


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

    # Use spatial binning for efficiency
    from scipy.spatial import cKDTree
    tree = cKDTree(centers)
    # Query slightly beyond d_hex to catch all 6 neighbors
    pairs = tree.query_pairs(r=d_hex * (1 + tol))

    for i, j in pairs:
        dist = np.linalg.norm(centers[i] - centers[j])
        if abs(dist - d_hex) < tol * d_hex:
            neighbors[i].append(j)
            neighbors[j].append(i)
            edge_list.append((i, j))

    return neighbors, edge_list


def assign_dem_to_hexagons(centers, polygons, dem, transform, nodata, rows, cols):
    """
    Assign mean DEM elevation to each hexagonal cell.
    Returns hex_elev (N,) array; NaN for cells with no valid DEM pixels.
    """
    N = len(polygons)
    shapes = [(p, i + 1) for i, p in enumerate(polygons)]
    hex_ids = rasterize(shapes, out_shape=(rows, cols),
                        transform=transform, fill=0, dtype=np.int32)

    dem_flat = dem.ravel().astype(np.float64)
    ids_flat = hex_ids.ravel()

    valid = (ids_flat > 0)
    if nodata is not None:
        valid &= (dem_flat != nodata)
    valid &= ~np.isnan(dem_flat)

    hex_sum = np.zeros(N + 1, dtype=np.float64)
    hex_cnt = np.zeros(N + 1, dtype=np.float64)
    np.add.at(hex_sum, ids_flat[valid], dem_flat[valid])
    np.add.at(hex_cnt, ids_flat[valid], 1)

    hex_elev = np.full(N, np.nan, dtype=np.float64)
    ok = hex_cnt[1:] > 0
    hex_elev[ok] = hex_sum[1:][ok] / hex_cnt[1:][ok]
    return hex_elev, hex_ids


# ─────────────────────────────────────────────────────────────────
# SECTION 2: Routing Engine
# ─────────────────────────────────────────────────────────────────

class HexalinkRoutingEngine:
    """
    Full six-node diffusive-kinematic routing on true hexagonal finite volumes.

    State vectors (length N = number of valid hexagonal cells):
        V_s    : surface water volume [m³]
        theta  : vadose zone volumetric moisture content [-]
        S_gw   : groundwater storage [m³]
        V_chan : channel storage [m³]

    Mass balance:
        Tracked explicitly at every timestep.
        Deficit from positivity enforcement is recorded, not silently discarded.
        Relative mass error reported at end of run (Eq 4.6).

    IMPORTANT — test coverage: run_closed_domain_test(), run_six_neighbor_test(),
    and run_infiltration_mass_balance_test() below are standalone, hand-rolled
    synthetic tests. None of them instantiate this class or call its methods
    (_route_surface, _D, run, etc.) — they reimplement a simplified version of the
    same logic independently. Passing all three tells you the *reference math* is
    self-consistent; it does NOT by itself verify that this class is free of bugs, as
    demonstrated by Fix #1 above, which crashed this class on first real use while
    all three bundled tests reported near-zero error. run_class_integration_test()
    (Test D) closes this specific gap with a fast (<1 s) smoke test that actually
    instantiates this class against a small synthetic DEM and checks it runs,
    conserves mass, and stays non-negative — but it is a smoke test, not a substitute
    for Tests A-C's specific numerical guarantees. A full run against a real DEM
    (see `__main__` below) remains the most complete check of actual production
    behavior.
    """

    def __init__(self, dem_path, params, hex_radius=150):
        print("=" * 65)
        print("  Hexalink v1.0  |  True Hexagonal Routing Engine")
        print("=" * 65)

        # ── Load DEM ──────────────────────────────────────────────
        with rasterio.open(dem_path) as src:
            dem_full = src.read(1).astype(np.float64)
            self.nodata = src.nodata
            self.transform = src.transform
            self.crs = src.crs
            self.raster_rows, self.raster_cols = src.shape
            self.bounds = src.bounds

        # ── Build hexagonal grid ───────────────────────────────────
        print(f"\n[1/4] Building hexagonal grid  (R = {hex_radius} m) ...")
        self.R = hex_radius
        centers, polygons, self.A_hex, self.d_hex, self.w_hex = \
            build_hex_grid(tuple(self.bounds), hex_radius)
        print(f"      Generated {len(centers):,} candidate hexagons")

        # ── Assign DEM elevations ──────────────────────────────────
        print("[2/4] Assigning DEM to hexagons ...")
        hex_elev_all, self.hex_ids_raster = assign_dem_to_hexagons(
            centers, polygons, dem_full,
            self.transform, self.nodata,
            self.raster_rows, self.raster_cols)

        # Keep only valid (non-NaN) hexagons
        self.valid_idx = np.where(~np.isnan(hex_elev_all))[0]
        self.N = len(self.valid_idx)
        self.z = hex_elev_all[self.valid_idx]          # bed elevation [m]
        self.centers_valid = centers[self.valid_idx]
        self.polygons_valid = [polygons[i] for i in self.valid_idx]

        # Map from original index → compact index
        full_to_compact = np.full(len(centers), -1, dtype=np.int64)
        full_to_compact[self.valid_idx] = np.arange(self.N)

        print(f"      Valid hexagonal cells: {self.N:,}")

        # ── Build adjacency (true 6-neighbor) ─────────────────────
        print("[3/4] Building 6-neighbor adjacency ...")
        neighbors_full, edge_list_full = build_adjacency(centers, self.d_hex)

        # Remap to compact indices; keep only edges where both endpoints are valid
        self.neighbors = [[] for _ in range(self.N)]
        self.edge_list = []          # (i_compact, j_compact)  i < j
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
        avg_neighbors = sum(len(nb) for nb in self.neighbors) / max(self.N, 1)
        print(f"      Internal edges: {n_edges:,}  |  Avg neighbors: {avg_neighbors:.2f}")

        # ── Identify outlet cell ───────────────────────────────────
        # Outlet = valid cell with the fewest neighbors (boundary) AND lowest elevation
        n_neighbors = np.array([len(nb) for nb in self.neighbors])
        boundary_mask = n_neighbors < 6
        boundary_idx = np.where(boundary_mask)[0]
        if len(boundary_idx) == 0:
            boundary_idx = np.arange(self.N)
        self.outlet = boundary_idx[np.argmin(self.z[boundary_idx])]
        print(f"      Outlet cell: {self.outlet}  "
              f"(z = {self.z[self.outlet]:.2f} m, "
              f"neighbors = {n_neighbors[self.outlet]})")

        # ── Parameters (Sec 3.3, 3.5, 3.6, 3.8) ──────────────────
        print("[4/4] Setting parameters ...")
        self.Ks       = params.get('Ks',       50.0) / 1000 / 3600   # m/s
        self.n_mann   = params.get('n_manning', 0.05)
        self.alpha_gw = params.get('alpha_gw',  0.005) / 3600        # 1/s
        self.ET_pot   = params.get('ET_pot',    5.0) / 1000 / 3600   # m/s
        self.theta_s  = params.get('theta_s',   0.40)
        self.theta_r  = params.get('theta_r',   0.08)
        self.vadose_depth = params.get('vadose_depth', 0.5)          # m

        # ── State vectors ──────────────────────────────────────────
        self.V_s    = np.zeros(self.N, dtype=np.float64)  # surface [m³]
        self.theta  = np.full(self.N, self.theta_r)        # vadose [-]
        self.S_gw   = np.zeros(self.N, dtype=np.float64)  # groundwater [m³]
        self.V_chan = np.zeros(self.N, dtype=np.float64)  # channel [m³]

        # ── Precompute friction slope ──────────────────────────────
        # For each cell, use mean absolute elevation difference to neighbors
        self.S_f = np.full(self.N, 1e-4, dtype=np.float64)
        for i in range(self.N):
            nb = self.neighbors[i]
            if nb:
                dz = np.abs(self.z[i] - self.z[nb])
                self.S_f[i] = max(np.mean(dz) / self.d_hex, 1e-4)

        print(f"\n  Grid ready. {self.N:,} cells | {n_edges:,} edges\n")

    # ── Diffusion coefficient ──────────────────────────────────────
    def _D(self, h, S_f):
        """D = (1/n) * h^(5/3) * S_f^(-1/2)   [m²/s]  (Eq 3.13)"""
        h_eff = np.maximum(h, 1e-6)
        return (1.0 / self.n_mann) * (h_eff ** (5/3)) / np.sqrt(S_f)

    # ── Surface routing with aggregate donor limiter ───────────────
    def _route_surface(self, dt):
        """
        Compute antisymmetric diffusive fluxes on all internal edges,
        then apply aggregate donor-based positivity limiter (Eqs 5.4-5.8).

        FIX #4 (outlet donor accounting): the outlet cell can lose water
        two ways in the same step — through internal edges to neighbors,
        and through the outlet boundary discharge. Previously these were
        limited independently, so the outlet cell could be overdrawn by
        the sum of both. The outlet's boundary discharge demand is now
        folded into V_req_out BEFORE lambda is computed, so a single
        donor-limiter scaling factor governs all outflow from that cell.

        Returns dV_s (N,): net volume change for each cell this timestep.
        Also returns Q_outlet: volume discharged at outlet this step.
        """
        h = self.V_s / self.A_hex                  # water depth [m]
        H = h + self.z                              # hydraulic head [m]

        # ── Raw flux on every edge (one per shared face) ──────────
        # Q_ij > 0 means flow from i to j
        Hi = H[self.edge_i]
        Hj = H[self.edge_j]
        Di = self._D(h[self.edge_i], self.S_f[self.edge_i])
        Dj = self._D(h[self.edge_j], self.S_f[self.edge_j])
        # FIX #5: harmonic mean (was arithmetic mean, contradicting the
        # docstring/comment). Harmonic mean correctly lets the lower-
        # conductivity face dominate, matching Sec 3.6.
        D_face = (2.0 * Di * Dj) / (Di + Dj + 1e-12)

        # Eq 3.12 adapted: Q_ij = D_face * (H_i - H_j) / d * w
        Q_raw = D_face * (Hi - Hj) / self.d_hex * self.w_hex  # [m³/s]
        Q_vol = Q_raw * dt                          # [m³] per timestep

        # ── Outlet boundary discharge demand (computed but NOT yet
        #    applied) — Manning discharge through the open outlet face ──
        h_out = max(self.V_s[self.outlet] / self.A_hex, 0.0)
        Q_out_rate = (1.0 / self.n_mann) * (h_out ** (5/3)) * \
                     np.sqrt(self.S_f[self.outlet]) * self.w_hex
        Q_outlet_demand = Q_out_rate * dt           # [m³], uncapped

        # ── Aggregate donor limiter (Eqs 5.4 – 5.8) ──────────────
        # Step 1: accumulate total outward demand per donor cell,
        # INCLUDING the outlet cell's own boundary discharge demand.
        V_req_out = np.zeros(self.N, dtype=np.float64)
        np.add.at(V_req_out, self.edge_i,  np.maximum( Q_vol, 0))
        np.add.at(V_req_out, self.edge_j,  np.maximum(-Q_vol, 0))
        V_req_out[self.outlet] += Q_outlet_demand

        # Step 2: available volume (surface only; sources already added this step)
        V_avail = np.maximum(self.V_s, 0.0)        # after precip added

        # Step 3: scaling factor λ_i  (Eq 5.6)
        eps = 1e-12
        lam = np.minimum(1.0, V_avail / (V_req_out + eps))
        lam = np.where(V_req_out < eps, 1.0, lam)  # no demand → no scaling

        # Step 4: scale outward fluxes by donor λ
        lam_i = lam[self.edge_i]
        lam_j = lam[self.edge_j]

        Q_final = np.where(
            Q_vol >= 0,
            Q_vol * lam_i,     # flow i→j: donor is i
           -(-Q_vol) * lam_j   # flow j→i: donor is j
        )

        # Step 5: antisymmetric assignment — one value, opposite signs (Eq 5.8)
        dV = np.zeros(self.N, dtype=np.float64)
        np.add.at(dV, self.edge_i, -Q_final)
        np.add.at(dV, self.edge_j,  Q_final)

        # ── Outlet discharge, scaled by the SAME lambda used for its
        #    internal edges, so total outlet-cell withdrawal (internal +
        #    boundary) never exceeds V_avail[outlet] ──
        Q_outlet = Q_outlet_demand * lam[self.outlet]
        dV[self.outlet] -= Q_outlet

        return dV, Q_outlet

    # ── Infiltration ───────────────────────────────────────────────
    def _infiltration(self, dt):
        """
        Saturation-excess (Eq 3.6). Returns I_vol [m³].

        FIX #3: previously I_vol was capped only by available surface
        water (V_s), then applied to theta and silently clipped to
        theta_s — any infiltration volume beyond remaining vadose
        capacity was subtracted from V_s but never actually stored,
        destroying mass. I_vol is now also capped by the remaining
        vadose storage capacity (theta_s - theta) BEFORE it is
        returned, so water that cannot physically infiltrate is left
        on the surface instead of disappearing.
        """
        h = self.V_s / self.A_hex
        I_rate = np.minimum(h / dt, self.Ks)
        I_vol = I_rate * self.A_hex * dt
        I_vol = np.minimum(I_vol, self.V_s)

        vadose_capacity_vol = np.maximum(
            (self.theta_s - self.theta) * self.A_hex * self.vadose_depth, 0.0)
        I_vol = np.minimum(I_vol, vadose_capacity_vol)
        return I_vol

    # ── Groundwater baseflow ───────────────────────────────────────
    def _baseflow(self, dt):
        """Linear reservoir (Eqs 3.10-3.11). Returns Q_base [m³]."""
        return np.minimum(self.alpha_gw * self.S_gw * dt, self.S_gw)

    # ── ET ─────────────────────────────────────────────────────────
    def _et(self, dt):
        """
        Potential ET with stress factor (Eqs 3.19-3.20).
        Returns (ET_s, ET_u, ET_g) volumes from surface, vadose, groundwater.
        Total ET never exceeds available water across all three stores.
        """
        ET_demand = self.ET_pot * self.A_hex * dt

        # Surface first
        ET_s = np.minimum(ET_demand, self.V_s)
        rem = ET_demand - ET_s

        # Vadose second
        V_vadose = (self.theta - self.theta_r) * self.A_hex * self.vadose_depth
        ET_u = np.minimum(rem, np.maximum(V_vadose, 0))
        rem -= ET_u

        # Groundwater last
        ET_g = np.minimum(rem, self.S_gw)

        return ET_s, ET_u, ET_g

    # ── Total storage accounting ────────────────────────────────────
    def _total_storage(self):
        """
        Sum of all storage compartments [m³]: surface, vadose, groundwater,
        channel. FIX #2: vadose moisture (theta) was previously omitted
        from M0/M_final, so any water sitting in the vadose zone was
        invisible to the mass-balance check — it looked like it had left
        the system even though it was still stored. Vadose volume above
        residual moisture, (theta - theta_r) * A_hex * vadose_depth, is
        now included.
        """
        V_vadose = (self.theta - self.theta_r) * self.A_hex * self.vadose_depth
        return (np.sum(self.V_s) + np.sum(self.S_gw) +
                np.sum(self.V_chan) + np.sum(V_vadose))

    # ── Main timestepping loop ─────────────────────────────────────
    def run(self, precip_rate_mmhr, duration_min, dt_sec,
            output_path="flood_depth_routing.tif",
            print_interval=5):
        """
        Execute explicit forward-Euler timestepping (Algorithm 1).
        All state updates are synchronous (all cells updated from time-n values).
        Mass deficit from positivity enforcement is tracked explicitly.
        """
        print(f"  Storm:    {precip_rate_mmhr} mm/hr  ×  {duration_min} min")
        print(f"  Timestep: {dt_sec} s")
        print("─" * 65)

        P_rate = (precip_rate_mmhr / 1000) / 3600          # m/s
        P_vol  = P_rate * self.A_hex * dt_sec              # m³/step/cell
        n_steps = int((duration_min * 60) / dt_sec)

        # ── Mass balance accumulators (Eq 4.6) ────────────────────
        M_in        = 0.0   # cumulative precip volume [m³]
        M_out       = 0.0   # cumulative outlet discharge [m³]
        M_ET        = 0.0   # cumulative ET [m³]
        M_deficit   = 0.0   # cumulative positivity-enforcement deficit [m³]
        M0 = self._total_storage()                          # initial storage (Fix #2: includes vadose)

        t_start = time.perf_counter()

        for step in range(1, n_steps + 1):

            # ── Step A: Precipitation (source) ────────────────────
            self.V_s += P_vol
            M_in += P_vol * self.N

            # ── Step B: Infiltration ──────────────────────────────
            # I_vol is now pre-capped by remaining vadose capacity (Fix #3),
            # so this clip is a numerical safety net only — it should never
            # actually truncate theta anymore (no more silent mass loss).
            I_vol = self._infiltration(dt_sec)
            self.V_s    -= I_vol
            # Transfer to vadose
            d_theta = I_vol / (self.A_hex * self.vadose_depth)
            self.theta   = np.clip(self.theta + d_theta,
                                   self.theta_r, self.theta_s)

            # ── Step C: Vadose → Groundwater recharge ─────────────
            excess_theta = np.maximum(self.theta - 0.9 * self.theta_s, 0)
            R_vol = excess_theta * self.A_hex * self.vadose_depth
            self.theta -= excess_theta
            self.S_gw  += R_vol

            # ── Step D: Groundwater baseflow ──────────────────────
            Q_base = self._baseflow(dt_sec)
            self.S_gw   -= Q_base
            self.V_chan += Q_base

            # ── Step E: Surface routing (aggregate limiter) ────────
            dV_surf, Q_outlet = self._route_surface(dt_sec)
            self.V_s    += dV_surf
            self.V_chan += Q_outlet * 0.0   # outlet exits domain
            M_out       += Q_outlet

            # ── Step F: Evapotranspiration ────────────────────────
            ET_s, ET_u, ET_g = self._et(dt_sec)
            self.V_s    -= ET_s
            dth = ET_u / (self.A_hex * self.vadose_depth)
            self.theta   = np.maximum(self.theta - dth, self.theta_r)
            self.S_gw   -= ET_g
            M_ET        += np.sum(ET_s + ET_u + ET_g)

            # ── Step G: Channel outlet drain ──────────────────────
            Q_chan_out = self.V_chan[self.outlet]  # drain outlet channel cell
            self.V_chan[self.outlet] = 0.0
            M_out += Q_chan_out

            # ── Step H: Positivity enforcement (TRACKED) ──────────
            # Any remaining negative volumes are deficits — recorded, not hidden
            neg_s   = np.minimum(self.V_s,   0.0)
            neg_gw  = np.minimum(self.S_gw,  0.0)
            neg_ch  = np.minimum(self.V_chan, 0.0)
            deficit_step = -(np.sum(neg_s) + np.sum(neg_gw) + np.sum(neg_ch))
            M_deficit += deficit_step

            self.V_s    = np.maximum(self.V_s,    0.0)
            self.S_gw   = np.maximum(self.S_gw,   0.0)
            self.V_chan  = np.maximum(self.V_chan,  0.0)
            self.theta   = np.clip(self.theta, self.theta_r, self.theta_s)

            # ── Progress report ───────────────────────────────────
            if step % max(1, int(print_interval * 60 / dt_sec)) == 0 \
               or step == n_steps:
                h_max = np.max(self.V_s) / self.A_hex
                h_mean = np.mean(self.V_s) / self.A_hex
                print(f"  t={step*dt_sec/60:5.1f} min | "
                      f"max depth: {h_max:.4f} m | "
                      f"mean depth: {h_mean:.4f} m | "
                      f"deficit: {M_deficit:.3e} m³")

        runtime = time.perf_counter() - t_start

        # ── Mass balance verification (Eq 4.6) ────────────────────
        M_final = self._total_storage()                     # Fix #2: includes vadose
        dM_computed = M_final - M0
        dM_expected = M_in - M_out - M_ET
        mass_error_abs = abs(dM_computed - dM_expected) - M_deficit
        rel_error = (mass_error_abs / max(M_in, 1e-9)) * 100

        print("\n" + "─" * 65)
        print("  MASS BALANCE SUMMARY")
        print("─" * 65)
        print(f"  Precipitation in   : {M_in:.4f} m³")
        print(f"  Outlet discharge   : {M_out:.4f} m³")
        print(f"  Evapotranspiration : {M_ET:.4f} m³")
        print(f"  ΔStorage (computed): {dM_computed:.4f} m³")
        print(f"  ΔStorage (expected): {dM_expected:.4f} m³")
        print(f"  Positivity deficit : {M_deficit:.4e} m³")
        print(f"  Residual error     : {mass_error_abs:.4e} m³")
        print(f"  Relative error     : {rel_error:.4e} %")
        print(f"  Runtime            : {runtime:.2f} s")
        print("─" * 65)

        # ── Export raster ──────────────────────────────────────────
        flood_depth_hex = self.V_s / self.A_hex           # depth per hex cell
        self._rasterize_output(flood_depth_hex, output_path)

        return flood_depth_hex, {
            'M_in': M_in, 'M_out': M_out, 'M_ET': M_ET,
            'M_deficit': M_deficit, 'rel_error_pct': rel_error,
            'runtime_s': runtime, 'N_cells': self.N,
            'N_edges': len(self.edge_list)
        }

    def _rasterize_output(self, hex_depths, path):
        """Convert per-hexagon depths back to raster GeoTIFF."""
        shapes = [(poly, float(d))
                  for poly, d in zip(self.polygons_valid, hex_depths)
                  if not np.isnan(d)]
        out = rasterize(shapes,
                        out_shape=(self.raster_rows, self.raster_cols),
                        transform=self.transform,
                        fill=np.nan, dtype=np.float32)
        with rasterio.open(
            path, 'w', driver='GTiff',
            height=self.raster_rows, width=self.raster_cols,
            count=1, dtype=np.float32,
            crs=self.crs, transform=self.transform, nodata=np.nan
        ) as dst:
            dst.write(out, 1)
        print(f"  Saved: {path}")


# ─────────────────────────────────────────────────────────────────
# SECTION 3: Verification Tests
# ─────────────────────────────────────────────────────────────────

def run_closed_domain_test(N_cells=100, hex_radius=150, dt=10, n_steps=360,
                           P_rate_mmhr=100):
    """
    Test A: Closed-domain mass conservation (Sec 5.3, Test 2).

    A flat synthetic domain receives uniform rainfall with no outlet,
    no infiltration, no ET. Expected result:
        ΔM = P_rate × A_total × duration  (to machine precision)
    """
    print("\n" + "=" * 55)
    print("  VERIFICATION TEST A: Closed-Domain Conservation")
    print("=" * 55)

    R = hex_radius
    A_hex = (3 * np.sqrt(3) / 2) * R**2
    d_hex = np.sqrt(3) * R
    w_hex = R

    # Build a synthetic flat grid: NxN hex centers on a flat plane
    side = int(np.ceil(np.sqrt(N_cells)))
    xs = np.arange(side) * 1.5 * R
    ys_even = np.arange(side) * d_hex
    ys_odd  = ys_even + d_hex / 2

    centers = []
    for col in range(side):
        ys = ys_odd if col % 2 else ys_even
        for row in range(side):
            centers.append((xs[col], ys[row]))
    centers = np.array(centers)
    N = len(centers)

    # Build adjacency
    from scipy.spatial import cKDTree
    tree = cKDTree(centers)
    pairs = tree.query_pairs(r=d_hex * 1.15)
    neighbors = [[] for _ in range(N)]
    edge_i, edge_j = [], []
    seen = set()
    for i, j in pairs:
        dist = np.linalg.norm(centers[i] - centers[j])
        if abs(dist - d_hex) < 0.15 * d_hex:
            neighbors[i].append(j)
            neighbors[j].append(i)
            key = (min(i,j), max(i,j))
            if key not in seen:
                seen.add(key)
                edge_i.append(min(i,j))
                edge_j.append(max(i,j))
    edge_i = np.array(edge_i)
    edge_j = np.array(edge_j)

    # Flat elevation — zero gradient, no flow expected
    z = np.zeros(N)
    S_f = np.full(N, 1e-4)
    n_mann = 0.05

    V_s = np.zeros(N)
    P_vol = (P_rate_mmhr / 1000 / 3600) * A_hex * dt
    M0 = 0.0
    M_in = 0.0

    for step in range(n_steps):
        V_s += P_vol
        M_in += P_vol * N

        # No outlet, no losses — pure conservation check
        # Routing still happens but on flat terrain all fluxes ≈ 0
        h = V_s / A_hex
        H = h + z
        Hi = H[edge_i]; Hj = H[edge_j]
        h_i = h[edge_i]; h_j = h[edge_j]
        D_i = (1/n_mann) * np.maximum(h_i,1e-6)**(5/3) / np.sqrt(S_f[edge_i])
        D_j = (1/n_mann) * np.maximum(h_j,1e-6)**(5/3) / np.sqrt(S_f[edge_j])
        D_f = 0.5*(D_i+D_j)
        Q_vol = D_f * (Hi - Hj) / d_hex * w_hex * dt

        # Aggregate limiter
        V_req = np.zeros(N)
        np.add.at(V_req, edge_i, np.maximum( Q_vol, 0))
        np.add.at(V_req, edge_j, np.maximum(-Q_vol, 0))
        eps = 1e-12
        lam = np.minimum(1.0, V_s / (V_req + eps))
        lam = np.where(V_req < eps, 1.0, lam)
        Q_f = np.where(Q_vol >= 0, Q_vol * lam[edge_i],
                                   -(-Q_vol) * lam[edge_j])
        dV = np.zeros(N)
        np.add.at(dV, edge_i, -Q_f)
        np.add.at(dV, edge_j,  Q_f)
        V_s += dV
        V_s = np.maximum(V_s, 0.0)

    M_final = np.sum(V_s)
    expected = M_in
    abs_err = abs(M_final - expected)
    rel_err = abs_err / max(expected, 1e-12) * 100

    print(f"  Cells         : {N}")
    print(f"  Steps         : {n_steps}  (dt = {dt} s)")
    print(f"  P rate        : {P_rate_mmhr} mm/hr")
    print(f"  Expected ΔM   : {expected:.6f} m³")
    print(f"  Computed ΔM   : {M_final:.6f} m³")
    print(f"  Absolute error: {abs_err:.4e} m³")
    print(f"  Relative error: {rel_err:.4e} %")
    status = "PASS ✓" if rel_err < 1e-8 else "FAIL ✗"
    print(f"  Result        : {status}")
    return rel_err


def run_six_neighbor_test(hex_radius=150, dt=10, n_steps=100):
    """
    Test B: Six-neighbor donor test (Sec 5.3, Test 4).

    One wet central cell surrounded by six dry cells.
    Verifies aggregate limiter prevents over-drainage.
    Expected: no negative depths, symmetric redistribution.
    """
    print("\n" + "=" * 55)
    print("  VERIFICATION TEST B: Six-Neighbor Donor")
    print("=" * 55)

    R = hex_radius
    A_hex = (3 * np.sqrt(3) / 2) * R**2
    d_hex = np.sqrt(3) * R
    w_hex = R

    # 7 cells: centre + 6 neighbors
    angles = np.linspace(0, 2*np.pi, 7)[:-1]
    centers = np.array(
        [[0.0, 0.0]] +
        [[d_hex * np.cos(a), d_hex * np.sin(a)] for a in angles]
    )
    N = 7

    # All 6 edges connect center (0) to neighbors (1-6)
    edge_i = np.zeros(6, dtype=int)
    edge_j = np.arange(1, 7, dtype=int)

    z = np.zeros(N)
    S_f = np.full(N, 1e-4)
    n_mann = 0.05

    # Start: only centre cell has water
    V_s = np.zeros(N)
    V_s[0] = 1000.0   # 1000 m³ in centre
    M0 = np.sum(V_s)

    min_depth_global = np.inf
    for step in range(n_steps):
        h = V_s / A_hex
        H = h + z
        Hi = H[edge_i]; Hj = H[edge_j]
        D_i = (1/n_mann) * np.maximum(h[edge_i],1e-6)**(5/3) / np.sqrt(S_f[edge_i])
        D_j = (1/n_mann) * np.maximum(h[edge_j],1e-6)**(5/3) / np.sqrt(S_f[edge_j])
        D_f = 0.5*(D_i+D_j)
        Q_vol = D_f * (Hi - Hj) / d_hex * w_hex * dt

        V_req = np.zeros(N)
        np.add.at(V_req, edge_i, np.maximum( Q_vol, 0))
        np.add.at(V_req, edge_j, np.maximum(-Q_vol, 0))
        eps = 1e-12
        lam = np.minimum(1.0, V_s / (V_req + eps))
        lam = np.where(V_req < eps, 1.0, lam)
        Q_f = np.where(Q_vol >= 0, Q_vol * lam[edge_i],
                                   -(-Q_vol) * lam[edge_j])
        dV = np.zeros(N)
        np.add.at(dV, edge_i, -Q_f)
        np.add.at(dV, edge_j,  Q_f)
        V_s += dV
        min_depth_global = min(min_depth_global, np.min(V_s / A_hex))
        V_s = np.maximum(V_s, 0.0)

    M_final = np.sum(V_s)
    abs_err = abs(M_final - M0)
    rel_err = abs_err / M0 * 100
    neighbor_depths = V_s[1:] / A_hex
    symmetry_cv = np.std(neighbor_depths) / (np.mean(neighbor_depths) + 1e-12)

    print(f"  Initial centre volume : {M0:.2f} m³")
    print(f"  Final total volume    : {M_final:.6f} m³")
    print(f"  Relative error        : {rel_err:.4e} %")
    print(f"  Minimum depth reached : {min_depth_global:.6f} m")
    print(f"  Neighbor depth CV     : {symmetry_cv:.4f}  (0=perfect symmetry)")
    neg_ok = "PASS ✓" if min_depth_global >= -1e-10 else "FAIL ✗"
    cons_ok = "PASS ✓" if rel_err < 1e-8 else "FAIL ✗"
    print(f"  Non-negativity        : {neg_ok}")
    print(f"  Conservation          : {cons_ok}")
    return rel_err, min_depth_global


def run_infiltration_mass_balance_test(N_cells=100, hex_radius=150, dt=10, n_steps=400,
                                        P_rate_mmhr=150, Ks_mmhr=50.0,
                                        theta_s=0.40, theta_r=0.08,
                                        vadose_depth=0.1):
    """
    Test C: Full-system infiltration mass balance (new — Sec 6, added
    alongside Fixes #2 and #3).

    Flat synthetic domain (no routing, no outlet, no ET, no groundwater),
    uniform rainfall exceeding the infiltration rate Ks so that every
    cell's vadose zone reaches saturation partway through the run. This
    specifically exercises the vadose-capacity cap introduced in Fix #3:
    once theta -> theta_s, further infiltration must be rejected and the
    corresponding rainfall must remain as surface storage rather than
    being subtracted from V_s and discarded.

    Expected: M_in (cumulative precipitation) = final surface storage +
    final vadose storage (above theta_r), to machine precision. Under the
    pre-fix code this test would FAIL once cells saturate, because excess
    infiltration was silently clipped away from theta without being
    credited back to V_s.
    """
    print("\n" + "=" * 55)
    print("  VERIFICATION TEST C: Infiltration Mass Balance")
    print("=" * 55)

    R = hex_radius
    A_hex = (3 * np.sqrt(3) / 2) * R**2

    N = N_cells
    Ks = Ks_mmhr / 1000 / 3600            # m/s
    P_rate = P_rate_mmhr / 1000 / 3600    # m/s
    P_vol = P_rate * A_hex * dt           # m³/step/cell

    V_s = np.zeros(N)
    theta = np.full(N, theta_r)

    M_in = 0.0
    for step in range(n_steps):
        # Precipitation (source)
        V_s += P_vol
        M_in += P_vol * N

        # Infiltration — mirrors the fixed HexalinkRoutingEngine._infiltration:
        # capped by BOTH available surface water AND remaining vadose capacity.
        h = V_s / A_hex
        I_rate = np.minimum(h / dt, Ks)
        I_vol = np.minimum(I_rate * A_hex * dt, V_s)
        vadose_capacity_vol = np.maximum((theta_s - theta) * A_hex * vadose_depth, 0.0)
        I_vol = np.minimum(I_vol, vadose_capacity_vol)

        V_s -= I_vol
        theta += I_vol / (A_hex * vadose_depth)   # no clip needed: already capped exactly

    V_vadose_final = (theta - theta_r) * A_hex * vadose_depth
    M_final = np.sum(V_s) + np.sum(V_vadose_final)
    abs_err = abs(M_final - M_in)
    rel_err = abs_err / max(M_in, 1e-12) * 100
    n_saturated = int(np.sum(theta >= theta_s - 1e-9))

    print(f"  Cells                  : {N}")
    print(f"  Steps                  : {n_steps}  (dt = {dt} s)")
    print(f"  P rate                 : {P_rate_mmhr} mm/hr  |  Ks = {Ks_mmhr} mm/hr")
    print(f"  Cells reaching theta_s : {n_saturated}/{N}  (exercises vadose-capacity cap)")
    print(f"  Cumulative P (M_in)    : {M_in:.6f} m³")
    print(f"  Final V_s + V_vadose   : {M_final:.6f} m³")
    print(f"  Absolute error         : {abs_err:.4e} m³")
    print(f"  Relative error         : {rel_err:.4e} %")
    status = "PASS ✓" if rel_err < 1e-8 else "FAIL ✗"
    print(f"  Result                 : {status}")
    return rel_err


def run_class_integration_test(dt=5.0, duration_min=30.0, hex_radius=150,
                                precip_rate_mmhr=80.0, tol_pct=1e-6):
    """
    Test D: Real-class integration smoke test (new — closes a test-coverage gap).

    Tests A, B, and C above each reimplement a simplified, standalone version of the
    routing/limiter logic and never instantiate HexalinkRoutingEngine or call its methods
    (_route_surface, _D, run, etc.). Passing all three verifies that the *reference math*
    is self-consistent; it does NOT verify that the actual class is free of bugs. (This is
    not hypothetical: a real indexing bug in HexalinkRoutingEngine._D() previously crashed
    the class on first real use while Tests A/B/C reported near-zero error throughout,
    because none of them ever exercised the buggy code path.)

    This test closes that gap directly: it builds a small synthetic sloped-ramp DEM
    on the fly, instantiates the real HexalinkRoutingEngine against it (via the same
    rasterio-based DEM-loading path used for every real basin run), executes a short
    simulation, and checks that the class (a) runs to completion without error,
    (b) conserves mass to within tolerance, (c) reports zero positivity deficit, and
    (d) produces no negative depths anywhere in the domain.
    """
    import tempfile
    import rasterio as _rasterio
    from rasterio.transform import from_origin as _from_origin

    print("\n" + "=" * 55)
    print("  VERIFICATION TEST D: Real-Class Integration Smoke Test")
    print("=" * 55)

    # Build a small synthetic sloped-ramp DEM (60x60 px @ 30 m = 1800 m x 1800 m),
    # elevation decreasing north->south with a small deterministic cross-slope texture
    # so flow is not perfectly one-dimensional. This gives one unambiguous outlet and
    # a real, multi-cell hexagonal grid with genuine internal edges.
    nrows, ncols, res = 60, 60, 30.0
    y = np.linspace(100.0, 40.0, nrows).reshape(-1, 1)
    dem = np.tile(y, (1, ncols)).astype(np.float32)
    xx, yy = np.meshgrid(np.arange(ncols), np.arange(nrows))
    dem += 0.5 * np.sin(xx / 5.0)
    nodata = -9999.0
    transform = _from_origin(500000.0, 5000000.0, res, res)

    tmp_dir = tempfile.mkdtemp()
    dem_path = os.path.join(tmp_dir, "synthetic_smoke_dem.tif")
    with _rasterio.open(
        dem_path, 'w', driver='GTiff', height=nrows, width=ncols, count=1,
        dtype='float32', crs='EPSG:32632', transform=transform, nodata=nodata,
    ) as dst:
        dst.write(dem, 1)

    params = {
        'n_manning': 0.05, 'Ks': 20.0, 'ET_pot': 2.0,
        'theta_s': 0.40, 'theta_r': 0.08, 'vadose_depth': 0.5, 'alpha_gw': 0.01,
    }

    try:
        engine = HexalinkRoutingEngine(dem_path, params, hex_radius=hex_radius)
        depths, stats = engine.run(precip_rate_mmhr=precip_rate_mmhr, duration_min=duration_min,
                                    dt_sec=dt, output_path=os.path.join(tmp_dir, "smoke_depth.tif"),
                                    print_interval=10**9)
        ran_ok = True
        error_msg = None
    except Exception as exc:
        ran_ok = False
        error_msg = repr(exc)
        depths, stats = None, None

    if not ran_ok:
        print(f"  Class instantiation/run : FAIL \u2717  ({error_msg})")
        print(f"  Result                  : FAIL \u2717")
        return float('inf')

    valid = depths[~np.isnan(depths)]
    has_negative = bool(np.any(valid < -1e-9))
    rel_err = abs(stats['rel_error_pct'])
    deficit = stats['M_deficit']

    print(f"  Cells / edges           : {stats['N_cells']} / {stats['N_edges']}")
    print(f"  Class ran to completion : PASS \u2713")
    print(f"  Relative mass error     : {rel_err:.4e} %")
    print(f"  Positivity deficit      : {deficit:.4e} m\u00b3")
    print(f"  Any negative depth      : {has_negative}")

    passed = ran_ok and (rel_err < tol_pct) and (deficit == 0.0) and (not has_negative)
    status = "PASS \u2713" if passed else "FAIL \u2717"
    print(f"  Result                  : {status}")
    return rel_err


# ─────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Run verification tests first (no DEM needed)
    err_a = run_closed_domain_test()
    err_b, min_d = run_six_neighbor_test()
    err_c = run_infiltration_mass_balance_test()
    err_d = run_class_integration_test()

    print("\n" + "=" * 55)
    print("  VERIFICATION SUMMARY")
    print("=" * 55)
    print(f"  Test A (closed domain)  : rel. error = {err_a:.4e} %")
    print(f"  Test B (6-neighbor)     : rel. error = {err_b[0] if isinstance(err_b,tuple) else err_b:.4e} %")
    print(f"  Test C (infiltration)   : rel. error = {err_c:.4e} %")
    print(f"  Test D (class smoke test): rel. error = {err_d:.4e} %")
    print("=" * 55)

    # Full basin run (requires DEM)
    dem_file = "reprojected elevation 2.tif"
    if os.path.exists(dem_file):
        params = {
            'Ks': 50.0, 'n_manning': 0.05, 'alpha_gw': 0.005,
            'ET_pot': 5.0, 'theta_s': 0.40, 'theta_r': 0.08,
            'vadose_depth': 0.5
        }
        engine = HexalinkRoutingEngine(dem_file, params, hex_radius=150)
        depths, stats = engine.run(
            precip_rate_mmhr=300,
            duration_min=60,
            dt_sec=10,
            output_path="flood_depth_routing.tif"
        )
        print(f"\n  Max surface depth : {np.max(depths):.4f} m")
        print(f"  Mean surface depth: {np.mean(depths[depths>0]):.4f} m")
    else:
        print(f"\n  (DEM not found at '{dem_file}' — skipping basin run)")
        print("  Verification tests above are DEM-independent.")
