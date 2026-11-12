"""Synthetic, local-only regression tests for the screening-pathway correction."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

from run_imo_scenarios import create_hexagonal_grid, run_imo_scenarios, screening_depth_by_hex


class ScreeningTests(unittest.TestCase):
    def test_area_weighted_routing_is_exact(self):
        depths, norm = screening_depth_by_hex(
            np.array([np.nan, 1, 2, 3, 4]),
            np.array([0, 1, 2, 3, 4]), 0.075,
        )
        weighted = np.dot(depths[1:], np.array([1, 2, 3, 4])) / 10
        self.assertAlmostEqual(weighted, 0.075, places=14)
        self.assertGreater(norm, 0)
        self.assertTrue(np.isnan(depths[0]))

    def test_zero_runoff_has_no_display_floor(self):
        depths, _ = screening_depth_by_hex(
            np.array([1., 3., 6.]), np.array([1., 10., 3.]), 0,
        )
        self.assertTrue(np.array_equal(depths, [0., 0., 0.]))

    def test_end_to_end_synthetic_dem_nodata(self):
        with tempfile.TemporaryDirectory() as t:
            folder = Path(t)
            path = folder / 'dem.tif'
            n = 120
            dem = np.linspace(1, 800, n*n, dtype=np.float64).reshape(n, n)
            dem[0:5, 0:10] = -9999.
            with rasterio.open(path, 'w', driver='GTiff', height=n, width=n,
                                count=1, dtype='float64', nodata=-9999.,
                                crs='EPSG:32632', transform=from_origin(400000., 900000., 30., 30.)) as dst:
                dst.write(dem, 1)
            result = run_imo_scenarios(path, output_dir=folder / 'result')
            self.assertEqual(len(result), 4)
            for i, r in enumerate(result):
                self.assertEqual(r['valid_pixels'], n*n - 50)
                self.assertLess(abs(r['volume_residual_m3']), 1e-5)
                self.assertLess(r['relative_volume_error_percent'], 1e-9)
                with rasterio.open(r['output']) as tif:
                    vals = tif.read(1)
                    self.assertEqual(tif.count, 1)
                    self.assertTrue(np.all(np.isnan(vals[0:5, 0:10])))
                    if i < 2:
                        self.assertTrue(np.all(vals[np.isfinite(vals)] == 0))
                    else:
                        self.assertTrue(np.all(vals[np.isfinite(vals)] >= 0))
            self.assertAlmostEqual(result[2]['area_weighted_mean_depth_m'], 0.075, places=12)
            self.assertAlmostEqual(result[3]['area_weighted_mean_depth_m'], 0.175, places=12)
            self.assertTrue((folder/'result'/'screening_metrics.csv').exists())

    def test_integer_dem_nodata_support(self):
        with tempfile.TemporaryDirectory() as t:
            folder = Path(t)
            path = folder / 'int16_dem.tif'
            values = np.arange(60 * 60, dtype=np.int16).reshape(60, 60)
            values[:3, :5] = -9999
            with rasterio.open(path, 'w', driver='GTiff', height=60, width=60,
                                count=1, dtype='int16', nodata=-9999,
                                crs='EPSG:32632', transform=from_origin(400000., 900000., 30., 30.)) as dst:
                dst.write(values, 1)
            results = run_imo_scenarios(path, output_dir=folder / 'output')
            self.assertEqual(results[2]['valid_pixels'], 60*60-15)
            self.assertAlmostEqual(results[3]['area_weighted_mean_depth_m'], 0.175, places=12)

    def test_grid_regular_hex_coverage(self):
        polys, _ = create_hexagonal_grid((0, 0, 6000, 6000))
        self.assertGreater(len(polys), 500)
        from rasterio.features import rasterize
        ids = rasterize(((p, i+1) for i,p in enumerate(polys)), out_shape=(200,200),
                        transform=from_origin(0,6000,30,30), fill=0, dtype='int32')
        self.assertFalse(np.any(ids == 0))


if __name__ == '__main__':
    unittest.main()
