"""Area-conservative St. John terrain-priority screening on the 300 m hex mesh.
The method ranks elevations and conserves prescribed effective-runoff volume;
its screening-equivalent depth is not simulated hydraulic inundation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
import rasterio
from rasterio.features import rasterize
from shapely.geometry import Polygon
from screening_core import screening_depth_by_hex, SCENARIOS


def sha256(p):
    d=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            d.update(b)
    return d.hexdigest()


def run(dem_path: Path, mesh_path: Path, output_dir: Path, radius: float=300,
        infiltration_mmhr: float=120, et_mmhr: float=5, duration_min: float=60):
    output_dir.mkdir(parents=True, exist_ok=True)
    if radius <= 0 or infiltration_mmhr < 0 or et_mmhr < 0 or duration_min <= 0:
        raise ValueError('Invalid physical inputs')
    mesh=np.load(mesh_path)
    for key in ('cell_id','x_m','y_m','active_area_m2'):
        if key not in mesh.files: raise ValueError(f'Missing mesh array {key}')
    ids=mesh['cell_id']
    xs=mesh['x_m'].astype(np.float64)
    ys=mesh['y_m'].astype(np.float64)
    area=mesh['active_area_m2'].astype(np.float64)
    if len(np.unique(ids))!=len(ids) or not np.array_equal(ids,np.arange(1,len(ids)+1)):
        raise ValueError('Expected contiguous 1-based cell IDs')
    with rasterio.open(dem_path) as src:
        if src.crs.to_epsg()!=32629: raise ValueError(f'Expected EPSG:32629, got {src.crs}')
        dem=src.read(1,masked=True).astype(np.float64).filled(np.nan)
        prof=src.profile.copy()
        tr=src.transform
        pxarea=abs(tr.a*tr.e-tr.b*tr.d)
        valid=np.isfinite(dem)
        if pxarea != 10000.0: raise ValueError(f'Expected 100 m pixels; got area {pxarea}')
        cell_area=valid.sum()*pxarea
    # Build hex geometry from the exact preexisting mesh centres.
    theta=np.arange(6)*np.pi/3
    offsets=np.column_stack([radius*np.cos(theta),radius*np.sin(theta)])
    polygons=(
        (Polygon(np.column_stack([x+offsets[:,0], y+offsets[:,1]])),int(cell_id))
        for x,y,cell_id in zip(xs,ys,ids)
    )
    labels=rasterize(polygons, out_shape=dem.shape,transform=tr,fill=0,
                     dtype=np.int32,all_touched=False)
    uncovered=int(np.count_nonzero(valid & (labels==0)))
    if uncovered: raise AssertionError(f'{uncovered} valid DEM pixels not in mesh hex footprint')
    count=np.bincount(labels[valid],minlength=len(ids)+1).astype(np.float64)
    summed=np.bincount(labels[valid],weights=dem[valid],minlength=len(ids)+1)
    observed_area=count[1:]*pxarea
    mismatch=observed_area-area
    bad=int(np.count_nonzero(mismatch))
    # Domain-area check is exact even if individual pixel-to-cell assignments differ.
    if not np.isclose(cell_area,area.sum(),rtol=0,atol=1e-3):
        raise AssertionError(f'Mesh/DEM domain support disagreement: {cell_area} vs {area.sum()}')
    elev=np.full(len(ids)+1,np.nan,dtype=np.float64)
    present=count>0
    elev[present]=summed[present]/count[present]
    # Screening weights are invariant across positive rain scenarios; calculate
    # from any arbitrary *unit* depth to retrieve categories; not a new formula.
    unit, norm=screening_depth_by_hex(elev,count,1.0)
    weights=np.full(len(ids)+1,np.nan,dtype=np.float64)
    weights[present]=unit[present]*norm
    wimg=np.full(dem.shape,np.nan,dtype=np.float64)
    wimg[valid]=weights[labels[valid]]
    if not set(np.unique(wimg[valid])).issubset({1.0,2.0,4.0}):
        raise AssertionError('Unexpected screening class')
    metrics=[]
    output_base=prof.copy()
    output_base.update(driver='GTiff',dtype='float64',count=1,nodata=np.nan,
                       compress='deflate',predictor=3,tiled=True)
    for name,rain in SCENARIOS.items():
        runoff_mm=max(0.0,(rain-infiltration_mmhr-et_mmhr)*duration_min/60.0)
        depth, normalizer=screening_depth_by_hex(elev,count,runoff_mm/1000.0)
        image=np.full(dem.shape,np.nan,dtype=np.float64)
        image[valid]=depth[labels[valid]]
        out=output_dir/f'StJohn_screening_{name}.tif'
        with rasterio.open(out,'w',**output_base) as ds:
            ds.write(image,1)
            ds.set_band_description(1,'Equivalent algebraic screening depth (m); not inundation')
            ds.update_tags(source='Hexalink area-conservative terrain screening',
                           hydraulic_depth='false',geometry='existing StJohn 300m hex mesh',
                           rainfall_intensity_mmhr=str(rain),
                           infiltration_rate_mmhr=str(infiltration_mmhr),
                           illustrative_et_rate_mmhr=str(et_mmhr))
        expected=runoff_mm/1000.0*cell_area
        actual=float(np.sum(image[valid],dtype=np.float64)*pxarea)
        error=actual-expected
        if not np.isclose(expected,actual,rtol=1e-12,atol=1e-6):
            raise AssertionError(f'Area balance failure: {name} {error}')
        # Truly independent volume readback of exported map
        with rasterio.open(out) as ds:
            check=ds.read(1,masked=True).filled(np.nan)
        check_vol=float(np.nansum(check)*pxarea)
        if not np.isclose(expected,check_vol,rtol=1e-12,atol=1e-6):
            raise AssertionError(f'GeoTIFF readback balance failure: {name}')
        met=dict(scenario=name,rainfall_mmhr=float(rain),runoff_mm=runoff_mm,
                 model_mesh_radius_m=radius,valid_pixels=int(valid.sum()),
                 area_km2=cell_area/1e6,mean_screening_depth_m=float(np.mean(image[valid])),
                 min_screening_depth_m=float(np.min(image[valid])),
                 max_screening_depth_m=float(np.max(image[valid])),
                 expected_volume_m3=expected,mapped_volume_m3=actual,
                 file_readback_volume_m3=check_vol,
                 volume_residual_m3=error,
                 rel_error_percent=(100*abs(error)/expected if expected else 0.0),
                 normalization=normalizer,filename=out.name)
        metrics.append(met)
        print(f'{name}: runoff {runoff_mm:.1f} mm | mapped {actual/1e6:.3f} million m3 | balance {error:.3e} m3',flush=True)
    # categorical 1/2/4 independent of scenario magnitude
    categorical=np.zeros(dem.shape,dtype=np.uint8)
    categorical[valid]=wimg[valid].astype(np.uint8)
    profcat=prof.copy()
    profcat.update(driver='GTiff',dtype='uint8',count=1,nodata=0,
                   compress='deflate',tiled=True)
    with rasterio.open(output_dir/'StJohn_screening_terrain_weight_class.tif','w',**profcat) as ds:
        ds.write(categorical,1)
        ds.set_band_description(1,'Global relative-elevation screening weight class 1,2,4')
        ds.update_tags(source='Hexalink terrain-priority 1/2/4 screening',
                       NOT_A_HYDRAULIC_FLOOD_MAP='true',
                       ranking='1/2/4 from absolute basin-wide hex mean elevation percentiles')
    fractions={str(v):float(np.mean(categorical[valid]==v)) for v in (1,2,4)}
    report={
      'status':'COMPLETED_ALGEBRAIC_SCREENING_NOT_OBSERVATIONAL_VALIDATION',
      'dem_sha256':sha256(dem_path),'mesh_sha256':sha256(mesh_path),
      'screening_source_sha256':sha256(Path(__file__).parent/'screening_core.py'),
      'crs_epsg':32629,'model_radius_m':radius,'DEM_pixel_resolution_m':100,
      'mesh_cells':int(len(ids)),'mesh_active_area_km2':float(area.sum()/1e6),
      'DEM_valid_pixel_area_km2':float(cell_area/1e6),
      'uncovered_DEM_pixels':uncovered,
      'cell_count_areas_mismatch':bad,'max_abs_cell_assignment_area_difference_m2':float(np.max(abs(mismatch))),
      'weight_class_area_fractions':fractions,
      'limitations':[
       'Equivalent screening depths are algebraic allocations, not inundation depths.',
       'Ranking reflects basin-wide relative elevations, not D8 flow accumulation or local relative topographic position.',
       'Observational comparison is a separate single-event exercise; no hydraulic inundation prediction is claimed.',
       'Prescribed illustrative infiltration and ET loss rates are not event-calibrated.',
      ],'scenario_metrics':metrics}
    with open(output_dir/'screening_metrics.csv','w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(metrics[0]));w.writeheader();w.writerows(metrics)
    (output_dir/'screening_qa.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--dem',type=Path,required=True)
    p.add_argument('--mesh',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    run(args.dem,args.mesh,args.out)
