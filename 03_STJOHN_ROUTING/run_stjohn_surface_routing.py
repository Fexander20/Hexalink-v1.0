#!/usr/bin/env python3
"""One synthetic, surface-only routing experiment on the St. John hex mesh.

The full nominal area of each hexagon is used for rainfall and storage depth.
The DEM catchment mask area is reported separately; boundary hexagons are not
polygon-clipped. This is a routing demonstration, not an observed event hindcast.
"""
import hashlib, json, csv, math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from routing_engine import HexalinkRoutingEngine
ROOT=Path(__file__).resolve().parent
RADIUS=300.0; DT=10.0; NSTEP=360; STORM_STEPS=180; RAIN_MMHR=300.0
AHEX=3*np.sqrt(3)/2*RADIUS**2
X_CENTER=484958.6373624129; Y_CENTER=802665.9038581452
PATCH_RADIUS_M=12000.0

def load_mesh():
    with np.load(ROOT/'hexalink_stjohn_config_300m.npz',allow_pickle=False) as f:
        x=f['x_m'].astype(float);y=f['y_m'].astype(float)
        z=f['elevation_merit_m'].astype(float)
        ei=f['edge_i'].astype(np.int64);ej=f['edge_j'].astype(np.int64)
        active=f['active_area_m2'].astype(float)
    n=len(x)
    assert n==68292 and len(ei)==203076 and len(ej)==203076
    count=np.bincount(ei,minlength=n)+np.bincount(ej,minlength=n)
    ad=np.abs(z[ei]-z[ej])
    sumdz=np.bincount(ei,weights=ad,minlength=n)+np.bincount(ej,weights=ad,minlength=n)
    sf=np.maximum(np.divide(sumdz,count,out=np.zeros(n),where=count>0)/(math.sqrt(3)*RADIUS),1e-4)
    boundary=np.flatnonzero(count<6)
    outlet=int(boundary[np.argmin(z[boundary])])
    patch=((x-X_CENTER)**2+(y-Y_CENTER)**2<=PATCH_RADIUS_M**2)
    assert patch.sum()==1945
    return x,y,z,ei,ej,sf,active,outlet,patch

def main():
    x,y,z,ei,ej,sf,active,outlet,patch=load_mesh()
    e=HexalinkRoutingEngine.__new__(HexalinkRoutingEngine)
    e.N=len(x);e.A_hex=AHEX;e.d_hex=math.sqrt(3)*RADIUS;e.w_hex=RADIUS
    e.z=z;e.S_f=sf;e.edge_i=ei;e.edge_j=ej
    e.outlet=outlet;e.n_mann=0.05;e.V_s=np.zeros(e.N,float)
    rain=np.where(patch,RAIN_MMHR/1000/3600*DT*AHEX,0.)
    total=0.; out=0.;check=[]
    for step in range(1,NSTEP+1):
        if step<=STORM_STEPS:
            e.V_s+=rain;total+=float(rain.sum())
        dv,q=e._route_surface(DT)
        e.V_s+=dv;out+=q
        if step in (STORM_STEPS,NSTEP):
            d=e.V_s/AHEX
            check.append({'elapsed_min':step*DT/60,'input_m3':total,'storage_m3':float(e.V_s.sum()),'outlet_m3':out,
            'mass_closure_fraction':float((e.V_s.sum()+out-total)/total),
            'wet_cells_beyond_footprint_gt_0p1mm':int(np.count_nonzero((d>1e-4)&~patch)),
            'max_depth_m':float(d.max())})
    d=e.V_s/AHEX
    result={'case':'Synthetic localized rainfall; nominal full hex geometry',
        'cells':len(d),'edges':len(ei),'radius_m':RADIUS,'nominal_area_m2':AHEX,
        'mesh_domain_nominal_area_km2':len(d)*AHEX/1e6,
        'independent_dem_mask_area_km2':float(active.sum()/1e6),
        'rainfall_patch_cells':int(patch.sum()),
        'rainfall_patch_nominal_area_km2':float(patch.sum()*AHEX/1e6),
        'rainfall_input_m3':total,'final_storage_m3':float(e.V_s.sum()),
        'outlet_m3':out,'mass_closure_fraction':float((e.V_s.sum()+out-total)/total),
        'wet_cells_outside_patch_gt_0p1mm':int(np.count_nonzero((d>1e-4)&~patch)),
        'max_depth_m':float(d.max()),'min_storage_m3':float(e.V_s.min()),
        'checkpoints':check,'engine_sha256':hashlib.sha256((ROOT.parent/'routing_engine.py').read_bytes()).hexdigest()}
    assert abs(result['mass_closure_fraction'])<1e-10
    assert result['min_storage_m3']>=-1e-7
    np.savez_compressed(ROOT/'stjohn_routing_depths.npz',x_m=x,y_m=y,depth_m=d.astype('float32'),rain_patch=patch.astype('uint8'))
    (ROOT/'stjohn_routing_results.json').write_text(json.dumps(result,indent=2))
    with (ROOT/'stjohn_routing_checkpoints.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(check[0]));w.writeheader();w.writerows(check)
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
