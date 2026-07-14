#!/usr/bin/env python3
"""Extract exact metrics from Hexalink outputs for Table 4 (Manuscript Sec 6.5)"""
import numpy as np, rasterio, os, json, glob
from datetime import datetime

def extract(tiff, name, method):
    print(f"\n📊 Processing: {name} | {tiff}")
    if not os.path.exists(tiff): print("❌ Not found"); return None
    with rasterio.open(tiff) as src:
        d=src.read(1).astype(np.float32); d[d==src.nodata]=np.nan
        mask=~np.isnan(d); v=d[mask]
        if len(v)==0: print("❌ No valid data"); return None
        area_m2=abs(src.res[0]*src.res[1])
        if src.crs and src.crs.is_projected: area_m2=abs(src.res[0]*src.res[1])
        else: lat=(src.bounds.bottom+src.bounds.top)/2; area_m2=abs(src.res[0]*111320*src.res[1]*111320*np.cos(np.radians(lat)))
        return {'basin':name,'method':method,'valid_px':int(mask.sum()),'coverage_%':float(mask.sum()/d.size*100),
                'min_m':float(v.min()),'max_m':float(v.max()),'mean_m':float(v.mean()),'std_m':float(v.std()),
                'area_km2':float(mask.sum()*area_m2/1e6),'timestamp':datetime.now().isoformat()}

if __name__=="__main__":
    files=[("imo_flood_depth.tif","Imo River Basin","Topographic accumulation")]
    res=[extract(f,n,m) for f,n,m in files if os.path.exists(f)]
    if res:
        with open('hexalink_metrics.json','w') as f: json.dump(res,f,indent=2)
        print("\n💾 Saved: hexalink_metrics.json")
        for r in res: print(f"{r['basin']:<15} Area:{r['area_km2']:<8.1f}km² | Cells:{r['valid_px']:<10,} | Max:{r['max_m']:.4f}m | Cov:{r['coverage_%']:.2f}%")