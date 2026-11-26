#!/usr/bin/env python3
"""Compare independently observed GFD event water against Hexalink screening classes.

GFD export coding: 0 observed dry, 1 observed flood, 255 unavailable/permanent water.
Valid event pixels only; compute at GFD export grid (250 m), never inflate the
sample count by upsampling satellite pixels to the 100 m screening map.
Screening is an ordinal 1/2/4 basin-wide elevation-class score, not depth.
"""
import argparse, csv, json, hashlib
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling


def score(labels, weights):
    y=np.asarray(labels, dtype=np.uint8)
    w=np.asarray(weights, dtype=np.uint8)
    if y.ndim!=1 or w.shape!=y.shape or not np.isin(y,[0,1]).all():
        raise ValueError('Invalid labels or scores')
    if len(y)<20 or y.sum()<10 or (y==0).sum()<10:
        raise ValueError('Fewer than 10 flood OR 10 dry observed pixels; not a defensible classification sample')
    if not np.isin(w,[1,2,4]).all():
        raise ValueError('Expected released screening weights 1, 2, or 4')
    from sklearn.metrics import roc_auc_score, average_precision_score
    d={'valid_observed_pixels':int(y.size),'observed_flood_pixels':int(y.sum()),
       'observed_dry_pixels':int((y==0).sum()),'observed_flood_fraction':float(y.mean()),
       'ROC_AUC':float(roc_auc_score(y,w)),
       'average_precision':float(average_precision_score(y,w)),
       'random_AP_reference':float(y.mean()),'classes': []}
    for val in (1,2,4):
        select=(w==val)
        d['classes'].append({'screening_class':val,'pixels':int(select.sum()),
         'fraction_of_valid_observation_support':float(select.mean()),
         'flooded_pixels':int(np.sum(select & (y==1))),
         'observed_flood_rate':float(y[select].mean()) if select.any() else None})
    for name, selected in (('high_class_only_4',w==4),('medium_or_high_2_and_4',w>=2)):
        tp=int((selected & (y==1)).sum()); fp=int((selected & (y==0)).sum())
        fn=int((~selected & (y==1)).sum()); tn=int((~selected & (y==0)).sum())
        d[name]={'actual_selected_area_fraction':float(selected.mean()),
                 'TP':tp,'FP':fp,'FN':fn,'TN':tn,
                 'recall':tp/(tp+fn),'precision':tp/(tp+fp) if tp+fp else None,
                 'false_positive_rate':fp/(fp+tn) if fp+tn else None}
    return d


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--observed',type=Path,required=True)
    p.add_argument('--screening',type=Path,default=Path('ScreeningReference/StJohn_screening_terrain_weight_class.tif'))
    p.add_argument('--output',type=Path,default=Path('Results/StJohn_GFD_validation.json'))
    p.add_argument('--event-reference',required=True,help='Real GFD event ID/date copied from Earth Engine')
    args=p.parse_args()
    if not args.observed.exists():raise SystemExit(f'Missing observed flood raster: {args.observed}')
    with rasterio.open(args.observed) as obs:
        if obs.crs is None:raise ValueError('Observed raster has no CRS')
        if obs.count!=1:raise ValueError('Observed product must be the single-band 0/1/255 export')
        if not obs.crs.is_projected:raise ValueError('Export should be projected in EPSG:32629 at 250 m')
        if max(obs.res)>300 or min(obs.res)<200:raise ValueError('Expected approx. 250 m GFD export resolution; inspect raster')
        observed=obs.read(1)
        transform,crs=obs.transform,obs.crs
        grid_shape=observed.shape
        valid_obs=np.isin(observed,[0,1])
        if not (observed==255).any():print('NOTE: No nodata/permanent water pixels coded 255 in this export')
        with rasterio.open(args.screening) as scr:
            if scr.crs is None:raise ValueError('Screening raster has no CRS')
            classes=np.zeros(grid_shape,dtype=np.uint8)
            reproject(source=rasterio.band(scr,1),destination=classes,
                      src_transform=scr.transform,src_crs=scr.crs,src_nodata=0,
                      dst_transform=transform,dst_crs=crs,dst_nodata=0,
                      resampling=Resampling.nearest)
            # Guard model boundary: almost all 250 m cell must fall inside 100 m domain
            valid_src=(scr.read(1)!=0).astype(np.uint8)
            coverage=np.zeros(grid_shape,dtype=np.float32)
            reproject(source=valid_src,destination=coverage,
                      src_transform=scr.transform,src_crs=scr.crs,src_nodata=None,
                      dst_transform=transform,dst_crs=crs,dst_nodata=None,
                      resampling=Resampling.average)
            domain= (coverage>=0.95) & (classes>0)
    n_domain=int(domain.sum());paired=domain & valid_obs
    if n_domain==0:raise ValueError('NO spatial overlap between event raster and model domain')
    n_pair=int(paired.sum()); f=int((observed[paired]==1).sum())
    context={'status':'EVALUATED','source':'Global Flood Database (Tellman et al. 2021), GEE GLOBAL_FLOOD_DB/MODIS_EVENTS/V1',
        'event_reference':args.event_reference, 'observation_legend':{'0':'clear non-water/dry','1':'observed flood, excluding permanent water','255':'unobserved or permanent water'},
        'geographic_area':'St John basin through Mt Finlay, Liberia',
        'evaluation_grid_crs':str(crs),'evaluation_grid_resolution_m':list(map(float, (abs(transform.a),abs(transform.e)))),
        'domain_eligible_250m_pixels':n_domain,'paired_eligible_pixels':n_pair,
        'paired_domain_fraction':n_pair/n_domain,
        'paired_flood_pixels':f,
        'source_screening':str(args.screening),'source_observed':str(args.observed),
        'interpretation':'Ordinal terrain-based flood-priority discrimination against mapped event flood (not water depth or dynamic routing). The screening classes are derived entirely from basin-wide elevation percentiles.',
        'limitations':['GFD satellite event is composite/maximum flood, possibly incomplete and cloud affected',
        '250 m satellite observations are not independent 100 m training samples',
        'Missing/permanent water pixels excluded, dry means satellite-observed dry within event footprint',
        'Event occurrence need not coincide with synthetic one-hour rainfall scenarios',
        'Current screening rank is elevation-percentile rank; elevation comparison is tautological']}
    if n_pair<20 or f<10 or n_pair-f<10:
        context['status']='INSUFFICIENT_OBSERVATION_SUPPORT'
        context['next_action']='Do not report classification skill; select a qualifying mapped event or alternate satellite dataset.'
    else:
        context.update(score((observed[paired]==1).astype(np.uint8),classes[paired]))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(context,indent=2),encoding='utf8')
    if context['status']=='EVALUATED':
        with args.output.with_suffix('.csv').open('w',newline='',encoding='utf8') as fout:
            writer=csv.DictWriter(fout,fieldnames=context['classes'][0].keys());writer.writeheader();writer.writerows(context['classes'])
    print('STATUS:',context['status'])
    print('VALID GFD PIXELS:',n_pair,'FLOOD:',f,'DOMAIN COVERAGE:',f'{100*n_pair/n_domain:.1f}%')
    if context['status']=='EVALUATED':print('AUC:',context['ROC_AUC'],'AP:',context['average_precision'])
    print('REPORT:',args.output)

if __name__=='__main__':main()
