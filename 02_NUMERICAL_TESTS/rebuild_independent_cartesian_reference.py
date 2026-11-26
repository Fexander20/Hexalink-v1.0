#!/usr/bin/env python3
"""Standalone reconstruction of the independent 2-m Cartesian slope reference.

This script uses square-grid fluxes with arithmetic-midpoint mobility,
NOT the hexagonal secant-mobility solver; the reference is therefore genuinely
separate in discretization and implementation. It regenerates only the
reference depth sampled on the fixed 10-m hexagonal cells (the raw reference
array included in the research package).

Usage: python rebuild_independent_cartesian_reference.py
The Cartesian computation is more expensive than the hexagonal checks.
"""
from pathlib import Path
from math import pi, sqrt
import numpy as np
from scipy.ndimage import map_coordinates

ROOT=Path(__file__).resolve().parent
R=10.; D=335.; DX=2.; DT=.008; STOP=32.
N_MANN=.05; SF=.01; K=1/(N_MANN*sqrt(SF)); P=5./3.; X0=-75.

def h0(x,y):
    return .5*np.maximum(1.-((x-X0)**2+y*y)/70.**2,0.)**.6

def make_hex_xy():
    Q=int(np.ceil(2*D/(1.5*R)))+3
    xy=[]
    for q in range(-Q,Q+1):
        for r in range(-Q,Q+1):
            x=1.5*R*q;y=sqrt(3)*R*(r+q/2)
            if x*x+y*y<D*D:xy.append((x,y))
    return np.asarray(xy)

def hex_avg(xy,func):
    vs=R*np.c_[np.cos(np.arange(6)*pi/3),np.sin(np.arange(6)*pi/3)]
    a=np.zeros(len(xy))
    for j in range(6):
        u,v=vs[j],vs[(j+1)%6]
        for p in (u/6+v/6,2*u/3+v/6,u/6+2*v/3):
            a+=func(xy[:,0]+p[0],xy[:,1]+p[1])/18
    return a

def main():
    extent=360.
    n=round(2*extent/DX)
    x=-extent+(np.arange(n)+.5)*DX
    X,Y=np.meshgrid(x,x,indexing='xy')
    z=-.025*X
    h=h0(X,Y)
    m0=h.sum()*DX*DX
    nsteps=round(STOP/DT)
    assert abs(nsteps*DT-STOP)<1e-8
    for _ in range(nsteps):
        head=z+h
        fx=K*(.5*(h[:,:-1]+h[:,1:]))**P*(head[:,:-1]-head[:,1:])/DX
        fy=K*(.5*(h[:-1,:]+h[1:,:]))**P*(head[:-1,:]-head[1:,:])/DX
        req=np.zeros_like(h)
        req[:,:-1]+=DT/DX*np.maximum(fx,0.); req[:,1:]+=DT/DX*np.maximum(-fx,0.)
        req[:-1,:]+=DT/DX*np.maximum(fy,0.); req[1:,:]+=DT/DX*np.maximum(-fy,0.)
        lam=np.minimum(1.,h/(req+1e-14)); lam=np.where(req<1e-14,1,lam)
        fx=np.where(fx>=0,fx*lam[:,:-1],fx*lam[:,1:])
        fy=np.where(fy>=0,fy*lam[:-1,:],fy*lam[1:,:])
        delta=np.zeros_like(h)
        delta[:,:-1]-=fx*DT/DX;delta[:,1:]+=fx*DT/DX
        delta[:-1,:]-=fy*DT/DX;delta[1:,:]+=fy*DT/DX
        h+=delta
        if h.min()<-1e-10:raise RuntimeError(f'negative Cartesian depth {h.min()}')
        h[h<0]=0
    xy=make_hex_xy()
    def interp(xi,yi):
        coords=np.array([(yi-x[0])/DX,(xi-x[0])/DX])
        return map_coordinates(h,coords,order=1,mode='constant',cval=0.)
    data=hex_avg(xy,interp)
    np.savez_compressed(ROOT/'cartesian_reference_regenerated.npz',xy=xy,reference=data)
    saved=np.load(ROOT/'independent_cartesian_reference_on_hex.npz')
    assert np.array_equal(saved['xy'],xy), 'hex point ordering mismatch'
    diff=np.max(np.abs(saved['reference']-data))
    print({'dx_m':DX,'dt_s':DT,'steps':nsteps,'initial_mass_m3':float(m0),
           'final_mass_m3':float(h.sum()*DX*DX),'maximum_difference_vs_reference_m':float(diff)})
    assert diff < 1e-9,'Cartesian reference does not reproduce'

if __name__=='__main__':main()
