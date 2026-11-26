#!/usr/bin/env python3
"""Numerical verification of the Hexalink finite-volume surface solver.
Includes a nonlinear-diffusion analytical solution, an independent sloping-
terrain Cartesian reference, timestep refinement and lake-at-rest equilibrium.
"""
from pathlib import Path
import json, csv, hashlib, sys, os, platform, importlib.metadata as md
from math import pi, sqrt
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from routing_engine import HexalinkRoutingEngine

ROOT=Path(__file__).resolve().parent
REF=ROOT/'independent_cartesian_reference_on_hex.npz'
sha=hashlib.sha256((ROOT.parent/'routing_engine.py').read_bytes()).hexdigest()
M=8./3.; P=1./(M-1)
N_MANN=.05
SF_FLAT=1e-4
SF_SLOPE=.01


def regular_hex_mesh(radius,domain_radius):
    qr=int(np.ceil(2*domain_radius/(1.5*radius)))+3
    centres=[]; lookup={}
    for q in range(-qr,qr+1):
        x=1.5*radius*q
        for r in range(-qr,qr+1):
            y=sqrt(3)*radius*(r+q/2)
            if (x*x+y*y <= domain_radius**2 if domain_radius==550. else x*x+y*y < domain_radius**2):
                lookup[(q,r)]=len(centres)
                centres.append((x,y))
    edges=[]
    for (q,r),i in lookup.items():
        for v in ((q+1,r),(q,r+1),(q+1,r-1)):
            if v in lookup: edges.append((i,lookup[v]))
    return np.asarray(centres,float),np.asarray(edges,np.int64)


def cell_averages(xy,R,fn):
    # Deterministic 18-point triangular quadrature for initial cell averages.
    vv=R*np.c_[np.cos(np.arange(6)*pi/3),np.sin(np.arange(6)*pi/3)]
    val=np.zeros(xy.shape[0])
    for i in range(6):
        a,b=vv[i],vv[(i+1)%6]
        for pt in (a/6+b/6,2*a/3+b/6,a/6+2*b/3):
            val += fn(xy[:,0]+pt[0],xy[:,1]+pt[1])/18
    return val


def create_engine(R,xy,edges,depth,z,Sf):
    e=HexalinkRoutingEngine.__new__(HexalinkRoutingEngine)
    e.N=len(xy);e.A_hex=3*sqrt(3)/2*R**2
    e.w_hex=R;e.d_hex=sqrt(3)*R
    e.edge_i=edges[:,0];e.edge_j=edges[:,1]
    e.z=np.asarray(z,dtype=np.float64)
    e.S_f=np.asarray(Sf,dtype=np.float64)
    e.outlet=int(np.argmax(xy[:,0]**2+xy[:,1]**2))
    e.n_mann=N_MANN
    e.V_s=np.asarray(depth,dtype=np.float64)*e.A_hex
    return e


def frozen_coefficient_tau(e):
    """Diagnostic A/sum(Tface), not a nonlinear/CFL stability theorem.
    Applies to explicit Euler with fixed Tface and no source/active outlet.
    """
    h=np.maximum(e.V_s/e.A_hex,0.)
    ei,ej=e.edge_i,e.edge_j
    hi,hj=h[ei],h[ej]
    K_i=1/(e.n_mann*np.sqrt(e.S_f[ei])); K_j=1/(e.n_mann*np.sqrt(e.S_f[ej]))
    Kface=2*K_i*K_j/(K_i+K_j)
    diff=hi-hj
    near=np.abs(diff)<1e-8*np.maximum(np.maximum(hi,hj),1e-5)
    mob=np.empty_like(diff)
    mob[near]=(.5*(hi[near]+hj[near]))**(M-1)
    normal=~near
    mob[normal]=(hi[normal]**M-hj[normal]**M)/(M*diff[normal])
    T=Kface*mob*(e.w_hex/e.d_hex) # m^3/s / m head; units m^2/s
    tsum=np.bincount(ei,weights=T,minlength=e.N)+np.bincount(ej,weights=T,minlength=e.N)
    return float(e.A_hex/np.max(tsum)), float(np.max(tsum)), int(np.argmax(tsum))


DIM=2; ALPHA=DIM/(DIM*(M-1)+2); BETA=1/(DIM*(M-1)+2)
T0=10.;TEND=40.;H0=.5
KFLAT=1/(N_MANN*sqrt(SF_FLAT))
R0=sqrt((2*KFLAT*DIM/(ALPHA*(M-1)))*H0**(M-1)*T0)
EXACT_MASS=pi*R0**2*H0/(P+1)


def exact_barenblatt(x,y,t):
    u=t/T0
    rr=R0*u**BETA
    return H0*u**(-ALPHA)*np.maximum(1-(x*x+y*y)/rr**2,0.)**P


def run_barenblatt():
    rows=[]; profiles=None
    for R,dt in [(40.,.2),(20.,.05),(10.,.02)]:
        xy,ed=regular_hex_mesh(R,550.)
        h0=cell_averages(xy,R,lambda x,y:exact_barenblatt(x,y,T0))
        eng=create_engine(R,xy,ed,h0,np.zeros(len(xy)),np.full(len(xy),SF_FLAT))
        init_mass=float(eng.V_s.sum())
        start_tau,_,_=frozen_coefficient_tau(eng)
        out_total=0.
        nsteps=round((TEND-T0)/dt)
        assert np.isclose(nsteps*dt,TEND-T0)
        for k in range(nsteps):
            dv,o=eng._route_surface(dt)
            eng.V_s+=dv
            out_total+=o
        final_tau,_,_=frozen_coefficient_tau(eng)
        exact=cell_averages(xy,R,lambda x,y:exact_barenblatt(x,y,TEND))
        actual=eng.V_s/eng.A_hex
        l1=float(np.abs(actual-exact).sum()*eng.A_hex/EXACT_MASS*100)
        row={'radius_m':R,'dt_s':dt,'cells':eng.N,'edges':len(ed),'steps':nsteps,
             'relative_L1_error_pct_exact_mass':l1,
             'initial_discretization_mass_bias_pct':(init_mass/EXACT_MASS-1)*100,
             'mass_closure_fraction':(eng.V_s.sum()+out_total-init_mass)/init_mass,
             'minimum_depth_m':float(actual.min()),'outlet_volume_m3':float(out_total),
             'frozen_coefficient_tau_initial_s':start_tau,'frozen_coefficient_tau_final_s':final_tau,
             'dt_over_initial_tau':dt/start_tau,'dt_over_final_tau':dt/final_tau}
        rows.append(row)
        print('BAR',row,flush=True)
        assert abs(row['mass_closure_fraction'])<1e-11
        if R==10: profiles={'xy':xy,'numerical':actual,'exact':exact,'initial':h0}
    assert rows[0]['relative_L1_error_pct_exact_mass']>rows[1]['relative_L1_error_pct_exact_mass']>rows[2]['relative_L1_error_pct_exact_mass'], 'Barenblatt mesh convergence absent'
    return rows,profiles

SLOPE=.025; R_PULSE=70.; PULSE_X=-75.; END_SLOPE=32.

def slope_depth0(x,y):
    rr2=((x-PULSE_X)**2+y*y)/R_PULSE**2
    return .5*np.maximum(1-rr2,0)**.6


def make_slope_engine(R):
    xy,ed=regular_hex_mesh(R,335.)
    depth=cell_averages(xy,R,slope_depth0)
    eng=create_engine(R,xy,ed,depth,-SLOPE*xy[:,0],np.full(len(xy),SF_SLOPE))
    # Match reference-simulation outlet: maximally distant centre from source pulse.
    eng.outlet=int(np.argmax((xy[:,0]-PULSE_X)**2+xy[:,1]**2))
    return eng,xy,ed


def run_slope_time():
    R=10.
    archive=np.load(REF)
    refxy=archive['xy']; reference=archive['reference']
    results=[]; final_by_dt={}; initial_mass=None; tau0=None
    for dt in [.08,.04,.02]:
        eng,xy,ed=make_slope_engine(R)
        assert xy.shape==refxy.shape and np.allclose(xy,refxy,rtol=0,atol=1e-12), 'Slope geometry mismatch stored reference'
        assert eng.N==len(reference)
        m0=float(eng.V_s.sum());initial_mass=m0
        tau,_,_=frozen_coefficient_tau(eng)
        if tau0 is None:tau0=tau
        nsteps=int(round(END_SLOPE/dt))
        assert np.isclose(nsteps*dt,END_SLOPE)
        outlet=0.
        for k in range(nsteps):
            dv,o=eng._route_surface(dt)
            eng.V_s+=dv;outlet+=o
        h=eng.V_s/eng.A_hex
        final_by_dt[dt]=h
        tau_final,_,_=frozen_coefficient_tau(eng)
        l1_ref=100*eng.A_hex*np.abs(h-reference).sum()/m0
        row={'dt_s':dt,'steps':nsteps,'cells':eng.N,'edges':len(ed),
             'mass_closure_fraction':float((eng.V_s.sum()+outlet-m0)/m0),
             'L1_difference_vs_independent_Cartesian_reference_pct_initial_mass':float(l1_ref),
             'com_x_m':float(np.dot(xy[:,0],h)/h.sum()),
             'wet_front_max_x_m_at_1mm':float(np.max(xy[h>.001,0])),
             'minimum_depth_m':float(h.min()),
             'frozen_coefficient_tau_initial_s':tau,'frozen_coefficient_tau_final_s':tau_final,
             'dt_over_initial_tau':dt/tau,'dt_over_final_tau':dt/tau_final}
        print('SLOPE',row,flush=True)
        assert abs(row['mass_closure_fraction'])<1e-11
        assert row['minimum_depth_m']>=-1e-9
        results.append(row)
    err1=100*eng.A_hex*np.abs(final_by_dt[.08]-final_by_dt[.04]).sum()/initial_mass
    err2=100*eng.A_hex*np.abs(final_by_dt[.04]-final_by_dt[.02]).sum()/initial_mass
    trend={'L1_0p08_vs_0p04_pct':float(err1),'L1_0p04_vs_0p02_pct':float(err2),'difference_ratio':float(err1/err2)}
    assert 1.8<trend['difference_ratio']<2.2,'Timestep halving trend not approximately 2'
    return results,trend,final_by_dt,xy,reference


def run_lake_at_rest():
    xy,ed=regular_hex_mesh(15.,335.)
    x,y=xy[:,0],xy[:,1]
    z=.015*x+.15*np.sin(x/65)*np.cos(y/90)
    h=np.maximum(0.-z,0)
    eng=create_engine(15.,xy,ed,h,z,np.full(len(xy),SF_SLOPE))
    nbr=[[] for _ in range(len(xy))]
    for i,j in ed:
        nbr[i].append(j);nbr[j].append(i)
    eng.S_f=np.array([max(np.mean(np.abs(z[i]-z[n]))/eng.d_hex,1e-4) if n else 1e-4 for i,n in enumerate(nbr)])
    eng.outlet=int(np.argmax(z))
    assert h[eng.outlet]==0
    m0=eng.V_s.sum(); initial=eng.V_s.copy()
    tau0,_,_=frozen_coefficient_tau(eng)
    for k in range(100):
        dv,o=eng._route_surface(.1)
        assert abs(o)<1e-12
        eng.V_s+=dv
    res={'cells':eng.N,'wet_cells':int((h>0).sum()),'dry_cells':int((h==0).sum()),
         'maximum_depth_change_m':float(np.abs(eng.V_s-initial).max()/eng.A_hex),
         'Sf_min':float(eng.S_f.min()),'Sf_max':float(eng.S_f.max()),
         'frozen_coefficient_tau_initial_s':tau0,
         'mass_residual_fraction':float((eng.V_s.sum()-m0)/m0),
         'minimum_depth_m':float((eng.V_s/eng.A_hex).min())}
    print('LAKE',res,flush=True)
    assert res['maximum_depth_change_m']<1e-11
    assert abs(res['mass_residual_fraction'])<1e-11
    return res


def main():
    flat,profile=run_barenblatt()
    slope,trend,fs,xy,ref=run_slope_time()
    lake=run_lake_at_rest()
    data={'test_identifier':'Hexalink numerical verification',
        'source_sha256':sha,
        'source_file':'../routing_engine.py',
        'analytical_reference':'2D Barenblatt--Pattle exact solution for h_t=(K/m)*Laplace(h^m); m=8/3, K=1/(n sqrt(Sf))',
        'analytical_reference_doi':'10.1007/s10915-024-02642-x',
        'analytical_parameters':{'m':M,'K':KFLAT,'initial_time_s':T0,'final_time_s':TEND,'central_depth_initial_m':H0,'front_radius_initial_m':R0,'initial_exact_mass_m3':EXACT_MASS},
        'flat_analytical':flat,'slope_numerical':slope,'timestep_refinement':trend,'lake_at_rest':lake,
        'stability_interpretation':{'frozen_coefficient_diagonal_timescale':'tau_i = A_hex / sum_j (Kface M_ij w_ij /d_ij), evaluated at selected states; valid as a fixed-coefficient explicit diagonal nonnegativity bound for the simple closed-source-free linearized update only',
             'positivity':'donor aggregate limiter controls instantaneous outward demand, with minor floating-point residuals; it is not a general nonlinear convergence/stability proof',
             'limitations':'No general CFL guarantee or unconditional stability for steep terrain, strongly nonlinear fronts, runoff boundary outflow or arbitrary dt.'},
        'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':md.version('scipy'),'pandas':md.version('pandas'),'rasterio':md.version('rasterio')},
        'provenance_note':'All hexagonal simulations call HexalinkRoutingEngine._route_surface. The sloping-terrain reference uses an independently discretized Cartesian finite-volume model; its rebuild script and reference arrays are included.'}
    (ROOT/'verification_results.json').write_text(json.dumps(data,indent=2))
    def export(name, rows):
        with (ROOT/name).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0].keys()));writer.writeheader();writer.writerows(rows)
    export('integrated_barenblatt.csv',flat)
    export('integrated_slope_timestep.csv',slope)
    np.savez_compressed(ROOT/'verification_profiles.npz',xy_flat=profile['xy'],numerical_flat=profile['numerical'],exact_flat=profile['exact'],xy_slope=xy,final_slope_dt008=fs[.08],final_slope_dt004=fs[.04],final_slope_dt002=fs[.02],cartesian_slope_reference=ref)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7.5,4.4),constrained_layout=True)
    ax.plot([r['radius_m'] for r in flat],[r['relative_L1_error_pct_exact_mass'] for r in flat],'o-',lw=2)
    ax.invert_xaxis();ax.grid(alpha=.25)
    ax.set(xlabel='Hexagon circumradius R (m), finer to right',ylabel='Relative L1 depth error (% of exact mass)',title='2D Barenblatt benchmark — integrated Hexalink engine')
    fig.savefig(ROOT/'integrated_barenblatt_convergence.png',dpi=175)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.5,4.6),constrained_layout=True)
    rr=np.sqrt(np.sum(profile['xy']**2,axis=1));bins=np.arange(0,460,15.)
    indices=np.digitize(rr,bins)-1
    rr0=(bins[1:]+bins[:-1])*.5
    for a,label in [(profile['numerical'],'Integrated Hexalink, R=10 m'),(profile['exact'],'Exact cell averages')]:
        means=np.array([np.mean(a[indices==i]) if (indices==i).any() else np.nan for i in range(len(rr0))])
        ax.plot(rr0,means,'o-' if 'Integrated' in label else '-',ms=3,label=label)
    ax.set(xlim=(0,350),xlabel='Radius (m)',ylabel='Surface depth (m)',title='2D nonlinear diffusion: analytical comparison')
    ax.grid(alpha=.25);ax.legend(frameon=False)
    fig.savefig(ROOT/'integrated_barenblatt_radial.png',dpi=175)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4),constrained_layout=True)
    dtlist=[r['dt_s'] for r in slope]
    l1=[r['L1_difference_vs_independent_Cartesian_reference_pct_initial_mass'] for r in slope]
    ax.plot(dtlist,l1,'o-',lw=2,label='vs independent 2-m Cartesian reference')
    ax.invert_xaxis();ax.grid(alpha=.25)
    ax.set(xlabel='Timestep (s), finer to right',ylabel='L1 depth difference (% of initial mass)',title='Integrated Hexalink — sloping terrain temporal sensitivity')
    fig.savefig(ROOT/'integrated_slope_timestep.png',dpi=175)
    plt.close(fig)
    print('VERIFICATION PASS',flush=True)

if __name__=='__main__':main()
