import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import explicit_ca_md as m

def make_small_setup(nca=1):
    x=m.initial_state(456,nca)
    _,h,z,p,_,_=m.sequence_features()
    return x, np.r_[np.tile(h,m.NCHAINS),np.zeros(nca)],np.r_[np.tile(z,m.NCHAINS),np.full(nca,2.)],np.r_[np.tile(p,m.NCHAINS),np.zeros(nca)]

def test_forces_antisymmetric_total():
    x,h,z,p=make_small_setup(2)
    f,u=m.forces(x,h,z,p,m.NPROT,m.NBEADS,2.)
    assert np.isfinite(u) and np.isfinite(f).all()
    assert np.max(np.abs(f.sum(axis=0)))<1e-8

def test_analytic_force_matches_numerical_energy_gradient():
    x,h,z,p=make_small_setup(2)
    f,u=m.forces(x,h,z,p,m.NPROT,m.NBEADS,2.)
    eps=1e-5
    for i,k in [(0,0),(4,1),(m.NPROT,2)]:
        xp=x.copy();xm=x.copy()
        xp[i,k]+=eps;xm[i,k]-=eps
        up=m.forces(xp,h,z,p,m.NPROT,m.NBEADS,2.)[1]
        um=m.forces(xm,h,z,p,m.NPROT,m.NBEADS,2.)[1]
        assert np.isclose(f[i,k],-(up-um)/(2*eps),rtol=1e-3,atol=3e-4), (i,k,f[i,k],-(up-um)/(2*eps))

def test_ion_controls():
    x,h,z,p=make_small_setup(0)
    assert m.calcium_metrics(x)==(0,0,0.)
    assert len(x)==m.NPROT
    a,b=make_small_setup(8)[0],make_small_setup(8)[0]
    assert np.allclose(a,b)

def test_saxs_debye_nonnegative():
    x,_,_,_=make_small_setup(1)
    from casein_v2 import independent_cluster_saxs
    y,_,_=independent_cluster_saxs(x[:m.NPROT],np.ones(m.NPROT),np.geomspace(.05,3.,15),m.NCHAINS,m.BOX)
    assert np.all(y>=-1e-10)
    assert np.all(np.isfinite(y))

def test_short_trajectory_finite():
    x,h,z,p=make_small_setup(2)
    traj,U=m.integrate(x,h,z,p,m.NPROT,m.NBEADS,2.,123,80,20)
    assert traj.shape==(5,m.NPROT+2,3)
    assert np.isfinite(traj).all() and np.isfinite(U).all()
