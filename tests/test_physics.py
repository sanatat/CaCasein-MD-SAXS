"""Smoke/physics checks; run pytest -q."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from cacasein_md_saxs import (debye_curve, pair_forces,initial_positions,
  cluster_scattering,contacts_clusters,integrate_baoab,radius_gyration)


def test_debye_zero_and_translation_invariance():
    x=np.array([[0.,0,0],[1.,0,0],[0,2.,0]])
    q=np.array([0.,.2,2.])
    curve=debye_curve(x,q)
    assert np.isclose(curve[0],9.)
    assert np.allclose(curve,debye_curve(x+np.array([20,5,-1]),q))
    assert np.all(curve>=-1e-9)


def test_force_matches_gradient():
    xyz=initial_positions(42)
    # perturb a single point so bond and interchain derivatives both matter
    xyz[17,0]+=0.08
    f,u=pair_forces(xyz,0.7,30.,6,100.)
    d=1.e-5
    xp=xyz.copy();xm=xyz.copy()
    xp[17,0]+=d;xm[17,0]-=d
    up=pair_forces(xp,.7,30.,6,100.)[1]
    um=pair_forces(xm,.7,30.,6,100.)[1]
    assert np.isclose(f[17,0],-(up-um)/(2*d),atol=1.e-3)


def test_cluster_dilute_limit_and_periodic_image():
    xyz=initial_positions(12)
    # PBC contact case: move two chains near opposite edges
    xyz[:6]=xyz[:6]-np.array([14.,14.,14.])
    xyz[6:12]=xyz[:6]+np.array([29.8,0,0])
    groups,shift=contacts_clusters(xyz)
    assert any({0,1}.issubset(g) for g in groups)
    q=np.array([0.,.2,1.])
    intens,hist,gr,sh=cluster_scattering(xyz,q)
    assert np.isclose(intens[0],6*sum(len(g)**2 for g in gr)/12)


def test_integrator_short_run_finite():
    t,x,u=integrate_baoab(initial_positions(1),.4,2,100,20,.004,1.,1.,30.,6)
    assert x.shape==(6,72,3)
    assert np.isfinite(x).all() and np.isfinite(u).all()
    assert np.isclose(t[-1],.4)
