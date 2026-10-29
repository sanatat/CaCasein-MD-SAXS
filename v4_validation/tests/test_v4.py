import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from casein_v4 import bulk_initial_state,form_factor_debye,scattering_unwrapped_clusters,NPROT,BOX

def test_bulk_initial_reproducible():
    x=bulk_initial_state(42,8)
    assert np.array_equal(x,bulk_initial_state(42,8))
    assert len(x)==NPROT+8 and np.all(np.isfinite(x))

def test_bulk_no_steric_overlaps():
    x=bulk_initial_state(42,16)
    for i in range(NPROT,len(x)):
        d=x[:i]-x[i];d-=BOX*np.rint(d/BOX)
        assert np.linalg.norm(d,axis=1).min()>1.15-1e-10

def test_debye_q0_is_square_of_sum_weights():
    xyz=np.array([[0.,0,0],[1.,0,0],[1.,2,0]])
    w=np.array([1.,2.,3.]);q=np.array([0.,.1,1.])
    out=form_factor_debye(xyz,w,q)
    assert out[0]==pytest.approx(36.)
    assert out[1]<out[0]

def test_debye_translation_rotation_invariance():
    x=np.array([[0.,0.,0.],[1.,0.,0.],[0.,2.,0.]])
    w=np.array([1.,1.,2.]);q=np.geomspace(.1,2.,12)
    a=form_factor_debye(x,w,q)
    y=x @ np.array([[0.,1.,0.],[-1.,0.,0.],[0.,0.,1.]]) + [11.,13.,17.]
    assert np.allclose(a,form_factor_debye(y,w,q))

def test_debye_permutation():
    rng=np.random.default_rng(55);x=rng.normal(size=(7,3));w=np.arange(1.,8.)
    p=[4,2,5,1,3,0,6];q=np.geomspace(.01,4.,20)
    assert np.allclose(form_factor_debye(x,w,q),form_factor_debye(x[p],w[p],q))
