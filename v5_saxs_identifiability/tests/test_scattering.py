from pathlib import Path
import sys
import numpy as np
import pytest
from scipy.spatial.distance import pdist
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from saxs_core import (point_debye,weighted_rg,guinier_rg,random_sphere_points,
                       sphere_form_factor,homometric_pair,pair_distribution,
                       posterior_from_saxs)
from analyze_experimental_saxs import analyze
from check_shape2sas_reference import compare


def test_debye_q0_is_one():
    rng=np.random.default_rng(3)
    assert np.isclose(point_debye(rng.normal(size=(14,3)),[0])[0],1,atol=1e-14)


def test_translation_invariance():
    x=np.random.default_rng(7).normal(size=(12,3))
    q=np.linspace(.01,2.,35)
    np.testing.assert_allclose(point_debye(x,q),point_debye(x+[100.,-21.,42.],q),atol=1e-13)


def test_rotation_invariance():
    x=np.random.default_rng(23).normal(size=(22,3))
    q=np.linspace(.0,3.,18)
    th=.31;rot=np.array([[np.cos(th),-np.sin(th),0],[np.sin(th),np.cos(th),0],[0,0,1]])
    np.testing.assert_allclose(point_debye(x,q),point_debye(x@rot.T,q),atol=1e-14)


def test_reordering_invariance():
    x=np.random.default_rng(5).normal(size=(8,3))
    w=np.array([1,2,3,4,1,1,2,2],float)
    perm=[6,2,1,5,7,3,4,0]
    np.testing.assert_allclose(point_debye(x,np.array([.3,1.,2.]),w),point_debye(x[perm],np.array([.3,1.,2.]),w[perm]))


def test_two_bead_formula():
    x=np.array([[0,0,0],[3,0,0]],float)
    q=np.array([0.,.2,1.,4.])
    expected=.5*(1+np.sinc(3*q/np.pi))
    np.testing.assert_allclose(point_debye(x,q),expected,atol=1e-13)


def test_homometry_exact():
    x,y=homometric_pair()
    np.testing.assert_array_equal(np.sort(pdist(x)),np.sort(pdist(y)))
    np.testing.assert_allclose(point_debye(x,np.geomspace(.001,5.,100)),
                               point_debye(y,np.geomspace(.001,5.,100)),atol=1e-14)


def test_homometry_not_congruent_in_1d():
    x,y=homometric_pair()
    xx=x[:,0];yy=y[:,0]
    assert not np.array_equal(xx,yy)
    assert not np.array_equal(np.sort(xx.max()-xx),yy)


def test_homometric_rg_and_pair_hist():
    x,y=homometric_pair()
    assert abs(weighted_rg(x)-weighted_rg(y))<1e-12
    assert np.array_equal(pair_distribution(x,bins=18,rmax=18)[1],pair_distribution(y,bins=18,rmax=18)[1])


def test_analytic_sphere_at_zero():
    q=np.array([0,1e-9,.1])
    assert sphere_form_factor(q,10)[0]==pytest.approx(1)
    assert sphere_form_factor(q,10)[1]==pytest.approx(1,abs=1e-12)


def test_discretized_sphere_vs_analytic():
    x=random_sphere_points(600,10,42)
    q=np.linspace(.01,.1,18)
    assert np.sqrt(np.mean((point_debye(x,q)-sphere_form_factor(q,10))**2))<.05


def test_guinier_agrees_with_geometry():
    x=random_sphere_points(500,7,12)
    q=np.linspace(.004,.055,22)
    rg,_,_=guinier_rg(q,point_debye(x,q),qmax=.055)
    assert abs(rg-weighted_rg(x))/weighted_rg(x)<.03


def test_identical_models_equal_posterior_under_uniform_prior():
    x,y=homometric_pair();q=np.linspace(.03,2.0,80)
    c={'A':point_debye(x,q),'B':point_debye(y,q)}
    z=c['A']+np.random.default_rng(7).normal(0,.04,len(q))
    p,chi=posterior_from_saxs(z,np.full(len(q),.04),c)
    assert p['A']==pytest.approx(.5,abs=1e-12)
    assert chi['A']==pytest.approx(chi['B'],abs=1e-10)


def test_nonuniform_prior_does_not_create_saxs_information():
    x,y=homometric_pair();q=np.linspace(.05,1.3,40)
    c={'A':point_debye(x,q),'B':point_debye(y,q)}
    p,_=posterior_from_saxs(c['A'],np.full(len(q),.01),c,{'A':.8,'B':.2})
    assert p['A']==pytest.approx(.8)


def test_negative_q_rejected():
    with pytest.raises(ValueError):point_debye([[0,0,0]],[-1.])


def test_zero_sum_weights_rejected():
    with pytest.raises(ValueError):point_debye([[0,0,0],[1,0,0]],[.2],weights=[1,-1])


def test_experimental_importer_on_synthetic_not_real(tmp_path):
    q=np.linspace(.007,.055,22);rg=4.
    inten=2*np.exp(-q*q*rg*rg/3.)
    f=tmp_path/'mock.csv'
    f.write_text('q,intensity,sigma\n'+''.join(f'{a},{b},0.01\n' for a,b in zip(q,inten)))
    result=analyze(f,'nm^-1',.055,'SYNTHETIC-TEST-ONLY')
    assert result['guinier_Rg_nm']==pytest.approx(rg,rel=1e-10)
    assert result['guinier_qRg_lte_1.3']


def test_shape2sas_external_file_parser_on_synthetic(tmp_path):
    q=np.linspace(.003,.055,40)
    I=5*sphere_form_factor(q,40.)
    f=tmp_path/'synthetic_Iq.dat'
    f.write_text('#synthetic\n'+''.join(f'{x} {y}\n' for x,y in zip(q,I)))
    result=compare(f,40.)
    assert result['rms_P_difference']<1e-12
