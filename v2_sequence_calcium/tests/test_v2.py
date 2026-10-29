import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from casein_v2 import (sequence_features,PHOSPHOSERINES,MATURE_BETA_CASEIN_A2,
                       total_forces,intensity_debye,init_bundle,NCHAINS,BOX,integrate)
from analyze_local_contacts import phospho_contacts


def test_sequence_phosphoserines():
    segs,h,z,p,w,n=sequence_features()
    assert len(MATURE_BETA_CASEIN_A2)==209
    assert all(MATURE_BETA_CASEIN_A2[i-1]=='S' for i in PHOSPHOSERINES)
    assert len(segs)==21 and sum(n)==5 and np.count_nonzero(p)==2
    assert np.isclose(sum(w),20.9)


def test_qzero_identity_and_translation():
    xyz=np.array([[0,0,0],[1,0,0],[0,1,0]],float)
    ws=np.array([1.,.7,1.1]);q=np.array([0.,.1,1.])
    curve=intensity_debye(xyz,ws,q)
    assert np.isclose(curve[0],ws.sum()**2)
    assert np.allclose(curve,intensity_debye(xyz+np.array([8,2,-4]),ws,q))


def test_gradients_with_calcium_and_screening():
    # Two nonbonded beads from two chains; compare analytic force to finite difference.
    xyz=np.array([[1.,1.,1.],[2.57,1.1,1.]],float)
    h=np.array([.5,.6]);z=np.array([-.5,-.2]);p=np.array([1.,.5])
    for ca in (0.,.7,1.4):
        force,energy=total_forces(xyz,h,z,p,2,ca,box=BOX)
        delta=1e-5
        left=xyz.copy();right=xyz.copy()
        left[0,0]-=delta;right[0,0]+=delta
        grad=(total_forces(right,h,z,p,2,ca,BOX)[1]-total_forces(left,h,z,p,2,ca,BOX)[1])/(2*delta)
        assert np.isclose(force[0,0],-grad,rtol=2e-5,atol=2e-5)
        assert np.allclose(force.sum(axis=0),0,atol=1e-9)


def test_ca_zero_when_no_sites():
    xyz=np.array([[1.,1.,1.],[2.57,1.1,1.]],float)
    h=np.array([.4,.7]);z=np.array([-.4,-.2]);p=np.array([0.,0.])
    f0,u0=total_forces(xyz,h,z,p,2,0.)
    f1,u1=total_forces(xyz,h,z,p,2,1.4)
    assert np.allclose(f0,f1) and u0==u1


def test_integrator_and_phos_contacts_smoke():
    _,h,z,p,_,_=sequence_features()
    xyz=init_bundle(123)
    traj=integrate(xyz,np.tile(h,NCHAINS),np.tile(z,NCHAINS),np.tile(p,NCHAINS),NCHAINS,.7,124,15,5)
    assert traj[1].shape==(4,NCHAINS*len(h),3)
    assert all(np.isfinite(a).all() for a in traj)
    cnt,dist=phospho_contacts(traj[1][-1],p)
    assert cnt>=0 and dist>0
