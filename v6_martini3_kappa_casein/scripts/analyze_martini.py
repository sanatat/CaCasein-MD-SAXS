"""Post-run CG backbone diagnostics; q-profile is only a geometric BB proxy.

Intentionally no quantitative SAXS comparison: coarse-grain beads require
appropriate electron contrast and hydration treatment for physical I(q).
"""
import argparse
import csv
from pathlib import Path
import numpy as np


def rg(coords):
    a=np.asarray(coords,dtype=float)
    if a.ndim!=2 or a.shape[1]!=3 or len(a)<2:
        raise ValueError('Expected N x 3 coordinates')
    return float(np.sqrt(np.mean(np.sum((a-a.mean(axis=0))**2,axis=1))))


def debye_proxy(coords,q):
    xyz=np.asarray(coords,dtype=float)
    q=np.asarray(q,dtype=float)
    diff=xyz[:,None,:]-xyz[None,:,:]
    dist=np.linalg.norm(diff,axis=-1)
    # Sum over all equivalent backbone point scatters, normalized by N squared.
    return np.sinc(q[:,None,None]*dist[None,:,:]/np.pi).sum(axis=(1,2))/(len(xyz)**2)


def analyze(tpr,xtc,out):
    import MDAnalysis as mda
    u=mda.Universe(tpr,xtc)
    backbone=u.select_atoms('name BB')
    if backbone.n_atoms < 20:
        raise RuntimeError('Fewer than 20 BB beads; is this a Martini protein topology?')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    q=np.linspace(0,3,121)
    rows=[]; profiles=[];frames=0
    for ts in u.trajectory:
        # Explicit guard: do not report values if the protein is fragmented by PBC.
        xyz=backbone.positions.astype(float)*0.1 # Angstrom -> nm
        if np.max(np.linalg.norm(np.diff(xyz,axis=0),axis=1))>1.6:
            raise RuntimeError('Protein may cross PBC or have broken backbone; unwrap/center trajectory before analysis')
        val=rg(xyz)
        rows.append([float(ts.time),val,float(np.linalg.norm(xyz[-1]-xyz[0]))])
        profiles.append(debye_proxy(xyz,q))
        frames+=1
    if not frames: raise RuntimeError('Empty trajectory')
    with (out/'rg_ree_by_frame.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['time_ps','Rg_nm','end_to_end_nm']);w.writerows(rows)
    with (out/'bb_geometry_debye_proxy.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['q_nm_inv','mean_normalized_I_BB_proxy']);w.writerows(zip(q,np.mean(profiles,axis=0)))
    print(f'Analyzed {frames} frames and {backbone.n_atoms} backbone beads; Rg {np.mean([r[1] for r in rows]):.3f} nm')
    print('NOTE: This is NOT experimental or quantitatively calibrated SAXS.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--tpr',required=True);p.add_argument('--xtc',required=True);p.add_argument('--out',default='results')
    a=p.parse_args();analyze(a.tpr,a.xtc,a.out)
