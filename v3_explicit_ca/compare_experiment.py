#!/usr/bin/env python3
"""EXPLORATORY comparison of experimental 2-column SAXS data and toy-model SAXS.

Not a validation: requires a user-justified physical length conversion; samples differ.
Input CSV/TXT has two columns q(nm^-1 or A^-1), I(q). Optional headers ignored.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read_exp(path,unit):
    pairs=[]
    for ln in Path(path).read_text(errors='replace').splitlines():
        parts=ln.replace(',',' ').replace(';',' ').split()
        if len(parts)<2:continue
        try:q,i=map(float,parts[:2])
        except ValueError:continue
        if q>0 and i>0 and np.isfinite(q+i):pairs.append((q,i))
    if len(pairs)<10:raise ValueError('At least 10 valid (q, I) rows required')
    z=np.array(sorted(pairs));
    if unit=='A^-1':z[:,0]*=10
    return z

def make_overlay(experimental_csv,model_csv,condition,replica,length_nm,out,q_unit='nm^-1'):
    if length_nm<=0:raise ValueError('Provide a POSITIVE independently justified nm per reduced length')
    exp=read_exp(experimental_csv,q_unit)
    with Path(model_csv).open() as f:
        rows=[r for r in csv.DictReader(f) if r['condition']==condition and r['replica']==str(replica)]
    if not rows:raise ValueError('No matching condition/replica in simulation CSV')
    model=np.array([(float(r['q_reduced']),float(r['I_shape_normalized_at_q_min'])) for r in rows]);
    model[:,0]/=length_nm
    # Independent scaling at smallest q; no ungrounded absolute SAXS match implied.
    exnorm=exp[:,1]/exp[0,1]
    overlap=(max(exp[0,0],model[0,0]),min(exp[-1,0],model[-1,0]))
    fig,ax=plt.subplots(figsize=(9,5))
    ax.loglog(exp[:,0],exnorm,lw=1.5,label='Experimental SAXS (from supplied file)')
    ax.loglog(model[:,0],model[:,1],lw=1.5,label=f'EXPLORATORY model: {condition} / replica {replica}')
    ax.set(xlabel='q (nm$^{-1}$; simulated scale ASSUMED)',ylabel='Separate I(q) normalizations',
           title='Visual inspection ONLY — not experimental validation or fitted force field')
    ax.legend();fig.tight_layout();fig.savefig(out,dpi=150);plt.close(fig)
    info={'reference':'EXPLORATORY ONLY','experimental_file':str(experimental_csv),'condition':condition,
          'assumed_nm_per_reduced_length':length_nm,'overlap_nm^-1':overlap,
          'overlap_exists':overlap[0]<overlap[1],
          'warning':'Experiment and computational sample may differ in composition, treatment, resolution, scale, solvent contrast.'}
    Path(out).with_suffix('.json').write_text(json.dumps(info,indent=2)+'\n')
    return info

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--experiment',required=True)
    p.add_argument('--model',default='results/saxs_curve_replicas.csv')
    p.add_argument('--condition',default='ca8_coord',choices=['no_ca','ca8_no_coord','ca8_coord','ca16_coord'])
    p.add_argument('--replica',type=int,default=1)
    p.add_argument('--q-unit',default='nm^-1',choices=['nm^-1','A^-1'])
    p.add_argument('--length-nm-per-reduced',type=float,required=True)
    p.add_argument('--out',default='EXPLORATORY_NOT_VALIDATION.png')
    a=p.parse_args()
    print(json.dumps(make_overlay(a.experiment,a.model,a.condition,a.replica,a.length_nm_per_reduced,a.out,a.q_unit),indent=2))
