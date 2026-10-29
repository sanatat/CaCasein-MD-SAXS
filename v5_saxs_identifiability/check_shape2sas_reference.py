#!/usr/bin/env python3
"""Optional comparison with externally GENERATED Shape2SAS sphere output.

The official Shape2SAS executable is NOT bundled, and was NOT run as part
of the supplied V5 results. Generate an official Iq_*.dat first.
"""
from __future__ import annotations
import argparse, json
import numpy as np
from saxs_core import sphere_form_factor


def load_two_columns(path):
    rows=[]
    with open(path, encoding='utf-8',errors='replace') as fh:
        for line in fh:
            clean=line.strip()
            if not clean or clean.startswith(('#',';')):continue
            columns=clean.replace(',',' ').split()
            try:
                q=float(columns[0]);y=float(columns[1])
            except (IndexError,ValueError):continue
            if q>0 and np.isfinite(q) and np.isfinite(y):rows.append((q,y))
    if len(rows)<8:raise ValueError('At least eight (q,I) numeric rows are needed')
    a=np.asarray(rows,float)
    ix=np.argsort(a[:,0]);return a[ix,0],a[ix,1]


def compare(path,radius_angstrom,qmax_times_radius=2.8):
    q,I=load_two_columns(path)
    keep=(q*radius_angstrom<=qmax_times_radius)&(I>0)
    if keep.sum()<8:raise ValueError('Insufficient values in low-to-mid q range')
    q=q[keep];I=I[keep]
    expected=sphere_form_factor(q,radius_angstrom)
    # Shape2SAS may output arbitrary intensity units, fit one global scale.
    scale=float(np.dot(expected,I)/np.dot(expected,expected))
    norm=I/scale
    return {'source_file':str(path),'radius_angstrom':float(radius_angstrom),
            'points_compared':int(len(q)),'q_max_angstrom_inverse':float(q.max()),
            'fitted_global_scale_factor':scale,
            'rms_P_difference':float(np.sqrt(np.mean((norm-expected)**2))),
            'max_absolute_P_difference':float(np.max(np.abs(norm-expected))),
            'source':'EXTERNAL file (not provided by this project); verify from Shape2SAS'}

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--shape2sas-iq',required=True)
    ap.add_argument('--sphere-radius-angstrom',type=float,required=True)
    args=ap.parse_args()
    print(json.dumps(compare(args.shape2sas_iq,args.sphere_radius_angstrom),indent=2))
