#!/usr/bin/env python3
"""Careful, standalone Guinier exploration for externally provided 3-column SAXS.

This tool does NOT compare uncalibrated CaCasein simulations to experiments.
A user must supply q-unit and material provenance, and choose a Guinier window.
No automatic baseline subtraction, aggregation correction or data fitting to casein.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from saxs_core import guinier_rg


def read_saxs_csv(path,delimiter=','):
    with open(path,encoding='utf-8-sig',newline='') as fh:
        reader=csv.DictReader(fh,delimiter=delimiter)
        if not reader.fieldnames:raise ValueError('Missing header')
        required={'q','intensity','sigma'}
        if not required.issubset(set(reader.fieldnames)):
            raise ValueError('CSV headers must be q,intensity,sigma')
        rows=[]
        for row in reader:
            try: vals=tuple(float(row[k]) for k in ('q','intensity','sigma'))
            except (ValueError,TypeError):continue
            if np.isfinite(vals).all() and vals[0]>0 and vals[1]>0 and vals[2]>0:rows.append(vals)
    if len(rows)<6:raise ValueError('Need at least 6 positive finite intensity rows')
    a=np.asarray(sorted(rows,key=lambda x:x[0]));return a[:,0],a[:,1],a[:,2]


def analyze(path,q_unit,guinier_qmax,source_label):
    if not source_label.strip():raise ValueError('Source label is required')
    q,I,s=read_saxs_csv(path)
    converted=q if q_unit=='nm^-1' else q*10.0
    cutoff=guinier_qmax if q_unit=='nm^-1' else guinier_qmax*10.
    rg_nm,I0,n=guinier_rg(converted,I,cutoff,s)
    diag=rg_nm*cutoff
    return {'source':source_label,'n_total_valid':len(q),'n_guinier_window':n,
            'guinier_Rg_nm':rg_nm,'guinier_I0':I0,'qmax_times_Rg':diag,
            'guinier_qRg_lte_1.3':bool(diag<=1.3),
            'warning':'Exploratory only: validate background, smearing, interactions, low-q linearity, source composition and uncertainty correlations. Not casein model validation.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--csv',required=True,type=Path)
    ap.add_argument('--q-unit',required=True,choices=['nm^-1','A^-1'])
    ap.add_argument('--guinier-qmax',required=True,type=float,help='same unit as the input q')
    ap.add_argument('--source',required=True,help='accurate dataset/source provenance')
    args=ap.parse_args()
    print(json.dumps(analyze(args.csv,args.q_unit,args.guinier_qmax,args.source),indent=2))
