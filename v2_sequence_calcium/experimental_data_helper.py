#!/usr/bin/env python3
"""Download external published casein SAXS archive or plot two-column data.

This helper DOES NOT claim experimental validation, since the 2026 Zenodo data
come from enzyme-driven gelation of whole skim milk, not isolated beta-casein
or a matched calcium titration. Explicit reduced->nm length conversion required.
"""
import argparse
import hashlib
from pathlib import Path
from urllib.request import urlopen
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ZENODO_URL='https://zenodo.org/records/22044394/files/substracted_merged_USAXS_SAXS.7z?download=1'
KNOWN_MD5='652981bb9f72efbf7949127b77393826'


def fetch(out,timeout=30):
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True)
    with urlopen(ZENODO_URL,timeout=timeout) as stream, out.open('wb') as f:
        while True:
            chunk=stream.read(1024*1024)
            if not chunk:break
            f.write(chunk)
    digest=hashlib.md5(out.read_bytes()).hexdigest()
    if digest!=KNOWN_MD5:
        out.unlink(missing_ok=True)
        raise RuntimeError(f'Download checksum mismatch: {digest}, expected {KNOWN_MD5}')
    return out


def read_numeric(path, q_unit='nm^-1'):
    # Handle tab, comma and whitespace, skip headers/comments.
    data=[]
    for ln in Path(path).read_text(errors='replace').splitlines():
        s=ln.strip().replace(',',' ').replace(';',' ')
        if not s or s.startswith(('#',';','%')):continue
        fields=s.split()
        if len(fields)<2:continue
        try:
            q,i=float(fields[0]),float(fields[1])
        except ValueError:continue
        if np.isfinite(q) and np.isfinite(i) and q>0 and i>0:data.append((q,i))
    if len(data)<10:raise ValueError('Need at least ten positive numeric q, I rows')
    arr=np.array(sorted(data))
    if q_unit=='A^-1':arr[:,0]*=10.0
    return arr


def plot_overlay(datafile,modelcsv,length_nm_per_reduced,condition,unit,out):
    if length_nm_per_reduced<=0:raise ValueError('length_nm_per_reduced must be positive')
    exp=read_numeric(datafile,unit)
    with Path(modelcsv).open() as f:
        rows=list(csv.DictReader(f))
    m=np.array([(float(r['q_reduced']),float(r['I_q_per_bead_scattering_weight']))
               for r in rows if r['condition']==condition and r['replica']=='1'])
    qmodel_nm=m[:,0]/length_nm_per_reduced
    # Per curve shape normalization (separately): NEVER pretend common absolute scale.
    exnorm=exp[:,1]/exp[0,1]
    modnorm=m[:,1]/m[0,1]
    fig,ax=plt.subplots(figsize=(8.8,5))
    ax.loglog(exp[:,0],exnorm,label='EXPERIMENT: enzyme-treated skim milk (not matched)',lw=1.7)
    ax.loglog(qmodel_nm,modnorm,label=f'SYNTHETIC: beta-casein toy model ({condition})',lw=1.7)
    ax.set(xlabel=r'$q$ (nm$^{-1}$; simulated q uses ASSUMED bead length)',
           ylabel='Each profile separately normalized at its first q',
           title='Exploratory overlay ONLY — not a validated structural fit')
    ax.legend();fig.tight_layout();fig.savefig(out,dpi=160);plt.close(fig)
    return {'experimental_q_range_nm_inv':[float(exp[0,0]),float(exp[-1,0])],
            'assumed_nm_per_reduced_length':length_nm_per_reduced,
            'overlap_q_range_exists':bool(max(exp[0,0],qmodel_nm[0])<min(exp[-1,0],qmodel_nm[-1]))}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    g=sub.add_parser('fetch');g.add_argument('--out',default='external/zenodo_casein_saxs_2026.7z')
    v=sub.add_parser('overlay');v.add_argument('--experiment',required=True)
    v.add_argument('--model',default='results/saxs_curve_replicas.csv')
    v.add_argument('--condition',choices=['off','low','high'],default='off')
    v.add_argument('--q-unit',choices=['nm^-1','A^-1'],default='nm^-1')
    v.add_argument('--length-nm-per-reduced',type=float,required=True,
                   help='Externally justified physical bead scale in nm. NO MODEL DEFAULT EXISTS.')
    v.add_argument('--out',default='exploratory_overlay_NOT_VALIDATION.png')
    args=p.parse_args()
    if args.command=='fetch':print(f'Downloaded verified archive: {fetch(args.out)}')
    else:print(plot_overlay(args.experiment,args.model,args.length_nm_per_reduced,args.condition,args.q_unit,args.out))

if __name__=='__main__':main()
