#!/usr/bin/env python3
"""Replica-level phosphate-site contact diagnostics for saved v2 trajectories."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from casein_v2 import sequence_features,NCHAINS,BOX,BRIDGE_LEVELS


def phospho_contacts(x,site_mask,cutoff=1.75,box=BOX):
    b=len(site_mask)
    ids=np.array([c*b+i for c in range(NCHAINS) for i,p in enumerate(site_mask) if p>0])
    ch=ids//b
    d=x[ids,None,:]-x[None,ids,:]
    d-=box*np.rint(d/box)
    dist=np.linalg.norm(d,axis=2)
    selected=np.triu(ch[:,None]!=ch[None,:],k=1)
    return int(np.count_nonzero(dist[selected]<cutoff)), float(np.mean(dist[selected]))


def analyze(result_dir, cutoff=1.75):
    root=Path(result_dir)
    _,_,_,sites,_,_=sequence_features()
    reps=[]
    for condition in BRIDGE_LEVELS:
        files=sorted(root.glob(f'{condition}_rep*_trajectory.npz'))
        if not files:raise FileNotFoundError(f'Missing {condition} trajectories in {root}')
        for number,path in enumerate(files,1):
            with np.load(path) as f:xs=f['xyz_reduced']
            late=xs[max(1,int(len(xs)*.55)):]
            contacts=[];distances=[]
            for x in late:
                c,d=phospho_contacts(x,sites,cutoff)
                contacts.append(c);distances.append(d)
            reps.append({'condition':condition,'replica':number,'mean_phosphosite_interchain_contacts':float(np.mean(contacts)),
                         'mean_phosphosite_interchain_pair_distance_reduced':float(np.mean(distances)),
                         'late_frames':len(late)})
    with (root/'phosphosite_contacts_by_replica.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(reps[0]));writer.writeheader();writer.writerows(reps)
    summary=[]
    for condition in BRIDGE_LEVELS:
        seq=np.array([row['mean_phosphosite_interchain_contacts'] for row in reps if row['condition']==condition])
        ds=np.array([row['mean_phosphosite_interchain_pair_distance_reduced'] for row in reps if row['condition']==condition])
        summary.append({'condition':condition,'n_independent_trajectories':len(seq),'phosphosite_contacts_mean':float(seq.mean()),
                        'phosphosite_contacts_replica_sd':float(seq.std(ddof=1)) if len(seq)>1 else None,
                        'phosphosite_mean_pair_distance_reduced':float(ds.mean()),
                        'phosphosite_pair_distance_replica_sd':float(ds.std(ddof=1)) if len(ds)>1 else None,
                        'contact_cutoff_reduced':cutoff})
    (root/'phosphosite_contact_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(7.8,4.7))
    names=[row['condition'] for row in summary]
    means=[row['phosphosite_contacts_mean'] for row in summary]
    sds=[row['phosphosite_contacts_replica_sd'] or 0 for row in summary]
    ax.bar(names,means,yerr=sds,alpha=.7,capsize=6)
    for n,cond in enumerate(names):
        for row in reps:
            if row['condition']==cond:
                ax.scatter(n+.04*(row['replica']-2),row['mean_phosphosite_interchain_contacts'],s=28,zorder=5)
    ax.set(xlabel='Implicit Ca-mediated bridge attraction (proxy)',ylabel=f'Between-chain phosphosite contact pairs (< {cutoff} reduced)',
           title='Localized phosphate-rich contacts: 3 independent replica means')
    fig.tight_layout();fig.savefig(root/'05_phosphosite_contact_replica_sd.png',dpi=160);plt.close(fig)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,default=Path('results'));p.add_argument('--cutoff',type=float,default=1.75)
    args=p.parse_args()
    for row in analyze(args.results,args.cutoff):print(row)
