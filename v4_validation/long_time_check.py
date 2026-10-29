#!/usr/bin/env python3
"""Length-sensitivity check; short reduced-unit trajectories, not proof of equilibrium."""
import argparse,sys,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from casein_v4 import bulk_initial_state,initial_state, integrate, calcium_metrics, sequence_features, NCHAINS,NBEADS,NPROT,DT, COORD_EPS, BOX

def run(out, steps=12000, stride=120, replicas=3, seed=20261029):
    out=Path(out);out.mkdir(exist_ok=True,parents=True)
    _,h,z,p,w,_=sequence_features()
    rows=[];blockrows=[]
    h0=np.r_[np.tile(h,NCHAINS),np.zeros(8)]
    z0=np.r_[np.tile(z,NCHAINS),np.full(8,2.)]
    p0=np.r_[np.tile(p,NCHAINS),np.zeros(8)]
    for mode in ('bulk','near'):
        for rep in range(replicas):
            sd=seed+rep*271
            x0=bulk_initial_state(sd,8) if mode=='bulk' else initial_state(sd,8)
            tr,U=integrate(x0,h0,z0,p0,NPROT,NBEADS,COORD_EPS,sd+31,steps,stride)
            vals=[]
            for frame,x in enumerate(tr):
                occ,br,frac=calcium_metrics(x)
                vals.append(occ)
                rows.append(dict(mode=mode,replica=rep+1,time_reduced=frame*stride*DT,
                                 near_ca_count=occ,bridge_ca_count=br))
            for block in range(3):
                segment=np.array(vals)[int(len(vals)*block/3):int(len(vals)*(block+1)/3)]
                blockrows.append(dict(mode=mode,replica=rep+1,block=block+1,
                                      calcium_contacts_mean=float(segment.mean())))
            print(f'{mode} rep{rep+1}: block means = '+str([round(r['calcium_contacts_mean'],3) for r in blockrows if r['mode']==mode and r['replica']==rep+1]),flush=True)
    for name,data in [('timeseries.csv',rows),('thirds_by_replica.csv',blockrows)]:
        with (out/name).open('w',newline='') as f:
            d=csv.DictWriter(f,fieldnames=list(data[0]));d.writeheader();d.writerows(data)
    plt.rcParams.update({'figure.dpi':140,'axes.spines.right':False,'axes.spines.top':False})
    fig,axs=plt.subplots(1,2,figsize=(12,4.5),sharey=True)
    for ax,mode in zip(axs,['bulk','near']):
        for rep in range(replicas):
            z=[r for r in rows if r['mode']==mode and r['replica']==rep+1]
            ax.plot([r['time_reduced'] for r in z], [r['near_ca_count'] for r in z],
                    label=f'Replica {rep+1}',lw=1.3)
        ax.set(title=f'8 Ca, {mode} start',xlabel='Reduced time (not physical ns)',ylabel='Ca within phospho cutoff')
        ax.legend(fontsize=8)
    fig.suptitle('Time-extension diagnostic: no reliable equilibrium sampling guaranteed')
    fig.tight_layout();fig.savefig(out/'05_long_time_sampling.png');plt.close(fig)
    summary={m:{'third'+str(b):{'replica_mean':float(np.mean([r['calcium_contacts_mean'] for r in blockrows if r['mode']==m and r['block']==b])),
                                'replica_sd':float(np.std([r['calcium_contacts_mean'] for r in blockrows if r['mode']==m and r['block']==b],ddof=1))}
                 for b in (1,2,3)}for m in ('bulk','near')}
    summary['steps']=steps;summary['replicas']=replicas;summary['note']='toy model; nonequilibrium initialization; reduced units only'
    (out/'long_time_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',default=HERE/'long_time_results')
    parser.add_argument('--steps',type=int,default=12000)
    parser.add_argument('--stride',type=int,default=120)
    parser.add_argument('--replicas',type=int,default=3)
    a=parser.parse_args()
    if a.steps%a.stride or a.replicas<2:parser.error('steps divisible by stride and >=2 replicas')
    print(json.dumps(run(a.out,a.steps,a.stride,a.replicas),indent=2))
